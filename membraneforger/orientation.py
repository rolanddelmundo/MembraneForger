#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Membrane orientation: one rigid transform, derived from the anchor chain(s), applied to the whole complex.
#//=============================================================
"""Orient the all-atom complex in the membrane frame from a user-selected anchor chain.

Convention of the MembraneForger orientation frame (the OPM/PPM convention, in A):
    membrane normal   -> +Z
    membrane midplane -> z = 0
    IN  (cytoplasmic, PPM "in")  -> negative z
    OUT (extracellular, PPM "out") -> positive z

The anchor chain(s) are oriented by a provider (an exact OPM reference entry, or a local PPM 3.0 run on the
anchor coordinates themselves); the rigid transform taking the original anchor onto its oriented copy is derived
once and applied to EVERY atom of the original complex. Nothing is ever rebuilt from a provider's coordinate file.
"""
import math
import os
import re
import shutil
import time
import urllib.error
import urllib.request
from collections import OrderedDict
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .alignment import align_sequences, kabsch, rmsd
from .config import AMINO, ONE_LETTER, Settings
from .runtools import log, run_command, sha256
from .structio import read_pdb, residues_in_order, source_element, write_pdb, xyz_nm

__all__ = ['ORIENTATION_MODES', 'NTERM_SIDES', 'OPM_ASSET_URL', 'PPM_MEMBRANE_CODES', 'OrientationRequest',
           'OrientationFailure', 'parse_residue_ranges', 'select_residues', 'embedded_segment', 'protein_chains', 'chain_sequence',
           'pdb_id_from_header', 'resolve_pdb_id',
           'select_anchor_chains', 'OPMReference', 'parse_opm_file', 'fetch_opm_reference', 'match_chains',
           'atom_pairs', 'rigid_fit', 'apply_transform', 'spanning_chains', 'nterm_side_from_reference',
           'orientation_from_opm', 'find_ppm_executable', 'LocalPPM', 'orientation_from_ppm',
           'validate_full_transform', 'orient_complex', 'check_orientation_preserved', 'rotation_angle_deg']

ORIENTATION_MODES = ("auto", "ppm", "opm", "none")
NTERM_SIDES = ("auto", "in", "out")
# OPM publishes every oriented entry (the same coordinates OPRLM serves) in this public bucket, keyed by lower-case PDB ID.
OPM_ASSET_URL = "https://opm-assets.storage.googleapis.com/pdb/{pdb_id}.pdb"
# PPM 3.0 membrane codes (ppm3_instructions). "" is PPM's "undefined membrane": a flat bilayer whose hydrophobic
# thickness is optimised, which is the only choice that assumes nothing about the protein's biological membrane.
PPM_MEMBRANE_CODES = {"", "PMm", "PMp", "PMf", "Erf", "ERm", "GOL", "LYS", "END", "VAC", "MOM", "MIM", "THp", "THb",
                      "GnO", "GnI", "GpI", "ARC", "LPC", "MPC", "OPC", "EPC", "MIC"}
BACKBONE = ("N", "CA", "C", "O")
PPM_INPUT_NAME = "anchor.pdb"  # PPM derives its output name from this: anchorout.pdb
# With --orientation none, --orient-residues names the transmembrane segment; protein CA outside it within this many
# A of the bilayer centre it gives (half the 3.0 nm hydrophobic belt of embedding.BELT_NM) are reported.
EMBEDDED_SEGMENT_HALF_THICKNESS_A = 15.0
# A residue selection submitted to PPM must hold at least one transmembrane helix (about 20 residues); fewer than
# this many residues per anchor chain is refused rather than oriented.
MIN_SEGMENT_RESIDUES = 10


@dataclass(frozen=True)
class OrientationRequest:
    """What the user asked for; the CLI fills it from the configuration header and the command line."""
    mode: str = "auto"              # auto | ppm | opm | none
    chains: tuple = ()              # anchor chain label(s); empty = decide automatically or fail
    nterm_side: str = "auto"        # auto | in | out: side of the N terminus of the first anchor chain
    pdb_id: str | None = None       # optional exact PDB ID for the OPM reference
    ppm_exe: str | None = None      # PPM 3.0 executable (default: $MEMBRANEFORGER_PPM, then immers/ppm3 on PATH)
    ppm_membrane: str = ""          # PPM membrane code; "" = undefined flat bilayer
    ppm_heteroatoms: bool = False   # submit the anchor chains' heteroatoms to PPM
    opm_file: Path | None = None    # an already downloaded OPM/OPRLM coordinate file (offline use)
    opm_cache: Path | None = None   # where downloaded OPM files are kept (default ~/.cache/membraneforger/opm)
    residues: tuple = ()            # ((first, last), ...) residue numbers of the anchor chain(s) submitted to PPM; () = all


class OrientationFailure(SystemExit):
    """Orientation could not be established safely; the message says what the user must supply."""


def parse_residue_ranges(text: str | None) -> tuple:
    """Turn "343-363" or "343-363,370-380" or "350" into ((343, 363), ...); the empty string selects everything."""
    ranges = []
    for part in (text or "").replace(" ", ",").split(","):
        if not part:
            continue
        found = re.fullmatch(r"(-?\d+)(?:-(-?\d+))?", part)
        if not found:
            raise OrientationFailure(f"--orient-residues {text!r}: {part!r} is not a residue number or a first-last range")
        first, last = int(found.group(1)), int(found.group(2) if found.group(2) is not None else found.group(1))
        if last < first:
            raise OrientationFailure(f"--orient-residues {text!r}: range {part} ends before it starts")
        ranges.append((first, last))
    return tuple(ranges)


def select_residues(anchors: OrderedDict, ranges: tuple) -> OrderedDict:
    """Keep the anchor residues whose numbers fall in the ranges; every anchor chain must keep a helix-sized segment."""
    if not ranges:
        return anchors
    inside = lambda resid: any(first <= resid <= last for first, last in ranges)
    selected = OrderedDict((label, [r for r in residues if inside(r[0])]) for label, residues in anchors.items())
    short = {label: len(residues) for label, residues in selected.items() if len(residues) < MIN_SEGMENT_RESIDUES}
    if short:
        spans = ", ".join(f"{first}-{last}" for first, last in ranges)
        raise OrientationFailure(f"--orient-residues {spans} selects too few residues of anchor chain(s) " + ", ".join(
            f"{label} ({n}; residues {anchors[label][0][0]}-{anchors[label][-1][0]} are present)" for label, n in short.items())
            + f"; at least {MIN_SEGMENT_RESIDUES} are needed. The numbers are those of the input file.")
    return selected


def protein_chains(atoms: list[dict]) -> OrderedDict:
    """Group the protein residues of a structure by chain label (chain ID, else segid, else A), in file order."""
    chains: OrderedDict[str, list[dict]] = OrderedDict()
    for a in atoms:
        if a["resname"] in AMINO:
            chains.setdefault(a["chain"] or a.get("segid") or "A", []).append(a)
    return OrderedDict((label, residues_in_order(members)) for label, members in chains.items())


def chain_sequence(residues: list) -> str:
    """One-letter sequence of a residue list from residues_in_order (unknown residues align as X)."""
    return "".join(ONE_LETTER.get(rn, "X") for _, rn, _ in residues)


def pdb_id_from_header(path: Path) -> str | None:
    """Return the PDB ID of a HEADER record (columns 63-66), the only authoritative ID a PDB file carries."""
    with path.open(errors="replace") as fh:
        for line in fh:
            if line.startswith("HEADER"):
                candidate = line[62:66].strip().upper()
                return candidate if re.fullmatch(r"[1-9][A-Z0-9]{3}", candidate) else None
            if line[:6] in ("ATOM  ", "HETATM"):
                break
    return None


def resolve_pdb_id(explicit: str | None, path: Path | None) -> tuple[str | None, str]:
    """Pick the PDB ID: --pdb-id first, then the HEADER record; never the file name."""
    if explicit:
        candidate = explicit.strip().upper()
        if not re.fullmatch(r"[1-9][A-Z0-9]{3}", candidate):
            raise OrientationFailure(f"--pdb-id {explicit!r} is not a PDB ID (four characters, starting with a digit)")
        return candidate, "command line"
    found = pdb_id_from_header(path) if path else None
    return found, "HEADER record" if found else "none"


def select_anchor_chains(chains: OrderedDict, requested: tuple, reference: "OPMReference | None" = None,
                         min_identity: float = 0.95) -> tuple[list[str], str]:
    """Choose the chain(s) that define the membrane orientation, or fail with what to pass instead."""
    labels = list(chains)
    if requested:
        missing = [c for c in requested if c not in chains]
        if missing:
            raise OrientationFailure(f"anchor chain(s) {missing} are not protein chains of the input; protein chains are "
                                     f"{', '.join(labels)}")
        return list(requested), "command line"
    if len(labels) == 1:
        return labels, "the only protein chain"
    if reference is not None:
        spanning = spanning_chains(reference)
        matched = {}
        for label in labels:
            best = match_chains(chains[label], reference.chains, min_identity)
            if best and best["chain"] in spanning:
                matched.setdefault(best["chain"], []).append(label)
        if len(matched) == 1 and len(next(iter(matched.values()))) == 1:
            label = next(iter(matched.values()))[0]
            return [label], f"the only chain that matches a membrane-spanning chain of OPM {reference.pdb_id}"
    raise OrientationFailure(
        f"Multiple protein chains were detected: {', '.join(labels)}.\n\n"
        "MembraneForger cannot safely determine which chain defines the membrane orientation.\n\n"
        f"Re-run with:\n    --orient-chain {labels[0]}\n\nor:\n    --orient-chains {labels[0]},{labels[1]}\n\n"
        "(choose the chain that actually spans or associates with the membrane)")


@dataclass
class OPMReference:
    """An OPM/OPRLM oriented coordinate file: membrane normal along z, midplane at z = 0, half thickness in A."""
    pdb_id: str | None
    path: Path
    source: str
    sha256: str
    half_thickness_a: float
    chains: OrderedDict  # protein residues per chain, DUM atoms excluded


def parse_opm_file(path: Path, pdb_id: str | None, source: str) -> OPMReference:
    """Read and validate an OPM coordinate file; the thickness REMARK and the dummy-atom planes must agree."""
    text = path.read_text(errors="replace")
    found = re.search(r"1/2 of bilayer thickness:\s*([-\d.]+)", text)
    atoms, _ = read_pdb(path)
    dummies = [a for a in atoms if a["resname"] == "DUM"]
    protein = [a for a in atoms if a["resname"] != "DUM"]
    if not found or not protein:
        raise OrientationFailure(f"{path} is not an OPM/OPRLM oriented coordinate file (no thickness REMARK or no atoms)")
    half = float(found.group(1))
    if half <= 0 or not math.isfinite(half):
        raise OrientationFailure(f"{path}: implausible half thickness {half}")
    if dummies:
        planes = sorted({round(a["z"], 1) for a in dummies})
        if planes != sorted({round(half, 1), round(-half, 1)}):
            raise OrientationFailure(f"{path}: dummy-atom planes {planes} do not match the half thickness {half}")
    header = [l for l in text.splitlines() if l.startswith("HEADER")]
    if pdb_id and header and header[0][62:66].strip().upper() not in ("", pdb_id):
        raise OrientationFailure(f"{path} is the OPM entry {header[0][62:66].strip()}, not {pdb_id}")
    if pdb_id is None and header and re.fullmatch(r"[1-9][A-Z0-9]{3}", header[0][62:66].strip().upper()):
        pdb_id = header[0][62:66].strip().upper()
    chains = protein_chains(protein)
    if not chains:
        raise OrientationFailure(f"{path}: no protein residues")
    return OPMReference(pdb_id, path, source, sha256(path), half, chains)


def fetch_opm_reference(pdb_id: str, cache: Path | None, out: Path | None = None) -> OPMReference | None:
    """Download (or reuse from the cache) the OPM entry of a PDB ID; None when OPM has no such entry."""
    cache = Path(cache or Path.home() / ".cache" / "membraneforger" / "opm")
    target = cache / f"{pdb_id.lower()}.pdb"
    url = OPM_ASSET_URL.format(pdb_id=pdb_id.lower())
    if not target.is_file() or target.stat().st_size == 0:
        cache.mkdir(parents=True, exist_ok=True)
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                data = response.read()
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                if out:
                    log(out, f"OPM has no entry for {pdb_id} ({url} -> 404)")
                return None
            raise OrientationFailure(f"OPM lookup of {pdb_id} failed: HTTP {exc.code} from {url}") from None
        except (urllib.error.URLError, OSError, ValueError) as exc:
            raise OrientationFailure(f"OPM lookup of {pdb_id} failed: {exc} ({url}); pass --opm-file with a downloaded "
                                     "OPM/OPRLM coordinate file to work offline") from None
        if b"1/2 of bilayer thickness" not in data:
            raise OrientationFailure(f"{url} did not return an OPM coordinate file")
        part = target.with_suffix(".part")
        part.write_bytes(data)
        os.replace(part, target)
        if out:
            log(out, f"OPM entry {pdb_id} downloaded from {url} ({len(data)} bytes)")
    elif out:
        log(out, f"OPM entry {pdb_id} read from the cache {target}")
    return parse_opm_file(target, pdb_id, url)


def match_chains(residues: list, candidates: OrderedDict, min_identity: float, min_coverage: float = 0.5) -> dict | None:
    """Find the candidate chain whose sequence matches a residue list best, if any passes the identity criteria."""
    seq, best = chain_sequence(residues), None
    for label, other in candidates.items():
        pairs = align_sequences(seq, chain_sequence(other))
        if not pairs:
            continue
        same = sum(seq[i] == chain_sequence(other)[j] for i, j in pairs)
        identity = same / len(pairs)
        coverage = len(pairs) / min(len(seq), len(other))
        if identity >= min_identity and coverage >= min_coverage and (best is None or same > best["same"]):
            best = {"chain": label, "pairs": pairs, "same": same, "identity": identity, "coverage": coverage,
                    "aligned": len(pairs)}
    return best


def atom_pairs(residues: list, others: list, pairs: list, names: tuple | None = BACKBONE, subset: set | None = None):
    """Collect matched atom coordinates (A) of aligned residues by atom name; names=None takes every heavy atom."""
    P, Q, labels = [], [], []
    for i, j in pairs:
        if subset is not None and j not in subset:
            continue
        mine = {a["atom"]: a for a in residues[i][2] if names is not None or source_element(a) != "H"}
        theirs = {a["atom"]: a for a in others[j][2]}
        for name in (names or mine):
            if name in mine and name in theirs:
                P.append((mine[name]["x"], mine[name]["y"], mine[name]["z"]))
                Q.append((theirs[name]["x"], theirs[name]["y"], theirs[name]["z"]))
                labels.append(f"{residues[i][1]}{residues[i][0]}:{name}")
    return np.array(P, dtype=float).reshape(-1, 3), np.array(Q, dtype=float).reshape(-1, 3), labels


def rigid_fit(P: np.ndarray, Q: np.ndarray, trim_floor_a: float | None = None, trim_factor: float = 2.5) -> dict:
    """Kabsch fit of P onto Q, optionally trimming outliers (deviation > max(factor x median, floor)) iteratively."""
    if len(P) < 3:
        raise OrientationFailure(f"only {len(P)} matched atoms; at least 3 are needed for a rigid fit")
    keep = np.ones(len(P), dtype=bool)
    R, t = kabsch(P, Q)
    error = np.linalg.norm(P @ R.T + t - Q, axis=1)
    if trim_floor_a is not None:
        for _ in range(10):
            trimmed = error <= max(trim_factor * float(np.median(error[keep])), trim_floor_a)
            if trimmed.sum() < 3 or (trimmed == keep).all():
                break
            keep = trimmed
            R, t = kabsch(P[keep], Q[keep])
            error = np.linalg.norm(P @ R.T + t - Q, axis=1)
    det = float(np.linalg.det(R))
    if abs(det - 1.0) > 1e-6:
        raise OrientationFailure(f"the rigid fit is not a proper rotation (determinant {det:.6f})")
    return {"R": R, "t": t, "error": error, "keep": keep, "rmsd": rmsd(error), "core_rmsd": rmsd(error[keep]),
            "core_fraction": float(keep.mean()), "max_deviation": float(error.max()), "det": det, "atoms": len(P)}


def apply_transform(atoms: list[dict], R: np.ndarray, t: np.ndarray) -> list[dict]:
    """Return copies of the atoms moved by x' = R x + t; every other field is kept verbatim."""
    moved = xyz_nm(atoms) @ R.T + t
    return [dict(a, x=float(x), y=float(y), z=float(z)) for a, (x, y, z) in zip(atoms, moved)]


def spanning_chains(reference: OPMReference) -> set[str]:
    """Chains of an OPM entry that have CA atoms on both sides of the hydrophobic slab (membrane-spanning)."""
    spanning = set()
    for label, residues in reference.chains.items():
        z = np.array([a["z"] for _, _, atoms in residues for a in atoms if a["atom"] == "CA"])
        if len(z) and (z < -reference.half_thickness_a).any() and (z > reference.half_thickness_a).any():
            spanning.add(label)
    return spanning


def nterm_side_from_reference(reference: OPMReference, chain: str) -> tuple[str | None, str]:
    """Read the N-terminal side of a reference chain from the z of its first CA (OPM: in < 0 < out)."""
    for _, name, atoms in reference.chains[chain]:
        ca = next((a for a in atoms if a["atom"] == "CA"), None)
        if ca is None:
            continue
        if abs(ca["z"]) <= reference.half_thickness_a:
            return None, (f"the first residue {name}{atoms[0]['resid']} of OPM chain {chain} lies inside the hydrophobic "
                          f"slab (z {ca['z']:.1f} A), so the N-terminal side cannot be read from the reference")
        side = "out" if ca["z"] > 0 else "in"
        return side, f"OPM {reference.pdb_id} chain {chain}: first residue {name}{atoms[0]['resid']} at z {ca['z']:+.1f} A"
    return None, f"OPM chain {chain} has no CA atoms"


def orientation_from_opm(anchors: OrderedDict, reference: OPMReference, settings: Settings) -> dict:
    """Derive the transform taking the user's anchor onto the OPM reference, fitted on the membrane-embedded core."""
    # The orientation of a membrane protein is set by its transmembrane core; soluble domains of the same entry may
    # have moved in the user's model. So the fit uses the backbone of reference residues inside the hydrophobic slab
    # and reports the whole-chain deviation separately. A reference whose embedded core does not superpose is refused.
    matches, P, Q, Pall, Qall = [], [], [], [], []
    for label, residues in anchors.items():
        best = match_chains(residues, reference.chains, settings.orient_min_identity)
        if best is None:
            raise OrientationFailure(f"anchor chain {label} has no chain of OPM {reference.pdb_id} with >= "
                                     f"{settings.orient_min_identity:.0%} sequence identity")
        other = reference.chains[best["chain"]]
        slab = {j for j, (_, _, atoms) in enumerate(other)
                if any(a["atom"] == "CA" and abs(a["z"]) <= reference.half_thickness_a for a in atoms)}
        p, q, _ = atom_pairs(residues, other, best["pairs"], BACKBONE, slab)
        pa, qa, _ = atom_pairs(residues, other, best["pairs"], BACKBONE)
        matches.append({"anchor_chain": label, "reference_chain": best["chain"], "aligned_residues": best["aligned"],
                        "identity": round(best["identity"], 4), "coverage": round(best["coverage"], 4),
                        "slab_residues_matched": len({j for _, j in best["pairs"] if j in slab}),
                        "slab_residues_in_reference": len(slab)})
        P.append(p), Q.append(q), Pall.append(pa), Qall.append(qa)
    P, Q, Pall, Qall = (np.vstack(x) for x in (P, Q, Pall, Qall))
    slab_residues = sum(m["slab_residues_matched"] for m in matches)
    if slab_residues < settings.orient_opm_min_slab_residues:
        raise OrientationFailure(f"only {slab_residues} anchor residues lie inside the hydrophobic slab of OPM "
                                 f"{reference.pdb_id} (at least {settings.orient_opm_min_slab_residues} needed)")
    fit = rigid_fit(P, Q, trim_floor_a=settings.orient_opm_trim_floor_a)
    whole = np.linalg.norm(Pall @ fit["R"].T + fit["t"] - Qall, axis=1)
    if fit["core_rmsd"] > settings.orient_opm_max_core_rmsd_a or fit["core_fraction"] < settings.orient_opm_min_core_fraction:
        raise OrientationFailure(f"the membrane-embedded core of the anchor does not superpose on OPM {reference.pdb_id}: "
                                 f"core RMSD {fit['core_rmsd']:.2f} A over {int(fit['keep'].sum())}/{len(P)} backbone atoms "
                                 f"(limits {settings.orient_opm_max_core_rmsd_a} A, {settings.orient_opm_min_core_fraction:.0%})")
    return {"R": fit["R"], "t": fit["t"], "matches": matches, "half_thickness_a": reference.half_thickness_a,
            "fit": {"fitted_on": "backbone N/CA/C/O of reference residues inside the hydrophobic slab",
                    "matched_atoms": int(len(P)), "core_atoms": int(fit["keep"].sum()), "core_fraction": round(fit["core_fraction"], 4),
                    "rmsd_A": round(fit["rmsd"], 4), "core_rmsd_A": round(fit["core_rmsd"], 4),
                    "max_deviation_A": round(fit["max_deviation"], 3), "determinant": round(fit["det"], 9),
                    "whole_chain_backbone_atoms": int(len(Pall)), "whole_chain_rmsd_A": round(rmsd(whole), 4),
                    "whole_chain_median_A": round(float(np.median(whole)), 4)}}


def find_ppm_executable(requested: str | None) -> Path:
    """Locate the PPM 3.0 executable and its residue library, refusing to guess when the requested one is absent."""
    candidates = [requested, os.environ.get("MEMBRANEFORGER_PPM")] + [shutil.which(n) for n in ("immers", "ppm3", "ppm")]
    chosen = next((c for c in candidates if c), None)
    if chosen is None:
        raise OrientationFailure("PPM 3.0 executable not found: pass --ppm-exe /path/to/immers, set MEMBRANEFORGER_PPM, or "
                                 "put the compiled PPM 3.0 program (immers) on PATH; its res.lib must sit next to it")
    exe = Path(shutil.which(chosen) or chosen).resolve()
    if not exe.is_file() or not os.access(exe, os.X_OK):
        raise OrientationFailure(f"PPM executable {chosen} is not an executable file")
    if not (exe.parent / "res.lib").is_file():
        raise OrientationFailure(f"PPM residue library res.lib is missing next to {exe}")
    return exe


class LocalPPM:
    """Run a locally installed PPM 3.0 (the standalone 'immers' program) on the anchor chain(s) only."""

    def __init__(self, exe: Path, membrane: str = "", heteroatoms: bool = False, timeout_s: int = 3600):
        """Remember the executable, the membrane code and whether the anchor heteroatoms are submitted."""
        if membrane not in PPM_MEMBRANE_CODES:
            raise OrientationFailure(f"unknown PPM membrane code {membrane!r}; known codes: "
                                     + ", ".join(sorted(c for c in PPM_MEMBRANE_CODES if c)) + " (empty = undefined bilayer)")
        self.exe, self.membrane, self.heteroatoms, self.timeout_s = exe, membrane, heteroatoms, timeout_s

    def input_text(self, chains: list[str], nterm: str) -> str:
        """PPM 3.0 multi-membrane input: one planar membrane, so the curved-membrane search is never enabled."""
        return (f"2\n{'yes' if self.heteroatoms else 'no'}\n{PPM_INPUT_NAME}\n1\n{self.membrane:<3s}\nplanar\n"
                f"{nterm:<3s}\n{','.join(chains)}\n")

    def run(self, atoms: list[dict], anchors: list[str], nterm: str, out: Path, work: Path, ranges: tuple = ()) -> dict:
        """Write the anchor (only the residues in ranges, when given), run PPM in isolation, and return its outputs."""
        if nterm not in ("in", "out"):
            raise OrientationFailure("PPM needs the N-terminal side (in or out) of the first anchor chain")
        workdir = work / "ppm"
        shutil.rmtree(workdir, ignore_errors=True)
        workdir.mkdir(parents=True)
        letters, relabel, submitted = "ABCDEFGHIJKLMNOPQRSTUVWXYZ", {}, []
        for label in anchors:
            relabel[label] = letters[len(relabel)]
        for a in atoms:
            label = a["chain"] or a.get("segid") or "A"
            if label in relabel and (a["resname"] in AMINO or self.heteroatoms) and (
                    not ranges or any(first <= a["resid"] <= last for first, last in ranges)):
                submitted.append(dict(a, chain=relabel[label], segid=""))
        if not submitted:
            raise OrientationFailure("no atoms to submit to PPM")
        write_pdb(submitted, workdir / PPM_INPUT_NAME)
        shutil.copy(self.exe.parent / "res.lib", workdir / "res.lib")
        inp = self.input_text([relabel[c] for c in anchors], nterm)
        (workdir / "ppm.inp").write_text(inp)
        started = time.time()
        transcript = run_command(out, [str(self.exe)], stdin=inp, cwd=workdir, timeout=self.timeout_s,
                                 produces=("anchorout.pdb", "datapar1"))
        if "wrong name" in transcript or "Too many" in transcript or "Place library" in transcript:
            raise OrientationFailure(f"PPM refused its input: {transcript.strip().splitlines()[-1]}")
        result = self.parse_outputs(workdir)
        residues = len({(a["chain"], a["resid"], a["resname"]) for a in submitted})
        result.update({"workdir": str(workdir), "seconds": round(time.time() - started, 1), "input": inp,
                       "chain_relabel": relabel, "atoms_submitted": len(submitted), "residues_submitted": residues,
                       "residue_ranges_submitted": [list(r) for r in ranges] or "all",
                       "heteroatoms_submitted": self.heteroatoms, "membrane_code": self.membrane or "(undefined)",
                       "curvature": "planar", "nterm_side": nterm, "executable": str(self.exe),
                       "executable_sha256": sha256(self.exe), "res_lib_sha256": sha256(self.exe.parent / "res.lib"),
                       "version": "PPM 3.0 standalone (immers); the program reports no version string",
                       "input_sha256": sha256(workdir / PPM_INPUT_NAME), "output_sha256": sha256(workdir / "anchorout.pdb"),
                       "stdout_tail": transcript.strip().splitlines()[-6:]})
        return result

    @staticmethod
    def parse_outputs(workdir: Path) -> dict:
        """Read anchorout.pdb, datapar1 and datasub1; anything missing, empty or malformed is a failure."""
        output = workdir / "anchorout.pdb"
        text = output.read_text(errors="replace")
        found = re.search(r"1/2 of bilayer thickness:\s*([-\d.]+)", text)
        atoms, _ = read_pdb(output)
        oriented = [a for a in atoms if a["resname"] != "DUM"]
        dummies = [a for a in atoms if a["resname"] == "DUM"]
        if not found or not oriented or not np.isfinite(xyz_nm(oriented)).all():
            raise OrientationFailure(f"PPM output {output} has no thickness REMARK, no atoms or non-finite coordinates")
        half = float(found.group(1))
        planes = sorted({round(a["z"], 1) for a in dummies})
        if not dummies or planes != sorted({round(half, 1), round(-half, 1)}):
            raise OrientationFailure(f"PPM output {output}: membrane dummy-atom planes {planes} do not match the half "
                                     f"thickness {half} (a curved or malformed result)")
        rows = [l for l in (workdir / "datapar1").read_text(errors="replace").splitlines() if l.strip()]
        row = next((l for l in rows if l.split(";")[0].strip() == PPM_INPUT_NAME), None)
        if row is None:
            raise OrientationFailure(f"PPM wrote no flat-membrane result for {PPM_INPUT_NAME} in datapar1 ({rows or 'empty'})")
        fields = [f.strip() for f in row.split(";")[1:] if f.strip()]
        try:
            thickness, thickness_sd, tilt, tilt_sd, energy = (float(f) for f in fields[:5])
        except ValueError:
            raise OrientationFailure(f"malformed datapar1 row: {row!r}") from None
        curved = (workdir / "datapar2").read_text(errors="replace").strip() if (workdir / "datapar2").is_file() else ""
        if PPM_INPUT_NAME in curved:
            raise OrientationFailure("PPM reported a curved membrane although a planar one was requested")
        segments = ""
        if (workdir / "datasub1").is_file():
            segments = next((l.strip() for l in (workdir / "datasub1").read_text(errors="replace").splitlines()
                             if l.split(";")[0].strip() == PPM_INPUT_NAME), "")
        return {"oriented": oriented, "half_thickness_a": half, "hydrophobic_thickness_a": thickness,
                "thickness_sd_a": thickness_sd, "tilt_deg": tilt, "tilt_sd_deg": tilt_sd,
                "transfer_energy_kcal_mol": energy, "transmembrane_segments": segments or None,
                "dummy_atoms_ignored": len(dummies)}


def orientation_from_ppm(anchors: OrderedDict, result: dict, settings: Settings) -> dict:
    """Derive the transform from the submitted anchor to PPM's copy of it; PPM only moves atoms, so the fit is exact."""
    oriented = protein_chains(result["oriented"])
    P, Q, matches = [], [], []
    for label, residues in anchors.items():
        best = match_chains(residues, oriented, 0.999)
        if best is None:
            raise OrientationFailure(f"PPM output has no chain with the sequence of anchor chain {label}")
        p, q, _ = atom_pairs(residues, oriented[best["chain"]], best["pairs"], names=None)
        matches.append({"anchor_chain": label, "ppm_chain": best["chain"], "aligned_residues": best["aligned"],
                        "identity": round(best["identity"], 4)})
        P.append(p), Q.append(q)
    P, Q = np.vstack(P), np.vstack(Q)
    heavy = sum(1 for residues in anchors.values() for _, _, atoms in residues for a in atoms if source_element(a) != "H")
    fraction = len(P) / heavy if heavy else 0.0
    if fraction < settings.orient_ppm_min_matched_fraction:
        raise OrientationFailure(f"only {len(P)} of {heavy} anchor heavy atoms were found in the PPM output "
                                 f"({fraction:.0%} < {settings.orient_ppm_min_matched_fraction:.0%})")
    fit = rigid_fit(P, Q)
    if fit["rmsd"] > settings.orient_ppm_max_fit_rmsd_a:
        raise OrientationFailure(f"PPM output is not a rigid copy of the submitted anchor: fit RMSD {fit['rmsd']:.3f} A over "
                                 f"{len(P)} atoms (limit {settings.orient_ppm_max_fit_rmsd_a} A, max deviation "
                                 f"{fit['max_deviation']:.2f} A)")
    return {"R": fit["R"], "t": fit["t"], "matches": matches, "half_thickness_a": result["half_thickness_a"],
            "fit": {"fitted_on": "every heavy atom of the anchor found in the PPM output", "matched_atoms": int(len(P)),
                    "matched_fraction": round(fraction, 4), "rmsd_A": round(fit["rmsd"], 5),
                    "max_deviation_A": round(fit["max_deviation"], 4), "determinant": round(fit["det"], 9)}}


def validate_full_transform(original: list[dict], oriented: list[dict], R: np.ndarray, t: np.ndarray,
                            pairs: int = 50000) -> dict:
    """Prove the whole complex received one proper rigid transform: distances invariant, every chain moved alike."""
    X, Y = xyz_nm(original), xyz_nm(oriented)
    if len(X) != len(Y) or any((a["atom"], a["resname"], a["resid"], a["chain"]) != (b["atom"], b["resname"], b["resid"], b["chain"])
                               for a, b in zip(original, oriented)):
        raise OrientationFailure("the oriented complex does not list the same atoms as the input")
    rng = np.random.default_rng(20260101)
    i, j = rng.integers(0, len(X), size=(2, min(pairs, len(X) * (len(X) - 1) // 2)))
    before, after = np.linalg.norm(X[i] - X[j], axis=1), np.linalg.norm(Y[i] - Y[j], axis=1)
    distance_change = float(np.abs(after - before).max())
    reproduced = float(np.abs(X @ R.T + t - Y).max())
    det = float(np.linalg.det(R))
    per_chain = {}
    for label in OrderedDict.fromkeys(a["chain"] or a.get("segid") or "A" for a in original):
        idx = [k for k, a in enumerate(original) if (a["chain"] or a.get("segid") or "A") == label]
        own = rigid_fit(X[idx], Y[idx]) if len(idx) >= 3 else None
        per_chain[label] = {"atoms": len(idx), "rmsd_to_global_transform_A": round(reproduced, 6) if own is None else
                            round(float(np.sqrt(((X[idx] @ R.T + t - Y[idx]) ** 2).sum(axis=1).mean())), 6),
                            "own_fit_rmsd_A": None if own is None else round(own["rmsd"], 6)}
    problems = []
    if distance_change > 1e-6:
        problems.append(f"pairwise distances changed by up to {distance_change:.2e} A")
    if reproduced > 1e-6:
        problems.append(f"the complex does not reproduce x' = R x + t (max {reproduced:.2e} A)")
    if abs(det - 1.0) > 1e-9:
        problems.append(f"rotation determinant {det:.9f}")
    if problems:
        raise OrientationFailure("full-complex transform validation failed: " + "; ".join(problems))
    return {"pairwise_distance_max_change_A": distance_change, "sampled_pairs": int(len(i)),
            "max_reproduction_error_A": reproduced, "determinant": det, "chains": per_chain}


def rotation_angle_deg(R: np.ndarray) -> float:
    """Angle of a rotation matrix in degrees."""
    return float(np.degrees(np.arccos(np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0))))


def frame_metrics(oriented: list[dict], anchors: list[str], half_thickness_a: float, ranges: tuple = ()) -> dict:
    """Describe where the oriented anchor sits in the membrane frame (and, with ranges, what lies outside the segment)."""
    ca = np.array([[a["x"], a["y"], a["z"]] for a in oriented
                   if a["atom"] == "CA" and (a["chain"] or a.get("segid") or "A") in anchors and a["resname"] in AMINO])
    inside = int((np.abs(ca[:, 2]) <= half_thickness_a).sum()) if len(ca) else 0
    whole = xyz_nm(oriented)
    return {"anchor_ca_atoms": int(len(ca)), "anchor_ca_inside_slab": inside,
            "anchor_ca_above_slab": int((ca[:, 2] > half_thickness_a).sum()) if len(ca) else 0,
            "anchor_ca_below_slab": int((ca[:, 2] < -half_thickness_a).sum()) if len(ca) else 0,
            "anchor_ca_centroid_A": [round(float(v), 3) for v in ca.mean(axis=0)] if len(ca) else None,
            "complex_z_range_A": [round(float(whole[:, 2].min()), 3), round(float(whole[:, 2].max()), 3)]} | (
            segment_metrics(oriented, anchors, half_thickness_a, ranges) if ranges else {})


def embedded_segment(aa_atoms: list[dict], request: OrientationRequest, ranges: tuple,
                     half_thickness_a: float = EMBEDDED_SEGMENT_HALF_THICKNESS_A) -> dict:
    """With --orientation none: the bilayer centre (A) that puts the named transmembrane segment in the middle of the slab.

    The input already has its membrane normal along z and is not rotated; --orient-residues then names the segment
    that sits in the bilayer, and the bilayer centre is the midpoint of that segment's CA z range. This replaces the
    hydrophobic-belt search over the whole protein (embedding.hydrophobic_belt), which an extracellular or
    cytoplasmic domain with hydrophobic surface can pull off the helix.
    """
    chains = protein_chains(aa_atoms)
    anchors, how = select_anchor_chains(chains, tuple(request.chains))
    selected = select_residues(OrderedDict((label, chains[label]) for label in anchors), ranges)
    z = [a["z"] for residues in selected.values() for _, _, atoms in residues for a in atoms if a["atom"] == "CA"]
    if not z:
        raise OrientationFailure(f"--orient-residues {', '.join(f'{f}-{l}' for f, l in ranges)} selects no CA atom")
    centre = 0.5 * (min(z) + max(z))
    shifted = [dict(a, z=a["z"] - centre) for a in aa_atoms]
    segment = segment_metrics(shifted, anchors, half_thickness_a, ranges)
    segment.update({"segment_ca_z_range_A": [round(min(z), 3), round(max(z), 3)], "half_thickness_A": half_thickness_a})
    return {"anchor_chains": anchors, "anchor_selection": how, "segment_bilayer_centre_A": round(centre, 3), "segment": segment}


def segment_metrics(oriented: list[dict], anchors: list[str], half_thickness_a: float, ranges: tuple) -> dict:
    """With a residue selection: the protein residues outside it whose CA lies in the hydrophobic slab.

    A single-pass protein oriented on its transmembrane helix should have everything else (ectodomain, linker,
    cytoplasmic tail, partner chains) outside the slab; residues listed here sit in the bilayer core."""
    in_segment = lambda a: (a["chain"] or a.get("segid") or "A") in anchors and any(f <= a["resid"] <= l for f, l in ranges)
    ca = [a for a in oriented if a["atom"] == "CA" and a["resname"] in AMINO]
    segment = [a for a in ca if in_segment(a)]
    buried = [a for a in ca if not in_segment(a) and abs(a["z"]) <= half_thickness_a]
    return {"segment_ranges": [list(r) for r in ranges], "segment_ca_atoms": len(segment),
            "segment_ca_inside_slab": sum(abs(a["z"]) <= half_thickness_a for a in segment),
            "outside_segment_ca_atoms": len(ca) - len(segment), "outside_segment_ca_inside_slab": len(buried),
            "outside_segment_residues_inside_slab": [f"{a['chain'] or a.get('segid') or 'A'}:{a['resname']}{a['resid']} "
                                                     f"z {a['z']:+.1f}" for a in buried]}


def nterm_check(oriented: list[dict], anchor: str, nterm: str | None, half_thickness_a: float, ranges: tuple = ()) -> dict:
    """Check that the first anchor residue (of the selection, with ranges) sits on the requested side when outside the slab."""
    residues = select_residues(protein_chains([a for a in oriented if (a["chain"] or a.get("segid") or "A") == anchor]),
                               ranges).get(anchor, [])
    ca = next((a for _, _, atoms in residues for a in atoms if a["atom"] == "CA"), None)
    if ca is None:
        return {"status": "no CA atom in the anchor"}
    report = {"first_residue": f"{residues[0][1]}{residues[0][0]}", "first_ca_z_A": round(ca["z"], 3), "requested": nterm}
    if nterm not in ("in", "out"):
        report["status"] = "not requested"
    elif abs(ca["z"]) <= half_thickness_a:
        report["status"] = "first residue inside the hydrophobic slab; sidedness not testable from it"
    else:
        observed = "out" if ca["z"] > 0 else "in"
        report["status"] = "PASS" if observed == nterm else f"FAIL: N terminus is {observed}"
        if observed != nterm:
            raise OrientationFailure(f"the oriented anchor has its N terminus {observed} (z {ca['z']:+.1f} A) although {nterm} "
                                     "was specified")
    return report


def orient_complex(aa_atoms: list[dict], aa_path: Path, request: OrientationRequest, settings: Settings,
                   out: Path, work: Path) -> dict:
    """Orient the complete complex from its anchor chain(s) and return the oriented atoms with the report."""
    if request.mode not in ORIENTATION_MODES:
        raise OrientationFailure(f"ORIENTATION must be one of {ORIENTATION_MODES}, not {request.mode!r}")
    if request.nterm_side not in NTERM_SIDES:
        raise OrientationFailure(f"NTERM_SIDE must be one of {NTERM_SIDES}, not {request.nterm_side!r}")
    ranges = tuple(request.residues)
    if ranges and (request.mode == "opm" or request.opm_file):
        raise OrientationFailure("--orient-residues selects the residues submitted to PPM; it does not apply to --orientation opm "
                                 "or --opm-file")
    if ranges and request.mode != "none" and request.nterm_side not in ("in", "out"):
        # An OPM entry gives the side of a chain's N terminus, which need not be the side of the selection's first residue.
        raise OrientationFailure("--orient-residues needs --nterm-side in or out: the side of the membrane on which the FIRST "
                                 "SELECTED residue lies (for a type I single-pass protein oriented on its helix: out)")
    report = {"enabled": request.mode != "none", "mode": request.mode, "provider": None, "convention":
              "membrane normal +z, midplane z = 0, IN (cytoplasmic) negative z, OUT positive z, coordinates in A",
              "input_file": str(aa_path), "input_sha256": sha256(aa_path), "request": {**asdict(request),
              "opm_file": str(request.opm_file) if request.opm_file else None,
              "opm_cache": str(request.opm_cache) if request.opm_cache else None}}
    if request.mode == "none":
        report.update({"anchor_chains": [], "status": "skipped: orientation disabled; the input coordinates are used as given"})
        if ranges:
            report.update(embedded_segment(aa_atoms, request, ranges))
            frame = report["segment"]
            buried = frame["outside_segment_residues_inside_slab"]
            log(out, f"embedded segment {', '.join(f'{f}-{l}' for f, l in ranges)} of chain(s) {', '.join(report['anchor_chains'])}: "
                     f"bilayer centre at z = {report['segment_bilayer_centre_A']:+.2f} A (midpoint of its CA z range "
                     f"{frame['segment_ca_z_range_A'][0]:+.1f} to {frame['segment_ca_z_range_A'][1]:+.1f} A); the input is not rotated; "
                     f"{len(buried)}/{frame['outside_segment_ca_atoms']} protein CA outside the segment lie within "
                     f"+-{frame['half_thickness_A']:.0f} A of that centre"
                     + (f" ({', '.join(buried[:12])}{', ...' if len(buried) > 12 else ''})" if buried else ""),
                "WARN" if buried else "INFO")
        return {"oriented": aa_atoms, "report": report, "R": np.eye(3), "t": np.zeros(3)}
    chains = protein_chains(aa_atoms)
    pdb_id, id_source = resolve_pdb_id(request.pdb_id, aa_path)
    report.update({"pdb_id": pdb_id, "pdb_id_source": id_source})
    reference = None
    if request.opm_file:
        reference = parse_opm_file(Path(request.opm_file), pdb_id, str(request.opm_file))
        log(out, f"OPM reference read from {request.opm_file} (sha256 {reference.sha256[:12]})")
    elif pdb_id and request.mode in ("auto", "opm", "ppm") and not ranges:  # a residue selection is oriented by PPM only
        try:
            reference = fetch_opm_reference(pdb_id, request.opm_cache, out)
        except OrientationFailure as exc:
            if request.mode == "opm":
                raise
            log(out, f"{exc}; continuing without an OPM reference", "WARN")
    if request.mode == "opm" and reference is None:
        raise OrientationFailure("ORIENTATION=opm needs an exact OPM/OPRLM entry: pass --pdb-id of an entry OPM holds, or "
                                 "--opm-file with its downloaded coordinate file" + (f" (no entry for {pdb_id})" if pdb_id else ""))
    anchors, how = select_anchor_chains(chains, tuple(request.chains), reference, settings.orient_min_identity)
    anchor_residues = OrderedDict((label, chains[label]) for label in anchors)
    segment = select_residues(anchor_residues, ranges)  # fails here, before anything expensive, on a bad selection
    report.update({"anchor_chains": anchors, "anchor_selection": how,
                   "anchor_residues": {label: len(chains[label]) for label in anchors},
                   "anchor_segment": {"ranges": [list(r) for r in ranges],
                                      "residues": {label: len(r) for label, r in segment.items()},
                                      "role": "only these residues are submitted to PPM; the whole complex follows them "
                                              "rigidly"} if ranges else None})
    log(out, f"orientation anchor chain(s) {','.join(anchors)} ({how}); PDB ID {pdb_id or 'none'} ({id_source})")
    # N-terminal side: the command line, else the exact reference, else unknown.
    nterm, nterm_source = (request.nterm_side, "command line") if request.nterm_side in ("in", "out") else (None, "unknown")
    if reference is not None:  # never with a residue selection: the reference is not loaded then
        best = match_chains(anchor_residues[anchors[0]], reference.chains, settings.orient_min_identity)
        if best is not None:
            side, why = nterm_side_from_reference(reference, best["chain"])
            if nterm is None and side is not None:
                nterm, nterm_source = side, why
            elif nterm is not None and side is not None and side != nterm:
                raise OrientationFailure(f"--nterm-side {nterm} contradicts the exact reference ({why}); check the chain or the side")
    report.update({"nterm_side": nterm, "nterm_side_source": nterm_source})
    derived, provider = None, None
    if reference is not None and request.mode in ("auto", "opm"):  # ORIENTATION=ppm uses the reference for the N-terminal side only
        try:
            derived = orientation_from_opm(anchor_residues, reference, settings)
            provider = "opm"
            report["reference"] = {"pdb_id": reference.pdb_id, "source": reference.source, "file": str(reference.path),
                                   "sha256": reference.sha256, "half_thickness_A": reference.half_thickness_a,
                                   "spanning_chains": sorted(spanning_chains(reference)), "chain_matches": derived["matches"],
                                   "role": "orientation reference only; no coordinate of the user's complex comes from it"}
            log(out, f"OPM {reference.pdb_id}: anchor core RMSD {derived['fit']['core_rmsd_A']} A over "
                     f"{derived['fit']['core_atoms']} backbone atoms; whole-chain backbone RMSD {derived['fit']['whole_chain_rmsd_A']} A")
        except OrientationFailure as exc:
            if request.mode == "opm":
                raise
            report["reference_rejected"] = str(exc)
            log(out, f"OPM reference not used: {exc}; orienting with PPM instead", "WARN")
    if derived is None:
        if nterm is None:
            raise OrientationFailure(
                "The biological membrane sidedness cannot be determined safely.\n\n"
                "PPM needs the side of the membrane on which the N terminus of the first anchor chain "
                f"({anchors[0]}) lies. Specify either:\n\n    --nterm-side in\n\nor:\n\n    --nterm-side out\n"
                + (f"\n({nterm_source})" if nterm_source != "unknown" else ""))
        exe = find_ppm_executable(request.ppm_exe)
        ppm = LocalPPM(exe, request.ppm_membrane, request.ppm_heteroatoms, settings.orient_ppm_timeout_s)
        log(out, f"running PPM 3.0 ({exe}) on {sum(len(r) for r in segment.values())} anchor residues"
                 + (f" (residues {', '.join(f'{a}-{b}' for a, b in ranges)})" if ranges else "") + ", "
                 f"membrane {request.ppm_membrane or 'undefined (flat bilayer)'}, N terminus {nterm}")
        result = ppm.run(aa_atoms, anchors, nterm, out, work, ranges)
        derived = orientation_from_ppm(segment, result, settings)
        provider = "ppm"
        report["ppm"] = {k: v for k, v in result.items() if k != "oriented"} | {"chain_matches": derived["matches"]}
        log(out, f"PPM: hydrophobic thickness {result['hydrophobic_thickness_a']} +- {result['thickness_sd_a']} A, tilt "
                 f"{result['tilt_deg']} +- {result['tilt_sd_deg']} deg, transfer energy {result['transfer_energy_kcal_mol']} "
                 f"kcal/mol, {result['seconds']} s; anchor fit RMSD {derived['fit']['rmsd_A']} A over {derived['fit']['matched_atoms']} atoms")
    R, t = derived["R"], derived["t"]
    oriented = apply_transform(aa_atoms, R, t)
    validation = validate_full_transform(aa_atoms, oriented, R, t)
    validation["frame"] = frame_metrics(oriented, anchors, derived["half_thickness_a"], ranges)
    validation["nterm"] = nterm_check(oriented, anchors[0], nterm, derived["half_thickness_a"], ranges)
    validation["status"] = "PASS"
    oriented_path = out / "oriented.pdb"
    write_pdb(oriented, oriented_path)
    report.update({"provider": provider, "fit": derived["fit"], "half_thickness_A": derived["half_thickness_a"],
                   "rotation_matrix": [[float(v) for v in row] for row in R], "translation_A": [float(v) for v in t],
                   "rotation_angle_deg": round(rotation_angle_deg(R), 3),
                   "membrane_normal": [0.0, 0.0, 1.0], "membrane_center_A": [0.0, 0.0, 0.0],
                   "validation": validation, "oriented_file": str(oriented_path), "oriented_sha256": sha256(oriented_path),
                   "status": "PASS"})
    log(out, f"orientation by {provider}: rotation {report['rotation_angle_deg']} deg, {len(oriented)} atoms moved as one body; "
             f"{validation['frame']['anchor_ca_inside_slab']}/{validation['frame']['anchor_ca_atoms']} anchor CA inside the "
             f"+-{derived['half_thickness_a']} A slab; pairwise distances changed by <= {validation['pairwise_distance_max_change_A']:.1e} A", "PASS")
    if ranges:
        frame = validation["frame"]
        buried = frame["outside_segment_residues_inside_slab"]
        log(out, f"anchor segment: {frame['segment_ca_inside_slab']}/{frame['segment_ca_atoms']} CA inside the slab; "
                 f"{len(buried)}/{frame['outside_segment_ca_atoms']} protein CA outside the segment lie inside it"
                 + (f" ({', '.join(buried[:12])}{', ...' if len(buried) > 12 else ''})" if buried else ""),
            "WARN" if buried else "INFO")
    return {"oriented": oriented, "report": report, "R": R, "t": t}


def check_orientation_preserved(oriented: list[dict], placed: list[dict], anchors: list[str], midplane_a: float) -> dict:
    """Prove the CG registration kept the orientation: only a rotation about z and a translation separate the two."""
    X, Y = xyz_nm(oriented), xyz_nm(placed)
    fit = rigid_fit(X, Y)
    R = fit["R"]
    normal = R @ np.array([0.0, 0.0, 1.0])
    tilt = float(np.degrees(np.arctan2(math.hypot(normal[0], normal[1]), normal[2])))  # well conditioned near 0
    anchor = [k for k, a in enumerate(oriented) if (a["chain"] or a.get("segid") or "A") in anchors and a["atom"] == "CA"]
    depth_before = float(X[anchor, 2].mean())
    depth_after = float(Y[anchor, 2].mean() - midplane_a)
    # Limits are numerical noise of a float64 fit over ~10^4 atoms, far below anything physical (1e-4 deg, 1e-6 A).
    if fit["rmsd"] > 1e-6 or tilt > 1e-4 or abs(depth_after - depth_before) > 1e-6:
        raise OrientationFailure(f"CG registration changed the orientation: normal tilted {tilt:.2e} deg, depth moved "
                                 f"{depth_after - depth_before:.2e} A, rigid residual {fit['rmsd']:.2e} A")
    return {"status": "PASS", "membrane_normal_after_registration": [round(float(v), 9) for v in normal],
            "normal_tilt_deg": tilt, "rotation_about_z_deg": round(float(np.degrees(np.arctan2(R[1, 0], R[0, 0]))), 4),
            "anchor_ca_depth_A": round(depth_before, 4), "anchor_ca_depth_after_registration_A": round(depth_after, 4),
            "rigid_residual_A": fit["rmsd"]}
