#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Input validation and output-directory protection.
#//=============================================================
"""Read and validate both inputs, and keep a build from overwriting inputs or reusing stale outputs."""
import json
import math
import shutil
import time
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .config import AMINO, GENERATED, SOLVENT
from .martini import MARTINI3_PROTEIN, MARTINI3_PROTEIN_ALIASES
from .runtools import LOG_NAME, sha256
from .structio import read_pdb, xyz_nm

__all__ = ['read_cg', 'read_all_atom', 'protect_inputs', 'SIMULATION_PATTERNS', 'existing_outputs', 'refuse_overwrite',
           'clear_stale_outputs', 'ARCHIVE_DIR', 'archive_previous_outputs', 'check_inputs_unchanged', 'frame_overlaps', 'align_frame_to_box']

# Files a GROMACS simulation started from a build leaves in its directory: checkpoints (mdrun -cpi continues from them, with
# -append or -noappend) and trajectories. The build itself runs steepest-descent minimizations only; it never reads, writes or
# deletes these files, but a rerun would replace the em.gro, topol.top and index files such a simulation depends on.
SIMULATION_PATTERNS = ("*.cpt", "*.xtc", "*.tng", "*.trr")

def read_cg(path: Path) -> tuple[list[dict], list[float]]:
    """Read a coarse-grained .gro or .pdb into beads (nm) and an orthorhombic box (nm)."""
    try:
        if path.suffix.lower() == ".gro":
            lines = path.read_text(errors="replace").splitlines()
            n = int(lines[1])
            if n <= 0 or len(lines) < n + 3:
                raise ValueError(f"header declares {n} atoms but the file has {max(len(lines) - 3, 0)} atom lines")
            atoms = [{"resid": int(l[0:5]), "resname": l[5:10].strip(), "atom": l[10:15].strip(), "chain": "",
                      "x": float(l[20:28]), "y": float(l[28:36]), "z": float(l[36:44])} for l in lines[2:2 + n]]
            cell = [float(v) for v in lines[2 + n].split()]
            if len(cell) not in (3, 9) or any(abs(v) > 1e-6 for v in cell[3:]):
                raise ValueError("box is missing or not orthorhombic")
            box = cell[:3]
        elif path.suffix.lower() == ".pdb":
            if sum(l.startswith("MODEL ") for l in path.read_text(errors="replace").splitlines()) > 1:
                raise ValueError("multiple MODEL records")
            raw, cryst1 = read_pdb(path)
            if not cryst1:
                raise ValueError("no CRYST1 box")
            cryst = cryst1.ljust(54)
            if not all(abs(float(cryst[i:i + 7]) - 90.0) < 0.02 for i in (33, 40, 47)):
                raise ValueError("box is not orthorhombic")
            box = [float(cryst[6:15]) / 10.0, float(cryst[15:24]) / 10.0, float(cryst[24:33]) / 10.0]
            atoms = [dict(a, x=a["x"] / 10.0, y=a["y"] / 10.0, z=a["z"] / 10.0) for a in raw if a["altloc"] in " A"]
        else:
            raise ValueError("expected a .pdb or .gro file")
        if not atoms:
            raise ValueError("no atoms")
        if not np.isfinite(xyz_nm(atoms)).all() or not all(math.isfinite(v) and v > 0 for v in box):
            raise ValueError("non-finite coordinates or a non-positive box")
    except (ValueError, IndexError) as exc:
        raise SystemExit(f"coarse-grained input {path.name}: {exc}") from None
    return atoms, box


def read_all_atom(path: Path, forcefield: Path) -> tuple[list[dict], list[str]]:
    """Read and validate the all-atom protein/ligand PDB."""
    # Returns its atoms (altloc A only) and notes for the log.
    text = path.read_text(errors="replace").splitlines()
    try:
        atoms, _ = read_pdb(path)
    except (ValueError, IndexError) as exc:
        raise SystemExit(f"all-atom input {path.name}: malformed ATOM/HETATM record ({exc})") from None
    if not atoms:
        raise SystemExit(f"all-atom input {path.name}: no ATOM/HETATM records")
    if sum(l.startswith("MODEL ") for l in text) > 1:
        raise SystemExit(f"all-atom input {path.name}: multiple MODEL records")
    if any(l[:6] in ("ATOM  ", "HETATM") and l[26:27].strip() for l in text):
        raise SystemExit(f"all-atom input {path.name}: residue insertion codes are not supported")
    if not np.isfinite(xyz_nm(atoms)).all():
        raise SystemExit(f"all-atom input {path.name}: non-finite coordinates")
    alternates = sorted({f"{a['chain']}:{a['resname']}{a['resid']}" for a in atoms if a["altloc"] not in " A"})
    atoms = [a for a in atoms if a["altloc"] in " A"]
    if not any(a["resname"] in AMINO for a in atoms):
        raise SystemExit(f"all-atom input {path.name}: no protein residues")
    seen = Counter((a["chain"], a["segid"], a["resid"], a["resname"], a["atom"]) for a in atoms)
    repeats = [k for k, n in seen.items() if n > 1]
    if repeats:
        raise SystemExit(f"all-atom input {path.name}: {len(repeats)} duplicate atoms, e.g. "
                         f"{repeats[0][0]}:{repeats[0][3]}{repeats[0][2]}:{repeats[0][4]}")
    others = sorted({a["resname"] for a in atoms if a["resname"] not in AMINO})
    bulk = [n for n in others if n in SOLVENT | {"HOH", "SOL", "WAT", "TIP", "NA", "CL", "K"}]
    if bulk:
        raise SystemExit(f"all-atom input {path.name}: contains solvent or bulk ions {bulk}; "
                         "supply only the protein/ligand complex")
    known = {p.stem for p in (forcefield / "toppar").glob("*.itp")} | ({"V6G"} if (forcefield / "V6G_gromacs").is_dir() else set())
    missing = [n for n in others if n not in known]
    if missing:
        raise SystemExit(f"all-atom input {path.name}: no topology in {forcefield} for residue(s) {missing}")
    return atoms, [f"{path.name}: {len(atoms)} atoms, {sum(a['resname'] in AMINO for a in atoms)} protein, "
                   f"ligands {', '.join(others) or 'none'}; dropped altloc B of {', '.join(alternates) or 'none'}"]


def protect_inputs(inputs: list[Path], out: Path, keep: Path | None = None, resources: tuple = ()) -> None:
    """Refuse an output directory in which the build would overwrite or delete an input or an installation resource."""
    # `keep` is the --membrane source: it may live in the output directory (or in its int/ folder, where an earlier
    # build kept it) only under the name membrane.pdb, the one generated file that route never writes.
    for path in inputs:
        if out not in path.parents:
            continue
        first = path.relative_to(out).parts[0]
        kept = path == keep and path.name == "membrane.pdb" and path.parent in (out, out / "int")
        if (first in GENERATED or first == LOG_NAME) and not kept:
            raise SystemExit(f"input {path} would be overwritten by the build; choose another --out")
    for resource in resources:
        if resource is not None and (out == resource or resource in out.parents or out in resource.parents):
            raise SystemExit(f"--out {out} overlaps the installation directory {resource}; choose another --out")


def existing_outputs(out: Path, keep: Path | None = None) -> dict:
    """What an output directory already holds: {"build": outputs of an earlier build, "simulation": checkpoints and trajectories}.

    "build" lists the config.GENERATED names present (the --membrane source `keep`, and int/ when it holds nothing else, do not
    count). "simulation" lists checkpoint and trajectory files that are not themselves GENERATED (the minimizations' em.cpt,
    em.trr and the like are), anywhere below the directory.
    """
    out = Path(out)
    build = []
    for name in GENERATED:
        target = out / name
        if target == keep or not (target.exists() or target.is_symlink()):
            continue
        if keep is not None and target.is_dir() and target in Path(keep).parents:
            if not any(child != keep for child in target.iterdir()):
                continue
        build.append(name)
    simulation = sorted({str(p.relative_to(out)) for pattern in SIMULATION_PATTERNS for p in out.rglob(pattern)
                         if p.is_file() and p.relative_to(out).parts[0] not in GENERATED and p.name not in GENERATED}) if out.is_dir() else []
    return {"build": build, "simulation": simulation}


def refuse_overwrite(out: Path, keep: Path | None = None, overwrite: bool = False, overwrite_restart: bool = False) -> None:
    """Refuse an output directory that holds an earlier build or a simulation's restart files unless the matching flag is given.

    A rerun deletes the earlier build's outputs (clear_stale_outputs), among them the em.gro, topol.top and index files a running
    or continued simulation reads. So an earlier build needs --overwrite, and checkpoints or trajectories need
    --overwrite-restart as well. The build never deletes those simulation files and never starts an mdrun with -cpi: GROMACS
    continuation (mdrun -cpi state.cpt, with -append or -noappend) in that directory is not affected by anything here.
    """
    found = existing_outputs(out, keep)
    if found["simulation"] and not (overwrite and overwrite_restart):
        shown = ", ".join(found["simulation"][:4]) + (f" and {len(found['simulation']) - 4} more" if len(found["simulation"]) > 4 else "")
        raise SystemExit(f"--out {out} holds simulation restart files ({shown}); a rerun would replace the em.gro, topol.top and index "
                         "files a simulation continued from them depends on. Choose another --out, or pass --overwrite and "
                         "--overwrite-restart to replace the build outputs (the simulation files themselves are never deleted)")
    if found["build"] and not overwrite:
        shown = ", ".join(found["build"][:5]) + (f" and {len(found['build']) - 5} more" if len(found["build"]) > 5 else "")
        raise SystemExit(f"--out {out} already holds the output of an earlier build ({shown}); choose another --out, or pass "
                         "--overwrite to replace it")


# --archive-previous: what is copied out of an earlier build before --overwrite replaces it. "reports" is the small set that explains
# the earlier outcome (tens of kB); "all" adds its structures (about 10 MB for a 200k-atom system). Files are looked up in the output
# directory and in int/, where a failed build leaves em.unverified.gro and its intermediates.
ARCHIVE_DIR = "previous_builds"
ARCHIVE_REPORTS = (LOG_NAME, "run_manifest.json", "audit.json", "backbone_check.json", "contacts.json", "ring_piercing.json",
                   "orientation_report.json", "aa_cg_mapping.tsv", "membrane_validation.json", "membrane_validation.md",
                   "membrane_validation_lipids.tsv")
ARCHIVE_STRUCTURES = ("em.gro", "em.unverified.gro", "topol.top", "index_ini.ndx")


def archive_previous_outputs(out: Path, mode: str, keep: Path | None = None) -> Path | None:
    """Copy the reports (mode "reports") or reports and structures (mode "all") of an earlier build to <out>/previous_builds/<time>_<status>/.

    Returns that folder, or None when `out` holds no earlier build. Call it after refuse_overwrite and before clear_stale_outputs: the
    folder is not a GENERATED name, so the rerun leaves it alone, and nothing is deleted here. A file that cannot be copied raises
    OSError, so the caller can stop before anything is replaced."""
    out = Path(out)
    if mode == "no" or not existing_outputs(out, keep)["build"]:
        return None
    status, stamp = "unknown", None
    manifest = out / "run_manifest.json"
    if manifest.is_file():
        try:
            record = json.loads(manifest.read_text())
            status = str(record.get("status", "unknown"))
            stamp = time.strftime("%Y%m%dT%H%M%S", time.strptime(str(record.get("started", ""))[:19], "%Y-%m-%dT%H:%M:%S"))
        except (ValueError, AttributeError):
            pass
    stamp = stamp or time.strftime("%Y%m%dT%H%M%S", time.localtime((manifest if manifest.is_file() else out).stat().st_mtime))
    folder, n = out / ARCHIVE_DIR / f"{stamp}_{status}", 1
    while folder.exists():
        n += 1
        folder = out / ARCHIVE_DIR / f"{stamp}_{status}_{n}"
    folder.mkdir(parents=True)
    names = ARCHIVE_REPORTS + (ARCHIVE_STRUCTURES if mode == "all" else ())
    copied = []
    for name in names:
        for source in (out / name, out / "int" / name):
            if source.is_file():
                shutil.copy2(source, folder / name)
                copied.append(name)
                break
    (folder / "ARCHIVED.json").write_text(json.dumps({"from": str(out), "mode": mode, "status": status, "files": copied,
                                                      "archived": time.strftime("%Y-%m-%dT%H:%M:%S%z")}, indent=1) + "\n")
    return folder


def clear_stale_outputs(out: Path, keep: Path | None = None) -> list[str]:
    """Delete every file a previous build could have left in the output directory; returns what was removed."""
    removed = []
    for name in GENERATED:
        target = out / name
        if target == keep or not (target.exists() or target.is_symlink()):
            continue
        if target.is_dir() and not target.is_symlink() and keep is not None and target in Path(keep).parents:
            for child in target.iterdir():  # the folder holding the --membrane source: clear everything else in it
                if child != keep:
                    shutil.rmtree(child) if child.is_dir() and not child.is_symlink() else child.unlink()
        elif target.is_dir() and not target.is_symlink():
            shutil.rmtree(target)
        else:
            target.unlink()
        removed.append(name)
    return removed


def check_inputs_unchanged(hashes: dict) -> None:
    """Fail if any input file no longer has the hash recorded at the start of the build."""
    for path, digest in hashes.items():
        if sha256(path) != digest:
            raise SystemExit(f"{path.name} changed during the build")


# A frame whose coordinates were rotated about z after the simulation (a rotational fit written without its box) is
# no longer periodic in the box it declares: wrapping it folds the corners of the rotated square onto its edges, so
# beads of different molecules land on top of each other there and the opposite corners stay empty. Every periodic
# operation downstream (slicing images, the seam, Voronoi cells, RDFs) would then be wrong near the edges. The test
# below counts bead pairs of DIFFERENT lipid residues closer than FRAME_OVERLAP_NM after wrapping (protein beads are
# left out: bonded side chains of neighbouring residues sit closer than that); an equilibrated Martini membrane has
# none (its closest inter-molecule pair is about 0.34 nm), a rotated one has thousands.
FRAME_OVERLAP_NM = 0.30
FRAME_SOLVENT = {"W", "SW", "TW", "ION", "NA", "CL", "CA", "K", "MG"}


def frame_overlaps(xyz: np.ndarray, owner: np.ndarray, box: np.ndarray, radius: float = FRAME_OVERLAP_NM) -> int:
    """Number of bead pairs of different residues closer than radius (3D) once the coordinates are wrapped into the box."""
    wrapped = np.mod(xyz, box)
    wrapped[wrapped >= box] = 0.0  # np.mod of a tiny negative coordinate rounds up to exactly the box edge, which cKDTree refuses
    pairs = cKDTree(wrapped, boxsize=box).query_pairs(radius, output_type="ndarray")
    return int((owner[pairs[:, 0]] != owner[pairs[:, 1]]).sum()) if len(pairs) else 0


def align_frame_to_box(atoms: list[dict], box: list[float]) -> tuple[list[dict], dict]:
    """Detect a frame rotated about z relative to its (square) box and rotate it back so that it is periodic again.

    The rotation angle is the one at which wrapping creates the fewest inter-residue bead overlaps (coarse 1-degree
    scan over the 90 degrees a square box allows, refined to 0.01 degree). All atoms, water and ions included, are
    rotated about the box centre. A frame that no rotation makes periodic is refused. Returns the atoms and a report.
    """
    cell, cell3 = np.array(box[:2], dtype=float), np.array(box[:3], dtype=float)
    kept = [i for i, a in enumerate(atoms) if a["resname"] not in FRAME_SOLVENT and a["atom"] not in FRAME_SOLVENT
            and MARTINI3_PROTEIN_ALIASES.get(a["resname"], a["resname"]) not in MARTINI3_PROTEIN]
    xyz = np.array([[atoms[i]["x"], atoms[i]["y"], atoms[i]["z"]] for i in kept])
    xy = xyz[:, :2]
    keys, owner, last = {}, np.zeros(len(kept), int), None
    for n, i in enumerate(kept):  # a residue is a run of equal (resid, resname), as in martini.cg_residues
        key = (atoms[i]["resid"], atoms[i]["resname"])
        if key != last:
            keys[key] = len(keys)
            last = key
        owner[n] = len(keys) - 1
    centre = 0.5 * cell
    rotated = lambda deg: (xy - centre) @ np.array([[math.cos(math.radians(deg)), math.sin(math.radians(deg))],
                                                   [-math.sin(math.radians(deg)), math.cos(math.radians(deg))]]) + centre
    score = lambda deg: frame_overlaps(np.column_stack([rotated(deg), xyz[:, 2]]), owner, cell3)
    as_read = score(0.0)
    square = bool(abs(cell[0] - cell[1]) < 0.01)
    best = (0.0, as_read)
    if square:
        for step, span in ((1.0, 90.0), (0.1, 1.0), (0.01, 0.1)):
            start = best[0] if span < 90.0 else 0.0
            for deg in np.arange(start - (span if span < 90.0 else 0.0), start + span, step):
                n = score(float(deg))
                if n < best[1]:
                    best = (round(float(deg), 2), n)
    angle, remaining = best
    tolerated = max(5, len(kept) // 2000)  # a handful of inherent near-contacts is data, thousands are the rotation
    if remaining > tolerated and remaining > 0.1 * as_read:
        raise SystemExit(f"the coarse-grained frame is not periodic in its {cell[0]:.3f} x {cell[1]:.3f} nm box: {as_read} bead pairs of "
                         f"different molecules overlap (< {FRAME_OVERLAP_NM} nm) once wrapped, and no rotation about z removes them "
                         f"({remaining} remain at {angle} degrees); its coordinates do not belong to the box line")
    if as_read <= tolerated or abs(angle) < 0.05 or remaining >= 0.5 * as_read:
        return atoms, {"rotation_about_z_deg": 0.0, "overlapping_pairs_as_read": as_read, "overlapping_pairs_after": as_read,
                       "beads_tested": len(kept), "square_box": square}
    c, s_ = math.cos(math.radians(angle)), math.sin(math.radians(angle))
    fixed = [dict(a, x=float(centre[0] + (a["x"] - centre[0]) * c - (a["y"] - centre[1]) * s_),   # the same rotation the scan scored
                  y=float(centre[1] + (a["x"] - centre[0]) * s_ + (a["y"] - centre[1]) * c)) for a in atoms]
    after = frame_overlaps(np.array([[fixed[i]["x"], fixed[i]["y"], fixed[i]["z"]] for i in kept]), owner, cell3)
    if after != remaining:
        raise SystemExit(f"frame alignment is inconsistent: {remaining} overlaps expected after rotating {angle} degrees, {after} found")
    return fixed, {"rotation_about_z_deg": angle, "overlapping_pairs_as_read": as_read, "overlapping_pairs_after": after,
                   "beads_tested": len(kept), "square_box": square}
