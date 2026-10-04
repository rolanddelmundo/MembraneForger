#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Membrane backmapping through the mstool worker.
#//=============================================================
"""Drive the mstool worker, verify what it returns, and assemble membrane.pdb."""
import json
import os
import shutil
from collections import Counter, OrderedDict
from pathlib import Path

import numpy as np

from .audit import molecule_rings
from .martini import gro_safe_aliases, structure_resname
from .mstool_worker import TESTED_MSTOOL
from .runtools import log, run_command
from .structio import source_element, write_gro, write_pdb, xyz_nm

__all__ = ['WORKER', 'read_mapping', 'read_backmapped_table', 'check_backmapped_inventory', 'ROCK_RESNAME',
           'BOND_RADII_A', 'ring_centres', 'rock_atoms', 'read_isomer_review', 'backmap_membrane',
           'assemble_membrane_pdb']

WORKER = Path(__file__).resolve().parent / "mstool_worker.py"


def read_mapping(out: Path, work: Path, python: str, data: Path) -> dict:
    """Ask the mstool interpreter for its version and residue mapping, and refuse an untested mstool."""
    job = {"mode": "describe", "data": str(data), "result": str(work / "mapping.json")}
    (work / "describe.json").write_text(json.dumps(job))
    run_command(out, [python, WORKER, work / "describe.json"], produces=("mapping.json",), cwd=work)
    described = json.loads((work / "mapping.json").read_text())
    if described["mstool_version"] not in TESTED_MSTOOL:
        raise SystemExit(f"mstool {described['mstool_version']} at {described['mstool_path']} is not a tested version "
                         f"{TESTED_MSTOOL}; pass --mstool-python for a supported install")
    return described


def read_backmapped_table(path: Path) -> list[list[dict]]:
    """Read the worker's atom table into per-residue atom lists (A)."""
    residues, last = [], None
    for line in path.read_text().splitlines()[1:]:
        chain, resid, resname, atom, x, y, z = line.split("\t")
        if (chain, resid, resname) != last:
            residues.append([])
            last = (chain, resid, resname)
        residues[-1].append({"atom": atom, "altloc": " ", "resname": resname, "x": float(x), "y": float(y), "z": float(z)})
    return residues


def check_backmapped_inventory(residues: list[list[dict]], membrane: list[dict], mapping: dict) -> None:
    """Fail if backmapping lost, added or truncated any membrane molecule."""
    expected = Counter(mol["aa"] for mol in membrane)
    found = Counter(res[0]["resname"] for res in residues)
    if found != expected:
        raise SystemExit("backmapping changed the membrane inventory: " + ", ".join(
            f"{name} {found.get(name, 0)}/{expected.get(name, 0)}" for name in sorted(set(found) | set(expected))
            if found.get(name, 0) != expected.get(name, 0)))
    for res in residues:
        if sorted(a["atom"] for a in res) != sorted(mapping[res[0]["resname"]]["atoms"]):
            raise SystemExit(f"backmapped {res[0]['resname']} is incomplete: {len(res)} atoms, "
                             f"expected {len(mapping[res[0]['resname']]['atoms'])}")
    if not np.isfinite(np.array([[a["x"], a["y"], a["z"]] for res in residues for a in res])).all():
        raise SystemExit("backmapping produced non-finite coordinates")


ROCK_RESNAME = "RCK"  # neutral name: mstool must treat the complex as an obstacle, never as a mapped residue


BOND_RADII_A = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07}  # covalent radii for the ring search


def ring_centres(placed: list[dict]) -> list[tuple]:
    """Find the centre of every 5- or 6-membered ring of the placed complex from heavy-atom distances."""
    centres, residues = [], OrderedDict()
    for a in placed:
        if source_element(a) in BOND_RADII_A:
            residues.setdefault((a["chain"], a["segid"], a["resid"], a["resname"]), []).append(a)
    for atoms in residues.values():
        if len(atoms) < 5:
            continue
        xyz = xyz_nm(atoms)
        radii = np.array([BOND_RADII_A[source_element(a)] for a in atoms])
        close = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2) <= 1.25 * (radii[:, None] + radii[None, :])
        bonds = [(i, j) for i in range(len(atoms)) for j in range(i + 1, len(atoms)) if close[i, j]]
        centres += [tuple(xyz[list(ring)].mean(axis=0)) for ring in molecule_rings(bonds)]
    return centres


def rock_atoms(placed: list[dict]) -> list[dict]:
    """Turn the placed complex into mstool's obstacle: neutral residue name, plus one point at every ring centre."""
    # Neutral name: mstool writes the rock back under its original residue names and then asserts that every mapped
    # residue (ALA, CYS, ...) has all the atoms of its chiral centres; an input without hydrogens would abort there.
    # Ring centres: the obstacle is made of atom-sized spheres, so the middle of an aromatic or proline ring is
    # open and a lipid tail can be built straight through it. A point there keeps lipids out of the ring.
    rock = [dict(a, resname=ROCK_RESNAME) for a in placed]
    last = max((a["resid"] for a in placed), default=0)
    for n, (x, y, z) in enumerate(ring_centres(placed), 1):
        rock.append({"atom": "CRC", "altloc": " ", "resname": ROCK_RESNAME, "chain": "", "resid": (last + n) % 10000,
                     "x": float(x), "y": float(y), "z": float(z), "segid": "RING", "elem": "C"})
    return rock


def read_isomer_review(path: Path, chirality: Path | None = None) -> dict:
    """Read mstool's isomer summary plus our own count of flipped centres at well-formed chirality definitions."""
    # "chiral" from mstool includes centres whose definition is malformed; "chiral_well_formed" does not, and
    # together with "cistrans" it is what the build refuses (a wrong epimer or double bond is permanent).
    counts = {}
    if chirality is not None and chirality.is_file():
        flips = json.loads(chirality.read_text())
        counts["chiral_well_formed"] = sum(n for centres in flips.values() for n in centres.values())
        counts["chiral_well_formed_detail"] = flips
    if path.is_file():
        summary = path.read_text().split("In summary, the number of residues with the flipped isomers:")[-1]
        for line in summary.splitlines():
            key, _, value = line.partition(":")
            if value.strip().isdigit():
                counts[key.strip()] = int(value)
    return counts


def backmap_membrane(membrane: list[dict], placed: list[dict], box: list[float], mapping: dict, out: Path, work: Path,
                     python: str, data: Path, threads: int, nsteps: int, seed: int = 1) -> tuple[list[list[dict]], dict]:
    """Backmap the membrane with mstool around the placed all-atom complex, which is held as a rigid rock."""
    alias = gro_safe_aliases([mol["aa"] for mol in membrane])
    beads = [dict(b, resid=n, resname=alias.get(mol["aa"], mol["aa"])) for n, mol in enumerate(membrane, 1) for b in mol["beads"]]
    write_gro(beads, box, work / "cg_membrane.gro", "membrane beads for backmapping")
    rock = rock_atoms(placed)
    write_pdb(rock, work / "rock.pdb")
    log(out, f"backmapping obstacle: {len(placed)} atoms of the placed complex and {len(rock) - len(placed)} ring-centre points")
    job = {"mode": "backmap", "data": str(data), "structure": "cg_membrane.gro", "rock": "rock.pdb",
           "workdir": "mstool", "nsteps": nsteps, "result": "membrane_aa.tsv", "seed": seed,
           "rename": {f":{short}": f":{name}" for name, short in alias.items()}}
    shutil.rmtree(work / "mstool", ignore_errors=True)  # a previous attempt's work must not be reused
    (work / "backmap.json").write_text(json.dumps(job))
    run_command(out, [python, WORKER, work / "backmap.json"],
                env=dict(os.environ, OPENMM_CPU_THREADS=str(threads), OMP_NUM_THREADS=str(threads)),
                produces=("membrane_aa.tsv", "membrane_aa.tsv.chirality.json"), cwd=work)
    residues = read_backmapped_table(work / "membrane_aa.tsv")
    check_backmapped_inventory(residues, membrane, mapping)
    return residues, read_isomer_review(work / "mstool" / "log.txt", work / "membrane_aa.tsv.chirality.json")


def assemble_membrane_pdb(placed: list[dict], lipids: list[list[dict]], box: list[float], path: Path) -> int:
    """Write the placed complex followed by the backmapped membrane as one PDB with the CG cell; returns lipid atoms."""
    for name in {res[0]["resname"] for res in lipids}:
        structure_resname(name)
    taken = {a["chain"] for a in placed}
    letters = [c for c in "LMNOQSTUVWXYZDEFGHIJKABCPR0123456789" if c not in taken]
    if len(lipids) > 9999 * len(letters):
        raise SystemExit(f"{len(lipids)} membrane residues exceed the PDB residue numbering available")
    membrane = [dict(a, chain=letters[n // 9999], resid=n % 9999 + 1, segid="MEMB") for n, res in enumerate(lipids) for a in res]
    write_pdb(placed + membrane, path,
              f"CRYST1{box[0] * 10:9.3f}{box[1] * 10:9.3f}{box[2] * 10:9.3f}  90.00  90.00  90.00 P 1           1")
    return len(membrane)
