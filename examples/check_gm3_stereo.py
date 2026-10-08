#!/usr/bin/env python3
#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Check the sugar stereocentres and ceramide trans bonds of every GM3 (GLPA) in a built system.
#//=============================================================
"""Report GM3 ganglioside sugar stereocentres and ceramide trans bonds with the wrong configuration in a .gro or .pdb.

    python examples/check_gm3_stereo.py em.gro            # or any frame of a trajectory, e.g. from gmx trjconv
    python examples/check_gm3_stereo.py membrane.pdb

Reads GLPA residues (the CHARMM36 topology name) by atom order, and GM3 residues (the backmapping name) by atom
name. Each centre is tested with the corrected definitions in backmap_data/map.dat (all 16 sugar stereocentres:
glucose C1-C5, galactose C1-C5, sialic acid C2 and C4-C8), and the ceramide C4=C5 double bond and amide with its
[ trans ] definitions. Neither a stereocentre nor a C=C bond can invert in a classical MD run, so the result for
the first frame holds for the whole trajectory (a cis amide can isomerize, but only on very long time scales).
"""
import sys
from collections import Counter
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from membraneforger.config import GM3_XML_TO_GLPA  # noqa: E402

# GM3 centre -> sugar position, for the report
POSITION = {"C1": "Glc C1", "C2": "Glc C2", "C3": "Glc C3", "C4": "Glc C4", "C5": "Glc C5",
            "C7": "Gal C1", "C8": "Gal C2", "C10": "Gal C3", "C11": "Gal C4", "C12": "Gal C5",
            "C15": "Neu5Ac C2", "C17": "Neu5Ac C4", "C20": "Neu5Ac C5", "C18": "Neu5Ac C6", "C19": "Neu5Ac C7",
            "C21": "Neu5Ac C8"}


TRANS = {("C3S", "C4S", "C5S", "C6S"): "ceramide C4=C5", ("C2S", "NF", "C1F", "C2F"): "ceramide amide"}


def section(name: str) -> list:
    """The rows of one [ name ] section of the GM3 block of backmap_data/map.dat."""
    block = (REPO / "backmap_data/map.dat").read_text().split("RESI GM3")[1].split("RESI ")[0]
    body = block.split(f"[ {name} ]")[1].split("[")[0] if f"[ {name} ]" in block else ""
    return [line.split() for line in body.splitlines() if line.strip()]


def definitions() -> list:
    """The GM3 sugar chirality definitions [target, centre, C, D, E] of backmap_data/map.dat."""
    return [d for d in section("chiral") if d[1] in POSITION]


def dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - b0 @ b1 * b1, b2 - b2 @ b1 * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


# A GLPA molecule (CHARMM-GUI GM3) is written as four consecutive residues: ceramide, Glc, Gal, Neu5Ac.
GLPA_PARTS = ("CER16", "BGLC", "BGAL", "ANE5")  # .gro truncates CER160 and ANE5AC to five characters


def read_atoms(path: Path) -> list:
    """(residue key, residue name, atom name, xyz in A) for every atom, in file order."""
    lines = path.read_text(errors="replace").splitlines()
    if path.suffix == ".gro":
        return [((int(l[0:5]), l[5:10].strip()), l[5:10].strip(), l[10:15].strip(),
                 np.array([float(l[20:28]), float(l[28:36]), float(l[36:44])]) * 10.0) for l in lines[2:2 + int(lines[1])]]
    return [((l[21], int(l[22:26]), l[17:21].strip()), l[17:21].strip(), l[12:16].strip(),
             np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
            for l in lines if l.startswith(("ATOM", "HETATM"))]


def read_residues(path: Path) -> dict:
    """(molecule number, first residue number, resname) -> list of (name, xyz in A), in file order: GM3/GLPA
    residues, or GLPA split into its CHARMM sugars (which restart their residue numbers in every molecule)."""
    atoms, molecules, k = read_atoms(path), {}, 0
    while k < len(atoms):
        key, resname = atoms[k][0], atoms[k][1]
        if resname in ("GLPA", "GM3"):
            run = [a for a in atoms[k:k + len(GM3_XML_TO_GLPA) + 1] if a[0] == key]
            molecules[(len(molecules) + 1, key[-2], resname)] = [(a[2], a[3]) for a in run]
            k += len(run)
        elif resname.startswith(GLPA_PARTS[0]):
            run = atoms[k:k + len(GM3_XML_TO_GLPA)]
            parts = [p for p in dict.fromkeys(a[1] for a in run)]
            if len(parts) != 4 or not all(n.startswith(part) for n, part in zip(parts, GLPA_PARTS)):
                raise SystemExit(f"{path}: the GLPA molecule starting at {resname} {key} is not CER160, BGLC, BGAL, ANE5AC")
            molecules[(len(molecules) + 1, key[-2], "GLPA")] = [(a[2], a[3]) for a in run]
            k += len(run)
        else:
            k += 1
    return molecules


def named(atoms: list, resname: str) -> dict:
    """GM3 names -> xyz. GLPA atoms are named by position (GLPA reuses C1..C6 in all three sugars)."""
    if resname == "GM3":
        return dict(atoms)
    if len(atoms) != len(GM3_XML_TO_GLPA):
        raise SystemExit(f"a GLPA residue has {len(atoms)} atoms, expected {len(GM3_XML_TO_GLPA)}")
    for (_, glpa), (name, _) in zip(GM3_XML_TO_GLPA, atoms):
        if name != glpa:
            raise SystemExit(f"GLPA atom order differs from the CHARMM36 topology: {name} where {glpa} is expected")
    return {gm3: xyz for (gm3, _), (_, xyz) in zip(GM3_XML_TO_GLPA, atoms)}


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    path = Path(sys.argv[1])
    residues = read_residues(path)
    if not residues:
        print(f"{path}: no GLPA or GM3 residues")
        return 0
    flipped, per_residue = Counter(), Counter()
    for key, atoms in residues.items():
        x = named(atoms, key[-1])
        for target, centre, c, d, e in definitions():
            t, o, p, q, r = (x[n] for n in (target, centre, c, d, e))
            if (t - o) @ np.cross(q - p, r - q) < 0:
                flipped[centre] += 1
                per_residue[key] += 1
        for d in map(tuple, section("trans")):
            if abs(dihedral(*(x[a] for a in d))) < 90:
                flipped[d] += 1
                per_residue[key] += 1
    n = len(residues)
    print(f"{path}: {n} GM3 molecules, 16 sugar stereocentres and {len(section('trans'))} ceramide trans bonds each")
    for centre, position in POSITION.items():
        print(f"  {position:10s} ({centre:3s}): {flipped[centre]:3d} wrong ({100 * flipped[centre] / n:5.1f}%)")
    for d in map(tuple, section("trans")):
        print(f"  {TRANS.get(d, '-'.join(d)):16s}: {flipped[d]:3d} cis   ({100 * flipped[d] / n:5.1f}%)")
    affected = len(per_residue)
    print(f"molecules with at least one wrong centre or bond: {affected}/{n}"
          + ("" if not affected else " -> molecule(s) " + ", ".join(f"{k[0]} (residue {k[1]})" for k in sorted(per_residue))))
    return 1 if affected else 0


if __name__ == "__main__":
    sys.exit(main())
