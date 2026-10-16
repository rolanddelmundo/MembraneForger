#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Geometric stereochemistry of GM3: the sugar stereocentres and ceramide trans bonds the topology restraints do not cover.
#//=============================================================
"""Check every GM3 (GLPA) in a structure for inverted sugar stereocentres and cis ceramide bonds.

The CHARMM36 GLPA topology restrains, under DIHRES, the ceramide C2/C3 centres, its C4=C5 double bond, the three
pyranose ring chairs and the sialic acid C7/C8 centres: 23 rows. None of them covers glucose C1-C5, galactose C1-C5
or sialic acid C2 and C4-C6, and a sugar carbon that is inverted can still sit in a correct chair, so an audit that
reads only those rows passes it. This module tests the 16 sugar stereocentres and the two ceramide trans bonds
geometrically, with the chirality definitions of backmap_data/map.dat (the ones backmapping restrains and reviews):
a centre is wrong when (target - centre) . ((D - C) x (E - D)) < 0 and a trans bond is wrong when its dihedral is
within 90 degrees of 0. Neither can change in a classical MD run, so one frame decides.

GLPA is read by atom order (its three sugars reuse the names C1..C6), GM3 by atom name; in a .gro the molecule is
four consecutive CHARMM residues (CER160, BGLC, BGAL, ANE5AC).
"""
from collections import Counter
from pathlib import Path

import numpy as np

from .config import GM3_XML_TO_GLPA

__all__ = ['GM3_POSITIONS', 'GM3_TRANS_NAMES', 'GLPA_PARTS', 'map_section', 'gm3_chiral_definitions', 'gm3_trans_definitions',
           'dihedral', 'read_structure_atoms', 'read_gm3_molecules', 'gm3_named', 'wrong_gm3_centres', 'check_gm3_stereo']

# GM3 centre name -> sugar position, for the report (the 16 sugar stereocentres).
GM3_POSITIONS = {"C1": "Glc C1", "C2": "Glc C2", "C3": "Glc C3", "C4": "Glc C4", "C5": "Glc C5",
                 "C7": "Gal C1", "C8": "Gal C2", "C10": "Gal C3", "C11": "Gal C4", "C12": "Gal C5",
                 "C15": "Neu5Ac C2", "C17": "Neu5Ac C4", "C20": "Neu5Ac C5", "C18": "Neu5Ac C6", "C19": "Neu5Ac C7",
                 "C21": "Neu5Ac C8"}
GM3_TRANS_NAMES = {("C3S", "C4S", "C5S", "C6S"): "ceramide C4=C5", ("C2S", "NF", "C1F", "C2F"): "ceramide amide"}
GLPA_PARTS = ("CER16", "BGLC", "BGAL", "ANE5")  # .gro truncates CER160 and ANE5AC to five characters


def map_section(data: Path, name: str) -> list:
    """The rows of one [ name ] section of the GM3 block of map.dat."""
    block = (data / "map.dat").read_text().split("RESI GM3")[1].split("RESI ")[0]
    body = block.split(f"[ {name} ]")[1].split("[")[0] if f"[ {name} ]" in block else ""
    return [line.split() for line in body.splitlines() if line.strip()]


def gm3_chiral_definitions(data: Path) -> list:
    """The sugar chirality definitions [target, centre, C, D, E] of map.dat, one per stereocentre of GM3_POSITIONS."""
    rows = [d for d in map_section(data, "chiral") if d[1] in GM3_POSITIONS]
    missing = sorted(set(GM3_POSITIONS) - {d[1] for d in rows})
    if missing:
        raise SystemExit(f"{data / 'map.dat'} has no chirality definition for GM3 centre(s) {missing}")
    return rows


def gm3_trans_definitions(data: Path) -> list:
    """The [ trans ] dihedral definitions of GM3 in map.dat."""
    return [tuple(d) for d in map_section(data, "trans")]


def dihedral(p0, p1, p2, p3) -> float:
    """Dihedral angle in degrees."""
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - b0 @ b1 * b1, b2 - b2 @ b1 * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def read_structure_atoms(path: Path) -> list:
    """(residue key, residue name, atom name, xyz in A) for every atom of a .gro or .pdb, in file order."""
    lines = path.read_text(errors="replace").splitlines()
    if path.suffix == ".gro":
        return [((int(l[0:5]), l[5:10].strip()), l[5:10].strip(), l[10:15].strip(),
                 np.array([float(l[20:28]), float(l[28:36]), float(l[36:44])]) * 10.0) for l in lines[2:2 + int(lines[1])]]
    return [((l[21], int(l[22:26]), l[17:21].strip()), l[17:21].strip(), l[12:16].strip(),
             np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
            for l in lines if l.startswith(("ATOM", "HETATM"))]


def read_gm3_molecules(path: Path) -> dict:
    """(molecule number, first residue number, resname) -> [(atom name, xyz in A)] for every GM3 / GLPA molecule."""
    atoms, molecules, k = read_structure_atoms(path), {}, 0
    while k < len(atoms):
        key, resname = atoms[k][0], atoms[k][1]
        if resname in ("GLPA", "GM3"):
            run = [a for a in atoms[k:k + len(GM3_XML_TO_GLPA) + 1] if a[0] == key]
            molecules[(len(molecules) + 1, key[-2], resname)] = [(a[2], a[3]) for a in run]
            k += len(run)
        elif resname.startswith(GLPA_PARTS[0]):
            run = atoms[k:k + len(GM3_XML_TO_GLPA)]
            parts = list(dict.fromkeys(a[1] for a in run))
            if len(parts) != 4 or not all(n.startswith(part) for n, part in zip(parts, GLPA_PARTS)):
                raise SystemExit(f"{path}: the GLPA molecule starting at {resname} {key} is not CER160, BGLC, BGAL, ANE5AC")
            molecules[(len(molecules) + 1, key[-2], "GLPA")] = [(a[2], a[3]) for a in run]
            k += len(run)
        else:
            k += 1
    return molecules


def gm3_named(atoms: list, resname: str) -> dict:
    """GM3 (map.dat) names -> xyz for one molecule; GLPA atoms are matched by position in the CHARMM36 topology order."""
    if resname == "GM3":
        return dict(atoms)
    if len(atoms) != len(GM3_XML_TO_GLPA):
        raise SystemExit(f"a GLPA residue has {len(atoms)} atoms, expected {len(GM3_XML_TO_GLPA)}")
    for (_, glpa), (name, _) in zip(GM3_XML_TO_GLPA, atoms):
        if name != glpa:
            raise SystemExit(f"GLPA atom order differs from the CHARMM36 topology: {name} where {glpa} is expected")
    return {gm3: xyz for (gm3, _), (_, xyz) in zip(GM3_XML_TO_GLPA, atoms)}


def wrong_gm3_centres(named: dict, chirals: list, trans: list) -> tuple[list, list]:
    """The inverted stereocentres (centre names) and cis trans-bonds (definitions) of one molecule given by GM3 names."""
    inverted = [centre for target, centre, c, d, e in chirals
                if (named[target] - named[centre]) @ np.cross(named[d] - named[c], named[e] - named[d]) < 0]
    cis = [d for d in trans if abs(dihedral(*(named[a] for a in d))) < 90]
    return inverted, cis


def check_gm3_stereo(path: Path, data: Path) -> tuple[dict, list]:
    """Audit every GM3 in a structure: counts per centre and bond, the affected molecules, and the failure messages."""
    molecules = read_gm3_molecules(path)
    chirals, trans = gm3_chiral_definitions(data), gm3_trans_definitions(data)
    report = {"molecules": len(molecules), "stereocentres_per_molecule": len(chirals), "trans_bonds_per_molecule": len(trans),
              "wrong_by_centre": {}, "cis_by_bond": {}, "affected_molecules": []}
    if not molecules:
        return report, []
    flipped, cis_count, affected = Counter(), Counter(), []
    for key, atoms in molecules.items():
        inverted, cis = wrong_gm3_centres(gm3_named(atoms, key[-1]), chirals, trans)
        flipped.update(inverted)
        cis_count.update(cis)
        if inverted or cis:
            affected.append({"molecule": key[0], "first_residue": key[1], "inverted": [GM3_POSITIONS[c] for c in inverted],
                             "cis": [GM3_TRANS_NAMES.get(d, "-".join(d)) for d in cis]})
    report["wrong_by_centre"] = {GM3_POSITIONS[c]: flipped.get(c, 0) for c in GM3_POSITIONS}
    report["cis_by_bond"] = {GM3_TRANS_NAMES.get(d, "-".join(d)): cis_count.get(d, 0) for d in trans}
    report["affected_molecules"] = affected
    fails = []
    if affected:
        worst = affected[0]
        fails.append(f"GM3: {len(affected)} of {len(molecules)} molecules have an inverted sugar stereocentre or a cis ceramide bond "
                     f"({sum(flipped.values())} centres, {sum(cis_count.values())} bonds), e.g. molecule {worst['molecule']} "
                     f"(residue {worst['first_residue']}): {', '.join(worst['inverted'] + worst['cis'])}; a configuration never "
                     "corrects itself in MD")
    return report, fails
