#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Geometric stereochemistry of GM3: the sugar stereocentres and ceramide trans bonds the topology restraints do not cover.
#//=============================================================
"""Check every GM3 (GLPA) in a structure for inverted sugar stereocentres and cis ceramide bonds.

The CHARMM36 GLPA topology restrains, under DIHRES, the ceramide C2/C3 centres, its C4=C5 double bond, its amide
C2S-NF-C1F-C2F, the three pyranose ring chairs and the sialic acid C7/C8 centres: 24 rows. None of them covers glucose C1-C5, galactose C1-C5
or sialic acid C2 and C4-C6, and a sugar carbon that is inverted can still sit in a correct chair, so an audit that
reads only those rows passes it. This module tests the 16 sugar stereocentres and the two ceramide trans bonds
geometrically, with the chirality definitions of membraneforger/backmap_data/map.dat (the ones backmapping restrains and reviews):
a centre is wrong when (target - centre) . ((D - C) x (E - D)) < 0 and a trans bond is wrong when its dihedral is
within 90 degrees of 0. Neither can change in a classical MD run, so one frame decides.

Every trans bond is also placed in one of three bands, reported per bond: trans (|dihedral| >= 150), twisted (between 30 and
150) and cis (<= 30). A bond is refused below 90 degrees (cis, or twisted past half-way); the 90-150 band is reported, not
refused, because a backmapped membrane starts a few tenths of a percent of its GM3 and PSM amides there. The GLPA amide
C2S-NF-C1F-C2F has a CHARMM36 torsion of zero force constant, so only its DIHRES row (6 4 56 58) holds it trans in emres.

Coordinates in a .gro or a PDB with a CRYST1 record may be wrapped into the cell atom by atom (mdrun's wrapping differs
between builds), so every molecule is first made whole: each atom is moved to the periodic image nearest the previous
atom of the molecule. A test on raw wrapped coordinates reads a split molecule as inverted - a false positive.

The six inositol carbons of SAPI25 (PIP2) are tested the same way at the end of this module (check_pi_stereo).

GLPA is read by atom order (its three sugars reuse the names C1..C6), GM3 by atom name; in a .gro the molecule is
four consecutive CHARMM residues (CER160, BGLC, BGAL, ANE5AC).
"""
import re
from collections import Counter
from pathlib import Path

import numpy as np

from .config import GM3_XML_TO_GLPA

__all__ = ['GM3_POSITIONS', 'GM3_TRANS_NAMES', 'GLPA_PARTS', 'TRANS_MIN_DEG', 'CIS_MAX_DEG', 'REFUSE_BELOW_DEG', 'trans_band',
           'map_section', 'gm3_chiral_definitions', 'gm3_trans_definitions', 'dihedral', 'read_structure_atoms', 'make_whole',
           'read_gm3_molecules', 'gm3_named', 'wrong_gm3_centres', 'check_gm3_stereo', 'PI_RESNAME', 'PI_RESNAMES', 'PI_POSITIONS',
           'pi_chiral_definitions', 'read_pi_molecules', 'check_pi_stereo']

# GM3 centre name -> sugar position, for the report (the 16 sugar stereocentres).
GM3_POSITIONS = {"C1": "Glc C1", "C2": "Glc C2", "C3": "Glc C3", "C4": "Glc C4", "C5": "Glc C5",
                 "C7": "Gal C1", "C8": "Gal C2", "C10": "Gal C3", "C11": "Gal C4", "C12": "Gal C5",
                 "C15": "Neu5Ac C2", "C17": "Neu5Ac C4", "C20": "Neu5Ac C5", "C18": "Neu5Ac C6", "C19": "Neu5Ac C7",
                 "C21": "Neu5Ac C8"}
GM3_TRANS_NAMES = {("C3S", "C4S", "C5S", "C6S"): "ceramide C4=C5", ("C2S", "NF", "C1F", "C2F"): "ceramide amide"}
GLPA_PARTS = ("CER16", "BGLC", "BGAL", "ANE5")  # .gro truncates CER160 and ANE5AC to five characters
TRANS_MIN_DEG, CIS_MAX_DEG, REFUSE_BELOW_DEG = 150.0, 30.0, 90.0  # trans / twisted / cis bands of |dihedral|; refused below the last


def trans_band(phi: float) -> str:
    """The band of a bond that should be trans: "trans" (|phi| >= 150), "cis" (<= 30) or "twisted" in between."""
    return "trans" if abs(phi) >= TRANS_MIN_DEG else "cis" if abs(phi) <= CIS_MAX_DEG else "twisted"


def map_section(data: Path, name: str, residue: str = "GM3") -> list:
    """The rows of one [ name ] section of the block of one residue (GM3 unless given) of map.dat."""
    text = (data / "map.dat").read_text()
    start = re.search(rf"(?m)^RESI[ \t]+{re.escape(residue)}[ \t]*\r?$", text)
    if start is None:
        raise SystemExit(f"{data / 'map.dat'} has no RESI {residue} block")
    block = text[start.end():].split("RESI ")[0]
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


def read_structure_atoms(path: Path) -> tuple[list, np.ndarray | None]:
    """(residue key, residue name, atom name, xyz in A) for every atom of a .gro or .pdb, in file order, and the cell.

    The cell is the rectangular box diagonal in A (the .gro box line or the PDB CRYST1 record), or None without one."""
    lines = path.read_text(errors="replace").splitlines()
    if path.suffix == ".gro":
        count = int(lines[1])
        box = np.array([float(x) for x in lines[2 + count].split()[:3]]) * 10.0
        return [((int(l[0:5]), l[5:10].strip()), l[5:10].strip(), l[10:15].strip(),
                 np.array([float(l[20:28]), float(l[28:36]), float(l[36:44])]) * 10.0) for l in lines[2:2 + count]], \
            box if (box > 0).all() else None
    cryst1 = next((l for l in lines if l.startswith("CRYST1")), None)
    box = np.array([float(cryst1[6:15]), float(cryst1[15:24]), float(cryst1[24:33])]) if cryst1 else None
    return [((l[21], int(l[22:26]), l[17:21].strip()), l[17:21].strip(), l[12:16].strip(),
             np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
            for l in lines if l.startswith(("ATOM", "HETATM"))], box if box is not None and (box > 0).all() else None


def make_whole(atoms: list, box: np.ndarray | None) -> list:
    """One molecule's [(atom name, xyz)] with every atom moved to the periodic image nearest the previous atom.

    Consecutive atoms of a lipid topology sit within a few bonds of each other, far under half the cell, so chaining
    the minimum image from atom to atom reassembles a molecule that mdrun wrapped across the boundary. Without a box
    the atoms are returned as they are."""
    if box is None:
        return atoms
    whole, previous = [], None
    for name, xyz in atoms:
        if previous is not None:
            xyz = xyz - box * np.round((xyz - previous) / box)
        whole.append((name, xyz))
        previous = xyz
    return whole


def read_gm3_molecules(path: Path) -> dict:
    """(molecule number, first residue number, resname) -> [(atom name, xyz in A)] for every GM3 / GLPA molecule, made whole."""
    (atoms, box), molecules, k = read_structure_atoms(path), {}, 0
    while k < len(atoms):
        key, resname = atoms[k][0], atoms[k][1]
        if resname in ("GLPA", "GM3"):
            run = [a for a in atoms[k:k + len(GM3_XML_TO_GLPA) + 1] if a[0] == key]
            molecules[(len(molecules) + 1, key[-2], resname)] = make_whole([(a[2], a[3]) for a in run], box)
            k += len(run)
        elif resname.startswith(GLPA_PARTS[0]):
            run = atoms[k:k + len(GM3_XML_TO_GLPA)]
            parts = list(dict.fromkeys(a[1] for a in run))
            if len(parts) != 4 or not all(n.startswith(part) for n, part in zip(parts, GLPA_PARTS)):
                raise SystemExit(f"{path}: the GLPA molecule starting at {resname} {key} is not CER160, BGLC, BGAL, ANE5AC")
            molecules[(len(molecules) + 1, key[-2], "GLPA")] = make_whole([(a[2], a[3]) for a in run], box)
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
    cis = [d for d in trans if abs(dihedral(*(named[a] for a in d))) < REFUSE_BELOW_DEG]  # cis, or twisted past half-way
    return inverted, cis


def check_gm3_stereo(path: Path, data: Path) -> tuple[dict, list]:
    """Audit every GM3 in a structure: counts per centre and bond, the affected molecules, and the failure messages."""
    molecules = read_gm3_molecules(path)
    chirals, trans = gm3_chiral_definitions(data), gm3_trans_definitions(data)
    report = {"molecules": len(molecules), "stereocentres_per_molecule": len(chirals), "trans_bonds_per_molecule": len(trans),
              "wrong_by_centre": {}, "cis_by_bond": {}, "bands_by_bond": {}, "affected_molecules": []}
    if not molecules:
        return report, []
    flipped, cis_count, affected, bands = Counter(), Counter(), [], {d: Counter() for d in trans}
    for key, atoms in molecules.items():
        named = gm3_named(atoms, key[-1])
        inverted, cis = wrong_gm3_centres(named, chirals, trans)
        for d in trans:
            bands[d][trans_band(dihedral(*(named[a] for a in d)))] += 1
        flipped.update(inverted)
        cis_count.update(cis)
        if inverted or cis:
            affected.append({"molecule": key[0], "first_residue": key[1], "inverted": [GM3_POSITIONS[c] for c in inverted],
                             "cis": [GM3_TRANS_NAMES.get(d, "-".join(d)) for d in cis]})
    report["wrong_by_centre"] = {GM3_POSITIONS[c]: flipped.get(c, 0) for c in GM3_POSITIONS}
    report["cis_by_bond"] = {GM3_TRANS_NAMES.get(d, "-".join(d)): cis_count.get(d, 0) for d in trans}
    report["bands_by_bond"] = {GM3_TRANS_NAMES.get(d, "-".join(d)): {b: bands[d].get(b, 0) for b in ("trans", "twisted", "cis")}
                               for d in trans}
    report["affected_molecules"] = affected
    fails = []
    if affected:
        worst = affected[0]
        fails.append(f"GM3: {len(affected)} of {len(molecules)} molecules have an inverted sugar stereocentre or a cis or half-twisted "
                     "(under 90 degrees) ceramide bond "
                     f"({sum(flipped.values())} centres, {sum(cis_count.values())} bonds), e.g. molecule {worst['molecule']} "
                     f"(residue {worst['first_residue']}): {', '.join(worst['inverted'] + worst['cis'])}; a configuration never "
                     "corrects itself in MD")
    return report, fails


# SAPI25 (PIP2) inositol ring: the six ring carbons carry the stereocentres. The CHARMM-GUI topology restrains only the ring
# chair (six +-60 rows, which an exchange of a hydrogen and a hydroxyl or phosphate at one carbon leaves satisfied), so an
# inverted inositol carbon passes an audit that reads the DIHRES rows. map.dat defines each centre like those of GM3.
PI_RESNAME = "SAPI25"
PI_RESNAMES = ("SAPI25", "SAPI2", "SAPI")  # as written in a topology, in a .gro (five characters) and in a PDB (four)
PI_POSITIONS = {"C11": "inositol C1", "C12": "inositol C2", "C13": "inositol C3", "C14": "inositol C4", "C15": "inositol C5",
                "C16": "inositol C6"}


def pi_chiral_definitions(data: Path) -> list:
    """The inositol chirality definitions [target, centre, C, D, E] of the SAPI25 block of map.dat, one per ring carbon."""
    rows = [d for d in map_section(data, "chiral", PI_RESNAME) if d[1] in PI_POSITIONS]
    missing = sorted(set(PI_POSITIONS) - {d[1] for d in rows})
    if missing:
        raise SystemExit(f"{data / 'map.dat'} has no chirality definition for SAPI25 centre(s) {missing}")
    return rows


def read_pi_molecules(path: Path) -> dict:
    """(molecule number, residue number) -> {atom name: xyz in A} for every SAPI25 residue of a .gro or .pdb, made whole."""
    atoms, box = read_structure_atoms(path)
    molecules, last = {}, None
    for key, resname, name, xyz in atoms:
        if resname not in PI_RESNAMES:
            last = None
            continue
        if key != last or name in molecules[next(reversed(molecules))]["names"]:  # a new residue (or the same number again)
            molecules[(len(molecules) + 1, key[-2])] = {"names": set(), "atoms": []}
        current = molecules[next(reversed(molecules))]
        current["names"].add(name)
        current["atoms"].append((name, xyz))
        last = key
    return {k: dict(make_whole(v["atoms"], box)) for k, v in molecules.items()}


def check_pi_stereo(path: Path, data: Path) -> tuple[dict, list]:
    """Audit every SAPI25 in a structure: inverted inositol stereocentres per centre, the affected molecules and the failures."""
    molecules = read_pi_molecules(path)
    chirals = pi_chiral_definitions(data)
    report = {"molecules": len(molecules), "stereocentres_per_molecule": len(chirals), "wrong_by_centre": {}, "affected_molecules": []}
    if not molecules:
        return report, []
    flipped, affected = Counter(), []
    for key, named in molecules.items():
        absent = [a for d in chirals for a in d if a not in named]
        if absent:
            raise SystemExit(f"{path}: SAPI25 molecule {key[0]} (residue {key[1]}) has no atom(s) {sorted(set(absent))}")
        inverted, _ = wrong_gm3_centres(named, chirals, [])
        flipped.update(inverted)
        if inverted:
            affected.append({"molecule": key[0], "first_residue": key[1], "inverted": [PI_POSITIONS[c] for c in inverted]})
    report["wrong_by_centre"] = {PI_POSITIONS[c]: flipped.get(c, 0) for c in PI_POSITIONS}
    report["affected_molecules"] = affected
    fails = []
    if affected:
        worst = affected[0]
        fails.append(f"SAPI25: {len(affected)} of {len(molecules)} molecules have an inverted inositol stereocentre "
                     f"({sum(flipped.values())} centres), e.g. molecule {worst['molecule']} (residue {worst['first_residue']}): "
                     f"{', '.join(worst['inverted'])}; a configuration never corrects itself in MD")
    return report, fails
