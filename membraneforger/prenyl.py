#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// E/Z configuration of the S-geranylgeranyl cysteine (CYSG / CYSM).
#//=============================================================
"""E/Z configuration of the geranylgeranyl group of CYSG (and CYSM, its methyl ester).

A geranylgeranyl group is all-E: at each of its three CH=C(CH3) double bonds (C2=C3, C7=C8, C12=C13) the chain carbons on the two
sides (C1/C5, C6/C10, C11/C15) lie trans, the dihedral about the bond is near 180 degrees. The torsion terms of the force field keep
a double bond planar but cannot turn a Z bond into an E one, so a Z group in the input stays Z through minimization and every later
simulation. The bundled GCGR-Gs frame (and the models derived from it) carried all three bonds in Z (dihedral 0 to 5 degrees).

review_prenyl reads the three dihedrals of every such residue; convert_to_e turns each Z bond into E by rotating everything beyond it,
as a rigid body, 180 degrees about the bond axis: bond lengths, bond angles and the planarity at the double bond are kept exactly, and
hydrogens move with their carbon. apply_prenyl_configuration is the pipeline policy (report / yes / fix / no).
"""
import numpy as np
from scipy.spatial import cKDTree

__all__ = ['PRENYL_RESIDUES', 'PRENYL_MODES', 'PRENYL_DOUBLE_BONDS', 'Z_BELOW_DEG', 'prenyl_dihedrals', 'review_prenyl', 'convert_to_e',
           'apply_prenyl_configuration']

PRENYL_RESIDUES = ("CYSG", "CYSM")
PRENYL_MODES = ("report", "yes", "fix", "no")
# (chain carbon before, the two double-bond carbons, chain carbon after); C4, C9 and C14 are the methyl branches, C17=C18 carries two methyls
PRENYL_DOUBLE_BONDS = (("C1", "C2", "C3", "C5"), ("C6", "C7", "C8", "C10"), ("C11", "C12", "C13", "C15"))
PRENYL_BONDS = (("SG", "C1"), ("C1", "C2"), ("C2", "C3"), ("C3", "C4"), ("C3", "C5"), ("C5", "C6"), ("C6", "C7"), ("C7", "C8"),
                ("C8", "C9"), ("C8", "C10"), ("C10", "C11"), ("C11", "C12"), ("C12", "C13"), ("C13", "C14"), ("C13", "C15"),
                ("C15", "C16"), ("C16", "C17"), ("C17", "C18"), ("C18", "C19"), ("C18", "C20"))
Z_BELOW_DEG = 90.0         # |dihedral| below this is Z, above it E (a planar bond is near 0 or 180)
H_BOND_A = 1.3             # a hydrogen belongs to the nearest carbon within this distance (A)
MIN_CLEARANCE_A = 2.0      # after a conversion no moved heavy atom may lie closer than this to a heavy atom of another residue (A)
HEAVY_NAMES = {n for bond in PRENYL_BONDS for n in bond}


def dihedral_deg(p0, p1, p2, p3) -> float:
    """Signed dihedral (degrees) of four points."""
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - b0 @ b1 * b1, b2 - b2 @ b1 * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def prenyl_dihedrals(coords: dict) -> list:
    """The three double-bond dihedrals (degrees) of one residue given {atom name: xyz}; None where an atom is missing."""
    out = []
    for names in PRENYL_DOUBLE_BONDS:
        out.append(dihedral_deg(*[np.asarray(coords[n], float) for n in names]) if all(n in coords for n in names) else None)
    return out


def prenyl_groups(atoms: list) -> list:
    """[(key, resname, [indices into atoms])] for every CYSG / CYSM residue, key = (chain, segid, resid)."""
    groups = {}
    for i, a in enumerate(atoms):
        if a["resname"] in PRENYL_RESIDUES:
            groups.setdefault((a["chain"], a.get("segid", ""), a["resid"], a["resname"]), []).append(i)
    return [((k[0], k[1], k[2]), k[3], idx) for k, idx in groups.items()]


def review_prenyl(atoms: list) -> list:
    """One record per CYSG / CYSM residue: its three dihedrals, the Z bonds, and any prenyl atom the residue lacks."""
    records = []
    for key, resname, idx in prenyl_groups(atoms):
        coords = {atoms[i]["atom"]: np.array([atoms[i]["x"], atoms[i]["y"], atoms[i]["z"]]) for i in idx}
        dih = prenyl_dihedrals(coords)
        z = [f"{b[1]}={b[2]}" for b, d in zip(PRENYL_DOUBLE_BONDS, dih) if d is not None and abs(d) < Z_BELOW_DEG]
        records.append({"residue": f"{key[0] or '-'}:{resname}{key[2]}", "dihedral_deg": [None if d is None else round(d, 1) for d in dih],
                        "z_bonds": z, "missing": sorted(HEAVY_NAMES - set(coords))})
    return records


def distal_atoms(coords: dict, first: str, second: str) -> set:
    """Names of the atoms beyond the bond first-second (second included): heavy atoms by the prenyl bond list, hydrogens with their carbon."""
    seen, stack = {second}, [second]
    while stack:
        x = stack.pop()
        for a, b in PRENYL_BONDS:
            for p, q in ((a, b), (b, a)):
                if p == x and q != first and q not in seen:
                    seen.add(q)
                    stack.append(q)
    heavy = [n for n in coords if n in HEAVY_NAMES]
    for n, xyz in coords.items():
        if n in HEAVY_NAMES:
            continue
        parent = min(heavy, key=lambda h: float(np.linalg.norm(coords[h] - xyz)))
        if float(np.linalg.norm(coords[parent] - xyz)) < H_BOND_A and parent in seen:
            seen.add(n)
    return seen


def convert_to_e(coords: dict) -> tuple:
    """Return (new coordinates, converted bond names) for one residue: every Z double bond rotated to E, from the proximal bond to the distal."""
    new = {n: np.asarray(v, float).copy() for n, v in coords.items()}
    converted = []
    for names in PRENYL_DOUBLE_BONDS:
        if not all(n in new for n in names):
            continue
        if abs(dihedral_deg(*[new[n] for n in names])) >= Z_BELOW_DEG:
            continue
        _, b, c, _ = names
        axis = new[c] - new[b]
        axis /= np.linalg.norm(axis)
        for n in distal_atoms(new, b, c) - {c}:
            r = new[n] - new[b]
            new[n] = new[b] + 2.0 * float(r @ axis) * axis - r      # rotation by 180 degrees about the bond axis
        converted.append(f"{b}={c}")
    return new, converted


def apply_prenyl_configuration(atoms: list, mode: str, source: str) -> tuple:
    """Policy for the all-atom input: returns (atoms, [(level, text)], record).

    report  list every Z double bond as a warning and go on (the default, as in earlier releases)
    yes     stop the build on a Z double bond
    fix     convert each Z double bond to E by the rigid rotation of convert_to_e, then check the moved atoms against the rest of the input
    no      do nothing
    """
    if mode not in PRENYL_MODES:
        raise SystemExit(f"--prenyl-configuration must be one of {PRENYL_MODES}, not {mode!r}")
    if mode == "no":
        return atoms, [], {"mode": mode}
    records = review_prenyl(atoms)
    bad = [r for r in records if r["z_bonds"]]
    notes = []
    for r in records:
        if r["missing"]:
            notes.append(("WARN", f"{source}: {r['residue']} lacks prenyl atom(s) {', '.join(r['missing'])}; its double bonds were not all reviewed"))
    record = {"mode": mode, "residues": len(records), "z_residues": [r["residue"] for r in bad], "reviewed": records}
    if not bad:
        if records:
            notes.append(("INFO", f"{source}: {len(records)} geranylgeranyl group(s), all double bonds E"))
        return atoms, notes, record
    listing = "; ".join(f"{r['residue']} {', '.join(r['z_bonds'])} (dihedral {r['dihedral_deg']})" for r in bad[:5])
    if mode == "yes":
        raise SystemExit(f"all-atom input {source}: geranylgeranyl double bond(s) in Z configuration (a natural prenyl group is all-E): {listing}; "
                         "rerun with --prenyl-configuration fix to convert them, or correct the input")
    if mode == "report":
        notes.append(("WARN", f"{source}: geranylgeranyl double bond(s) in Z configuration (a natural prenyl group is all-E and minimization cannot "
                              f"change it): {listing}; --prenyl-configuration fix converts them, yes stops the build"))
        return atoms, notes, record
    conversions = []
    for key, resname, idx in prenyl_groups(atoms):
        coords = {atoms[i]["atom"]: np.array([atoms[i]["x"], atoms[i]["y"], atoms[i]["z"]]) for i in idx}
        new, converted = convert_to_e(coords)
        if not converted:
            continue
        own = set(idx)
        others = [j for j in range(len(atoms)) if j not in own and not atoms[j]["atom"].startswith(("H", "Q"))
                  and atoms[j]["resname"] not in ("TIP3", "SOL")]
        if others:
            tree = cKDTree(np.array([[atoms[j]["x"], atoms[j]["y"], atoms[j]["z"]] for j in others]))
            moved = [i for i in idx if atoms[i]["atom"] in HEAVY_NAMES and not np.allclose(new[atoms[i]["atom"]], coords[atoms[i]["atom"]])]
            close = min((float(tree.query(new[atoms[i]["atom"]])[0]) for i in moved), default=99.0)
            if close < MIN_CLEARANCE_A:
                raise SystemExit(f"all-atom input {source}: converting {resname}{key[2]} to E leaves a prenyl atom {close:.2f} A "
                                 f"from another residue (limit {MIN_CLEARANCE_A} A); re-place the group in the input")
        for i in idx:
            x, y, z = new[atoms[i]["atom"]]
            atoms[i]["x"], atoms[i]["y"], atoms[i]["z"] = float(x), float(y), float(z)
        conversions.append(f"{key[0] or '-'}:{resname}{key[2]} {', '.join(converted)}")
    record["converted"] = conversions
    notes.append(("WARN", f"{source}: geranylgeranyl double bond(s) converted from Z to E by a rigid 180 degree rotation of the distal group "
                          f"(bond lengths and angles unchanged): {'; '.join(conversions)}"))
    return atoms, notes, record
