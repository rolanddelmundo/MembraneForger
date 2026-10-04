#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Embedding an all-atom complex into a membrane simulated with another protein.
#//=============================================================
"""Place an all-atom complex into a pre-equilibrated Martini 3 membrane, trim the patch, and edit its lipids."""
import math
from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

from .config import AMINO, MIN_Z_PAD_TOTAL_NM, SLAB_Z_PAD_NM
from .structio import element, wrap, xyz_nm

__all__ = ['LIPID_NAMES', 'LIPID_ALIASES', 'CONVERTIBLE', 'BELT_NM', 'OVERLAP_NM', 'HYDROPHOBIC', 'CHARGED',
           'EXPOSURE_RADIUS_NM', 'BURIED_ABOVE', 'lipid_name', 'hydrophobic_belt', 'periodic_mean', 'embed_complex',
           'trim_membrane', 'edit_lipids', 'check_box_z']

# User-facing Martini 3 lipid names (INSANE spelling) and how they are spelled inside the classifier.
LIPID_NAMES = ("CHOL", "POPC", "DOPC", "POPE", "DOPE", "POPS", "DOPS", "PSM", "DPG3", "SAP6")
LIPID_ALIASES = {"DPG3": "GM3", "GM3": "GM3", "PIP2": "SAP6", "SAPI25": "SAP6", "CHL1": "CHOL", "DPSM": "PSM"}
# Lipids with the same twelve-bead Martini 3 layout: one can be turned into another by renaming its beads.
CONVERTIBLE = ("POPC", "DOPC", "POPE", "DOPE", "POPS", "DOPS", "PSM")
BELT_NM = 3.0  # hydrophobic thickness of the bilayer: the window that locates the transmembrane region
OVERLAP_NM = 0.40  # a lipid with a bead this close to a protein heavy atom is removed before backmapping
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "CYS", "GLY", "PRO"}
CHARGED = {"ASP", "GLU", "LYS", "ARG"}
# A residue is on the surface when at most this many protein heavy atoms lie within this radius of its side chain.
# On the example receptor-G protein complex surface side chains have 40-60 such neighbours and buried ones 80-110.
EXPOSURE_RADIUS_NM, BURIED_ABOVE = 0.8, 60


def lipid_name(name: str) -> str:
    """Return the classifier's spelling of a user-given lipid name, refusing names that are not supported."""
    upper = name.strip().upper()
    canonical = LIPID_ALIASES.get(upper, upper)
    if canonical not in {LIPID_ALIASES.get(n, n) for n in LIPID_NAMES}:
        raise SystemExit(f"unknown lipid {name}; choose from {', '.join(LIPID_NAMES)}")
    return canonical


def hydrophobic_belt(aa_atoms: list[dict], width_nm: float = BELT_NM, step_nm: float = 0.05,
                     exposure_radius_nm: float = EXPOSURE_RADIUS_NM, buried_above: int = BURIED_ABOVE) -> dict:
    """Find the z (A) of the slab of width `width_nm` whose exposed residues are the most hydrophobic."""
    # The all-atom input must already have its membrane normal along z; this only decides where along z the
    # bilayer centre lies, which is where the membrane's midplane is put. Only surface residues count (a side chain
    # with at most `buried_above` protein heavy atoms within `exposure_radius_nm` of its centroid): the surface of a
    # transmembrane segment is hydrophobic, the surface of a soluble domain is not, whatever their cores hold.
    heavy = [a for a in aa_atoms if a["resname"] in AMINO and element(a["atom"]) != "H"]
    residues = {}
    for a in heavy:
        residues.setdefault((a["chain"], a.get("segid", ""), a["resid"], a["resname"]), []).append(a)
    if len(residues) < 10:
        raise SystemExit("the all-atom input has fewer than 10 protein residues; cannot locate its transmembrane region")
    xyz = xyz_nm(heavy) / 10.0
    tree = cKDTree(xyz)
    centroids, kind = [], []
    for key, atoms in residues.items():
        side = [a for a in atoms if a["atom"] not in ("N", "CA", "C", "O", "OXT")] or [a for a in atoms if a["atom"] == "CA"] or atoms
        centroids.append(xyz_nm(side).mean(axis=0) / 10.0)
        kind.append(1 if key[3] in HYDROPHOBIC else -1 if key[3] in CHARGED else 0)
    centroids, kind = np.array(centroids), np.array(kind)
    neighbours = np.array([len(hits) for hits in tree.query_ball_point(centroids, exposure_radius_nm)])
    exposed = neighbours <= buried_above
    z = centroids[:, 2]
    order = np.argsort(z)
    z, weight = z[order], np.where(exposed[order], kind[order], 0)
    cumulative = np.concatenate([[0], np.cumsum(weight)])
    centres = np.arange(z.min() + width_nm / 2.0, z.max() - width_nm / 2.0 + step_nm, step_nm)
    if not len(centres):
        centres = np.array([0.5 * (z.min() + z.max())])
    low = np.searchsorted(z, centres - width_nm / 2.0, side="left")
    high = np.searchsorted(z, centres + width_nm / 2.0, side="right")
    scores = cumulative[high] - cumulative[low]
    best = int(np.argmax(scores))
    inside = weight[low[best]:high[best]]
    return {"z_a": round(float(centres[best]) * 10.0, 2), "width_nm": width_nm, "score": int(scores[best]),
            "hydrophobic": int((inside == 1).sum()), "charged": int((inside == -1).sum()),
            "residues": int(high[best] - low[best]), "exposed_residues": int(exposed.sum()),
            "exposure": {"radius_nm": exposure_radius_nm, "buried_above": buried_above}}


def periodic_mean(values: np.ndarray, length: float) -> float:
    """Mean of coordinates on a periodic axis of the given length, inside [0, length)."""
    angle = 2.0 * math.pi * values / length
    return (math.atan2(np.sin(angle).mean(), np.cos(angle).mean()) / (2.0 * math.pi) * length) % length


def embed_complex(aa_atoms: list[dict], cg_protein: list[list[dict]], membrane: list[dict], box: list[float],
                  slab: tuple[float, float], bilayer_z_a: float | None = None, overlap_nm: float = OVERLAP_NM) -> dict:
    """Put the complex where the frame's own protein was, with its hydrophobic belt on the midplane, and clear overlaps."""
    # Returns the identity rotation and the translation (A) applied to the complex, the placed atoms, the membrane
    # molecules that remain, and notes/metrics for the log and manifest. The complex is not rotated.
    cell = np.array(box)
    midplane = 0.5 * (slab[0] + slab[1])
    belt = hydrophobic_belt(aa_atoms) if bilayer_z_a is None else {"z_a": float(bilayer_z_a), "given": True}
    centre_z = belt["z_a"] / 10.0
    bb = np.array([[b["x"], b["y"], b["z"]] for res in cg_protein for b in res if b["atom"] == "BB"])
    inside = bb[(bb[:, 2] >= slab[0]) & (bb[:, 2] <= slab[1])] if len(bb) else bb
    if len(inside):
        target_xy = np.array([periodic_mean(inside[:, 0], cell[0]), periodic_mean(inside[:, 1], cell[1])])
        where = f"the centre of the frame's own transmembrane protein ({len(inside)} BB beads in the bilayer)"
    else:
        target_xy = 0.5 * cell[:2]
        where = "the centre of the membrane patch (the frame has no protein in the bilayer)"
    heavy = np.array([element(a["atom"]) != "H" for a in aa_atoms])
    xyz = xyz_nm(aa_atoms) / 10.0
    tm = heavy & (np.abs(xyz[:, 2] - centre_z) <= BELT_NM / 2.0)
    if tm.sum() < 10:
        raise SystemExit(f"fewer than 10 heavy atoms within {BELT_NM / 2.0:.1f} nm of the bilayer centre z = "
                         f"{belt['z_a']:.1f} A; orient the input with the membrane normal along z or pass --bilayer-z")
    shift_nm = np.array([target_xy[0] - xyz[tm, 0].mean(), target_xy[1] - xyz[tm, 1].mean(), midplane - centre_z])
    t = shift_nm * 10.0
    placed = [dict(a, x=a["x"] + t[0], y=a["y"] + t[1], z=a["z"] + t[2]) for a in aa_atoms]
    tree = cKDTree(wrap(xyz[heavy] + shift_nm, box), boxsize=box)
    kept, removed = [], Counter()
    for mol in membrane:
        beads = wrap(np.array([[b["x"], b["y"], b["z"]] for b in mol["beads"]]), box)
        if any(tree.query_ball_point(beads, overlap_nm)):
            removed[mol["cg"]] += 1
        else:
            kept.append(mol)
    if not kept:
        raise SystemExit("every membrane molecule overlaps the placed complex; the input is not oriented along z")
    notes = [f"bilayer centre of the all-atom input at z = {belt['z_a']:.1f} A "
             + ("(given)" if belt.get("given") else f"(hydrophobic belt: {belt['hydrophobic']} hydrophobic and "
                                                    f"{belt['charged']} charged residues in {belt['width_nm']:.1f} nm)"),
             f"complex moved by ({t[0]:+.1f}, {t[1]:+.1f}, {t[2]:+.1f}) A onto {where}, midplane z = {midplane:.2f} nm; "
             "it is not rotated: orient the input with the membrane normal along z before the build",
             f"{sum(removed.values())} membrane molecule(s) within {overlap_nm:.2f} nm of the complex removed: "
             + (", ".join(f"{name} {n}" for name, n in sorted(removed.items())) or "none")]
    metrics = {"mode": "embed", "bilayer_centre_input_A": belt["z_a"], "hydrophobic_belt": belt,
               "translation_A": [round(float(v), 3) for v in t], "target_xy_nm": [round(float(v), 3) for v in target_xy],
               "midplane_nm": round(midplane, 3), "overlap_cutoff_nm": overlap_nm, "removed_lipids": dict(removed),
               "membrane_molecules_kept": len(kept)}
    return {"R": np.eye(3), "t": t, "placed": placed, "membrane": kept, "removed": removed, "notes": notes, "metrics": metrics}


def trim_membrane(membrane: list[dict], placed: list[dict], box: list[float], xy_nm: tuple[float, float]) -> dict:
    """Cut the membrane patch down to xy_nm around the placed complex and move everything into the new cell."""
    # Molecules whose centre falls outside the new cell are dropped; the complex ends up centred in x and y.
    (x, y), cell = xy_nm, np.array(box)
    if x > cell[0] + 1e-6 or y > cell[1] + 1e-6:
        raise SystemExit(f"requested box {x * 10:.0f} x {y * 10:.0f} A exceeds the membrane patch "
                         f"{cell[0] * 10:.0f} x {cell[1] * 10:.0f} A; the box can only be trimmed, not enlarged")
    if x <= 0 or y <= 0:
        raise SystemExit("box edges must be positive")
    heavy = [a for a in placed if element(a["atom"]) != "H"]
    centre = xyz_nm(heavy).mean(axis=0)[:2] / 10.0
    shift = np.array([0.5 * x - centre[0], 0.5 * y - centre[1], 0.0])
    kept, removed = [], Counter()
    for mol in membrane:
        xyz = np.array([[b["x"], b["y"], b["z"]] for b in mol["beads"]]) + shift
        middle = xyz.mean(axis=0)
        image = np.array([-cell[0] * math.floor(middle[0] / cell[0]), -cell[1] * math.floor(middle[1] / cell[1]), 0.0])
        xyz, middle = xyz + image, middle + image
        if middle[0] < x and middle[1] < y:
            for bead, (bx, by, bz) in zip(mol["beads"], xyz):
                bead["x"], bead["y"], bead["z"] = float(bx), float(by), float(bz)
            kept.append(mol)
        else:
            removed[mol["cg"]] += 1
    if not kept:
        raise SystemExit("no membrane molecule is left inside the requested box")
    moved = [dict(a, x=a["x"] + shift[0] * 10.0, y=a["y"] + shift[1] * 10.0) for a in placed]
    note = (f"membrane trimmed to {x:.2f} x {y:.2f} nm around the complex: {sum(removed.values())} molecule(s) outside "
            f"removed ({', '.join(f'{n} {k}' for n, k in sorted(removed.items())) or 'none'}); lipids across the new "
            "periodic seam were not equilibrated together, so equilibrate the system before production")
    return {"membrane": kept, "placed": moved, "box": [x, y, float(cell[2])], "removed": removed, "note": note,
            "shift_A": [round(float(v) * 10.0, 3) for v in shift]}


def edit_lipids(membrane: list[dict], delete: list[str], add: str | None, mapping: dict) -> tuple[list[dict], dict]:
    """Remove every molecule of the named lipids or, with `add`, turn them into that lipid by renaming their beads."""
    names = {lipid_name(n) for n in delete}
    target = lipid_name(add) if add else None
    if target and not names:
        raise SystemExit("--addlipid replaces the lipids removed by --dellipid; say which lipid to replace with --dellipid")
    if target and (not {*names, target} <= set(CONVERTIBLE)):
        raise SystemExit(f"cannot turn {', '.join(sorted(names))} into {target}: only {', '.join(CONVERTIBLE)} share a "
                         "Martini 3 bead layout; cholesterol, DPG3 and SAP6 can be removed but not replaced")
    if target and target in names:
        raise SystemExit(f"--addlipid {target} is also named by --dellipid")
    if target and target not in mapping:
        raise SystemExit(f"the installed mapping has no residue {target}")
    kept, removed, converted = [], Counter(), Counter()
    for mol in membrane:
        if mol["cg"] not in names:
            kept.append(mol)
        elif target:
            position = {name: i for i, name in enumerate(mapping[mol["aa"]]["beads"])}
            beads = [dict(b, atom=mapping[target]["beads"][position[b["atom"]]], resname=target) for b in mol["beads"]]
            kept.append(dict(mol, cg=target, aa=target, beads=beads, dropped=[]))
            converted[mol["cg"]] += 1
        else:
            removed[mol["cg"]] += 1
    if not kept:
        raise SystemExit("every membrane molecule was removed by --dellipid")
    missing = sorted(n for n in names if not removed[n] and not converted[n])
    if missing:
        raise SystemExit(f"the membrane has no {', '.join(missing)} to remove")
    return kept, {"deleted": dict(removed), "converted_to_" + target if target else "converted": dict(converted)}


def check_box_z(placed: list[dict], slab: tuple[float, float], requested_z_nm: float) -> float:
    """Refuse a requested box z that cannot hold the placed complex, before the slow backmapping is paid for."""
    # Same rule as rebox_system (water padding above and below the complex, measured from the bilayer centre),
    # estimated from the coarse-grained bilayer; rebox_system re-checks on the all-atom membrane. Returns the need (nm).
    z = np.array([a["z"] for a in placed if element(a["atom"]) != "H"]) / 10.0
    centre = 0.5 * (slab[0] + slab[1])
    need = max(2.0 * (max(z.max() - centre, centre - z.min()) + SLAB_Z_PAD_NM),
               max(z.max(), slab[1]) - min(z.min(), slab[0]) + MIN_Z_PAD_TOTAL_NM)
    if requested_z_nm < need - 1e-6:
        raise SystemExit(f"requested box z {requested_z_nm * 10:.0f} A is below the {need * 10:.0f} A this complex needs "
                         f"({SLAB_Z_PAD_NM} nm of water above and below it); checked before backmapping")
    return need
