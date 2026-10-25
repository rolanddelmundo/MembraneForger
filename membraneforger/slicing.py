#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// BOX = auto: slice the coarse-grained membrane to the complex plus a buffer, before backmapping.
#//=============================================================
"""Cut a smaller periodic sample out of the equilibrated coarse-grained membrane, lipid by lipid, under PBC.

The membrane is cropped BEFORE the expensive backmapping, so mstool only backmaps the lipids that are kept.

    window (per axis) = [min, max] of the solute + buffer on each side (or the requested --box width, centred)
    an axis whose window is as wide as the periodic cell is not cropped (the original periodicity is kept)
    new cell = window width (cropped axes) or the original cell width; z is left to the box stage

Selection is by ONE representative site per lipid, its headgroup anchor (lipids.LIPID_SLICE_ANCHORS: PO4, ROH, the
GM3 sugar centroid; otherwise the centre of geometry of the whole residue):

    1. the lipid is made whole under the periodic boundaries of the SOURCE cell;
    2. the periodic image of its anchor that falls into the half-open window  lower <= a < lower + width  is found.
       The window is never wider than the cell, so exactly one image qualifies: a lipid is selected at most once and
       a lipid exactly on a boundary is counted on one side only. This is the decision one gets by tiling the
       membrane 3 x 3 in x and y and cutting the window out of the tiling, without building the tiling;
    3. if the anchor is in the window the COMPLETE residue is kept, moved as one rigid body by the same shift;
    4. the window's lower corner becomes the origin, so every anchor lies in [0, width) of the new cell while tail
       beads may reach a little beyond it, as in any periodic system (mstool runs with pbc=True and GROMACS wraps).

The earlier rule, "keep a lipid only if ALL its beads lie inside the window", is kept for the diagnostics only. It
deleted every lipid with one tail bead across the window edge while the new cell kept the full window area, so a
cut membrane lost 15-30 % of its lipids per unit area (about 50 A^2 per lipid in the embedded frame became 56-59
A^2 in the cut, with the two leaflets depleted by different amounts inside one shared box), a defect that neither
minimization nor NPT box relaxation can repair. The report counts the lipids the old rule would have rejected
although their anchor was inside, so that the cause stays visible.

Crop position. The window is centred on the complex, but it may be shifted by a fraction of the buffer on a cropped
axis without the complex leaving its margin. A coarse deterministic search over such offsets picks the one whose
lipid counts best reproduce the reference leaflet density and composition (score = per-leaflet |% deviation of the
lipid count from the count the reference APL predicts for the accessible area| + composition distance in mole
percent + a penalty per hard seam overlap); ties go to the smallest offset, so the result is reproducible.

Seam. Cropping makes opposite window edges periodic neighbours that never were in the simulation: the tails of a
lipid kept on one edge now lie where the tails of the discarded lipids beyond the other edge used to be. The lateral
density is right on average, but the two edges do not interlock, so some bead pairs of different lipids end up closer
than any pair in the equilibrated frame. Deleting a lipid of every such pair (the earlier remedy) removes 5-10 % of
the membrane and recreates the packing defect; instead the seam is relaxed in place at the coarse-grained level: a
short steepest descent with a soft repulsion between beads of DIFFERENT lipids that the new periodicity brought
closer than seam_min_bead_nm, with the same repulsion from the placed complex, while every bead is held near its
position and every lipid keeps its shape through intramolecular distance restraints. Contacts that were already close
in the source frame are data and are never acted on. Only a lipid still in a hard-core overlap afterwards
(seam_hard_core_nm) is removed, and that is reported. The displacements are a few tenths of a nanometre, confined to
the seam, far smaller than the perturbation the seam itself is; the cut edges were never equilibrated together, so
the built system must still be equilibrated before production.

Orientation determines the reference frame and slicing operates in it; nothing here knows about PPM or OPM.
"""
from collections import Counter

import numpy as np
from scipy.spatial import Voronoi, cKDTree

from .config import Settings
from .lipids import anchor_bead, anchor_xyz, leaflet_of
from .structio import element, xyz_nm

__all__ = ['xy_periodic_check', 'periodic_image_gap', 'voronoi_shells', 'solute_extent', 'crop_windows',
           'image_molecule', 'make_whole', 'lipid_anchor', 'anchor_image_shift',
           'all_beads_inside', 'select_lipids', 'created_pairs', 'pair_keys', 'offset_slack', 'relax_seam',
           'remove_hard_core_overlaps',
           'composition_distance', 'score_offset', 'choose_offset', 'leaflet_counts', 'selection_diagnostics',
           'slice_membrane_cg']

def solute_extent(placed: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    """Lowest and highest x, y of every placed solute atom (protein and ligands), in nm."""
    xy = np.array([[a["x"], a["y"]] for a in placed], dtype=float) / 10.0
    return xy.min(axis=0), xy.max(axis=0)


def crop_windows(low: np.ndarray, high: np.ndarray, cell: np.ndarray, buffer_nm: float,
                 requested: tuple | None = None, min_buffer_nm: float = 0.0) -> dict:
    """Per-axis window: solute extent + buffer, or the requested width centred on the solute; full cell if it fits."""
    centre, windows = 0.5 * (low + high), []
    for axis in range(2):
        extent = float(high[axis] - low[axis])
        if requested is not None:
            width = float(requested[axis])
            if width > cell[axis] + 0.01:
                raise SystemExit(f"BOX {'xy'[axis]} = {width:.3f} nm is larger than the {cell[axis]:.3f} nm of the "
                                 "coarse-grained cell; a box cannot be larger than the membrane it is cut from")
            if width < cell[axis] - 0.01 and width < extent + 2.0 * min_buffer_nm:  # the margin rule applies only to a cut axis
                raise SystemExit(f"BOX {'xy'[axis]} = {width:.3f} nm leaves less than {min_buffer_nm} nm of membrane on each side "
                                 f"of the complex ({extent:.3f} nm wide); at least {extent + 2.0 * min_buffer_nm:.3f} nm is needed")
        else:
            width = extent + 2.0 * buffer_nm
        cropped = width < cell[axis] - 0.01
        windows.append({"axis": "xy"[axis], "extent_nm": extent, "width_nm": float(width if cropped else cell[axis]),
                        "cropped": bool(cropped), "lower_nm": float(centre[axis] - (width if cropped else cell[axis]) / 2.0),
                        "centre_nm": float(centre[axis])})
    return {"axes": windows, "lower": np.array([w["lower_nm"] for w in windows]),
            "size": np.array([w["width_nm"] for w in windows]), "centre": centre}


def image_molecule(beads_xy: np.ndarray, centre: np.ndarray, cell_xy: np.ndarray) -> np.ndarray:
    """Whole-cell shift (nm, x and y) that puts the molecule's centroid in the periodic image nearest the window centre."""
    delta = beads_xy.mean(axis=0) - centre
    return -cell_xy * np.round(delta / cell_xy)


def make_whole(xyz: np.ndarray, cell: np.ndarray, reference: int = 0) -> np.ndarray:
    """Unwrap one molecule under periodic boundaries so every bead sits in the image nearest the reference bead."""
    delta = xyz - xyz[reference]
    return xyz[reference] + delta - cell * np.round(delta / cell)


def lipid_anchor(mol: dict, cell: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Whole-molecule coordinates (nm) and the anchor position of one lipid (lipids.LIPID_SLICE_ANCHORS)."""
    names = [b["atom"] for b in mol["beads"]]
    xyz = make_whole(np.array([[b["x"], b["y"], b["z"]] for b in mol["beads"]], dtype=float), cell)
    return xyz, anchor_xyz(xyz, names, anchor_bead(mol["cg"]))


def anchor_image_shift(anchor_xy: np.ndarray, lower: np.ndarray, size: np.ndarray, cell_xy: np.ndarray) -> np.ndarray | None:
    """Whole-cell shift that puts the anchor into the half-open window [lower, lower + size) in x and y, or None."""
    # size <= cell on every axis, so at most one image can qualify: no lipid is ever selected twice.
    shift = -cell_xy * np.floor((anchor_xy - lower) / cell_xy)  # the image with lower <= a < lower + cell
    moved = anchor_xy + shift
    return shift if bool(np.all(moved < lower + size - 1e-9)) else None


def all_beads_inside(xy: np.ndarray, lower: np.ndarray, size: np.ndarray, cropped: list[bool]) -> bool:
    """The earlier selection rule, kept for the diagnostics: every bead inside the window on every cropped axis."""
    return all((not cropped[a]) or (xy[:, a].min() >= lower[a] and xy[:, a].max() <= lower[a] + size[a]) for a in range(2))


def select_lipids(whole: list[tuple], lower: np.ndarray, size: np.ndarray, cell: np.ndarray, cropped: list[bool],
                  centre: np.ndarray, midplane_nm: float) -> tuple[list, list]:
    """Apply the anchor rule to every (molecule, whole xyz, anchor) record; returns kept (mol, xyz, index) and the decisions."""
    kept, decisions = [], []
    for index, (mol, xyz, anchor) in enumerate(whole):
        shift = anchor_image_shift(anchor[:2], lower, size, cell)
        old_xy = xyz[:, :2] + image_molecule(xyz[:, :2], centre, cell)
        decisions.append({"index": index, "cg": mol["cg"], "leaflet": leaflet_of(float(anchor[2]), midplane_nm),
                          "anchor_inside": shift is not None, "all_beads_inside": all_beads_inside(old_xy, lower, size, cropped)})
        if shift is None:
            continue
        moved = xyz.copy()
        moved[:, :2] = xyz[:, :2] + shift - lower  # the window's lower corner becomes the origin; the residue stays whole
        kept.append((mol, moved, index))
    return kept, decisions


def inherent_pairs(points: np.ndarray, owner: np.ndarray, size: np.ndarray, cell_z: float, cropped: list[bool], radius: float) -> set:
    """Inter-molecule bead pairs within radius under the SOURCE periodicity (not periodic across a cropped axis)."""
    old_box = np.array([1e4 if cropped[0] else size[0], 1e4 if cropped[1] else size[1], cell_z])
    shifted = points - points.min(axis=0)
    pairs = cKDTree(np.mod(shifted, old_box), boxsize=old_box).query_pairs(radius, output_type="ndarray")
    return {(min(int(a), int(b)), max(int(a), int(b))) for a, b in pairs if owner[a] != owner[b]}


def created_pairs(points: np.ndarray, owner: np.ndarray, box: np.ndarray, radius: float, inherent: set) -> np.ndarray:
    """Inter-molecule bead pairs within radius in the periodic cell, excluding those already close in the source frame."""
    pairs = cKDTree(np.mod(points, box), boxsize=box).query_pairs(radius, output_type="ndarray")
    if not len(pairs):
        return np.zeros((0, 2), int)
    a, b = np.minimum(pairs[:, 0], pairs[:, 1]), np.maximum(pairs[:, 0], pairs[:, 1])
    if isinstance(inherent, np.ndarray):  # sorted pair keys a * n + b (pair_keys): vectorized, for large systems
        known = np.isin(a * len(points) + b, inherent, assume_unique=False)
    else:
        known = np.array([(int(i), int(j)) in inherent for i, j in zip(a, b)], dtype=bool)
    keep = (owner[a] != owner[b]) & ~known
    return np.column_stack([a[keep], b[keep]])


def pair_keys(pairs: set, n: int) -> np.ndarray:
    """A set of bead index pairs (i < j) as sorted integer keys i * n + j, the form created_pairs looks up fastest."""
    return np.sort(np.array([i * n + j for i, j in pairs], dtype=np.int64))


def relax_seam(xyz: list[np.ndarray], size: np.ndarray, cell_z: float, cropped: list[bool], protein_nm: np.ndarray,
               contact_nm: float, steps: int = 400, repulsion_nm: float = 0.40, max_step_nm: float = 0.02,
               protein_repulsion_nm: float | None = None, protein_clearance_nm: float | None = None) -> dict:
    """Relax the new periodic seam in place: push apart beads of different lipids that only the cut brought together.

    Steepest descent on: a soft repulsion (force k (r0 - r), potential k (r0 - r)^2 / 2, k = 1) for inter-lipid bead
    pairs closer than repulsion_nm that were NOT that close in the source frame, the same repulsion from the placed
    complex's heavy atoms (which do not move), weak position restraints holding every bead near where it was, and distance
    restraints between all beads of one lipid that keep its shape. It stops when no such pair is closer than contact_nm (the seam contact
    threshold) or
    after `steps` steps. Returns the moved coordinates and what moved.

    The same relaxation makes room for a complex placed into an uncut membrane (embedding.embed_complex): with
    cropped = [False, False] no inter-lipid pair is "created" at the start, the complex repels beads within
    protein_repulsion_nm (default repulsion_nm), and with protein_clearance_nm set the relaxation also continues
    until no bead is closer than that to a heavy atom of the complex.
    """
    owner = np.concatenate([np.full(len(x), k) for k, x in enumerate(xyz)])
    start = np.vstack(xyz).astype(float)
    new_box = np.array([size[0], size[1], cell_z])
    inherent = pair_keys(inherent_pairs(start, owner, size, cell_z, cropped, repulsion_nm), len(start))
    offsets = np.concatenate([[0], np.cumsum([len(x) for x in xyz])])
    intra = np.array([(offsets[k] + i, offsets[k] + j) for k, x in enumerate(xyz) for i in range(len(x)) for j in range(i + 1, len(x))])
    d0 = np.linalg.norm(start[intra[:, 1]] - start[intra[:, 0]], axis=1)
    protein = cKDTree(np.mod(protein_nm, new_box), boxsize=new_box) if len(protein_nm) else None
    protein_range = repulsion_nm if protein_repulsion_nm is None else protein_repulsion_nm
    k_rep, k_intra, k_pos = 1.0, 1.0, 0.05
    minimum_image = lambda d: d - new_box * np.round(d / new_box)
    scatter = lambda index, values: np.stack([np.bincount(index, values[:, k], len(points)) for k in range(3)], axis=1)
    skin = 0.1  # neighbour list: pairs within repulsion_nm + skin, rebuilt once a bead has moved more than skin / 2
    # A lipid at its starting coordinates has no shape or position force, so only lipids that something has pushed enter the
    # intra-lipid restraints; once active a lipid stays active.
    lipid_of, active, active_intra = owner[intra[:, 0]], np.zeros(len(xyz), bool), np.zeros(len(intra), bool)
    points, remaining, step, crowded = start.copy(), np.zeros((0, 2), int), 0, False
    listed, candidates = points.copy(), created_pairs(points, owner, new_box, repulsion_nm + skin, inherent)
    while step < steps:
        step += 1
        if np.linalg.norm(points - listed, axis=1).max() > 0.5 * skin:
            listed, candidates = points.copy(), created_pairs(points, owner, new_box, repulsion_nm + skin, inherent)
        delta = minimum_image(points[candidates[:, 1]] - points[candidates[:, 0]]) if len(candidates) else np.zeros((0, 3))
        dist = np.linalg.norm(delta, axis=1)
        close = dist <= repulsion_nm
        pairs, delta, dist = candidates[close], delta[close], dist[close]
        remaining = pairs[dist < contact_nm]
        near, index = (protein.query(np.mod(points, new_box), distance_upper_bound=protein_range)
                       if protein is not None else (np.full(len(points), np.inf), None))
        crowded = protein_clearance_nm is not None and bool((near < protein_clearance_nm).any())
        if not len(remaining) and not crowded:
            break
        force = np.zeros_like(points)
        unit = delta / np.maximum(dist, 1e-9)[:, None]
        push = (k_rep * (repulsion_nm - dist))[:, None] * unit
        force += scatter(pairs[:, 1], push) - scatter(pairs[:, 0], push)
        if protein is not None:
            hit = np.isfinite(near)
            away = minimum_image(np.mod(points[hit], new_box) - protein.data[index[hit]])
            force[hit] += (k_rep * (protein_range - near[hit]))[:, None] * away / np.maximum(near[hit], 1e-9)[:, None]
        pushed = np.bincount(owner, force.any(axis=1), len(xyz)) > 0
        if (pushed & ~active).any():
            active |= pushed
            active_intra = active[lipid_of]
            first, second, rest = intra[active_intra, 0], intra[active_intra, 1], d0[active_intra]
        if active_intra.any():
            bond = points[second] - points[first]
            length = np.linalg.norm(bond, axis=1)
            pull = (k_intra * (length - rest))[:, None] * bond / np.maximum(length, 1e-9)[:, None]
            force += scatter(first, pull) - scatter(second, pull)
        force -= k_pos * (points - start)
        largest = float(np.linalg.norm(force, axis=1).max())
        points += force * (min(max_step_nm / largest, 0.1) if largest > 0 else 0.0)
    moved = np.linalg.norm(points - start, axis=1)
    per_lipid = np.array([moved[offsets[k]:offsets[k + 1]].max() for k in range(len(xyz))])
    pairs = created_pairs(points, owner, new_box, contact_nm, inherent)
    gaps = np.linalg.norm(minimum_image(points[pairs[:, 1]] - points[pairs[:, 0]]), axis=1) if len(pairs) else np.zeros(0)
    shape = np.abs(np.linalg.norm(points[intra[:, 1]] - points[intra[:, 0]], axis=1) - d0)
    clearance = float(protein.query(np.mod(points, new_box))[0].min()) if protein is not None else None
    return {"xyz": [points[offsets[k]:offsets[k + 1]] for k in range(len(xyz))], "steps": step,
            "converged": bool(len(remaining) == 0 and not crowded), "pairs_created_by_new_periodicity": int(len(pairs)),
            "closest_created_pair_nm": round(float(gaps.min()), 4) if len(gaps) else None,
            "lipids_moved": int((per_lipid > 1e-6).sum()), "max_bead_displacement_nm": round(float(moved.max()), 4),
            "mean_bead_displacement_nm": round(float(moved.mean()), 4), "bead_displacements_nm": moved,
            "max_intramolecular_distance_change_nm": round(float(shape.max()), 4) if len(shape) else 0.0,
            "pairs_already_close_in_source_frame": len(inherent), "remaining_pairs": pairs, "remaining_pair_distances_nm": gaps,
            "owner": owner, "closest_bead_to_complex_nm": round(clearance, 4) if clearance is not None else None}


def remove_hard_core_overlaps(kept: list, relaxed: dict, hard_core_nm: float) -> tuple[list, Counter]:
    """Drop the few lipids still in a hard-core overlap after the seam relaxation (the one with most overlaps first)."""
    pairs = relaxed["remaining_pairs"][relaxed["remaining_pair_distances_nm"] < hard_core_nm] if len(relaxed["remaining_pairs"]) else []
    pairs = [(int(relaxed["owner"][a]), int(relaxed["owner"][b])) for a, b in pairs]
    alive, removed = set(range(len(kept))), Counter()
    while pairs:
        counts = Counter(i for pair in pairs for i in pair)
        worst = max(counts, key=lambda i: (counts[i], i))
        removed[kept[worst][0]["cg"]] += 1
        alive.discard(worst)
        pairs = [p for p in pairs if worst not in p]
    return [kept[i] for i in sorted(alive)], removed


def composition_distance(reference: Counter, sample: Counter) -> float:
    """Half the sum of absolute mole-fraction differences between two compositions (0 = identical, 1 = disjoint)."""
    n_ref, n = sum(reference.values()), sum(sample.values())
    if not n_ref or not n:
        return 1.0
    return 0.5 * sum(abs(reference[k] / n_ref - sample[k] / n) for k in set(reference) | set(sample))


def score_offset(decisions: list[dict], reference: dict, accessible_nm2: dict, hard_overlaps: int) -> dict:
    """Score one crop position from its lipid decisions: density deviation per leaflet + composition + seam penalty."""
    # reference: {"lower": {"apl_nm2": ..., "composition": Counter}, "upper": {...}} measured on the embedded membrane;
    # accessible_nm2: the lipid-accessible area of the new cell per leaflet (window area minus the protein footprint).
    kept = [d for d in decisions if d["anchor_inside"]]
    score, detail = 0.0, {}
    for leaflet in ("lower", "upper"):
        rows = [d for d in kept if d["leaflet"] == leaflet]
        target = accessible_nm2[leaflet] / reference[leaflet]["apl_nm2"] if reference[leaflet]["apl_nm2"] else 0.0
        density = 100.0 * (len(rows) - target) / target if target else 0.0
        distance = composition_distance(reference[leaflet]["composition"], Counter(d["cg"] for d in rows))
        detail[leaflet] = {"lipids": len(rows), "expected_lipids": round(float(target), 1),
                           "density_deviation_percent": round(float(density), 2), "composition_distance": round(float(distance), 4)}
        score += abs(density) + 100.0 * distance
    score += 1.0 * hard_overlaps
    return {"score": round(float(score), 3), "hard_core_overlaps": int(hard_overlaps), **detail}


def offset_slack(settings: Settings) -> float:
    """How far (nm) choose_offset may shift a cropped window on each side: the search range, never past the minimum margin."""
    if not settings.slice_optimize_offset:
        return 0.0
    return min(settings.slice_offset_search_nm, max(settings.box_xy_buffer_nm - settings.box_xy_min_buffer_nm, 0.0))


def choose_offset(whole: list[tuple], windows: dict, cell: np.ndarray, cell_z: float, midplane_nm: float, reference: dict | None,
                  protein_area_nm2: dict, settings: Settings) -> dict:
    """Try shifting the window on the cropped axes by a coarse grid of offsets and keep the best-scoring position."""
    cropped = [w["cropped"] for w in windows["axes"]]
    lower0, size, centre = windows["lower"], windows["size"], windows["centre"]
    area = float(size[0] * size[1])
    accessible = {leaflet: max(area - protein_area_nm2.get(leaflet, 0.0), 1e-6) for leaflet in ("lower", "upper")}
    slack = offset_slack(settings)
    steps = int(round(slack / settings.slice_offset_step_nm)) if settings.slice_offset_step_nm > 0 else 0
    grid = [k * settings.slice_offset_step_nm for k in range(-steps, steps + 1)]
    candidates = sorted({(dx if cropped[0] else 0.0, dy if cropped[1] else 0.0) for dx in grid for dy in grid},
                        key=lambda o: (abs(o[0]) + abs(o[1]), o))  # smallest offset first: a tie keeps the centred window
    if reference is None or not settings.slice_optimize_offset or not any(cropped):
        candidates = [(0.0, 0.0)]
    best, tried = None, []
    for offset in candidates:
        lower = lower0 + np.array(offset)
        kept, decisions = select_lipids(whole, lower, size, cell, cropped, centre + np.array(offset), midplane_nm)
        hard = 0
        if reference is not None and kept and any(cropped):
            xyz = [x for _, x, _ in kept]
            owner = np.concatenate([np.full(len(x), k) for k, x in enumerate(xyz)])
            pts = np.vstack(xyz)
            hard = len(created_pairs(pts, owner, np.array([size[0], size[1], cell_z]), settings.seam_hard_core_nm,
                                     inherent_pairs(pts, owner, size, cell_z, cropped, settings.seam_hard_core_nm)))
        scored = score_offset(decisions, reference, accessible, hard) if reference is not None else {"score": 0.0}
        tried.append({"offset_nm": [round(float(v), 3) for v in offset], **scored})
        if best is None or scored["score"] < best["score"] - 1e-9:
            best = {"offset": np.array(offset), "lower": lower, "kept": kept, "decisions": decisions, "score": scored["score"], "detail": scored}
    return {**best, "tried": tried, "unshifted": tried[0], "accessible_nm2": accessible}


def leaflet_counts(membrane: list[dict], midplane_nm: float) -> dict:
    """Lipid counts per species in the lower and upper leaflet (headgroup anchor below or above the midplane)."""
    lower, upper = Counter(), Counter()
    for mol in membrane:
        xyz = np.array([[b["x"], b["y"], b["z"]] for b in mol["beads"]], dtype=float)
        z = anchor_xyz(xyz, [b["atom"] for b in mol["beads"]], anchor_bead(mol["cg"]))[2]
        (lower if leaflet_of(float(z), midplane_nm) == "lower" else upper)[mol["cg"]] += 1
    return {"lower": dict(sorted(lower.items())), "upper": dict(sorted(upper.items()))}


def selection_diagnostics(decisions: list[dict]) -> dict:
    """Compare the anchor rule with the earlier all-beads-inside rule, by leaflet and species."""
    new = [d for d in decisions if d["anchor_inside"]]
    old = [d for d in decisions if d["all_beads_inside"]]
    only_tails_out = [d for d in new if not d["all_beads_inside"]]
    by = lambda rows, key: dict(sorted(Counter(d[key] for d in rows).items()))
    return {"criterion": "headgroup anchor inside the half-open window; the complete residue is kept",
            "kept_by_anchor_rule": len(new), "kept_by_all_beads_inside_rule": len(old), "difference": len(new) - len(old),
            "rejected_by_old_rule_with_anchor_inside": len(only_tails_out),
            "rejected_by_old_rule_with_anchor_inside_by_species": by(only_tails_out, "cg"),
            "rejected_by_old_rule_with_anchor_inside_by_leaflet": by(only_tails_out, "leaflet"),
            "kept_by_leaflet": {"anchor_rule": by(new, "leaflet"), "all_beads_inside_rule": by(old, "leaflet")},
            "kept_by_species": {"anchor_rule": by(new, "cg"), "all_beads_inside_rule": by(old, "cg")},
            "note": "the old rule deleted every lipid with one bead across the window edge but kept the full window "
                    "area, which under-packed the cut membrane; the anchor rule keeps the lipid population of the source"}


def periodic_image_gap(points_xy: np.ndarray, cell: np.ndarray) -> float:
    """Closest distance (nm) between a set of points and its own periodic images in x and y."""
    tree, nearest = cKDTree(points_xy[:, :2]), np.inf
    for sx, sy in ((1, 0), (0, 1), (1, 1), (1, -1)):
        nearest = min(nearest, float(tree.query(points_xy[:, :2] + np.array([sx, sy]) * cell)[0].min()))
    return nearest


def voronoi_shells(protein_xy: np.ndarray, lipid_xy: np.ndarray, cell: np.ndarray) -> dict:
    """Lipid shells around the protein of one leaflet by Voronoi adjacency, periodic in x and y, with no shell width."""
    # The protein atoms and the lipid anchors of the leaflet are generators of a 2D Voronoi tessellation of the cell tiled
    # 3 x 3. First-shell lipids have a cell that shares an edge with a protein atom's cell; shell k+1 are the lipids
    # adjacent to shell k not yet reached, counted from the central protein through the tiling. The number of lipid rows
    # between the protein and its own image is the smallest shell number of a lipid (in any tile) that is adjacent to an
    # image's protein atom; each protein owns half of those rows, so shells_per_side is that number halved and rounded down. This
    # follows the idea of defining shells by molecular adjacency instead of a fixed width; it is not a re-implementation of
    # any published code (those use further criteria, e.g. a boundary-proximity cut-off, that are not applied here).
    n_lipid, n_protein = len(lipid_xy), len(protein_xy)
    if n_lipid < 3 or n_protein == 0:
        return {"lipids_per_shell": [0, 0, 0], "lipid_rows_between_images": None, "shells_per_side": 0}
    wrapped = np.vstack([lipid_xy, protein_xy]) % cell
    shifts = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
    tiled = np.vstack([wrapped + cell * s for s in shifts])
    total, central = len(wrapped), shifts.index((0, 0))

    def block(g: int) -> int:  # which tile a generator belongs to
        return g // total

    def is_protein(g: int) -> bool:  # local index g % total is in `wrapped`
        return g % total >= n_lipid

    neighbours = [[] for _ in range(len(tiled))]
    for a, b in Voronoi(tiled).ridge_points:
        neighbours[a].append(b)
        neighbours[b].append(a)
    depth, frontier = {}, []
    for g in range(central * total + n_lipid, (central + 1) * total):  # protein atoms of the central tile
        for h in neighbours[g]:
            if not is_protein(h) and h not in depth:
                depth[h] = 1
                frontier.append(h)
    while frontier:
        following = []
        for g in frontier:
            for h in neighbours[g]:
                if not is_protein(h) and h not in depth:
                    depth[h] = depth[g] + 1
                    following.append(h)
        frontier = following
    in_central = [g for g in depth if block(g) == central]
    per_shell = [sum(1 for g in in_central if depth[g] == k) for k in (1, 2, 3)]
    rows = [depth[g] for g in depth if any(is_protein(h) and block(h) != central for h in neighbours[g])]
    between = min(rows) if rows else None
    return {"lipids_per_shell": per_shell, "lipid_rows_between_images": between,
            "shells_per_side": 0 if between is None else between // 2}


def xy_periodic_check(complex_a: list[dict], sliced: list[dict], cell: np.ndarray, midplane_nm: float, settings: Settings,
                      cropped: tuple = ()) -> dict:
    """Lipid shells around the transmembrane region and the complex's distance to its own periodic images in x and y."""
    # Reported, not enforced: a violation is listed in `flags` and the caller logs it; no dimension is changed here.
    # `cropped` says per axis whether the membrane was cut; a flag on an axis that was not cut means the source cell
    # itself is too small, which a larger buffer cannot fix, and the message says so.
    # estimated_shell_equivalents_per_side is a geometric proxy: half the closest approach of the TMD to its image through
    # the membrane (each of the two proteins owns half of the lipids between them) in units of lipid_shell_nm. The
    # shells counted from molecular adjacency (voronoi_shells, per leaflet, from the lipid coordinates) are the measure
    # that supports a statement that two lipid shells surround the receptor.
    heavy = np.array([[a["x"], a["y"], a["z"]] for a in complex_a if element(a["atom"]) != "H"], dtype=float) / 10.0
    tmd = heavy[np.abs(heavy[:, 2] - midplane_nm) <= settings.tmd_half_nm]
    if not len(tmd):
        return {"checked": False, "reason": "no complex heavy atom within tmd_half_nm of the bilayer midplane", "flags": []}
    tmd_gap, complex_gap = periodic_image_gap(tmd, cell), periodic_image_gap(heavy, cell)
    anchors = np.array([anchor_xyz(m["xyz"], [b["atom"] for b in m["beads"]], anchor_bead(m["cg"])) for m in sliced])
    distance = cKDTree(tmd[:, :2] % cell, boxsize=cell).query(anchors[:, :2] % cell)[0]
    width = settings.lipid_shell_nm
    counts = [int(((distance >= i * width) & (distance < (i + 1) * width)).sum()) for i in range(3)]
    per_side = tmd_gap / 2.0 / width
    upper = anchors[:, 2] >= midplane_nm
    voronoi = {}
    for name, mask, inside in (("upper", upper, (heavy[:, 2] >= midplane_nm) & (heavy[:, 2] <= midplane_nm + settings.tmd_half_nm)),
                               ("lower", ~upper, (heavy[:, 2] < midplane_nm) & (heavy[:, 2] >= midplane_nm - settings.tmd_half_nm))):
        voronoi[name] = voronoi_shells(heavy[inside][:, :2], anchors[mask][:, :2], cell)
    least = min(v["shells_per_side"] for v in voronoi.values())
    flags = []
    remedy = ("the source membrane cell is not cut on every axis and is too small for the complex: use a larger "
              "coarse-grained membrane (the bundled frames are 18 nm wide)" if not all(cropped) and len(cropped) else "raise --xy-buffer")
    if per_side < settings.min_lipid_shells:
        flags.append(f"an estimated {per_side:.1f} lipid shell equivalents ({0.5 * tmd_gap:.2f} nm) per side between the transmembrane "
                     f"region and its image; at least {settings.min_lipid_shells:g} recommended (Feng et al. 2023): {remedy}")
    if least < settings.min_lipid_shells:
        flags.append(f"Voronoi adjacency gives {voronoi['upper']['shells_per_side']} (upper) and {voronoi['lower']['shells_per_side']} "
                     f"(lower) lipid shells per side between the transmembrane region and its image; at least "
                     f"{settings.min_lipid_shells:g} recommended: {remedy}")
    if complex_gap < settings.min_image_clearance_nm:
        flags.append(f"the complex is {complex_gap:.2f} nm from its own periodic image in x/y, below the "
                     f"{settings.min_image_clearance_nm:g} nm cutoff: {remedy}")
    return {"checked": True, "tmd_atoms": int(len(tmd)), "lipid_shell_nm": width, "tmd_image_gap_nm": round(tmd_gap, 3),
            "shells_between_images": round(tmd_gap / width, 2), "estimated_shell_equivalents_per_side": round(per_side, 2),
            "lipids_per_shell": counts, "voronoi_shells": voronoi, "voronoi_shells_per_side": least,
            "complex_image_clearance_nm": round(complex_gap, 3),
            "min_lipid_shells": settings.min_lipid_shells, "min_image_clearance_nm": settings.min_image_clearance_nm, "flags": flags}


def slice_membrane_cg(membrane: list[dict], placed: list[dict], box: list[float], midplane_nm: float, settings: Settings,
                      requested: tuple | None = None, reference: dict | None = None) -> dict:
    """Crop the CG membrane to the placed complex plus the buffer; returns the new membrane, complex, box and a report.

    `reference` (optional) is the packing of the membrane being cut, {"lower"/"upper": {"apl_nm2", "composition",
    "protein_area_nm2"}} as packing.reference_packing gives it; with it the crop position is optimized.
    """
    cell = np.array(box[:2], dtype=float)
    cell3 = np.array(box, dtype=float)
    low, high = solute_extent(placed)
    windows = crop_windows(low, high, cell, settings.box_xy_buffer_nm, requested, settings.box_xy_min_buffer_nm)
    size = windows["size"]
    cropped = [w["cropped"] for w in windows["axes"]]
    before = Counter(m["cg"] for m in membrane)
    whole = [(mol, *lipid_anchor(mol, cell3)) for mol in membrane]
    protein_area = {leaflet: reference[leaflet].get("protein_area_nm2", 0.0) for leaflet in ("lower", "upper")} if reference else {}
    chosen = choose_offset(whole, windows, cell, float(box[2]), midplane_nm, reference, protein_area, settings)
    lower, kept, decisions = chosen["lower"], chosen["kept"], chosen["decisions"]
    outside = sum(not d["anchor_inside"] for d in decisions)
    moved_complex = [dict(a, x=a["x"] - 10.0 * lower[0], y=a["y"] - 10.0 * lower[1]) for a in placed]
    removed_seam, relaxed = Counter(), None
    if any(cropped) and kept:
        protein_nm = xyz_nm([a for a in moved_complex if element(a["atom"]) != "H"]) / 10.0
        relaxed = relax_seam([x for _, x, _ in kept], size, float(box[2]), cropped, protein_nm, settings.seam_min_bead_nm,
                             settings.seam_relax_steps)
        kept = [(mol, xyz, index) for (mol, _, index), xyz in zip(kept, relaxed["xyz"])]
        kept, removed_seam = remove_hard_core_overlaps(kept, relaxed, settings.seam_hard_core_nm)
    sliced = []
    for mol, xyz, _ in kept:
        beads = [dict(b, x=float(p[0]), y=float(p[1]), z=float(p[2])) for b, p in zip(mol["beads"], xyz)]
        sliced.append(dict(mol, beads=beads, xyz=xyz.copy()))
    after = Counter(m["cg"] for m in sliced)
    leaflets = leaflet_counts(sliced, midplane_nm)
    phosphate = [m for m in sliced if any(b["atom"] in ("PO4", "P") for b in m["beads"])]
    if (not leaflets["lower"] or not leaflets["upper"] or len(sliced) < settings.slice_min_lipids
            or len(phosphate) < settings.slice_min_lipids):
        raise SystemExit(f"slicing the membrane to {size[0]:.2f} x {size[1]:.2f} nm leaves {len(sliced)} lipids "
                         f"({len(phosphate)} phospholipids; lower leaflet {sum(leaflets['lower'].values())}, upper "
                         f"{sum(leaflets['upper'].values())}); at least {settings.slice_min_lipids} are needed in a two-leaflet "
                         "bilayer. Use a larger buffer (--xy-buffer) or a larger coarse-grained cell")
    new_box = [float(size[0]), float(size[1]), float(box[2])]
    margins = {"low_x_nm": round(float(low[0] - lower[0]), 3), "high_x_nm": round(float(lower[0] + size[0] - high[0]), 3),
               "low_y_nm": round(float(low[1] - lower[1]), 3), "high_y_nm": round(float(lower[1] + size[1] - high[1]), 3)}
    seam = {"threshold_nm": settings.seam_min_bead_nm, "hard_core_nm": settings.seam_hard_core_nm,
            "pairs_created_by_new_periodicity": 0, "pairs_already_close_in_source_frame": 0,
            "relaxation": None, "lipids_removed": sum(removed_seam.values()), "removed_by_species": dict(sorted(removed_seam.items())),
            "note": "seam contacts are relaxed in place at the coarse-grained level; only hard-core overlaps that survive the "
                    "relaxation cost a lipid; pairs already close in the source frame are data and are left alone"}
    if relaxed is not None:
        seam.update(pairs_created_by_new_periodicity=relaxed["pairs_created_by_new_periodicity"],
                    pairs_already_close_in_source_frame=relaxed["pairs_already_close_in_source_frame"],
                    relaxation={k: relaxed[k] for k in ("steps", "converged", "closest_created_pair_nm", "lipids_moved",
                                                        "max_bead_displacement_nm", "mean_bead_displacement_nm",
                                                        "max_intramolecular_distance_change_nm")})
    report = {"mode": "user" if requested is not None else "auto",
              "algorithm": "complex extent + buffer in x and y; each lipid made whole under the source periodicity and kept, "
                           "complete, when the periodic image of its headgroup anchor (lipids.LIPID_SLICE_ANCHORS) falls in "
                           "the half-open window; crop position chosen by lipid density and composition; seam relaxed in place",
              "buffer_nm": settings.box_xy_buffer_nm, "requested_xy_nm": list(requested[:2]) if requested is not None else None,
              "source_cell_nm": [float(v) for v in box], "new_cell_nm": new_box, "axes": windows["axes"],
              "cropped": any(cropped), "area_fraction_kept": round(float(size[0] * size[1] / (cell[0] * cell[1])), 4),
              "solute_extent_nm": [[round(float(v), 3) for v in low], [round(float(v), 3) for v in high]],
              "membrane_margin_around_complex_nm": margins,
              "crop_offset": {"offset_nm": [round(float(v), 3) for v in chosen["offset"]],
                              "optimized": bool(reference is not None and settings.slice_optimize_offset and any(cropped)),
                              "candidates": len(chosen["tried"]), "chosen": chosen["detail"], "unshifted": chosen["unshifted"],
                              "accessible_area_nm2": {k: round(float(v), 3) for k, v in chosen["accessible_nm2"].items()}},
              "lipids_before": sum(before.values()), "lipids_after": len(sliced),
              "lipids_outside_window": outside, "composition_before": dict(sorted(before.items())),
              "composition_after": dict(sorted(after.items())), "leaflets_before": leaflet_counts(membrane, midplane_nm),
              "leaflets_after": leaflets, "bilayer_midplane_nm": round(midplane_nm, 4),
              "selection": selection_diagnostics(decisions), "seam": seam,
              "xy_check": xy_periodic_check(moved_complex, sliced, np.array(new_box[:2]), midplane_nm, settings, tuple(cropped))}
    return {"membrane": sliced, "placed": moved_complex, "box": new_box, "composition": after, "report": report,
            "kept_indices": [index for _, _, index in kept]}
