#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// BOX = auto: slice the coarse-grained membrane to the complex plus a buffer, before backmapping.
#//=============================================================
"""Crop the coarse-grained membrane in x and y around the placed complex, keeping whole molecules only.

This is the slicing algorithm of the earlier MembraneForger workflow (archive/scripts/aa_stage4.py, slice_membrane,
buffer 1.0 nm), moved in front of the expensive step. The earlier version cropped the already backmapped all-atom
system; here the coarse-grained molecules are cropped, so mstool only backmaps the lipids that are kept.

    window (per axis) = [min, max] of the solute + buffer on each side
    a lipid is kept if ALL of its beads lie inside the window (its periodic image nearest the window centre)
    an axis whose window is as wide as the periodic cell is not cropped (the original periodicity is kept)
    new cell = window width (cropped axes) or the original cell width; z is left to the box stage

Cropping makes opposite window edges periodic neighbours that never were in the simulation. Bead pairs of different
lipids that are closer than the seam threshold ONLY through that new periodicity are seam clashes: they are
removed (the lipid with most seam clashes first) and reported. Contacts that were already close in the source frame
belong to the data, are not touched, and are counted separately.

Orientation determines the reference frame and slicing operates in it; nothing here knows about PPM or OPM.
"""
from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

from .config import Settings

__all__ = ['solute_extent', 'crop_windows', 'image_molecule', 'seam_clashes', 'leaflet_counts', 'slice_membrane_cg']


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


def seam_clashes(xyz: list[np.ndarray], size: np.ndarray, cell_z: float, threshold_nm: float,
                 cropped: list[bool]) -> tuple[list, int]:
    """Bead pairs of different molecules closer than threshold only through the new periodic cell, plus the old count."""
    # "Already close" is judged with the periodicity of the SOURCE frame: periodic in z and in every uncropped axis,
    # not periodic across a cropped axis. A pair that is close in the new cell but not in that source arrangement was
    # created by the new periodic boundary.
    owner = np.concatenate([np.full(len(x), k) for k, x in enumerate(xyz)])
    points = np.vstack(xyz)
    points = points - points.min(axis=0)
    new_box = np.array([size[0], size[1], cell_z])
    old_box = np.array([1e4 if cropped[0] else size[0], 1e4 if cropped[1] else size[1], cell_z])
    periodic = cKDTree(np.mod(points, new_box), boxsize=new_box).query_pairs(threshold_nm, output_type="ndarray")
    source = cKDTree(np.mod(points, old_box), boxsize=old_box).query_pairs(threshold_nm, output_type="ndarray")
    inherent = {(min(int(a), int(b)), max(int(a), int(b))) for a, b in source if owner[a] != owner[b]}
    created = [(int(owner[a]), int(owner[b])) for a, b in periodic
               if owner[a] != owner[b] and (min(int(a), int(b)), max(int(a), int(b))) not in inherent]
    return created, len(inherent)


def leaflet_counts(membrane: list[dict], midplane_nm: float) -> dict:
    """Lipid counts per species in the lower and upper leaflet (centroid below or above the midplane)."""
    lower, upper = Counter(), Counter()
    for mol in membrane:
        z = float(np.mean([b["z"] for b in mol["beads"]]))
        (lower if z < midplane_nm else upper)[mol["cg"]] += 1
    return {"lower": dict(sorted(lower.items())), "upper": dict(sorted(upper.items()))}


def slice_membrane_cg(membrane: list[dict], placed: list[dict], box: list[float], midplane_nm: float, settings: Settings,
                      requested: tuple | None = None) -> dict:
    """Crop the CG membrane to the placed complex plus the buffer; returns the new membrane, complex, box and a report."""
    cell = np.array(box[:2], dtype=float)
    low, high = solute_extent(placed)
    windows = crop_windows(low, high, cell, settings.box_xy_buffer_nm, requested, settings.box_xy_min_buffer_nm)
    lower, size, centre = windows["lower"], windows["size"], windows["centre"]
    cropped = [w["cropped"] for w in windows["axes"]]
    before = Counter(m["cg"] for m in membrane)
    kept, outside = [], 0
    for mol in membrane:
        xyz = np.array([[b["x"], b["y"], b["z"]] for b in mol["beads"]], dtype=float)
        shift = image_molecule(xyz[:, :2], centre, cell)
        xy = xyz[:, :2] + shift
        inside = all((not cropped[a]) or (xy[:, a].min() >= lower[a] and xy[:, a].max() <= lower[a] + size[a]) for a in range(2))
        if not inside:
            outside += 1
            continue
        moved = xyz.copy()
        moved[:, :2] = xy - lower  # the window's lower corner becomes the origin
        kept.append((mol, moved))
    if not any(cropped):
        kept = [(mol, np.array([[b["x"], b["y"], b["z"]] for b in mol["beads"]], dtype=float) + np.array([-lower[0], -lower[1], 0.0]))
                for mol in membrane]
        # nothing is cropped: the molecules keep their own images; only the origin moves with the complex
    removed_seam = Counter()
    created, inherent = [], 0
    if any(cropped) and kept:
        created, inherent = seam_clashes([x for _, x in kept], size, float(box[2]), settings.seam_min_bead_nm, cropped)
        pairs = list(created)
        alive = set(range(len(kept)))
        while pairs:
            counts = Counter(i for pair in pairs for i in pair)
            worst = max(counts, key=lambda i: (counts[i], i))
            removed_seam[kept[worst][0]["cg"]] += 1
            alive.discard(worst)
            pairs = [p for p in pairs if worst not in p]
        kept = [kept[i] for i in sorted(alive)]
    sliced = []
    for mol, xyz in kept:
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
    moved_complex = [dict(a, x=a["x"] - 10.0 * lower[0], y=a["y"] - 10.0 * lower[1]) for a in placed]
    new_box = [float(size[0]), float(size[1]), float(box[2])]
    margins = {"low_x_nm": round(float(low[0] - lower[0]), 3), "high_x_nm": round(float(lower[0] + size[0] - high[0]), 3),
               "low_y_nm": round(float(low[1] - lower[1]), 3), "high_y_nm": round(float(lower[1] + size[1] - high[1]), 3)}
    report = {"mode": "user" if requested is not None else "auto", "algorithm": "complex extent + buffer in x and y, whole molecules only "
              "(archive/scripts/aa_stage4.py slice_membrane), applied to the coarse-grained membrane before backmapping",
              "buffer_nm": settings.box_xy_buffer_nm, "requested_xy_nm": list(requested[:2]) if requested is not None else None,
              "source_cell_nm": [float(v) for v in box], "new_cell_nm": new_box, "axes": windows["axes"],
              "cropped": any(cropped), "area_fraction_kept": round(float(size[0] * size[1] / (cell[0] * cell[1])), 4),
              "solute_extent_nm": [[round(float(v), 3) for v in low], [round(float(v), 3) for v in high]],
              "membrane_margin_around_complex_nm": margins,
              "lipids_before": sum(before.values()), "lipids_after": len(sliced),
              "lipids_outside_window": outside, "composition_before": dict(sorted(before.items())),
              "composition_after": dict(sorted(after.items())), "leaflets_after": leaflets,
              "bilayer_midplane_nm": round(midplane_nm, 4),
              "seam": {"threshold_nm": settings.seam_min_bead_nm, "pairs_created_by_new_periodicity": len(created),
                       "lipids_removed": sum(removed_seam.values()), "removed_by_species": dict(sorted(removed_seam.items())),
                       "pairs_already_close_in_source_frame": inherent,
                       "note": "pairs already close in the source frame are data, not seam artefacts, and are left alone"}}
    return {"membrane": sliced, "placed": moved_complex, "box": new_box, "composition": after, "report": report}
