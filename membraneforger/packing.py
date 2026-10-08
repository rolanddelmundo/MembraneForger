#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Leaflet packing: periodic Voronoi area per lipid, global and per-species APL, composition.
#//=============================================================
"""Measure how a membrane is packed, leaflet by leaflet, at either resolution, and compare two such measurements.

Two quantities that answer different questions are kept apart throughout:

    global leaflet APL      APL_global = (A_xy - A_protein) / N_lipids      overall packing of a leaflet
    local, per-lipid area   A_i = area of lipid i's cell in a 2D periodic Voronoi tessellation of its leaflet
    species APL             APL_species = mean(A_i for lipids i of that species)

The tessellation uses one headgroup anchor per lipid (lipids.LIPID_SLICE_ANCHORS / LIPID_AA_ANCHORS) and the heavy
atoms of the placed complex that lie in the leaflet's z range as additional generators; the protein cells are
discarded, so the protein footprint is never assigned to a lipid and  sum(A_i) = A_xy - A_protein  exactly. The
global leaflet APL is therefore the mean of the per-lipid areas, and the protein footprint is A_xy - sum(A_i).
"Available area divided by the number of lipids of species X" is NOT an APL of X in a mixed bilayer and is not
computed. A second tessellation of the phospholipid headgroups alone (sterols left out: their anchor sits below
the phosphate plane) gives the area per lipid head in the headgroup plane.

Every function works in nm and reports areas in both nm^2 and A^2. A "stage" is one membrane at one point of the
build: {"name", "resolution", "box_nm", "midplane_nm", "lipids": [per-lipid records], "protein_nm": heavy atoms}.
"""
import math
from collections import Counter

import numpy as np
from scipy.spatial import Voronoi

from .lipids import anchor_bead, anchor_xyz, cg_name, is_sterol, leaflet_of
from .structio import element, xyz_nm

__all__ = ['LEAFLETS', 'cg_stage', 'aa_stage', 'lipid_midplane', 'periodic_voronoi_areas', 'leaflet_slab',
           'tessellate_leaflet', 'measure_packing', 'reference_packing', 'species_table', 'composition_table',
           'expected_composition_distance', 'compare_packing', 'per_lipid_rows', 'region_reference', 'window_distribution']

LEAFLETS = ("lower", "upper")
NM2_TO_A2 = 100.0


def lipid_midplane(z: np.ndarray) -> float:
    """Bilayer midplane from anchor z values: halfway between the medians of the two halves of their extent."""
    middle = 0.5 * float(z.min() + z.max())
    lower, upper = z[z < middle], z[z >= middle]
    if not len(lower) or not len(upper):
        return middle
    return 0.5 * (float(np.median(lower)) + float(np.median(upper)))


def cg_stage(name: str, membrane: list[dict], placed: list[dict], box: list[float], midplane_nm: float) -> dict:
    """One analysis stage from coarse-grained membrane molecules (whole, nm) and the placed complex (A)."""
    lipids = []
    for k, mol in enumerate(membrane):
        xyz = np.array([[b["x"], b["y"], b["z"]] for b in mol["beads"]], dtype=float)
        anchor = anchor_xyz(xyz, [b["atom"] for b in mol["beads"]], anchor_bead(mol["cg"]))
        lipids.append({"index": k, "resid": k + 1, "resname": mol["cg"], "leaflet": leaflet_of(float(anchor[2]), midplane_nm),
                       "x_nm": float(anchor[0]), "y_nm": float(anchor[1]), "z_nm": float(anchor[2]),
                       "z_range_nm": (float(xyz[:, 2].min()), float(xyz[:, 2].max()))})
    protein = xyz_nm([a for a in placed if element(a["atom"]) != "H"]) / 10.0 if placed else np.zeros((0, 3))
    return {"name": name, "resolution": "cg", "box_nm": [float(v) for v in box], "midplane_nm": float(midplane_nm),
            "lipids": lipids, "protein_nm": protein}


def aa_stage(name: str, residues: list[list[dict]], protein_nm: np.ndarray, box_nm: list[float], scale: float = 1.0) -> dict:
    """One analysis stage from all-atom lipid residues (atom dicts with x, y, z in units of 1/scale nm) and the protein."""
    # scale = 0.1 for PDB coordinates in A; protein_nm holds the heavy atoms of the complex (nm).
    rows = []
    for k, atoms in enumerate(residues):
        heavy = [a for a in atoms if element(a["atom"]) != "H"]
        xyz = xyz_nm(heavy) * scale
        anchor = anchor_xyz(xyz, [a["atom"] for a in heavy], anchor_bead(atoms[0]["resname"], atomistic=True))
        rows.append({"index": k, "resid": int(atoms[0].get("resid", k + 1)), "resname": cg_name(atoms[0]["resname"]),
                     "x_nm": float(anchor[0]), "y_nm": float(anchor[1]), "z_nm": float(anchor[2]),
                     "z_range_nm": (float(xyz[:, 2].min()), float(xyz[:, 2].max()))})
    midplane = lipid_midplane(np.array([r["z_nm"] for r in rows])) if rows else 0.0
    for r in rows:
        r["leaflet"] = leaflet_of(r["z_nm"], midplane)
    return {"name": name, "resolution": "aa", "box_nm": [float(v) for v in box_nm], "midplane_nm": float(midplane),
            "lipids": rows, "protein_nm": np.asarray(protein_nm, dtype=float).reshape(-1, 3)}


def periodic_voronoi_areas(points: np.ndarray, box_xy: np.ndarray, extra: np.ndarray | None = None) -> np.ndarray:
    """Area of the periodic 2D Voronoi cell of every point, with `extra` generators competing for area but not reported."""
    # The cell is tiled 3 x 3 so that every central-image generator has all its neighbours; the cell of a generator
    # is the convex polygon of its region's vertices, whose area is the shoelace sum after sorting them by angle.
    generators = np.vstack([points] + ([extra] if extra is not None and len(extra) else []))
    wrapped = np.mod(generators, box_xy)
    tiled = np.vstack([wrapped + box_xy * (i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)])
    central = 4 * len(wrapped)  # the (0, 0) image block
    vor = Voronoi(tiled)
    areas = np.zeros(len(points))
    for k in range(len(points)):
        region = vor.regions[vor.point_region[central + k]]
        if -1 in region or len(region) < 3:
            raise SystemExit("periodic Voronoi tessellation failed: an open cell in the central image (box too small for its lipids?)")
        v = vor.vertices[region] - wrapped[k]
        v = v[np.argsort(np.arctan2(v[:, 1], v[:, 0]))]
        areas[k] = 0.5 * abs(float(np.dot(v[:, 0], np.roll(v[:, 1], -1)) - np.dot(v[:, 1], np.roll(v[:, 0], -1))))
    return areas


HEAD_PAD_NM = 0.2  # the headgroup plane extends this far beyond the farthest anchor (NC3/NH3 beads, choline atoms)


def leaflet_slab(stage: dict, leaflet: str) -> tuple[float, float]:
    """z range (nm) in which protein atoms compete with the lipids of one leaflet: midplane to the headgroup plane."""
    # Bounded by the farthest ANCHOR, not the farthest bead: a single tail bead dipping below the heads would pull
    # in protein atoms that lie in the water below the membrane (a G protein, say) and let them steal area from
    # lipids whose heads sit above them.
    rows = [r for r in stage["lipids"] if r["leaflet"] == leaflet]
    mid = stage["midplane_nm"]
    if leaflet == "lower":
        return (min(r["z_nm"] for r in rows) - HEAD_PAD_NM, mid)
    return (mid, max(r["z_nm"] for r in rows) + HEAD_PAD_NM)


def tessellate_leaflet(stage: dict, leaflet: str) -> dict:
    """Per-lipid Voronoi areas of one leaflet (all lipids, and phospholipid heads alone), with the protein as competitor."""
    rows = [r for r in stage["lipids"] if r["leaflet"] == leaflet]
    if len(rows) < 3:
        raise SystemExit(f"the {leaflet} leaflet has {len(rows)} lipids; at least 3 are needed for a tessellation")
    box_xy = np.array(stage["box_nm"][:2], dtype=float)
    lo, hi = leaflet_slab(stage, leaflet)
    protein = stage["protein_nm"]
    inside = protein[(protein[:, 2] >= lo) & (protein[:, 2] <= hi)][:, :2] if len(protein) else np.zeros((0, 2))
    points = np.array([[r["x_nm"], r["y_nm"]] for r in rows])
    areas = periodic_voronoi_areas(points, box_xy, inside)
    heads = [i for i, r in enumerate(rows) if not is_sterol(r["resname"])]
    head_areas = periodic_voronoi_areas(points[heads], box_xy, inside) if len(heads) >= 3 else np.full(len(heads), np.nan)
    head_area = dict(zip(heads, head_areas))
    for i, r in enumerate(rows):
        r["voronoi_area_nm2"] = float(areas[i])
        r["head_area_nm2"] = float(head_area[i]) if i in head_area else None
    area = float(box_xy[0] * box_xy[1])
    return {"lipids": len(rows), "area_nm2": area, "protein_atoms_in_slab": int(len(inside)),
            "protein_area_nm2": area - float(areas.sum()), "accessible_area_nm2": float(areas.sum()),
            "apl_nm2": float(areas.mean()), "apl_A2": float(areas.mean()) * NM2_TO_A2,
            "head_apl_A2": float(np.nanmean(head_areas)) * NM2_TO_A2 if len(heads) >= 3 else None, "heads": len(heads),
            "slab_nm": [round(lo, 4), round(hi, 4)], "composition": dict(sorted(Counter(r["resname"] for r in rows).items()))}


def measure_packing(stage: dict) -> dict:
    """Tessellate both leaflets of a stage; fills the per-lipid areas in place and returns the leaflet summaries."""
    return {leaflet: tessellate_leaflet(stage, leaflet) for leaflet in LEAFLETS}


def reference_packing(stage: dict) -> dict:
    """The compact reference a slicer needs from an embedded membrane: APL, composition and protein area per leaflet."""
    measured = measure_packing(stage)
    return {leaflet: {"apl_nm2": m["apl_nm2"], "composition": Counter(m["composition"]), "protein_area_nm2": m["protein_area_nm2"]}
            for leaflet, m in measured.items()}


def quartiles(values: np.ndarray) -> dict:
    """Mean, median, standard deviation and interquartile range of a sample (A^2 when given nm^2 x 100)."""
    if not len(values):
        return {"mean": None, "median": None, "sd": None, "iqr": None}
    q1, q3 = np.percentile(values, [25, 75])
    return {"mean": round(float(values.mean()), 2), "median": round(float(np.median(values)), 2),
            "sd": round(float(values.std(ddof=1)), 2) if len(values) > 1 else 0.0, "iqr": round(float(q3 - q1), 2)}


def species_table(stage: dict, reference: dict | None = None) -> list[dict]:
    """Per-leaflet, per-species local APL statistics (A^2) from the per-lipid Voronoi areas, with N shown."""
    rows = []
    for leaflet in LEAFLETS:
        lipids = [r for r in stage["lipids"] if r["leaflet"] == leaflet and "voronoi_area_nm2" in r]
        total = len(lipids)
        for species in sorted({r["resname"] for r in lipids}):
            areas = np.array([r["voronoi_area_nm2"] for r in lipids if r["resname"] == species]) * NM2_TO_A2
            heads = np.array([r["head_area_nm2"] for r in lipids if r["resname"] == species and r.get("head_area_nm2") is not None]) * NM2_TO_A2
            row = {"stage": stage["name"], "leaflet": leaflet, "lipid": species, "n": len(areas), "mole_fraction": round(len(areas) / total, 4),
                   **{f"apl_{k}_A2": v for k, v in quartiles(areas).items()},
                   "head_apl_mean_A2": round(float(heads.mean()), 2) if len(heads) else None,
                   "well_sampled": len(areas) >= 10}
            if reference is not None:
                ref = next((q for q in reference if q["leaflet"] == leaflet and q["lipid"] == species), None)
                row["apl_delta_vs_reference_A2"] = round(row["apl_mean_A2"] - ref["apl_mean_A2"], 2) if ref else None
                change = 100.0 * (row["apl_mean_A2"] - ref["apl_mean_A2"]) / ref["apl_mean_A2"] if ref and ref["apl_mean_A2"] else None
                row["apl_delta_vs_reference_percent"] = round(change, 2) if change is not None else None
                row["mole_fraction_delta_points"] = round(100.0 * (row["mole_fraction"] - ref["mole_fraction"]), 2) if ref else None
            rows.append(row)
    return rows


def expected_composition_distance(reference: Counter, n: int) -> float:
    """Expected composition distance of a random sample of n lipids drawn from the reference mole fractions."""
    # Multinomial sampling: E|x_s - x| ~ sqrt(2/pi) sqrt(x (1 - x) / n) per species; the distance is half their sum.
    total = sum(reference.values())
    if not total or not n:
        return 0.0
    return 0.5 * sum(math.sqrt(2.0 / math.pi) * math.sqrt(x * (1.0 - x) / n) for x in (c / total for c in reference.values()))


def composition_table(reference: dict, sample: dict) -> dict:
    """Per-leaflet species counts and mole fractions of the reference and the sample, with the distance between them."""
    out = {}
    for leaflet in LEAFLETS:
        ref, new = Counter(reference[leaflet]["composition"]), Counter(sample[leaflet]["composition"])
        n_ref, n_new = sum(ref.values()), sum(new.values())
        rows = []
        for species in sorted(set(ref) | set(new)):
            x_ref, x_new = (ref[species] / n_ref if n_ref else 0.0), (new[species] / n_new if n_new else 0.0)
            sigma = math.sqrt(x_ref * (1 - x_ref) / n_new) if n_new else 0.0  # sampling noise of a random crop of n_new lipids
            rows.append({"lipid": species, "n_reference": ref[species], "n_sample": new[species],
                         "x_reference": round(x_ref, 4), "x_sample": round(x_new, 4), "delta_points": round(100.0 * (x_new - x_ref), 2),
                         "sampling_sd_points": round(100.0 * sigma, 2), "flag": abs(x_new - x_ref) > 3.0 * sigma + 1e-9})
        distance = 0.5 * sum(abs(r["x_sample"] - r["x_reference"]) for r in rows)
        expected = expected_composition_distance(ref, n_new)
        out[leaflet] = {"n_reference": n_ref, "n_sample": n_new, "species": rows, "distance": round(distance, 4),
                        "expected_distance_random_crop": round(expected, 4),
                        "within_finite_crop_variability": distance <= 3.0 * expected + 1e-9,
                        "flagged_species": [r["lipid"] for r in rows if r["flag"]]}
    return out


def compare_packing(reference: dict, sample: dict, tolerance_percent: float, reference_measure: dict | None = None,
                    sample_measure: dict | None = None) -> dict:
    """Global leaflet APL of a sample against its reference: absolute and relative change, pass/fail per leaflet."""
    ref = reference_measure or measure_packing(reference)
    new = sample_measure or measure_packing(sample)
    leaflets = {}
    for leaflet in LEAFLETS:
        a, b = ref[leaflet], new[leaflet]
        delta = b["apl_A2"] - a["apl_A2"]
        percent = 100.0 * delta / a["apl_A2"] if a["apl_A2"] else float("nan")
        expected = b["accessible_area_nm2"] / a["apl_nm2"] if a["apl_nm2"] else float("nan")
        leaflets[leaflet] = {"reference_lipids": a["lipids"], "lipids": b["lipids"],
                             "reference_area_A2": round(a["area_nm2"] * NM2_TO_A2, 1), "area_A2": round(b["area_nm2"] * NM2_TO_A2, 1),
                             "reference_protein_area_A2": round(a["protein_area_nm2"] * NM2_TO_A2, 1),
                             "protein_area_A2": round(b["protein_area_nm2"] * NM2_TO_A2, 1),
                             "reference_accessible_area_A2": round(a["accessible_area_nm2"] * NM2_TO_A2, 1),
                             "accessible_area_A2": round(b["accessible_area_nm2"] * NM2_TO_A2, 1),
                             "reference_apl_A2": round(a["apl_A2"], 2), "apl_A2": round(b["apl_A2"], 2),
                             "reference_head_apl_A2": round(a["head_apl_A2"], 2) if a["head_apl_A2"] else None,
                             "head_apl_A2": round(b["head_apl_A2"], 2) if b["head_apl_A2"] else None,
                             "delta_apl_A2": round(delta, 2), "delta_apl_percent": round(percent, 2),
                             "expected_lipids_from_reference_apl": round(expected, 1),
                             "lipid_count_deviation_percent": round(100.0 * (b["lipids"] - expected) / expected, 2) if expected else None,
                             "pass": abs(percent) <= tolerance_percent}
    return {"tolerance_percent": tolerance_percent, "leaflets": leaflets, "pass": all(v["pass"] for v in leaflets.values()),
            "composition": composition_table(ref, new)}


def per_lipid_rows(stage: dict) -> list[dict]:
    """Flat per-lipid records (stage, leaflet, resid, resname, x, y, local area) for the TSV output."""
    return [{"stage": stage["name"], "leaflet": r["leaflet"], "resid": r["resid"], "resname": r["resname"],
             "x_A": round(10.0 * r["x_nm"], 3), "y_A": round(10.0 * r["y_nm"], 3), "z_A": round(10.0 * r["z_nm"], 3),
             "voronoi_area_A2": round(r["voronoi_area_nm2"] * NM2_TO_A2, 3) if "voronoi_area_nm2" in r else None,
             "head_area_A2": round(r["head_area_nm2"] * NM2_TO_A2, 3) if r.get("head_area_nm2") is not None else None}
            for r in stage["lipids"]]


def region_reference(stage: dict, indices: list[int]) -> dict:
    """Leaflet APL of the reference membrane restricted to the lipids a crop kept, measured in the reference's own cell.

    The crop is a sample of the parent membrane; its lipids' Voronoi areas in the FULL periodic parent tell what the
    sampled region was packed like before the cut. A cut that keeps the parent's packing reproduces these values
    (up to the seam); the difference between this region mean and the whole-cell mean is sampling heterogeneity of a
    finite crop, not a construction error, and is reported separately as representativeness.
    """
    chosen = set(indices)
    out = {}
    for leaflet in LEAFLETS:
        rows = [r for r in stage["lipids"] if r["leaflet"] == leaflet and "voronoi_area_nm2" in r]
        inside = [r for r in rows if r["index"] in chosen]
        areas = np.array([r["voronoi_area_nm2"] for r in inside])
        whole = np.array([r["voronoi_area_nm2"] for r in rows])
        out[leaflet] = {"lipids": len(inside), "apl_nm2": float(areas.mean()) if len(areas) else float("nan"),
                        "apl_A2": round(float(areas.mean()) * NM2_TO_A2, 2) if len(areas) else None,
                        "whole_cell_apl_A2": round(float(whole.mean()) * NM2_TO_A2, 2) if len(whole) else None,
                        "composition": dict(sorted(Counter(r["resname"] for r in inside).items()))}
        if len(areas) and len(whole):
            out[leaflet]["region_vs_whole_cell_percent"] = round(100.0 * (areas.mean() - whole.mean()) / whole.mean(), 2)
    return out


def window_distribution(stage: dict, size_nm: np.ndarray, cropped: list[bool], samples: int = 12) -> dict:
    """Leaflet APL and composition of equal-size windows at many positions over the periodic embedded membrane.

    Every window is a crop the slicer could have taken (anchor rule, half-open interval), and its lipids keep their
    Voronoi areas of the uncut membrane, so the spread of the window means is the finite-window variability of this
    membrane at this crop size: what a correct crop may differ from the whole cell by, before any construction error.
    Returns per leaflet the arrays of window APL (A^2) and composition distance to the whole cell, plus the whole-cell
    values, for percentile grading (qc.classify_percentile).
    """
    cell = np.array(stage["box_nm"][:2], dtype=float)
    size = np.array(size_nm, dtype=float)
    offsets = [np.linspace(0.0, cell[a], samples, endpoint=False) if cropped[a] else np.array([0.0]) for a in range(2)]
    out = {}
    for leaflet in LEAFLETS:
        rows = [r for r in stage["lipids"] if r["leaflet"] == leaflet and "voronoi_area_nm2" in r]
        xy = np.array([[r["x_nm"], r["y_nm"]] for r in rows])
        areas = np.array([r["voronoi_area_nm2"] for r in rows]) * NM2_TO_A2
        names = np.array([r["resname"] for r in rows])
        whole = Counter(names.tolist())
        apl, distance, counts = [], [], []
        for ox in offsets[0]:
            for oy in offsets[1]:
                shifted = np.mod(xy - np.array([ox, oy]), cell)
                inside = np.all(shifted < size - 1e-9, axis=1)
                if inside.sum() < 3:
                    continue
                apl.append(float(areas[inside].mean()))
                distance.append(0.5 * sum(abs(whole[k] / len(rows) - Counter(names[inside].tolist())[k] / inside.sum()) for k in whole))
                counts.append(int(inside.sum()))
        out[leaflet] = {"window_apl_A2": np.array(apl), "window_composition_distance": np.array(distance), "window_lipids": np.array(counts),
                        "whole_cell_apl_A2": float(areas.mean()) if len(areas) else float("nan"), "windows": len(apl)}
    return out
