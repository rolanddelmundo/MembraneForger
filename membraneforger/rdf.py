#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Lateral (2D) radial distribution functions of lipid headgroup anchors, per leaflet.
#//=============================================================
"""Does a cut, a backmapping or a minimization keep the lateral organization of the membrane?

The in-plane pair distribution of the headgroup anchors (lipids.LIPID_SLICE_ANCHORS / LIPID_AA_ANCHORS) of one
leaflet, g(r), is computed with the periodic minimum image of that stage's own x, y cell:

    g(r) = n_pairs(r, r + dr) / (N_pairs_total * 2 pi r dr / A)

for all anchors of a leaflet, and for the same-species pairs of every species with enough molecules. Each curve is
summarized by its first peak (position, height), first minimum, and the coordination number up to that minimum.

The decisive comparison is the embedded CG membrane against the sliced CG membrane: same representation, same force
field, so the curves must agree within the counting noise of the smaller sample. Comparisons across resolutions
(Martini anchors against atomistic P/O3/N atoms) are reported but judged loosely: the atomistic force field has its
own excluded-volume distances, so the first peak moving after backmapping and minimization is expected physics, not
a construction error. Distances are reported in A.
"""
import numpy as np
from scipy.spatial import cKDTree

__all__ = ['lateral_rdf', 'rdf_features', 'leaflet_rdfs', 'region_rdfs', 'window_rdf_distribution', 'compare_rdf', 'rdf_comparison']

MIN_PAIR_SPECIES = 20   # same-species curves are computed for species with at least this many molecules in the leaflet
DR_NM = 0.1             # 1 A bins: a leaflet of 150-250 anchors gives about 10 pairs per bin in the first shell
FIRST_SHELL_NM = (0.3, 1.2)  # where the first neighbour shell of headgroup anchors lies at either resolution
R_CAP_NM = 2.5


def lateral_rdf(points: np.ndarray, box_xy: np.ndarray, dr: float = DR_NM, r_max: float | None = None,
                partners: np.ndarray | None = None, area: float | None = None) -> dict:
    """2D periodic RDF of one set of x, y points (or between two sets) in a rectangular cell; r in nm, g dimensionless.

    `area` is the area the points can occupy (the leaflet's lipid-accessible area, cell minus protein footprint);
    with it g(r) tends to 1 at large r whatever fraction of the cell the protein covers, so curves of cells with
    different protein fractions (the whole embedded membrane and its slice) are comparable. Default: the cell area.
    """
    box_xy = np.asarray(box_xy, dtype=float)
    r_max = min(0.5 * float(box_xy.min()) - 1e-6, R_CAP_NM) if r_max is None else r_max
    edges = np.arange(0.0, r_max + dr, dr)
    area = float(box_xy[0] * box_xy[1]) if area is None else float(area)
    a = cKDTree(np.mod(points, box_xy), boxsize=box_xy)
    if partners is None:
        pairs = a.query_pairs(r_max, output_type="ndarray")
        d = np.linalg.norm(minimum_image(points[pairs[:, 0]] - points[pairs[:, 1]], box_xy), axis=1) if len(pairs) else np.zeros(0)
        n_pairs = 0.5 * len(points) * (len(points) - 1)
        density = (len(points) - 1) / area
    else:
        b = cKDTree(np.mod(partners, box_xy), boxsize=box_xy)
        matrix = a.sparse_distance_matrix(b, r_max, output_type="coo_matrix")
        d = matrix.data
        n_pairs = float(len(points) * len(partners))
        density = len(partners) / area
    counts, _ = np.histogram(d, bins=edges)
    shells = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    g = counts / (n_pairs * shells / area) if n_pairs else np.zeros(len(counts))
    r = 0.5 * (edges[1:] + edges[:-1])
    return {"r_nm": r, "g": g, "counts": counts, "n_points": int(len(points)), "density_nm2": float(density), "dr_nm": dr,
            "expected_counts_per_unit_g": n_pairs * shells / area}


def minimum_image(delta: np.ndarray, box_xy: np.ndarray) -> np.ndarray:
    """Minimum-image x, y separations."""
    return delta - box_xy * np.round(delta / box_xy)


def rdf_features(curve: dict, smooth: int = 5, shell: tuple = FIRST_SHELL_NM) -> dict:
    """First-shell peak and minimum of one g(r), robust to counting noise; positions in A.

    The curve is smoothed with a `smooth`-bin running mean; the first peak is the maximum of the smoothed curve within
    the first-shell range and the first minimum the minimum of the smoothed curve in the 0.8 nm beyond the peak. The
    first-shell centroid, the g - 1 weighted mean position of the excess density between the shell start and the
    first minimum, is reported as well: it moves with the whole shell, not with one noisy bin.
    """
    r, g = curve["r_nm"], curve["g"]
    s = np.convolve(g, np.ones(smooth) / smooth, mode="same")
    window = np.nonzero((r >= shell[0]) & (r <= shell[1]))[0]
    if not len(window):
        return {"first_peak_A": None, "first_peak_height": None, "first_minimum_A": None, "first_minimum_value": None,
                "first_shell_centroid_A": None, "coordination_number": None, "peak_resolved": False}
    peak = int(window[np.argmax(s[window])])
    beyond = np.nonzero((r > r[peak]) & (r <= r[peak] + 0.8))[0]
    trough = int(beyond[np.argmin(s[beyond])]) if len(beyond) else None
    coordination = centroid = None
    if trough is not None:
        dr = curve["dr_nm"]
        shells = np.pi * ((r[:trough + 1] + 0.5 * dr) ** 2 - (r[:trough + 1] - 0.5 * dr) ** 2)
        coordination = float((curve["density_nm2"] * g[:trough + 1] * shells).sum())
        inside = (r >= shell[0]) & (r <= r[trough])
        excess = np.clip(g[inside] - 1.0, 0.0, None)
        centroid = float((r[inside] * excess).sum() / excess.sum()) if excess.sum() > 0 else None
    return {"first_peak_A": round(10.0 * float(r[peak]), 2), "first_peak_height": round(float(s[peak]), 3),
            "first_minimum_A": round(10.0 * float(r[trough]), 2) if trough is not None else None,
            "first_minimum_value": round(float(s[trough]), 3) if trough is not None else None,
            "first_shell_centroid_A": round(10.0 * centroid, 2) if centroid is not None else None,
            "coordination_number": round(coordination, 2) if coordination is not None else None,
            "peak_resolved": bool(s[peak] >= 1.05)}  # a shell weaker than 5 % above uniform has no position worth grading


def leaflet_rdfs(stage: dict, min_species: int = MIN_PAIR_SPECIES, accessible_nm2: dict | None = None) -> dict:
    """All-anchor and same-species RDFs of each leaflet of a stage: {leaflet: {label: curve + features}}."""
    # accessible_nm2: {leaflet: lipid-accessible area} from the tessellation (packing.measure_packing), for the normalization.
    out = {}
    box_xy = np.array(stage["box_nm"][:2], dtype=float)
    for leaflet in ("lower", "upper"):
        rows = [r for r in stage["lipids"] if r["leaflet"] == leaflet]
        points = np.array([[r["x_nm"], r["y_nm"]] for r in rows])
        area = (accessible_nm2 or {}).get(leaflet)
        curves = {}
        if len(points) >= 3:
            curves["all"] = {**lateral_rdf(points, box_xy, area=area), "well_sampled": len(points) >= min_species}
        for species in sorted({r["resname"] for r in rows}):
            sub = np.array([[r["x_nm"], r["y_nm"]] for r in rows if r["resname"] == species])
            if len(sub) >= min_species:
                curves[f"{species}-{species}"] = {**lateral_rdf(sub, box_xy, area=area), "well_sampled": True}
        for curve in curves.values():
            curve["features"] = rdf_features(curve)
        out[leaflet] = curves
    return out


def region_rdfs(stage: dict, indices: list[int], r_max_nm: float, accessible_nm2: dict | None = None,
                min_species: int = MIN_PAIR_SPECIES) -> dict:
    """RDFs of a subset of a stage's lipids (those a crop kept) measured in the stage's own, uncut periodic cell.

    The same lipids with the same neighbours, before the cut: the construction reference for a slice's RDF. r_max is
    the slice's (half its smaller edge), and the accessible area is the slice's, so the two curves share one scale.
    """
    chosen = set(indices)
    out = {}
    box_xy = np.array(stage["box_nm"][:2], dtype=float)
    for leaflet in ("lower", "upper"):
        rows = [r for r in stage["lipids"] if r["leaflet"] == leaflet and r["index"] in chosen]
        points = np.array([[r["x_nm"], r["y_nm"]] for r in rows])
        area = (accessible_nm2 or {}).get(leaflet)
        curves = {}
        if len(points) >= 3:
            curves["all"] = {**lateral_rdf(points, box_xy, r_max=r_max_nm, area=area), "well_sampled": len(points) >= min_species}
        for species in sorted({r["resname"] for r in rows}):
            sub = np.array([[r["x_nm"], r["y_nm"]] for r in rows if r["resname"] == species])
            if len(sub) >= min_species:
                curves[f"{species}-{species}"] = {**lateral_rdf(sub, box_xy, r_max=r_max_nm, area=area), "well_sampled": True}
        for curve in curves.values():
            curve["features"] = rdf_features(curve)
        out[leaflet] = curves
    return out


def window_rdf_distribution(stage: dict, size_nm, cropped: list[bool], r_max_nm: float, samples: int = 12) -> dict:
    """RMS deviation of the all-anchor g(r) of equal-size windows of a stage from the stage's whole-cell curve, per leaflet.

    Every window is a crop the slicer could have taken (anchor rule, half-open interval), measured in the uncut cell
    with the lipid density of the window itself, so the spread of these deviations is the finite-window variability
    of the lateral order at this crop size. A slice's deviation from the whole cell is graded against it (central
    95 % PASS, 95-99 % WARNING, beyond FAIL): "within the parent-window distribution".
    """
    cell = np.array(stage["box_nm"][:2], dtype=float)
    size = np.array(size_nm, dtype=float)
    offsets = [np.linspace(0.0, cell[a], samples, endpoint=False) if cropped[a] else np.array([0.0]) for a in range(2)]
    out = {}
    for leaflet in ("lower", "upper"):
        rows = [r for r in stage["lipids"] if r["leaflet"] == leaflet]
        xy = np.array([[r["x_nm"], r["y_nm"]] for r in rows])
        if len(xy) < 3:
            out[leaflet] = {"rms": np.zeros(0), "windows": 0}
            continue
        whole = lateral_rdf(xy, cell, r_max=r_max_nm)
        deviations = []
        for ox in offsets[0]:
            for oy in offsets[1]:
                inside = np.all(np.mod(xy - np.array([ox, oy]), cell) < size - 1e-9, axis=1)
                if inside.sum() < 3:
                    continue
                window = lateral_rdf(xy[inside], cell, r_max=r_max_nm, area=float(size[0] * size[1]))
                n = min(len(window["g"]), len(whole["g"]))
                r = whole["r_nm"][:n]
                band = (r >= 0.3) & (r <= 2.0)
                deviations.append(float(np.sqrt(np.mean((window["g"][:n][band] - whole["g"][:n][band]) ** 2))))
        out[leaflet] = {"rms": np.array(deviations), "windows": len(deviations)}
    return out


def compare_rdf(reference: dict, sample: dict, r_window_nm: tuple = (0.3, 2.0)) -> dict:
    """Difference between two curves on their common r grid: peak shift, RMS difference, and that RMS in noise units."""
    n = min(len(reference["g"]), len(sample["g"]))
    r = reference["r_nm"][:n]
    inside = (r >= r_window_nm[0]) & (r <= r_window_nm[1])
    diff = sample["g"][:n][inside] - reference["g"][:n][inside]
    rms = float(np.sqrt(np.mean(diff ** 2))) if inside.any() else float("nan")
    # Counting noise of each curve: var(g) = g / expected counts at g = 1 (Poisson), summed over both curves.
    noise = np.zeros(n)
    for curve in (reference, sample):
        expected = curve["expected_counts_per_unit_g"][:n]
        noise += np.where(expected > 0, np.maximum(reference["g"][:n], 1e-3) / np.maximum(expected, 1e-9), 0.0)
    noise_rms = float(np.sqrt(np.mean(noise[inside]))) if inside.any() else float("nan")
    a, b = reference["features"], sample["features"]
    shift = (b["first_peak_A"] - a["first_peak_A"]) if a["first_peak_A"] is not None and b["first_peak_A"] is not None else None
    return {"first_peak_shift_A": round(shift, 2) if shift is not None else None,
            "first_peak_height_change": round(b["first_peak_height"] - a["first_peak_height"], 3)
            if a["first_peak_height"] is not None and b["first_peak_height"] is not None else None,
            "rms_difference": round(rms, 4), "expected_rms_noise": round(noise_rms, 4),
            "rms_in_noise_units": round(rms / noise_rms, 2) if noise_rms else None,
            "window_A": [10.0 * r_window_nm[0], 10.0 * r_window_nm[1]]}


def rdf_comparison(reference: dict, sample: dict, max_peak_shift_a: float, max_noise_units: float, warning_shift_a: float = 0.5) -> dict:
    """Compare every leaflet curve of a sample with the reference's; the 'all' curve of each leaflet decides the status.

    First-shell peak shift: PASS within warning_shift_a, WARNING within max_peak_shift_a, FAIL beyond; the RMS
    difference of the curves must stay within max_noise_units of their combined counting noise (else FAIL).
    """
    out = {"pass": True, "status": "PASS", "leaflets": {}}
    for leaflet in ("lower", "upper"):
        rows = {}
        for label, curve in sample.get(leaflet, {}).items():
            if label not in reference.get(leaflet, {}):
                continue
            comparison = compare_rdf(reference[leaflet][label], curve)
            comparison["well_sampled"] = bool(curve["well_sampled"] and reference[leaflet][label]["well_sampled"])
            shift = comparison["first_peak_shift_A"]
            noise = comparison["rms_in_noise_units"]
            resolved = reference[leaflet][label]["features"].get("peak_resolved") and curve["features"].get("peak_resolved")
            if not resolved:
                shift = 0.0  # no resolved first shell on both curves: only the curve difference is graded
                comparison["note"] = "first-shell peak too weak to locate on both curves; graded on the RMS difference alone"
            status = "INSUFFICIENT SAMPLING" if (shift is None or noise is None or not comparison["well_sampled"]) else \
                "FAIL" if (abs(shift) > max_peak_shift_a or noise > max_noise_units) else \
                "WARNING" if abs(shift) > warning_shift_a else "PASS"
            comparison["status"], comparison["pass"] = status, status != "FAIL"
            rows[label] = comparison
            if label == "all":
                out["status"] = max(out["status"], status, key=("PASS", "WARNING", "INSUFFICIENT SAMPLING", "FAIL").index)
                out["pass"] = out["pass"] and status != "FAIL"
        out["leaflets"][leaflet] = rows
    out["criteria"] = {"first_peak_shift_A": {"pass": warning_shift_a, "warning": max_peak_shift_a}, "max_rms_in_noise_units": max_noise_units,
                       "decided_by": "the all-anchor curve of each leaflet; species curves are reported"}
    return out
