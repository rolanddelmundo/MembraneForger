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

from .packing import leaflet_slab

__all__ = ['lateral_rdf', 'rdf_features', 'leaflet_rdfs', 'region_rdfs', 'window_rdf_distribution', 'compare_rdf', 'rdf_comparison',
           'protein_density_profiles', 'open_area_grid', 'shell_areas', 'lateral_rdf_shell', 'leaflet_rdfs_shell',
           'random_reference_rdfs', 'compare_normalizations']

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
    of the lateral order at this crop size. A slice's deviation from the whole cell is graded against it one-sided
    (qc.classify_percentile: PASS below the 95th percentile, WARNING below the 99th, FAIL above; a deviation smaller
    than the windows' is never a defect): "within the parent-window distribution".
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


# ------------------------------------------------------------------------------------- around the protein, per species
PROFILE_R_MAX_NM = 3.0      # lateral distance from the protein covered by the profiles
PROFILE_DR_NM = 0.1         # 1 A shells
PROFILE_GRID_NM = 0.05      # grid that measures the membrane area in each shell around an irregular protein
PROFILE_ACCESSIBLE_NM = 0.3  # an anchor never sits closer than this to a protein heavy atom: nearer area is not lipid area
PROFILE_FIRST_SHELL_NM = 1.0  # the annulus of lipids in contact with the protein (summarized as one enrichment number)
PROFILE_HEAD_BAND_NM = 0.6  # protein atoms within this of the leaflet's median anchor z: its cross-section at headgroup depth
PROFILE_MIN_SPECIES = 10    # species with fewer molecules in a leaflet get no curve; below MIN_PAIR_SPECIES, "few"


def protein_density_profiles(stage: dict, r_max_nm: float = PROFILE_R_MAX_NM, dr_nm: float = PROFILE_DR_NM,
                             grid_nm: float = PROFILE_GRID_NM, band_nm: float | None = PROFILE_HEAD_BAND_NM) -> dict:
    """How each lipid species is distributed around the protein, leaflet by leaflet: density relative to its leaflet mean.

    For every lipid of a leaflet, the in-plane (x, y; periodic) distance from its headgroup anchor to the nearest protein
    heavy atom of the protein's cross-section at headgroup depth (within band_nm of the leaflet's median anchor z) is
    binned in dr_nm shells. The cross-section is taken at the anchors' own depth because a tilted helix shifts laterally
    by up to 2 nm across a leaflet: measured against all protein atoms of the leaflet's slab (band_nm=None,
    packing.leaflet_slab), a headgroup above the helix's buried end looks close to a protein that is not beside it, and
    even the total lipid density reads 1.3-1.4 at the protein (a single-pass protein tilted 53 degrees; 0.9-1.1 at
    headgroup depth). Each shell's membrane area is
    measured on a grid_nm grid (the protein is not a cylinder), so a shell's density is lipids / area. That density is
    divided by the species' mean density over the leaflet's lipid-accessible area (farther than PROFILE_ACCESSIBLE_NM
    from any protein atom): 1 = as common as anywhere in the leaflet, > 1 enriched next to the protein, < 1 depleted.
    Shells of one frame hold zero or one molecule of a species, so the cumulative ratio is reported too: the species'
    density within R of the protein over its leaflet mean (a depletion-enrichment index), whose noise falls as R grows
    and which tends to 1 far from the protein. The first-shell value is the cumulative ratio at PROFILE_FIRST_SHELL_NM.
    A leaflet without protein atoms at that depth has no profile.
    """
    cell = np.array(stage["box_nm"][:2], dtype=float)
    edges = np.arange(0.0, r_max_nm + dr_nm / 2, dr_nm)
    nx, ny = max(1, int(round(cell[0] / grid_nm))), max(1, int(round(cell[1] / grid_nm)))
    gx, gy = np.meshgrid((np.arange(nx) + 0.5) * cell[0] / nx, (np.arange(ny) + 0.5) * cell[1] / ny, indexing="ij")
    grid = np.column_stack([gx.ravel(), gy.ravel()])
    cell_area = float(cell[0] * cell[1]) / len(grid)
    out = {}
    protein = np.asarray(stage["protein_nm"], dtype=float).reshape(-1, 3)
    for leaflet in ("lower", "upper"):
        rows = [r for r in stage["lipids"] if r["leaflet"] == leaflet]
        if not rows:
            out[leaflet] = None
            continue
        if band_nm is None:
            low, high = leaflet_slab(stage, leaflet)
        else:  # the protein's cross-section in the headgroup layer, where the anchors are
            plane = float(np.median([r["z_nm"] for r in rows]))
            low, high = plane - band_nm, plane + band_nm
        inside = protein[(protein[:, 2] >= low) & (protein[:, 2] <= high)]
        if not len(inside):
            out[leaflet] = None
            continue
        tree = cKDTree(np.mod(inside[:, :2], cell), boxsize=cell)
        grid_d = tree.query(grid)[0]
        shell_area = np.histogram(grid_d, edges)[0] * cell_area
        open_d = grid_d[grid_d >= PROFILE_ACCESSIBLE_NM]
        accessible = float(len(open_d) * cell_area)
        within_area = np.array([(open_d < r).sum() for r in edges[1:]], dtype=float) * cell_area  # lipid area within each R
        first_area = float((open_d < PROFILE_FIRST_SHELL_NM).sum() * cell_area)
        anchors = np.mod(np.array([[r["x_nm"], r["y_nm"]] for r in rows]), cell)
        lipid_d = tree.query(anchors)[0]
        names = np.array([r["resname"] for r in rows])
        curves = {}
        for species in ["all"] + sorted(set(names)):
            mine = lipid_d if species == "all" else lipid_d[names == species]
            if species != "all" and len(mine) < PROFILE_MIN_SPECIES:
                continue
            bulk = len(mine) / accessible if accessible > 0 else float("nan")
            counts = np.histogram(mine, edges)[0]
            within = np.array([(mine < r).sum() for r in edges[1:]], dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(shell_area > 0, counts / shell_area / bulk, np.nan)
                cumulative = np.where(within_area > 0, within / within_area / bulk, np.nan)
            first = int((mine < PROFILE_FIRST_SHELL_NM).sum())
            curves[species] = {"r_nm": 0.5 * (edges[1:] + edges[:-1]), "density_ratio": ratio, "counts": counts, "n": len(mine),
                               "R_nm": edges[1:].copy(), "cumulative_ratio": cumulative,
                               "first_shell_n": first,
                               "first_shell_ratio": round(first / first_area / bulk, 3) if first_area > 0 and bulk > 0 else None,
                               "well_sampled": len(mine) >= MIN_PAIR_SPECIES}
        out[leaflet] = {"protein_atoms_in_slab": len(inside), "accessible_area_nm2": round(accessible, 3),
                        "first_shell_area_nm2": round(first_area, 3), "curves": curves}
    return out


# ----------------------------------------------------------------- accessible-shell-area normalization (a separate analysis)
# lateral_rdf divides the pair count of each shell by the area 2 pi r dr of the whole annulus times the leaflet's mean
# density N / A_acc (A_acc: the cell minus the protein footprint). Near a protein that overstates the room a lipid has: part
# of its annulus lies inside the protein, where no anchor can be, so even lipids placed at random in the open area give
# g(r) < 1 out to the distance at which most annuli stop touching the protein, and the deficit grows with the protein's
# share of the cell. The accessible-shell normalization below counts, for every anchor i, only the area of its annulus that
# is open to lipids,
#
#     E_k = sum_i rho A_i(r_k, r_k + dr),   rho = (N - 1) / A_acc  (M / A_acc for partners),
#     g_shell(r_k) = n_k / E_k,
#
# with A_i measured on the same periodic grid that gives A_acc, so g_shell = 1 exactly for uncorrelated anchors in the open
# area (random_reference_rdfs checks this against a geometry-matched random placement). It is reported next to the
# existing g(r), never in place of it; no gate reads it.


def open_area_grid(protein_xy: np.ndarray, box_xy: np.ndarray, grid_nm: float = 0.05, accessible_nm: float = 0.3) -> dict:
    """The periodic cell on a grid_nm grid, split into area open to lipid anchors and area closed by the protein cross-section.

    A grid cell is open when its centre is at least accessible_nm from every protein heavy atom of the cross-section
    (protein_xy, x and y in nm; none: the whole cell is open). Returns the open cell centres, the area of one cell and the
    open area A_acc = (open cells) x (cell area).
    """
    box_xy = np.asarray(box_xy, dtype=float)
    nx, ny = max(1, int(round(box_xy[0] / grid_nm))), max(1, int(round(box_xy[1] / grid_nm)))
    gx, gy = np.meshgrid((np.arange(nx) + 0.5) * box_xy[0] / nx, (np.arange(ny) + 0.5) * box_xy[1] / ny, indexing="ij")
    grid = np.column_stack([gx.ravel(), gy.ravel()])
    cell_area = float(box_xy[0] * box_xy[1]) / len(grid)
    protein_xy = np.asarray(protein_xy, dtype=float).reshape(-1, 2)
    if len(protein_xy):
        distance = cKDTree(np.mod(protein_xy, box_xy), boxsize=box_xy).query(grid)[0]
        grid = grid[distance >= accessible_nm]
    return {"open_xy": grid, "cell_area_nm2": cell_area, "accessible_nm2": float(len(grid) * cell_area), "box_xy": box_xy,
            "protein_atoms": int(len(protein_xy)), "grid_nm": grid_nm, "accessible_nm": accessible_nm}


def shell_areas(points: np.ndarray, open_grid: dict, edges: np.ndarray) -> np.ndarray:
    """Area (nm^2) of each anchor's annulus [edges[k], edges[k+1]) that is open to lipids: an (N, bins) array, periodic cell."""
    box_xy, open_xy, cell_area = open_grid["box_xy"], open_grid["open_xy"], open_grid["cell_area_nm2"]
    points = np.asarray(points, dtype=float).reshape(-1, 2)
    out = np.zeros((len(points), len(edges) - 1))
    for i, point in enumerate(points):
        delta = open_xy - point
        delta -= box_xy * np.round(delta / box_xy)
        out[i] = np.histogram(np.hypot(delta[:, 0], delta[:, 1]), bins=edges)[0] * cell_area
    return out


def lateral_rdf_shell(points: np.ndarray, open_grid: dict, dr: float = DR_NM, r_max: float | None = None,
                      partners: np.ndarray | None = None) -> dict:
    """g(r) of x, y points (or between two sets) normalized by each anchor's open annulus area; same keys as lateral_rdf.

    The expected counts per unit g are E_k = sum_i rho A_i(k) (halved for one set, whose pairs are counted once), so a curve
    of uncorrelated anchors in the open area is 1 at every r whatever the protein's share of the cell.
    """
    box_xy = open_grid["box_xy"]
    r_max = min(0.5 * float(box_xy.min()) - 1e-6, R_CAP_NM) if r_max is None else r_max
    edges = np.arange(0.0, r_max + dr, dr)
    points = np.asarray(points, dtype=float).reshape(-1, 2)
    area = open_grid["accessible_nm2"]
    a = cKDTree(np.mod(points, box_xy), boxsize=box_xy)
    if partners is None:
        pairs = a.query_pairs(r_max, output_type="ndarray")
        d = np.linalg.norm(minimum_image(points[pairs[:, 0]] - points[pairs[:, 1]], box_xy), axis=1) if len(pairs) else np.zeros(0)
        density, weight = (len(points) - 1) / area, 0.5
    else:
        b = cKDTree(np.mod(partners, box_xy), boxsize=box_xy)
        d = a.sparse_distance_matrix(b, r_max, output_type="coo_matrix").data
        density, weight = len(partners) / area, 1.0
    counts, _ = np.histogram(d, bins=edges)
    expected = weight * density * shell_areas(points, open_grid, edges).sum(axis=0)
    g = np.where(expected > 0, counts / np.maximum(expected, 1e-12), 0.0)
    r = 0.5 * (edges[1:] + edges[:-1])
    return {"r_nm": r, "g": g, "counts": counts, "n_points": int(len(points)), "density_nm2": float(density), "dr_nm": dr,
            "expected_counts_per_unit_g": expected, "normalization": "accessible shell area"}


def leaflet_rdfs_shell(stage: dict, band_nm: float = PROFILE_HEAD_BAND_NM, min_species: int = MIN_PAIR_SPECIES,
                       grid_nm: float = 0.05, accessible_nm: float = 0.3, r_max_nm: float | None = None) -> dict:
    """The all-anchor and same-species accessible-shell RDFs of each leaflet of a stage: {leaflet: {label: curve + features}}.

    The protein's cross-section is taken at the leaflet's headgroup depth (within band_nm of the median anchor z), as in
    protein_density_profiles. A leaflet without protein atoms at that depth is normalized by its whole cell.
    """
    box_xy = np.array(stage["box_nm"][:2], dtype=float)
    protein = np.asarray(stage["protein_nm"], dtype=float).reshape(-1, 3)
    out = {}
    for leaflet in ("lower", "upper"):
        rows = [r for r in stage["lipids"] if r["leaflet"] == leaflet]
        points = np.array([[r["x_nm"], r["y_nm"]] for r in rows])
        curves = {}
        if len(points) >= 3:
            plane = float(np.median([r["z_nm"] for r in rows]))
            section = protein[(protein[:, 2] >= plane - band_nm) & (protein[:, 2] <= plane + band_nm)][:, :2]
            grid = open_area_grid(section, box_xy, grid_nm, accessible_nm)
            curves["all"] = {**lateral_rdf_shell(points, grid, r_max=r_max_nm), "well_sampled": len(points) >= min_species,
                             "accessible_nm2": grid["accessible_nm2"]}
            for species in sorted({r["resname"] for r in rows}):
                sub = np.array([[r["x_nm"], r["y_nm"]] for r in rows if r["resname"] == species])
                if len(sub) >= min_species:
                    curves[f"{species}-{species}"] = {**lateral_rdf_shell(sub, grid, r_max=r_max_nm), "well_sampled": True,
                                                      "accessible_nm2": grid["accessible_nm2"]}
            for curve in curves.values():
                curve["features"] = rdf_features(curve)
        out[leaflet] = curves
    return out


def random_reference_rdfs(open_grid: dict, n_points: int, replicates: int = 200, hard_core_nm: float = 0.0, seed: int = 0,
                          dr: float = DR_NM, r_max: float | None = None) -> dict:
    """g(r) of n_points random anchors in the open area of a grid, with both normalizations, averaged over replicates.

    The reference matches the geometry of a real leaflet (cell, protein cross-section, number of anchors) and carries no
    lipid-lipid structure: positions are uniform over the open cells (hard_core_nm = 0: an ideal gas) or placed by random
    sequential addition with a minimum separation hard_core_nm. Returns the mean and the standard deviation over replicates
    of the existing g(r) (area = the open area) and of the accessible-shell g(r), on one r grid. The shell curve has to be
    1 beyond the hard core (an ideal gas: everywhere) within its replicate scatter; the existing one is not.
    """
    rng = np.random.default_rng(seed)
    box_xy, open_xy, cell = open_grid["box_xy"], open_grid["open_xy"], open_grid["grid_nm"]
    raw, shell = [], []
    for _ in range(replicates):
        chosen, attempts = [], 0
        while len(chosen) < n_points:
            attempts += 1
            if attempts > 2000 * n_points:
                raise SystemExit(f"cannot place {n_points} anchors {hard_core_nm} nm apart in {open_grid['accessible_nm2']:.1f} nm^2: "
                                 "the coverage is beyond random sequential addition (about 0.55)")
            candidate = open_xy[rng.integers(len(open_xy))] + rng.uniform(-0.5 * cell, 0.5 * cell, size=2)
            if hard_core_nm > 0 and chosen:
                delta = np.array(chosen) - candidate
                delta -= box_xy * np.round(delta / box_xy)
                if (np.hypot(delta[:, 0], delta[:, 1]) < hard_core_nm).any():
                    continue
            chosen.append(candidate)
        points = np.array(chosen)
        raw.append(lateral_rdf(points, box_xy, dr=dr, r_max=r_max, area=open_grid["accessible_nm2"]))
        shell.append(lateral_rdf_shell(points, open_grid, dr=dr, r_max=r_max))
    stack = lambda curves: np.array([c["g"] for c in curves])
    return {"r_nm": raw[0]["r_nm"], "replicates": replicates, "n_points": n_points, "hard_core_nm": hard_core_nm,
            "accessible_nm2": open_grid["accessible_nm2"], "protein_atoms": open_grid["protein_atoms"],
            "raw_mean": stack(raw).mean(axis=0), "raw_sd": stack(raw).std(axis=0, ddof=1),
            "shell_mean": stack(shell).mean(axis=0), "shell_sd": stack(shell).std(axis=0, ddof=1)}


def compare_normalizations(raw_curve: dict, shell_curve: dict, band_nm: tuple = (1.0, 2.5)) -> dict:
    """Mean g(r) of one leaflet curve over band_nm with both normalizations, and the offset of each from 1 (percent)."""
    inside = (raw_curve["r_nm"] >= band_nm[0]) & (raw_curve["r_nm"] <= band_nm[1])
    raw, shell = float(raw_curve["g"][inside].mean()), float(shell_curve["g"][inside].mean())
    return {"band_nm": list(band_nm), "mean_g_existing": round(raw, 4), "mean_g_shell": round(shell, 4),
            "offset_existing_percent": round(100 * (raw - 1), 2), "offset_shell_percent": round(100 * (shell - 1), 2)}
