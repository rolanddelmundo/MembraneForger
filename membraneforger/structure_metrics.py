#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Structural metrics of one membrane stage: thickness, protein tilt and depth, z profiles, hydration, interdigitation.
#//=============================================================
"""Does a stage of the build still look like the membrane it was cut from, and is the protein still where it was?

Everything here is a pure function of coordinates (nm internally; A where a key says so) that works alike on Martini
beads and CHARMM36 atoms, because every quantity is defined on one headgroup anchor per lipid (lipids.anchor_bead)
and on heavy-atom coordinates:

    thickness          median anchor z of the upper leaflet minus that of the lower leaflet (sterols left out), and a
                       periodic XY map of the same difference per cell
    protein tilt       angle between z and the principal axis of the atoms within a slab around the midplane
    insertion depth    mean z of those transmembrane atoms relative to the midplane (and the centre-of-mass z)
    z profiles         number density of water, headgroups, tails, sterol, protein and ions along z
    core hydration     water number density in the hydrophobic core relative to bulk water, and whether a connected
                       column of water crosses the core (a transmembrane water path)
    interdigitation    overlap integral of the upper- and lower-leaflet tail density profiles
    asymmetry          per-species counts and mole fractions per leaflet

Graded quantities return qc.metric records classified with qc.classify against (pass, warning) ceilings that the
CALLER supplies, because the same metric has different tolerances at the coarse-grained slice gate and after
backmapping. A "stage" is the record of packing.py: {"name", "resolution", "box_nm", "midplane_nm", "lipids", "protein_nm"}.
"""
import math
import re
from collections import Counter

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from .lipids import LIPID_AA_ANCHORS, anchor_bead, anchor_xyz, cg_name, is_sterol, leaflet_of
from .qc import classify, metric
from .structio import element

__all__ = ['TAIL_CARBON', 'CERAMIDE_CARBON', 'SOLVENT_RESNAMES', 'WATER_OXYGEN', 'bilayer_thickness', 'thickness_change',
           'protein_orientation', 'orientation_change', 'z_density_profiles', 'aa_groups_from_atoms', 'core_hydration',
           'tail_interdigitation', 'interdigitation_change', 'leaflet_asymmetry', 'aa_leaflet_tail_split', 'is_tail_carbon']

TAIL_CARBON = re.compile(r"^C[23]\d+$")      # CHARMM36 acyl chain carbons C22..C218 (sn-2) and C32..C318 (sn-1)
CERAMIDE_CARBON = re.compile(r"^C\d+[SF]$")  # PSM / GLPA ceramide chain carbons (sphingosine S, fatty acyl F)
SOLVENT_RESNAMES = frozenset({"TIP3", "SOD", "CLA", "POT"})
WATER_OXYGEN = ("TIP3", "OH2")
MIN_HEADS_PER_LEAFLET = 5   # fewer non-sterol anchors than this and the thickness falls back to all lipids
WATER_LINK_NM = 0.35        # hydrogen-bond distance between water oxygens: links a connected water column
PATH_BIN_NM = 0.2
NM_TO_A = 10.0


def _xyz(a) -> np.ndarray:
    """Any coordinate input as an (n, 3) float array."""
    return np.asarray(a, dtype=float).reshape(-1, 3)


def _grid_cell(xy: np.ndarray, box_xy: np.ndarray, shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    """Periodic grid cell (i, j) of every x, y point."""
    frac = np.mod(xy, box_xy) / box_xy
    i = np.minimum((frac[:, 0] * shape[0]).astype(int), shape[0] - 1)
    j = np.minimum((frac[:, 1] * shape[1]).astype(int), shape[1] - 1)
    return i, j


def _cell_means(xyz: np.ndarray, box_xy: np.ndarray, shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    """Mean z per periodic XY cell and the count per cell."""
    total, count = np.zeros(shape), np.zeros(shape)
    if len(xyz):
        i, j = _grid_cell(xyz[:, :2], box_xy, shape)
        np.add.at(total, (i, j), xyz[:, 2])
        np.add.at(count, (i, j), 1.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(count > 0, total / count, np.nan), count


def bilayer_thickness(stage: dict, bins_nm: float = 1.0) -> dict:
    """Headgroup-to-headgroup thickness (median upper anchor z - median lower anchor z) and its periodic XY map."""
    # Sterols are left out: their anchor (ROH / O3) sits ~0.5 nm below the phosphate plane and would thin the
    # apparent headgroup plane in a cholesterol-rich leaflet. The local map uses cell means, which are noisier than
    # the global medians, so a cell holding only a tilted lipid or two is read with its count in mind.
    box_xy = np.array(stage["box_nm"][:2], dtype=float)
    anchors = {}
    for leaflet in ("lower", "upper"):
        rows = [r for r in stage["lipids"] if r["leaflet"] == leaflet]
        heads = [r for r in rows if not is_sterol(r["resname"])]
        chosen = heads if len(heads) >= MIN_HEADS_PER_LEAFLET else rows
        anchors[leaflet] = np.array([[r["x_nm"], r["y_nm"], r["z_nm"]] for r in chosen], dtype=float).reshape(-1, 3)
    planes = {k: float(np.median(v[:, 2])) if len(v) else float("nan") for k, v in anchors.items()}
    thickness = planes["upper"] - planes["lower"]
    shape = (max(1, int(round(box_xy[0] / bins_nm))), max(1, int(round(box_xy[1] / bins_nm))))
    upper, n_up = _cell_means(anchors["upper"], box_xy, shape)
    lower, n_lo = _cell_means(anchors["lower"], box_xy, shape)
    local = np.where((n_up > 0) & (n_lo > 0), upper - lower, np.nan)
    finite = local[np.isfinite(local)]
    return {"thickness_nm": thickness, "thickness_A": thickness * NM_TO_A, "leaflet_planes_nm": planes,
            "n_per_leaflet": {k: int(len(v)) for k, v in anchors.items()}, "local_map": local, "bins_nm": bins_nm,
            "cell_nm": [float(box_xy[0] / shape[0]), float(box_xy[1] / shape[1])],
            "local_sd_nm": float(finite.std(ddof=1)) if len(finite) > 1 else float("nan"),
            "local_min_nm": float(finite.min()) if len(finite) else float("nan"),
            "local_max_nm": float(finite.max()) if len(finite) else float("nan"),
            "cells_filled": int(len(finite)), "cells": int(local.size)}


def thickness_change(reference: dict, sample: dict, pass_percent: float, warning_percent: float, stage: str = "sample") -> dict:
    """Percent change of the bilayer thickness of a sample against its reference, graded with the caller's ceilings."""
    # Both arguments are bilayer_thickness() results. The local-map spread is reported alongside, not graded: it
    # mixes undulation with the lateral heterogeneity of a mixed membrane.
    ref, new = reference["thickness_nm"], sample["thickness_nm"]
    percent = 100.0 * (new - ref) / ref if ref else float("nan")
    status = classify(percent, pass_percent, warning_percent)
    records = [metric("bilayer thickness", stage, new * NM_TO_A, "A", reference=ref * NM_TO_A, deviation=percent,
                      n=sum(sample["n_per_leaflet"].values()), status=status,
                      reason=f"{percent:+.2f} % vs reference (PASS <= {pass_percent} %, WARNING <= {warning_percent} %)"),
               metric("local thickness spread", stage, sample["local_sd_nm"] * NM_TO_A, "A", reference=reference["local_sd_nm"] * NM_TO_A,
                      n=sample["cells_filled"], status="NOT RUN", reason="standard deviation of the local thickness map; reported, not graded")]
    return {"reference_A": ref * NM_TO_A, "thickness_A": new * NM_TO_A, "delta_A": (new - ref) * NM_TO_A, "percent": percent,
            "status": status, "records": records}


def protein_orientation(protein_heavy_nm, midplane_nm: float, ca_mask=None, slab_half_nm: float = 1.5,
                        headgroup_planes_nm: tuple[float, float] | None = None, reference: dict | None = None) -> dict:
    """Tilt (deg) of the protein's transmembrane principal axis from z, its insertion depth (nm) and inversion vs a reference."""
    # The axis is the leading eigenvector of the covariance of the atoms within slab_half_nm of the midplane: only
    # the membrane-spanning part decides, so a G protein or an antibody hanging off one side does not tilt the
    # answer. An eigenvector has no sign; it is pointed towards the side the rest of the protein lies on (the
    # extramembrane mass tells top from bottom), so that an upside-down protein flips the sign relative to the
    # reference while tilt alone (arccos |axis_z|) cannot see the inversion.
    xyz = _xyz(protein_heavy_nm)
    sites = xyz[np.asarray(ca_mask, dtype=bool)] if ca_mask is not None else xyz
    dz = sites[:, 2] - midplane_nm
    slab = sites[np.abs(dz) < slab_half_nm]
    com_z = float(xyz[:, 2].mean() - midplane_nm) if len(xyz) else float("nan")
    if len(slab) < 3:
        axis, tilt, depth = None, float("nan"), float("nan")
    else:
        values, vectors = np.linalg.eigh(np.cov(slab.T))
        axis = vectors[:, int(np.argmax(values))]
        axis = axis if axis[2] >= 0 else -axis
        tilt = math.degrees(math.acos(min(1.0, abs(float(axis[2])))))
        depth = float(slab[:, 2].mean() - midplane_nm)
    polarity = 0 if not len(slab) or not np.isfinite(com_z) else int(np.sign(com_z - depth)) if abs(com_z - depth) > 1e-6 else 0
    signed_axis = None if axis is None else (axis * (polarity if polarity else 1)).tolist()
    inverted = None
    if reference is not None and polarity and reference.get("polarity"):
        inverted = bool(polarity != reference["polarity"])
    mismatch = None
    if headgroup_planes_nm is not None and len(slab):
        lo, hi = sorted(float(p) for p in headgroup_planes_nm)
        mismatch = float(np.mean((slab[:, 2] < lo) | (slab[:, 2] > hi)))
    embedded = None
    if headgroup_planes_nm is not None and len(xyz):
        lo, hi = sorted(float(p) for p in headgroup_planes_nm)
        embedded = float(np.mean((xyz[:, 2] >= lo) & (xyz[:, 2] <= hi)))
    return {"tilt_deg": tilt, "axis": None if axis is None else axis.tolist(), "signed_axis": signed_axis, "polarity": polarity,
            "depth_nm": depth, "depth_A": depth * NM_TO_A, "com_z_nm": com_z, "com_z_A": com_z * NM_TO_A,
            "n_slab": int(len(slab)), "n_sites": int(len(sites)), "slab_half_nm": slab_half_nm, "inverted": inverted,
            "hydrophobic_mismatch_fraction": mismatch, "embedded_fraction": embedded}


def orientation_change(reference: dict, sample: dict, tilt_pass_deg: float, tilt_warning_deg: float, depth_pass_A: float,
                       depth_warning_A: float, stage: str = "sample") -> dict:
    """Tilt change (deg) and depth change (A) of a protein_orientation() result against its reference, graded."""
    d_tilt = sample["tilt_deg"] - reference["tilt_deg"]
    d_depth = (sample["depth_nm"] - reference["depth_nm"]) * NM_TO_A
    d_com = (sample["com_z_nm"] - reference["com_z_nm"]) * NM_TO_A
    inverted = sample.get("inverted")
    if inverted is None and sample.get("polarity") and reference.get("polarity"):
        inverted = bool(sample["polarity"] != reference["polarity"])
    records = [metric("protein tilt", stage, sample["tilt_deg"], "deg", reference=reference["tilt_deg"], deviation=d_tilt, n=sample["n_slab"],
                      status=classify(d_tilt, tilt_pass_deg, tilt_warning_deg),
                      reason=f"{d_tilt:+.2f} deg vs reference (PASS <= {tilt_pass_deg}, WARNING <= {tilt_warning_deg})"),
               metric("protein insertion depth", stage, sample["depth_A"], "A", reference=reference["depth_A"], deviation=d_depth, n=sample["n_slab"],
                      status=classify(d_depth, depth_pass_A, depth_warning_A),
                      reason=f"{d_depth:+.2f} A vs reference (PASS <= {depth_pass_A}, WARNING <= {depth_warning_A})"),
               metric("protein centre of mass z", stage, sample["com_z_A"], "A", reference=reference["com_z_A"], deviation=d_com,
                      status="NOT RUN", reason="whole-protein centre relative to the midplane; reported, not graded"),
               metric("protein inversion", stage, bool(inverted) if inverted is not None else None, "flag",
                      status="FAIL" if inverted else ("PASS" if inverted is False else "INSUFFICIENT SAMPLING"),
                      reason="transmembrane axis points to the other side of the membrane than in the reference" if inverted
                      else "same side as the reference" if inverted is False else "no extramembrane mass to tell top from bottom")]
    return {"tilt_change_deg": d_tilt, "depth_change_A": d_depth, "com_change_A": d_com, "inverted": inverted, "records": records}


def z_density_profiles(groups: dict, box_nm, bin_nm: float = 0.1) -> dict:
    """Number density (nm^-3) along z of every labelled group of coordinates, after wrapping z into [0, box_z)."""
    # The bin width is box_z / n so that the bins tile the cell exactly; sum(rho) * A_xy * bin = N for every group.
    box = np.asarray(box_nm, dtype=float)
    area = float(box[0] * box[1])
    n_bins = max(1, int(math.ceil(box[2] / bin_nm)))
    edges = np.linspace(0.0, box[2], n_bins + 1)
    width = float(box[2] / n_bins)
    density, counts, n = {}, {}, {}
    for label, xyz in groups.items():
        pts = _xyz(xyz)
        z = np.mod(pts[:, 2], box[2]) if len(pts) else np.zeros(0)
        c, _ = np.histogram(z, bins=edges)
        counts[label], n[label] = c, int(len(pts))
        density[label] = c / (area * width)
    return {"z_nm": 0.5 * (edges[1:] + edges[:-1]), "bin_nm": width, "area_nm2": area, "density_nm3": density, "counts": counts, "n": n}


def is_tail_carbon(name: str) -> bool:
    """True for a CHARMM36 acyl or ceramide chain carbon name."""
    return bool(TAIL_CARBON.match(name) or CERAMIDE_CARBON.match(name))


def aa_groups_from_atoms(atoms_nm: list[dict], lipid_resnames, anchor_atom_by_resname: dict | None = None,
                         solvent=SOLVENT_RESNAMES) -> dict:
    """Water, ions, cholesterol, headgroup, tail and protein coordinate groups from all-atom records (nm)."""
    # Headgroups are the anchor atoms (P, or the GLPA amide N); tails the chain carbons; cholesterol all of its heavy
    # atoms (its ring carbons C20..C27 would otherwise be mistaken for chain carbons, so sterols are split off first).
    anchors = LIPID_AA_ANCHORS if anchor_atom_by_resname is None else anchor_atom_by_resname
    lipids = set(lipid_resnames)
    out = {k: [] for k in ("water", "ions", "cholesterol", "headgroups", "tails", "protein")}
    for a in atoms_nm:
        rn, name, p = a["resname"], a["atom"], (a["x"], a["y"], a["z"])
        if rn == WATER_OXYGEN[0]:
            if name == WATER_OXYGEN[1]:
                out["water"].append(p)
        elif rn in solvent:
            out["ions"].append(p)
        elif rn in lipids:
            if element(name) == "H":
                continue
            if is_sterol(cg_name(rn)):
                out["cholesterol"].append(p)
            elif is_tail_carbon(name):
                out["tails"].append(p)
            else:
                wanted = anchors.get(rn)
                wanted = (wanted,) if isinstance(wanted, str) else tuple(wanted or ())
                if name in wanted:
                    out["headgroups"].append(p)
        elif element(name) != "H":
            out["protein"].append(p)
    return {k: np.array(v, dtype=float).reshape(-1, 3) for k, v in out.items()}


def _components(points: np.ndarray, box: np.ndarray, link_nm: float) -> np.ndarray:
    """Connected-component label of every point when points closer than link_nm (periodic) are linked."""
    if not len(points):
        return np.zeros(0, dtype=int)
    tree = cKDTree(np.mod(points, box), boxsize=box)
    pairs = tree.query_pairs(link_nm, output_type="ndarray")
    n = len(points)
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    return connected_components(graph, directed=False)[1]


def core_hydration(water_nm, protein_nm, midplane_nm: float, core_half_nm: float, box_nm, bulk_region: tuple[float, float] | None = None,
                   protein_exclusion_nm: float = 0.5, thickness_nm: float | None = None, stage: str = "sample") -> dict:
    """Water in the hydrophobic core relative to bulk water (percent, graded 1 / 5 %) and the transmembrane water path test."""
    # Waters within protein_exclusion_nm of a protein heavy atom are "protein-associated" (a pore, a vestibule, a
    # crevice) and left out of the core count: they are the protein's hydration, not a defect of the membrane.
    # The core density is divided by the density of bulk water (|z - midplane| beyond the headgroups + 1 nm, or the
    # given absolute z range), so the answer does not depend on the water model's density. One frame cannot tell a
    # transient water from a stable pore; the path test therefore grades only PASS (no spanning column) or FAIL.
    box = np.asarray(box_nm, dtype=float)
    water, protein = _xyz(water_nm), _xyz(protein_nm)
    area = float(box[0] * box[1])
    dz = water[:, 2] - midplane_nm
    dz = dz - box[2] * np.round(dz / box[2])  # minimum image: the membrane may straddle the z boundary of the cell
    core = np.abs(dz) < core_half_nm
    associated = np.zeros(len(water), dtype=bool)
    if len(protein) and len(water):
        tree = cKDTree(np.mod(protein, box), boxsize=box)
        near = tree.query_ball_point(np.mod(water, box), protein_exclusion_nm, return_length=True)
        associated = np.asarray(near) > 0
    free_core = core & ~associated
    if bulk_region is None:
        half = 0.5 * thickness_nm if thickness_nm else core_half_nm + 1.0
        bulk = np.abs(dz) > half + 1.0
        bulk_height = max(box[2] - 2.0 * (half + 1.0), 0.0)
    else:
        lo, hi = sorted(float(v) for v in bulk_region)
        bulk = (water[:, 2] >= lo) & (water[:, 2] < hi)
        bulk_height = hi - lo
    core_density = float(free_core.sum() / (area * 2.0 * core_half_nm))
    bulk_density = float(bulk.sum() / (area * bulk_height)) if bulk_height > 0 else float("nan")
    ratio = 100.0 * core_density / bulk_density if bulk_density and np.isfinite(bulk_density) else float("nan")
    # Water path: every 0.2 nm bin across the core occupied AND one hydrogen-bonded component reaching both faces.
    edges = np.arange(-core_half_nm, core_half_nm + 1e-9, PATH_BIN_NM)
    edges = np.append(edges, core_half_nm) if edges[-1] < core_half_nm - 1e-9 else edges
    occupied, _ = np.histogram(dz[free_core], bins=edges)
    wide = (np.abs(dz) < core_half_nm + 0.3) & ~associated
    labels = _components(water[wide], box, WATER_LINK_NM)
    sizes, spanning = [], False
    for k in np.unique(labels):
        member = dz[wide][labels == k]
        sizes.append(int(len(member)))
        spanning = spanning or (member.min() <= -core_half_nm and member.max() >= core_half_nm)
    path = bool(np.all(occupied > 0)) and spanning
    status = classify(ratio, 1.0, 5.0)
    records = [metric("core water", stage, ratio, "% of bulk density", n=int(free_core.sum()), status=status,
                      reason=f"{int(free_core.sum())} free water(s) in |z - midplane| < {core_half_nm} nm, "
                             f"{int((core & associated).sum())} protein-associated excluded"),
               metric("transmembrane water path", stage, path, "flag", n=int(wide.sum()), status="FAIL" if path else "PASS",
                      reason=("a hydrogen-bonded water column spans the core" if path else "no water column spans the core")
                      + "; one frame cannot distinguish a transient from a persistent path")]
    return {"n_core": int(core.sum()), "n_core_free": int(free_core.sum()), "n_protein_associated": int((core & associated).sum()),
            "n_bulk": int(bulk.sum()), "core_density_nm3": core_density, "bulk_density_nm3": bulk_density, "ratio_percent": ratio,
            "status": status, "transmembrane_water_path": path, "component_sizes": sorted(sizes, reverse=True),
            "bins_occupied": int((occupied > 0).sum()), "bins": int(len(occupied)), "records": records}


def tail_interdigitation(tail_upper_nm, tail_lower_nm, midplane_nm: float, bin_nm: float = 0.1) -> dict:
    """Overlap of the upper- and lower-leaflet tail z profiles: int min(rho_u, rho_l) dz / int 0.5 (rho_u + rho_l) dz."""
    # 0 = the two leaflets' tails never share a z slab, 1 = identical profiles (full interdigitation). The profiles
    # are linear number densities (atoms per nm) relative to the midplane; a bin width cancels in the ratio.
    up, lo = _xyz(tail_upper_nm)[:, 2] - midplane_nm, _xyz(tail_lower_nm)[:, 2] - midplane_nm
    both = np.concatenate([up, lo])
    if not len(up) or not len(lo):
        return {"overlap": float("nan"), "z_nm": np.zeros(0), "rho_upper": np.zeros(0), "rho_lower": np.zeros(0),
                "n_upper": int(len(up)), "n_lower": int(len(lo))}
    start = math.floor(both.min() / bin_nm) * bin_nm
    edges = np.arange(start, both.max() + bin_nm, bin_nm)
    rho_u = np.histogram(up, bins=edges)[0] / bin_nm
    rho_l = np.histogram(lo, bins=edges)[0] / bin_nm
    mean = 0.5 * (rho_u + rho_l)
    overlap = float(np.minimum(rho_u, rho_l).sum() / mean.sum()) if mean.sum() else float("nan")
    return {"overlap": overlap, "z_nm": 0.5 * (edges[1:] + edges[:-1]), "rho_upper": rho_u, "rho_lower": rho_l,
            "n_upper": int(len(up)), "n_lower": int(len(lo)), "bin_nm": bin_nm}


def interdigitation_change(reference_value: float | None, value: float, stage: str = "sample") -> dict:
    """Relative change (percent) of the tail overlap against a reference, graded 10 / 20 %; REFERENCE NEEDED without one."""
    if reference_value is None or not np.isfinite(reference_value) or not reference_value:
        record = metric("tail interdigitation", stage, value, "overlap", status="REFERENCE NEEDED",
                        reason="no matched reference overlap; the value is reported, not graded")
        return {"relative_percent": None, "status": "REFERENCE NEEDED", "record": record}
    percent = 100.0 * (value - reference_value) / reference_value
    status = classify(percent, 10.0, 20.0)
    record = metric("tail interdigitation", stage, value, "overlap", reference=reference_value, deviation=percent, status=status,
                    reason=f"{percent:+.1f} % vs reference (PASS <= 10 %, WARNING <= 20 %)")
    return {"relative_percent": percent, "status": status, "record": record}


def leaflet_asymmetry(stage: dict) -> dict:
    """Per-leaflet species counts and mole fractions, and the upper/lower count ratio (overall and per species)."""
    counts = {leaflet: Counter(r["resname"] for r in stage["lipids"] if r["leaflet"] == leaflet) for leaflet in ("lower", "upper")}
    totals = {k: sum(v.values()) for k, v in counts.items()}
    species = sorted(set(counts["lower"]) | set(counts["upper"]))
    ratio = lambda up, lo: up / lo if lo else float("inf") if up else float("nan")
    return {"n": totals, "ratio_upper_over_lower": ratio(totals["upper"], totals["lower"]),
            "species": {s: {"n_upper": counts["upper"][s], "n_lower": counts["lower"][s],
                            "x_upper": counts["upper"][s] / totals["upper"] if totals["upper"] else 0.0,
                            "x_lower": counts["lower"][s] / totals["lower"] if totals["lower"] else 0.0,
                            "ratio_upper_over_lower": ratio(counts["upper"][s], counts["lower"][s])} for s in species}}


def aa_leaflet_tail_split(residues: list[list[dict]], midplane_nm: float) -> tuple[np.ndarray, np.ndarray]:
    """Tail carbon coordinates of whole all-atom lipid residues (nm) split by the leaflet their ANCHOR sits in."""
    # Assignment by the headgroup anchor, never by the tail atom's own z: a tail reaching past the midplane is the
    # interdigitation being measured, not a lipid of the other leaflet.
    upper, lower = [], []
    for atoms in residues:
        heavy = [a for a in atoms if element(a["atom"]) != "H"]
        if not heavy or is_sterol(cg_name(heavy[0]["resname"])):
            continue
        xyz = np.array([[a["x"], a["y"], a["z"]] for a in heavy], dtype=float)
        names = [a["atom"] for a in heavy]
        z = float(anchor_xyz(xyz, names, anchor_bead(heavy[0]["resname"], atomistic=True))[2])
        tails = [xyz[i] for i, n in enumerate(names) if is_tail_carbon(n)]
        (upper if leaflet_of(z, midplane_nm) == "upper" else lower).extend(tails)
    return np.array(upper, dtype=float).reshape(-1, 3), np.array(lower, dtype=float).reshape(-1, 3)
