#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Membrane validation across the build: packing, composition and lateral order of every stage, and the slice gate.
#//=============================================================
"""Measure the membrane at every stage of a build and compare each with the embedded membrane the cut was taken from.

Stages (packing.cg_stage / packing.aa_stage records) are measured with the same tessellation and RDF code at both
resolutions: the coarse-grained frame, the embedded membrane (the reference: the lipids that were actually cut),
the slice, the backmapped membrane (membrane.pdb) and the minimized system (em.gro). The slice gate runs right after
slicing, before anything is backmapped:

    construction check   leaflet APL of the slice against the SAME lipids' areas in the uncut embedded membrane
                         (packing.region_reference): what the cut did to the packing of the region it sampled.
                         Fails the build (Settings.apl_validate) beyond Settings.apl_slice_tolerance_percent.
    representativeness   leaflet APL of the slice against the whole embedded cell, and the lipid count against the
                         count the reference APL predicts for the accessible area: how typical the sampled region is
                         of the parent membrane. Reported and warned about, never fatal: a finite crop around a
                         protein samples the membrane next to that protein, not the cell average.
    composition          per-leaflet species counts and mole fractions against the reference, judged against the
                         multinomial noise of a random crop of that many lipids.
    lateral order        all-anchor headgroup RDF per leaflet against the reference, within counting noise.
    integrity            every kept lipid complete, every anchor inside the cell, no duplicate, no hard-core overlap
                         left at the seam, and no bead closer to the complex than before the cut.

Outputs: membrane_validation.json (everything), membrane_validation.md (the tables), membrane_validation_lipids.tsv
(one row per lipid and stage with its Voronoi area) and, when matplotlib is importable, four PNG figures.
"""
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .audit import itp_atoms, topology_includes, topology_molecules
from .config import AMINO, SOLVENT, Settings
from .lipids import LIPID_AA_ANCHORS, display_name
from .packing import (
    LEAFLETS,
    aa_stage,
    cg_stage,
    compare_packing,
    composition_table,
    measure_packing,
    per_lipid_rows,
    reference_packing,
    region_reference,
    species_table,
    window_distribution,
)
from .qc import classify, classify_percentile, metric, overall_status, records_table
from .rdf import PROFILE_FIRST_SHELL_NM, leaflet_rdfs, protein_density_profiles, rdf_comparison, region_rdfs, window_rdf_distribution
from .runtools import log
from .structio import element, read_gro, read_pdb, residues_in_order, xyz_nm
from .structure_metrics import (
    aa_groups_from_atoms,
    aa_leaflet_tail_split,
    bilayer_thickness,
    core_hydration,
    interdigitation_change,
    orientation_change,
    protein_orientation,
    tail_interdigitation,
    thickness_change,
)

__all__ = ['STAGE_LABELS', 'MembraneValidation', 'stage_from_membrane_pdb', 'stage_from_gro', 'check_integrity',
           'stage_summary', 'write_outputs', 'render_figures']

STAGE_LABELS = {"cg_frame": "CG frame", "embed": "Embedded CG membrane (reference)", "slice": "Sliced CG membrane",
                "backmap": "Backmapped all-atom membrane", "minimize": "Minimized all-atom system",
                "equilibrated": "Equilibrated (trajectory)", "final": "Final production window"}
AA_LIPIDS = set(LIPID_AA_ANCHORS) | {"POPC", "DOPC", "POPE", "DOPE", "POPS", "DOPS", "PSM"}
# Figure colours: the first categorical slots of the validated reference palette, assigned by stage, never cycled.
# Hard guardrails of the gates (docs/membraneforger_tutorial.md 5.3 and 6): (PASS ceiling, WARNING ceiling).
CG_GATE = {"thickness_percent": (3.0, 5.0), "tilt_deg": (1.0, 3.0), "depth_A": (0.5, 1.5)}
AA_GATE = {"thickness_percent": (10.0, 15.0), "tilt_deg": (2.0, 5.0), "depth_A": (1.0, 2.5)}
CORE_HALF_NM = 0.8   # the hydrophobic core: within this of the midplane no water belongs (about the C4-C14 tail region)
STAGE_COLOURS = {"cg_frame": "#52514e", "embed": "#2a78d6", "slice": "#eb6834", "backmap": "#1baf7a", "minimize": "#eda100",
                 "old_rule": "#e34948"}
LEAFLET_HATCH = {"upper": None, "lower": "///"}


def stage_from_membrane_pdb(path: Path, name: str = "backmap") -> dict:
    """Analysis stage from membrane.pdb (A): the placed complex followed by the backmapped lipids (segid MEMB)."""
    atoms, cryst1 = read_pdb(path)
    cryst = (cryst1 or "").ljust(54)
    box = [float(cryst[6:15]) / 10.0, float(cryst[15:24]) / 10.0, float(cryst[24:33]) / 10.0]
    lipids = [a for a in atoms if a.get("segid") == "MEMB"]
    solute = [a for a in atoms if a.get("segid") != "MEMB" and element(a["atom"]) != "H"]
    residues = [[dict(a, x=a["x"] / 10.0, y=a["y"] / 10.0, z=a["z"] / 10.0) for a in res] for _, _, res in residues_in_order_by_chain(lipids)]
    stage = aa_stage(name, residues, xyz_nm(solute) / 10.0, box)
    stage["tails"] = aa_leaflet_tail_split(residues, stage["midplane_nm"])
    stage["water_nm"] = np.zeros((0, 3))  # membrane.pdb holds no water yet
    return stage


def residues_in_order_by_chain(atoms: list[dict]) -> list[tuple]:
    """Group lipid atoms into residues by chain and residue number, keeping file order."""
    seen, order = {}, []
    for a in atoms:
        key = (a.get("chain", ""), a["resid"], a["resname"])
        if key not in seen:
            seen[key] = []
            order.append(key)
        seen[key].append(a)
    return [(k[1], k[2], seen[k]) for k in order]


# How a lipid's residue name can appear in a .gro: the format keeps 5 characters (SAPI25 -> SAPI2), and the GM3 topology
# (GLPA.itp) splits each molecule into four residues, CER160 BGLC BGAL ANE5AC, written CER16 BGLC BGAL ANE5A.
GRO_RESNAME_ALIASES = {"SAPI2": "SAPI25", "SAPI": "SAPI25"}
GLPA_PARTS = ("CER16", "BGLC", "BGAL", "ANE5A")
# Water and ions as residue or moleculetype names (CHARMM TIP3/SOD/CLA/POT, GROMACS SOL/NA/CL/K).
SOLVENT_NAMES = frozenset(SOLVENT) | {"SOL", "HOH", "WAT", "NA", "CL", "K"}


def gro_molecules(atoms: list[dict], topology: Path) -> list[tuple[str, list[int]]] | None:
    """Split a .gro of the built system into molecules with the [ molecules ] table of its topology and the .itp files.

    Each molecule is (moleculetype, atom indices). This is exact whatever the .gro does to residue names (truncated to
    5 characters, GM3 split into CHARMM sugar residues). None when the topology does not describe these atoms.
    """
    try:
        sizes = {mol: len(entries) for mol, entries in itp_atoms(topology_includes(topology, set())).items()}
        table = topology_molecules(topology)
    except (OSError, SystemExit, ValueError):
        return None
    if any(mol not in sizes for mol, _ in table) or sum(sizes[mol] * n for mol, n in table) != len(atoms):
        return None
    molecules, start = [], 0
    for mol, n in table:
        for _ in range(n):
            molecules.append((mol, list(range(start, start + sizes[mol]))))
            start += sizes[mol]
    return molecules


def gro_molecules_by_name(atoms: list[dict]) -> list[tuple[str, list[int]]]:
    """Fallback without a topology: molecules from residue runs, with the .gro spellings mapped back to the lipid names.

    A residue is a run of equal (resid, resname) (numbers wrap at 99999); SAPI2 is SAPI25, and a GM3 is one run of its
    four parts CER16 BGLC BGAL ANE5A, a new molecule starting at each CER16.
    """
    molecules, last = [], None
    for i, a in enumerate(atoms):
        name = GRO_RESNAME_ALIASES.get(a["resname"], a["resname"])
        if name in GLPA_PARTS:
            if name == GLPA_PARTS[0] and last != (a["resid"], name) or not molecules or molecules[-1][0] != "GLPA":
                molecules.append(("GLPA", []))
            molecules[-1][1].append(i)
            last = (a["resid"], name)
            continue
        if (a["resid"], name) != last:
            molecules.append((name, []))
            last = (a["resid"], name)
        molecules[-1][1].append(i)
    return molecules


def stage_from_gro(path: Path, name: str = "minimize", topology: Path | None = None) -> dict:
    """Analysis stage from a GROMACS .gro of the built system (nm).

    Molecules come from the topology (`topology`, else topol.top next to the .gro) when it describes the file, else from
    residue names (gro_molecules_by_name). Lipids are the molecules whose type is an all-atom lipid; protein and ligands
    are only protein residues and the topology's other non-solvent molecules. Anything else is reported in
    stage["unrecognized"] and kept out of both, instead of being counted as protein.
    """
    atoms, box = read_gro(path)
    top = topology if topology is not None else path.parent / "topol.top"
    molecules = gro_molecules(atoms, top) if top.is_file() else None
    source = f"topology {top.name}" if molecules is not None else "residue names"
    if molecules is None:
        molecules = gro_molecules_by_name(atoms)
    residues, solute, others, unrecognized = [], [], [], Counter()
    for mol, members in molecules:
        if mol in AA_LIPIDS:
            residues.append([dict(atoms[i], resname=mol) for i in members])
        elif mol in SOLVENT_NAMES:
            others += [atoms[i] for i in members]
        elif source.startswith("topology") or all(atoms[i]["resname"] in AMINO for i in members):
            solute += [atoms[i] for i in members]
            others += [atoms[i] for i in members]
        else:
            unrecognized[mol] += 1
    solute = [a for a in solute if element(a["atom"]) != "H"]
    cell = np.array(box[:3])
    for res in residues:  # make each lipid whole: mdrun may write atoms in different periodic images
        xyz = xyz_nm(res)
        xyz = xyz[0] + (xyz - xyz[0]) - cell * np.round((xyz - xyz[0]) / cell)
        for a, (x, y, z) in zip(res, xyz):
            a["x"], a["y"], a["z"] = float(x), float(y), float(z)
    stage = aa_stage(name, residues, xyz_nm(solute), box[:3])
    stage["residues"] = residues
    stage["molecule_source"] = source
    stage["unrecognized"] = dict(unrecognized)
    stage["tails"] = aa_leaflet_tail_split(residues, stage["midplane_nm"])
    groups = aa_groups_from_atoms([a for res in residues for a in res] + others, AA_LIPIDS)
    stage["water_nm"] = groups["water"]
    stage["z_groups"] = groups
    return stage


def check_integrity(cut: dict, membrane: list[dict], placed: list[dict], box: list[float]) -> dict:
    """Molecular and periodic integrity of a slice: complete residues, anchors in the cell, no duplicates, no new clashes."""
    beads = {m["cg"]: len(m["beads"]) for m in membrane}
    complete = all(len(m["beads"]) == beads[m["cg"]] for m in cut["membrane"])
    duplicates = len(cut["kept_indices"]) - len(set(cut["kept_indices"]))
    size = np.array(cut["box"][:2])
    anchors = np.array([[r["x_nm"], r["y_nm"]] for r in cg_stage("check", cut["membrane"], cut["placed"], cut["box"], 0.0)["lipids"]])
    relaxed = cut["report"]["seam"]["relaxation"] or {}
    slack = relaxed.get("max_bead_displacement_nm", 0.0) + 1e-9  # the seam relaxation may nudge an edge anchor past the cell edge
    inside = bool(np.all(anchors >= -slack) and np.all(anchors < size + slack))
    heavy = lambda atoms: xyz_nm([a for a in atoms if element(a["atom"]) != "H"]) / 10.0
    before_pts = np.vstack([np.array([[b["x"], b["y"], b["z"]] for b in m["beads"]]) for m in membrane])
    after_pts = np.vstack([m["xyz"] for m in cut["membrane"]])
    cell_before, cell_after = np.array(box[:3]), np.array(cut["box"][:3])
    before = float(cKDTree(np.mod(heavy(placed), cell_before), boxsize=cell_before).query(np.mod(before_pts, cell_before))[0].min())
    after = float(cKDTree(np.mod(heavy(cut["placed"]), cell_after), boxsize=cell_after).query(np.mod(after_pts, cell_after))[0].min())
    seam = cut["report"]["seam"]
    hard_free = seam["relaxation"] is None or seam["relaxation"]["closest_created_pair_nm"] is None or \
        seam["relaxation"]["closest_created_pair_nm"] >= seam["hard_core_nm"]
    return {"complete_residues": complete, "duplicate_lipids": duplicates, "anchors_inside_cell": inside,
            "closest_bead_to_complex_before_nm": round(before, 4), "closest_bead_to_complex_after_nm": round(after, 4),
            "no_new_protein_lipid_clash": after >= min(before, 0.40) - 0.02, "no_hard_core_seam_overlap": bool(hard_free),
            "pass": bool(complete and not duplicates and inside and hard_free and after >= min(before, 0.40) - 0.02)}


def window_summary(windows: dict) -> dict:
    """Count and spread of the equal-size windows a leaflet was graded against."""
    if not windows["windows"]:
        return {"count": 0}
    return {"count": windows["windows"],
            "apl_A2_p2_5_p97_5": [round(float(v), 2) for v in np.percentile(windows["window_apl_A2"], [2.5, 97.5])],
            "composition_distance_p95": round(float(np.percentile(windows["window_composition_distance"], 95)), 4)}


def stage_summary(stage: dict, measure: dict, rdfs: dict, reference: dict | None = None) -> dict:
    """JSON-able summary of one measured stage, with the comparison against the reference stage when given."""
    summary = {"name": stage["name"], "label": STAGE_LABELS.get(stage["name"], stage["name"]), "resolution": stage["resolution"],
               "box_nm": [round(v, 4) for v in stage["box_nm"]], "midplane_nm": round(stage["midplane_nm"], 4), "leaflets": {}}
    for leaflet in LEAFLETS:
        m = measure[leaflet]
        summary["leaflets"][leaflet] = {"lipids": m["lipids"], "area_A2": round(m["area_nm2"] * 100, 1),
                                        "protein_area_A2": round(m["protein_area_nm2"] * 100, 1),
                                        "accessible_area_A2": round(m["accessible_area_nm2"] * 100, 1),
                                        "apl_A2": round(m["apl_A2"], 2), "head_apl_A2": round(m["head_apl_A2"], 2) if m["head_apl_A2"] else None,
                                        "composition": m["composition"],
                                        "rdf": {label: {"n": c["n_points"], "well_sampled": c["well_sampled"], **c["features"]}
                                                for label, c in rdfs.get(leaflet, {}).items()}}
    summary["species"] = species_table(stage, reference["species"] if reference else None)
    if reference is not None:
        comparison = compare_packing(None, None, reference["tolerance_percent"], reference["measure"], measure)
        summary["vs_reference"] = {"apl": {l: {k: v for k, v in comparison["leaflets"][l].items() if k in
                                               ("delta_apl_A2", "delta_apl_percent", "expected_lipids_from_reference_apl",
                                                "lipid_count_deviation_percent", "lipids")} for l in LEAFLETS},
                                   "composition": {l: {k: v for k, v in comparison["composition"][l].items() if k in
                                                       ("distance", "expected_distance_random_crop", "within_finite_crop_variability",
                                                        "flagged_species")} for l in LEAFLETS},
                                   "rdf": rdf_comparison(reference["rdfs"], rdfs, reference["settings"].rdf_max_peak_shift_a,
                                                         reference["settings"].rdf_max_noise_units, reference["settings"].rdf_peak_warning_a)}
    return summary


class MembraneValidation:
    """Collects the measured stages of one build and writes the validation outputs."""

    def __init__(self, out: Path, settings: Settings):
        """Start with no stages; the embedded membrane becomes the reference when it is added."""
        self.out, self.settings = out, settings
        self.stages, self.summaries, self.reference, self.slice_check, self.cut = {}, [], None, None, None

    def add_cg(self, name: str, membrane: list[dict], placed: list[dict], box: list[float], midplane_nm: float,
               protein_nm: np.ndarray | None = None) -> dict:
        """Measure a coarse-grained stage (optionally with CG protein beads as the protein instead of the complex)."""
        stage = cg_stage(name, membrane, placed, box, midplane_nm)
        if protein_nm is not None:
            stage["protein_nm"] = np.asarray(protein_nm, dtype=float).reshape(-1, 3)
        return self.add(stage)

    def add(self, stage: dict) -> dict:
        """Measure any stage: tessellation, RDFs, thickness, protein orientation, all-atom core metrics; 'embed' sets the reference."""
        measure = measure_packing(stage)
        rdfs = leaflet_rdfs(stage, accessible_nm2={l: measure[l]["accessible_area_nm2"] for l in LEAFLETS})
        thickness = bilayer_thickness(stage)
        planes = (thickness["leaflet_planes_nm"]["lower"], thickness["leaflet_planes_nm"]["upper"])
        orientation = protein_orientation(stage["protein_nm"], stage["midplane_nm"], headgroup_planes_nm=planes,
                                          reference=self.reference["orientation"] if self.reference else None) if len(stage["protein_nm"]) else None
        if stage["name"] == "embed":
            self.reference = {"stage": stage, "measure": measure, "rdfs": rdfs, "species": species_table(stage),
                              "tolerance_percent": self.settings.apl_slice_tolerance_percent, "settings": self.settings,
                              "packing": reference_packing(stage), "thickness": thickness, "orientation": orientation}
        summary = stage_summary(stage, measure, rdfs, self.reference if stage["name"] != "embed" else None)
        summary["thickness"] = {k: thickness[k] for k in ("thickness_A", "local_sd_nm", "cells_filled", "n_per_leaflet")}
        summary["protein"] = {k: orientation[k] for k in ("tilt_deg", "depth_A", "com_z_A", "n_slab", "embedded_fraction")} if orientation else None
        density = protein_density_profiles(stage) if len(stage["protein_nm"]) else {}
        summary["protein_density"] = density_summary(density)
        records = self.structure_records(stage, thickness, orientation)
        if "molecule_source" in stage:  # an all-atom .gro stage: say how its molecules were identified, and what was not
            unknown = stage.get("unrecognized") or {}
            records.append(metric("molecules identified", stage["name"], len(stage["lipids"]), "lipids", n=len(stage["lipids"]),
                                  status="WARNING" if unknown else "PASS",
                                  reason=f"from {stage['molecule_source']}"
                                         + ("; not lipid, protein or solvent, left out of every measurement: "
                                            + ", ".join(f"{k} x{v}" for k, v in sorted(unknown.items())) if unknown else "")))
            summary["molecule_source"], summary["unrecognized"] = stage["molecule_source"], dict(unknown)
        summary["records"] = records
        self.stages[stage["name"]] = {"stage": stage, "measure": measure, "rdfs": rdfs, "thickness": thickness, "orientation": orientation,
                                      "protein_density": density,
                                      "records": records, **{k: stage[k] for k in ("interdigitation", "core_hydration") if k in stage}}
        self.summaries = [s for s in self.summaries if s["name"] != stage["name"]] + [summary]
        return summary

    def structure_records(self, stage: dict, thickness: dict, orientation: dict | None) -> list[dict]:
        """Thickness, protein orientation and (all-atom) interdigitation and core hydration records of one stage."""
        name = stage["name"]
        if self.reference is None or name == "embed" or name == "cg_frame":
            return [metric("bilayer thickness", name, thickness["thickness_A"], "A", n=sum(thickness["n_per_leaflet"].values()),
                           status="NOT RUN", reason="reference stage" if name == "embed" else "before embedding; reported")]
        gate = CG_GATE if stage["resolution"] == "cg" else AA_GATE
        records = thickness_change(self.reference["thickness"], thickness, *gate["thickness_percent"], stage=name)["records"]
        if orientation is not None and self.reference["orientation"] is not None:
            records += orientation_change(self.reference["orientation"], orientation, *gate["tilt_deg"], *gate["depth_A"], stage=name)["records"]
        if stage["resolution"] == "aa" and "tails" in stage:
            inter = tail_interdigitation(stage["tails"][0], stage["tails"][1], stage["midplane_nm"])
            first = next((self.stages[n] for n in ("backmap",) if n in self.stages and n != name), None)
            reference_value = first["interdigitation"]["overlap"] if first and "interdigitation" in first else None
            stage["interdigitation"] = inter
            change = interdigitation_change(reference_value, inter["overlap"], stage=name)
            change["record"]["reason"] += "" if reference_value is not None else " (first all-atom stage: no matched reference yet)"
            records.append(change["record"])
            if len(stage.get("water_nm", ())):
                hydration = core_hydration(stage["water_nm"], stage["protein_nm"], stage["midplane_nm"], CORE_HALF_NM, stage["box_nm"],
                                           thickness_nm=thickness["thickness_nm"], stage=name)
                records += hydration["records"]
                stage["core_hydration"] = {k: hydration[k] for k in ("n_core", "n_core_free", "n_protein_associated", "ratio_percent",
                                                                     "transmembrane_water_path", "status")}
        return records

    def slice_reference(self) -> dict | None:
        """The compact per-leaflet reference the slicer uses to place the window (None before the embed stage)."""
        return self.reference["packing"] if self.reference else None

    def validate_slice(self, cut: dict, membrane: list[dict], placed: list[dict], box: list[float]) -> dict:
        """Run the slice gate on a cut (after add('slice')); returns the record and raises SystemExit when it fails."""
        ref, cur, settings = self.reference, self.stages["slice"], self.settings
        summary = next(s for s in self.summaries if s["name"] == "slice")
        region = region_reference(ref["stage"], cut["kept_indices"])
        cropped = [w["cropped"] for w in cut["report"]["axes"]]
        windows = window_distribution(ref["stage"], cut["box"][:2], cropped, settings.slice_window_samples)
        warning, tolerance = settings.apl_slice_warning_percent, settings.apl_slice_tolerance_percent
        leaflets, records = {}, []
        for leaflet in LEAFLETS:
            apl = cur["measure"][leaflet]["apl_A2"]
            change = 100.0 * (apl - region[leaflet]["apl_A2"]) / region[leaflet]["apl_A2"] if region[leaflet]["apl_A2"] else float("nan")
            whole = summary["vs_reference"]["apl"][leaflet]
            status = classify(change, warning, tolerance)
            repr_status, repr_percentile = classify_percentile(apl, windows[leaflet]["window_apl_A2"])
            if repr_status == "INSUFFICIENT SAMPLING":
                repr_status = "NOT RUN"  # too few equal-size windows (an uncut axis): representativeness is reported, not graded
            comp = summary["vs_reference"]["composition"][leaflet]
            comp_status, comp_percentile = classify_percentile(comp["distance"], windows[leaflet]["window_composition_distance"], one_sided=True)
            comp_status = "PASS" if comp_status == "PASS" or comp["within_finite_crop_variability"] and comp_status != "FAIL" else comp_status
            n = cur["measure"][leaflet]["lipids"]
            leaflets[leaflet] = {"lipids": n, "region_lipids_in_reference": region[leaflet]["lipids"],
                                 "apl_A2": round(apl, 2), "region_reference_apl_A2": region[leaflet]["apl_A2"],
                                 "whole_cell_reference_apl_A2": region[leaflet]["whole_cell_apl_A2"],
                                 "construction_delta_percent": round(change, 2), "construction_status": status,
                                 "construction_pass": status != "FAIL",
                                 "representativeness_delta_percent": whole["delta_apl_percent"], "representativeness_status": repr_status,
                                 "representativeness_percentile": repr_percentile,
                                 "region_vs_whole_cell_percent": region[leaflet].get("region_vs_whole_cell_percent"),
                                 "expected_lipids_from_reference_apl": whole["expected_lipids_from_reference_apl"],
                                 "lipid_count_deviation_percent": whole["lipid_count_deviation_percent"],
                                 "composition_status": comp_status, "composition_percentile": comp_percentile,
                                 "equivalent_windows": window_summary(windows[leaflet])}
            records.append(metric("leaflet APL, same lipids uncut", "slice", apl, "A^2", leaflet, region[leaflet]["apl_A2"], None, n,
                                  round(change, 2), status, f"{change:+.2f} % (PASS <= {warning} %, WARNING <= {tolerance} %)"))
            records.append(metric("leaflet APL vs whole embedded cell", "slice", apl, "A^2", leaflet, region[leaflet]["whole_cell_apl_A2"],
                                  None, n, whole["delta_apl_percent"], repr_status if repr_status != "FAIL" else "WARNING",
                                  f"{whole['delta_apl_percent']:+.2f} %; percentile {repr_percentile} of {windows[leaflet]['windows']} "
                                  "equal-size windows (representativeness, never fatal"
                                  + ("; fewer than 20 windows, not graded)" if repr_status == "NOT RUN" else ")")))
            records.append(metric("composition distance", "slice", comp["distance"], "mole fraction", leaflet, 0.0, None, n, comp["distance"],
                                  comp_status, f"random crop of this size: {comp['expected_distance_random_crop']}; percentile "
                                  f"{comp_percentile} of equal-size windows"
                                  + (f"; flagged {comp['flagged_species']}" if comp["flagged_species"] else "")))
        region_curves = region_rdfs(ref["stage"], cut["kept_indices"], min(0.5 * min(cut["box"][:2]) - 1e-6, 2.5),
                                    {l: cur["measure"][l]["accessible_area_nm2"] for l in LEAFLETS})
        rdf = rdf_comparison(region_curves, cur["rdfs"], settings.rdf_max_peak_shift_a, settings.rdf_max_noise_units, settings.rdf_peak_warning_a)
        rdf["reference"] = "the same lipids in the uncut embedded membrane (construction reference)"
        rdf["whole_cell"] = summary["vs_reference"]["rdf"]
        r_max = min(0.5 * min(cut["box"][:2]) - 1e-6, 2.5)
        rdf_windows = window_rdf_distribution(ref["stage"], cut["box"][:2], cropped, r_max, settings.slice_window_samples)
        for leaflet in LEAFLETS:
            rdf_all = rdf["leaflets"][leaflet].get("all")
            if rdf_all:
                n = cur["measure"][leaflet]["lipids"]
                records.append(metric("headgroup RDF first-shell peak shift", "slice", rdf_all["first_peak_shift_A"], "A", leaflet, 0.0, None, n,
                                      rdf_all["first_peak_shift_A"], rdf_all["status"] if settings.rdf_validate else "NOT RUN",
                                      f"vs the same lipids uncut: RMS difference {rdf_all['rms_difference']} = {rdf_all['rms_in_noise_units']} x "
                                      f"counting noise" + (f"; {rdf_all['note']}" if rdf_all.get("note") else "")))
            whole_all = rdf["whole_cell"]["leaflets"][leaflet].get("all")
            if whole_all and rdf_windows[leaflet]["windows"]:
                status, percentile = classify_percentile(whole_all["rms_difference"], rdf_windows[leaflet]["rms"], one_sided=True)
                status = {"INSUFFICIENT SAMPLING": "NOT RUN"}.get(status, status)
                rdf["whole_cell"]["leaflets"][leaflet]["all"]["window_percentile"] = percentile
                rdf["whole_cell"]["leaflets"][leaflet]["all"]["window_status"] = status
                records.append(metric("headgroup RDF deviation vs parent windows", "slice", whole_all["rms_difference"], "g(r) RMS", leaflet,
                                      round(float(np.median(rdf_windows[leaflet]["rms"])), 4), None, rdf_windows[leaflet]["windows"], percentile,
                                      status if settings.rdf_validate else "NOT RUN",
                                      f"RMS deviation of the slice from the whole-cell g(r) at percentile {percentile} of "
                                      f"{rdf_windows[leaflet]['windows']} equal-size windows (one-sided: PASS below the 95th percentile, "
                                      "WARNING below the 99th, FAIL above; a smaller deviation than the windows' is never graded down)"))
        integrity = check_integrity(cut, membrane, placed, box)
        records.append(metric("slice integrity", "slice", int(integrity["pass"]), "bool", None, 1, None, len(cut["membrane"]), None,
                              "PASS" if integrity["pass"] else "FAIL",
                              f"complete residues {integrity['complete_residues']}, duplicates {integrity['duplicate_lipids']}, anchors inside "
                              f"{integrity['anchors_inside_cell']}, hard-core seam overlaps left {not integrity['no_hard_core_seam_overlap']}"))
        records.append(metric("closest lipid bead to the complex", "slice", integrity["closest_bead_to_complex_after_nm"], "nm", None,
                              integrity["closest_bead_to_complex_before_nm"], None, None,
                              round(integrity["closest_bead_to_complex_after_nm"] - integrity["closest_bead_to_complex_before_nm"], 4),
                              "PASS" if integrity["no_new_protein_lipid_clash"] else "FAIL",
                              "no lipid bead moved closer to the complex than before the cut"))
        composition = summary["vs_reference"]["composition"]
        verdict = {"integrity": integrity["pass"], "upper_apl": leaflets["upper"]["construction_pass"],
                   "lower_apl": leaflets["lower"]["construction_pass"],
                   "composition": all(leaflets[l]["composition_status"] != "FAIL" for l in LEAFLETS),
                   "rdf": (rdf["pass"] and all(rdf["whole_cell"]["leaflets"][l].get("all", {}).get("window_status", "PASS") != "FAIL"
                                               for l in LEAFLETS)) or not settings.rdf_validate,
                   "no_new_protein_lipid_clash": integrity["no_new_protein_lipid_clash"]}
        overall = overall_status(records)
        self.slice_check = {"reference": "embedded CG membrane (after embed_complex and lipid edits), the lipids that were cut",
                            "warning_percent": warning, "tolerance_percent": tolerance, "leaflets": leaflets, "integrity": integrity,
                            "composition": composition, "rdf": rdf, "records": records, "verdict": verdict, "pass": all(verdict.values()),
                            "status": overall["status"], "overall": overall,
                            "apl_validate": settings.apl_validate, "rdf_validate": settings.rdf_validate}
        self.log_slice_check(leaflets, composition, rdf, integrity)
        if settings.apl_validate and not self.slice_check["pass"]:
            failed = [l for l in LEAFLETS if not leaflets[l]["construction_pass"]]
            raise SystemExit("the sliced membrane does not preserve the embedded membrane: "
                             + "; ".join(f"{l} leaflet APL {leaflets[l]['region_reference_apl_A2']} -> {leaflets[l]['apl_A2']} A^2 "
                                         f"({leaflets[l]['construction_delta_percent']:+.1f} %, tolerance {tolerance} %)" for l in failed)
                             + "; ".join([""] + [f"{k} FAIL" for k, v in verdict.items() if not v and k not in ("upper_apl", "lower_apl")])
                             + ". See membrane_validation.md; --apl-tolerance raises the limit, --apl-validate no only reports it")
        return self.slice_check

    def log_slice_check(self, leaflets: dict, composition: dict, rdf: dict, integrity: dict) -> None:
        """Write the slice gate to the build log, one line per leaflet plus composition, RDF and integrity."""
        for l in LEAFLETS:
            v = leaflets[l]
            log(self.out, f"slice {l} leaflet: {v['lipids']} lipids, APL {v['apl_A2']} A^2 vs {v['region_reference_apl_A2']} A^2 for the same "
                          f"lipids in the embedded membrane ({v['construction_delta_percent']:+.2f} %: {v['construction_status']}, "
                          f"PASS <= {self.settings.apl_slice_warning_percent} %, FAIL > {self.settings.apl_slice_tolerance_percent} %); whole "
                          f"embedded cell {v['whole_cell_reference_apl_A2']} A^2 ({v['representativeness_delta_percent']:+.2f} %, "
                          f"{v['expected_lipids_from_reference_apl']} lipids expected; {v['representativeness_status']} against "
                          f"{v['equivalent_windows']['count']} equal-size windows)",
                "PASS" if v["construction_status"] == "PASS" else "WARN")
            c = composition[l]
            log(self.out, f"slice {l} leaflet composition distance {c['distance']} (random crop of this size: {c['expected_distance_random_crop']}; "
                          f"{v['composition_status']} against equal-size windows)"
                          + (f"; flagged species {c['flagged_species']}" if c["flagged_species"] else ""),
                "INFO" if v["composition_status"] == "PASS" else "WARN")
            r = rdf["leaflets"][l].get("all")
            if r and self.settings.rdf_validate:
                log(self.out, f"slice {l} leaflet headgroup RDF vs the same lipids uncut: first-shell peak shift {r['first_peak_shift_A']} A, "
                              f"RMS difference {r['rms_difference']} = {r['rms_in_noise_units']} x counting noise: {r['status']}",
                    "INFO" if r["pass"] else "WARN")
        log(self.out, f"slice integrity: complete residues {integrity['complete_residues']}, duplicates {integrity['duplicate_lipids']}, anchors "
                      f"inside {integrity['anchors_inside_cell']}, closest bead-complex {integrity['closest_bead_to_complex_before_nm']} -> "
                      f"{integrity['closest_bead_to_complex_after_nm']} nm, hard-core seam overlaps left: "
                      f"{not integrity['no_hard_core_seam_overlap']}",
            "PASS" if integrity["pass"] else "WARN")

    def old_rule_row(self, cut: dict) -> dict | None:
        """What the all-beads-inside rule would have given in the same window: lipid counts and the APL they imply."""
        if not self.reference or "slice" not in self.stages:
            return None
        selection = cut["report"]["selection"]
        measure = self.stages["slice"]["measure"]
        row = {"name": "old_rule", "label": "Old all-beads-inside slice (same window, estimated)", "leaflets": {}}
        for l in LEAFLETS:
            n_old = selection["kept_by_leaflet"]["all_beads_inside_rule"].get(l, 0)
            accessible = measure[l]["accessible_area_nm2"] * 100
            apl = accessible / n_old if n_old else None
            ref = self.reference["measure"][l]["apl_A2"]
            row["leaflets"][l] = {"lipids": n_old, "apl_A2": round(apl, 2) if apl else None,
                                  "delta_apl_percent": round(100 * (apl - ref) / ref, 2) if apl else None,
                                  "lost_with_anchor_inside": selection["rejected_by_old_rule_with_anchor_inside_by_leaflet"].get(l, 0)}
        return row

    def all_records(self) -> list[dict]:
        """Every metric record of the build so far, in stage order (the slice gate's first)."""
        records = list(self.slice_check["records"]) if self.slice_check else []
        for name in record_stage_order(self.stages):
            records += [r for r in self.stages[name].get("records", []) if name != "slice" or r["metric"] not in {q["metric"] for q in records}]
        return records

    def record(self, cut: dict | None = None) -> dict:
        """Everything for run_manifest.json and membrane_validation.json."""
        records = self.all_records()
        graded = [r for r in records if r["stage"] in ("slice", "backmap", "minimize")]
        return {"reference_stage": "embed", "stages": self.summaries, "slice_check": self.slice_check, "records": records,
                "stage_matrix": stage_matrix(self.stages), "overall": overall_status(graded) if graded else None,
                "gates": {"1 CG slicing fidelity": overall_status([r for r in records if r["stage"] == "slice"]) if self.slice_check else None,
                          "2 backmapping integrity": overall_status([r for r in records if r["stage"] in ("backmap", "minimize")])
                          if any(n in self.stages for n in ("backmap", "minimize")) else None,
                          "3 equilibration convergence": None, "4 equilibrium physical properties": None},
                "old_rule_estimate": self.old_rule_row(cut) if cut else None,
                "slice_selection": cut["report"]["selection"] if cut else None, "crop_offset": cut["report"]["crop_offset"] if cut else None,
                "seam": cut["report"]["seam"] if cut else None}

    def write(self, cut: dict | None = None) -> list[Path]:
        """Write membrane_validation.json/.md, the per-lipid TSV and the figures; returns the files written."""
        return write_outputs(self.out, self.record(cut), self.stages)


# ----------------------------------------------------------------------------------------------------------- outputs

def fmt(value, digits=1) -> str:
    """Table cell: a number to `digits` decimals, or a dash."""
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "-"
    return f"{value:.{digits}f}" if isinstance(value, float) else str(value)


def signed(value, digits=1, unit="") -> str:
    """Table cell for a difference: signed number with a unit, or a dash."""
    return "-" if value is None else f"{value:+.{digits}f}{unit}"


def stage_table(record: dict) -> str:
    """The stage comparison table in Markdown."""
    head = ("| Stage / algorithm | Upper lipids | Lower lipids | Upper APL (Å²) | Lower APL (Å²) | ΔAPL upper vs embed | "
            "ΔAPL lower vs embed | Composition deviation (upper / lower) | RDF deviation (upper / lower, noise units) | Pass/Fail |\n"
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|\n")
    rows = []
    check = record.get("slice_check") or {}
    for s in record["stages"]:
        up, lo = s["leaflets"]["upper"], s["leaflets"]["lower"]
        if s["name"] == "embed":
            delta = ("reference", "reference", "reference", "reference", "reference")
        elif "vs_reference" in s:
            v = s["vs_reference"]
            rdf_source = check["rdf"] if (s["name"] == "slice" and check) else v["rdf"]
            rdf_u, rdf_l = rdf_source["leaflets"]["upper"].get("all", {}), rdf_source["leaflets"]["lower"].get("all", {})
            delta = (signed(v["apl"]["upper"]["delta_apl_percent"], 1, " %"), signed(v["apl"]["lower"]["delta_apl_percent"], 1, " %"),
                     f"{fmt(v['composition']['upper']['distance'], 3)} / {fmt(v['composition']['lower']['distance'], 3)}",
                     f"{fmt(rdf_u.get('rms_in_noise_units'))} / {fmt(rdf_l.get('rms_in_noise_units'))}",
                     ("PASS" if check.get("pass") else "FAIL") if s["name"] == "slice" and check else
                     ("cross-resolution, reported" if s["resolution"] == "aa" else "-"))
        else:
            delta = ("-", "-", "-", "-", "-")
        rows.append(f"| {s['label']} | {up['lipids']} | {lo['lipids']} | {fmt(up['apl_A2'])} | {fmt(lo['apl_A2'])} | " + " | ".join(delta) + " |")
        if s["name"] == "slice" and record.get("old_rule_estimate"):
            o = record["old_rule_estimate"]
            rows.append(f"| {o['label']} | {o['leaflets']['upper']['lipids']} | {o['leaflets']['lower']['lipids']} | "
                        f"{fmt(o['leaflets']['upper']['apl_A2'])} | {fmt(o['leaflets']['lower']['apl_A2'])} | "
                        f"{signed(o['leaflets']['upper']['delta_apl_percent'], 1, ' %')} | "
                        f"{signed(o['leaflets']['lower']['delta_apl_percent'], 1, ' %')} | - | - | FAIL (would have) |")
    return head + "\n".join(rows) + "\n"


def species_markdown(record: dict) -> str:
    """The per-stage, per-leaflet, per-species table in Markdown."""
    head = ("| Stage | Leaflet | Lipid | N | Mean Voronoi APL (Å²) | Median (Å²) | SD (Å²) | Head-plane APL (Å²) | Mole fraction | "
            "ΔAPL vs embed | Δ composition (mol %) |\n|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
    rows = []
    for s in record["stages"]:
        for r in s["species"]:
            rows.append(f"| {s['label']} | {r['leaflet']} | {display_name(r['lipid'])} | {r['n']}{'' if r['well_sampled'] else ' (few)'} | "
                        f"{fmt(r['apl_mean_A2'])} | {fmt(r['apl_median_A2'])} | {fmt(r['apl_sd_A2'])} | {fmt(r['head_apl_mean_A2'])} | "
                        f"{fmt(r['mole_fraction'], 3)} | {signed(r.get('apl_delta_vs_reference_percent'), 1, ' %')} | "
                        f"{signed(r.get('mole_fraction_delta_points'), 1)} |")
    return head + "\n".join(rows) + "\n"


def algorithm_table(record: dict) -> str:
    """Old slicer against new slicer, from the measured values of this build."""
    check, old, sel = record.get("slice_check") or {}, record.get("old_rule_estimate"), record.get("slice_selection") or {}
    seam = record.get("seam") or {}
    if not (check and old):
        return "(the slice was not validated in this build)\n"
    L = check["leaflets"]
    comp = check["composition"]
    rdf = check["rdf"]["leaflets"]
    rows = [("Selection criterion", "every bead of the lipid inside the window", "headgroup anchor inside the half-open window"),
            ("Retains the whole lipid", "only when no bead crosses the edge", "yes, always"),
            ("PBC-aware", "periodic image nearest the window centre, then a hard cut", "lipid made whole, anchor imaged into the window, "
                                                                                     "one image only (half-open), seam relaxed in place"),
            ("Lipids kept (upper / lower)", f"{old['leaflets']['upper']['lipids']} / {old['leaflets']['lower']['lipids']}",
             f"{L['upper']['lipids']} / {L['lower']['lipids']}"),
            ("Boundary lipids lost with the anchor inside", f"{sel.get('rejected_by_old_rule_with_anchor_inside', '-')}",
             f"{seam.get('lipids_removed', 0)} (hard-core seam overlaps only)"),
            ("Upper ΔAPL vs embedded cell", signed(old["leaflets"]["upper"]["delta_apl_percent"], 1, " %"),
             signed(L["upper"]["representativeness_delta_percent"], 1, " %")),
            ("Lower ΔAPL vs embedded cell", signed(old["leaflets"]["lower"]["delta_apl_percent"], 1, " %"),
             signed(L["lower"]["representativeness_delta_percent"], 1, " %")),
            ("Upper ΔAPL vs the same lipids uncut", "-", signed(L["upper"]["construction_delta_percent"], 1, " %")),
            ("Lower ΔAPL vs the same lipids uncut", "-", signed(L["lower"]["construction_delta_percent"], 1, " %")),
            ("Composition distance (upper / lower)", "-", f"{comp['upper']['distance']} / {comp['lower']['distance']}"),
            ("CG headgroup RDF, RMS difference in noise units (upper / lower)", "-",
             f"{fmt(rdf['upper'].get('all', {}).get('rms_in_noise_units'))} / {fmt(rdf['lower'].get('all', {}).get('rms_in_noise_units'))}"),
            ("Preserves the parent lipid density", "no (the lipid count drops, the area does not)", "PASS" if check["pass"] else "FAIL")]
    return "| Property | Old slicer | New slicer |\n|---|---|---|\n" + "\n".join(f"| {a} | {b} | {c} |" for a, b, c in rows) + "\n"


def write_markdown(record: dict) -> str:
    """membrane_validation.md: the three tables and the slice verdict."""
    check = record.get("slice_check")
    overall = record.get("overall")
    text = ["# Membrane validation\n",
            f"Overall status: **{overall['status'] if overall else 'INSUFFICIENT SAMPLING'}**"
            + ("" if overall else " (no graded stage yet)") + "\n",
            "## Stage matrix\n", matrix_markdown(record.get("stage_matrix", [])),
            "Legend: ✓ ran in this build, o runs at that stage when the data exist (equilibration and production need "
            "`python -m membraneforger.equilibration`), blank not applicable.\n",
            "All areas are Å² per lipid from a periodic 2D Voronoi tessellation of one headgroup anchor per lipid, with the heavy atoms of "
            "the complex in the leaflet's headgroup-to-midplane slab as competing generators whose cells are discarded. The leaflet APL is "
            "the mean cell area (= (A_xy - A_protein) / N). Species APL is the mean cell area of that species, never total area over its "
            "count. 'Head-plane APL' tessellates the phosphate/amide anchors alone (sterols left out). N is shown for every species; "
            "values from fewer than 10 molecules are marked '(few)'.\n",
            "## 1. Construction validation (gate 1: CG slicing fidelity)\n", "### Stages\n", stage_table(record)]
    if check:
        text += ["### Slice gate\n",
                 f"Reference: {check['reference']}. Construction check: leaflet APL of the slice against the Voronoi areas of the same "
                 f"lipids in the uncut embedded membrane, tolerance {check['tolerance_percent']} %. Representativeness: the slice against "
                 "the whole embedded cell (reported, not a failure criterion).\n",
                 "| Leaflet | Lipids | Same lipids uncut (Å²) | Slice (Å²) | Construction Δ | Status | Whole embedded cell (Å²) | "
                 "Representativeness Δ | Percentile among equal windows | Lipids expected from reference APL | Composition |\n"
                 "|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---|"]
        for l in LEAFLETS:
            v = check["leaflets"][l]
            text.append(f"| {l} | {v['lipids']} | {fmt(v['region_reference_apl_A2'])} | {fmt(v['apl_A2'])} | "
                        f"{signed(v['construction_delta_percent'], 2, ' %')} | {v['construction_status']} | "
                        f"{fmt(v['whole_cell_reference_apl_A2'])} | "
                        f"{signed(v['representativeness_delta_percent'], 2, ' %')} | {fmt(v['representativeness_percentile'], 3)} | "
                        f"{fmt(v['expected_lipids_from_reference_apl'])} | {v['composition_status']} |")
        verdict = check["verdict"]
        text += ["", "Verdict: " + ", ".join(f"{k} {'PASS' if v else 'FAIL'}" for k, v in verdict.items())
                 + f" -> **{check['status']}**\n", "### Metric records\n",
                 records_table([r for r in record.get("records", check["records"]) if r["stage"] in ("slice", "cg_frame", "embed")])]
        text += ["",
                 "### Old versus new slicing algorithm\n", algorithm_table(record)]
        if record.get("crop_offset"):
            co = record["crop_offset"]
            text += [f"Crop position: offset {co['offset_nm']} nm chosen from {co['candidates']} candidate(s) (score {co['chosen'].get('score')}; "
                     f"centred window score {co['unshifted'].get('score')}).\n"]
        if record.get("seam") and record["seam"].get("relaxation"):
            r = record["seam"]["relaxation"]
            text += [f"Seam: {record['seam']['pairs_already_close_in_source_frame']} bead pairs were already within "
                     f"{record['seam']['threshold_nm']} nm in the source frame (left alone); the relaxation moved {r['lipids_moved']} lipids "
                     f"by at most {r['max_bead_displacement_nm']} nm (mean {r['mean_bead_displacement_nm']} nm) in {r['steps']} steps "
                     f"({'converged' if r['converged'] else 'not converged'}; closest remaining created pair {r['closest_created_pair_nm']} nm); "
                     f"{record['seam']['lipids_removed']} lipid(s) removed for hard-core overlaps.\n"]
    aa = [r for r in record.get("records", []) if r["stage"] in ("backmap", "minimize")]
    text += ["## 2. Backmapping validation (gate 2: backmapping integrity)\n",
             (records_table(aa) if aa else "No all-atom stage has been measured yet (membrane.pdb and em.gro are measured as they appear).\n"),
             "Ring threading and lipid-solute clashes are checked by the build itself (ring_piercing.json, the lipid scan, the audit); "
             "cross-resolution APL and RDF are reported in the stage table, not graded.\n",
             "## 3. Equilibration convergence (gate 3)\n",
             "Not part of the build: run `python -m membraneforger.equilibration` on the equilibration trajectory (box series, frames) to "
             "grade area/volume/APL/thickness convergence, core hydration, protein tilt and depth, and tail order; the results are "
             "appended here.\n",
             "## 4. Equilibrium physical validation (gate 4)\n",
             "Lateral diffusion, area compressibility K_A and leaflet tension need a restraint-free production window; same tool, "
             "same report. Without matched references these are REFERENCE NEEDED, without enough sampling INSUFFICIENT SAMPLING.\n",
             "## Species\n", species_markdown(record),
             "## Lipids around the protein\n", density_markdown(record),
             "## Reading the tables\n",
             "- The construction check is the slicing test: the cut keeps the same lipids with the same neighbours except at the new seam, so "
             "its APL must match theirs in the uncut membrane within a few percent.\n"
             "- Representativeness tells how typical the sampled region is of the whole embedded cell; a crop around a protein samples the "
             "membrane next to that protein.\n"
             "- Across resolutions the anchors change (Martini PO4/ROH to atomistic P/O3), so the all-atom rows are shown for continuity: "
             "backmapping keeps the lipid count and the cell, hence the leaflet APL; minimization relaxes short contacts and may move the "
             "RDF's first peak toward the atomistic contact distance without changing the leaflet mean.\n"]
    return "\n".join(text)


# Species colours of the density-profile figure: one fixed categorical slot per species, so a species keeps its colour in
# every stage. Ten species exceed the eight slots, but no leaflet holds more than eight: the species found in both
# leaflets take slots 1-5, the upper-only (PSM, GM3) and lower-only (POPS, DOPS, SAP6) ones share slots 6-8.
DENSITY_FIGURE_FROM_A = 6.0  # the cumulative curves start here: nearer, the lipid area within R is a fraction of one lipid's
CATEGORICAL = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
SPECIES_SLOT = {"CHOL": 0, "POPC": 1, "DOPC": 2, "POPE": 3, "DOPE": 4, "PSM": 5, "GM3": 6, "POPS": 5, "DOPS": 6, "SAP6": 7}


def density_summary(density: dict) -> dict:
    """The JSON form of protein_density_profiles: per leaflet and species, N, first-shell N and ratio, and the curve."""
    out = {}
    for leaflet, block in density.items():
        if not block:
            out[leaflet] = None
            continue
        out[leaflet] = {"protein_atoms_in_slab": block["protein_atoms_in_slab"], "first_shell_nm": PROFILE_FIRST_SHELL_NM,
                        "species": {name: {"n": c["n"], "first_shell_n": c["first_shell_n"], "first_shell_ratio": c["first_shell_ratio"],
                                           "well_sampled": c["well_sampled"], "r_A": [round(10 * float(v), 2) for v in c["r_nm"]],
                                           "density_ratio": [None if not np.isfinite(v) else round(float(v), 3) for v in c["density_ratio"]],
                                           "R_A": [round(10 * float(v), 2) for v in c["R_nm"]],
                                           "cumulative_ratio": [None if not np.isfinite(v) else round(float(v), 3)
                                                                for v in c["cumulative_ratio"]]}
                                    for name, c in block["curves"].items()}}
    return out


def density_stages(record: dict) -> list[dict]:
    """Stages with a density profile around the complex (the CG frame's profile is around the frame's own receptor)."""
    return [s for s in record["stages"] if s.get("name") != "cg_frame" and any((s.get("protein_density") or {}).values())]


def density_markdown(record: dict) -> str:
    """First-shell enrichment of each species around the protein, per leaflet and stage (the table view of the figure)."""
    stages = density_stages(record)
    if not stages:
        return "No stage has protein atoms inside the bilayer.\n"
    text = [f"Density of each species within {10 * PROFILE_FIRST_SHELL_NM:.0f} A of the protein (in-plane, headgroup anchor to the "
            "nearest protein heavy atom in that leaflet), divided by the species' mean density over the leaflet: 1 = as common as "
            "anywhere, > 1 enriched next to the protein, < 1 depleted; +- is the counting (Poisson) uncertainty, ratio / sqrt(n) for "
            "n molecules in the shell, so a ratio from one or two molecules is not a finding. N = molecules of the species in the "
            "leaflet / within the first shell. Curves: membrane_validation_protein_density.png.\n"]
    for leaflet in ("upper", "lower"):
        species = sorted({sp for s in stages for sp in ((s["protein_density"].get(leaflet) or {}).get("species") or {})},
                         key=lambda sp: (sp != "all", sp))
        text += [f"\n**{leaflet} leaflet**\n", "| Species | " + " | ".join(s["label"] for s in stages) + " |",
                 "|---|" + "---:|" * len(stages)]
        for sp in species:
            cells = []
            for s in stages:
                c = ((s["protein_density"].get(leaflet) or {}).get("species") or {}).get(sp)
                if not c:
                    cells.append("-")
                    continue
                ratio, shell = c["first_shell_ratio"], c["first_shell_n"]
                error = f" ± {ratio / math.sqrt(shell):.2f}" if ratio is not None and shell > 0 else ""
                cells.append(f"{fmt(ratio, 2)}{error} (N {c['n']}/{shell})" + ("" if c["well_sampled"] else " few"))
            text.append(f"| {'all lipids' if sp == 'all' else display_name(sp)} | " + " | ".join(cells) + " |")
    return "\n".join(text) + "\n"


def density_figure(out: Path, record: dict, plt) -> Path | None:
    """Per-species density around the protein: rows = leaflets, columns = stages; one fixed colour per species."""
    stages = density_stages(record)
    if not stages:
        return None
    fig, axes = plt.subplots(2, len(stages), figsize=(3.1 * len(stages) + 1.6, 6.2), sharex=True, sharey="row", squeeze=False)
    for row, leaflet in enumerate(("upper", "lower")):
        for col, s in enumerate(stages):
            ax = axes[row][col]
            block = (s["protein_density"].get(leaflet) or {}).get("species") or {}
            ax.axhline(1.0, color="#c3c2b7", linewidth=0.8)
            for name in sorted(block, key=lambda sp: (sp == "all", SPECIES_SLOT.get(sp, 99))):
                c = block[name]
                r = np.array(c["R_A"])
                y = np.array([np.nan if v is None else v for v in c["cumulative_ratio"]], dtype=float)
                y[r < DENSITY_FIGURE_FROM_A] = np.nan  # the innermost radii hold a fraction of a lipid's area
                if name == "all":
                    ax.plot(r, y, color="#52514e", linewidth=1.4, linestyle="--", label="all lipids")
                else:
                    ax.plot(r, y, color=CATEGORICAL[SPECIES_SLOT.get(name, 7)], linewidth=1.4,
                            label=f"{display_name(name)} (N={c['n']})" + ("" if c["well_sampled"] else ", few"))
            if row == 0:
                ax.set_title(s["label"], fontsize=9)
            if col == 0:
                ax.set_ylabel(f"{leaflet} leaflet\ndensity within R / leaflet mean")
            if row == 1:
                ax.set_xlabel("R, distance from the protein (Å)")
            ax.axvline(10 * PROFILE_FIRST_SHELL_NM, color="#c3c2b7", linewidth=0.8, linestyle=":")
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            ax.tick_params(labelsize=7)
        handles, labels = axes[row][-1].get_legend_handles_labels()
        axes[row][-1].legend(handles, labels, fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    fig.suptitle("Lipids around the protein: density of each species within R of the protein, over its leaflet mean\n"
                 "(headgroup anchors; in-plane distance to the protein at headgroup depth; dotted line = first shell)", fontsize=8.5)
    fig.tight_layout()
    path = out / "membrane_validation_protein_density.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def write_outputs(out: Path, record: dict, stages: dict) -> list[Path]:
    """Write the JSON, Markdown, TSV and figures of the validation."""
    written = []
    (out / "membrane_validation.json").write_text(json.dumps(record, indent=2, default=str) + "\n")
    (out / "membrane_validation.md").write_text(write_markdown(record))
    written += [out / "membrane_validation.json", out / "membrane_validation.md"]
    rows = [r for name in record_stage_order(stages) for r in per_lipid_rows(stages[name]["stage"])]
    if rows:
        keys = list(rows[0])
        (out / "membrane_validation_lipids.tsv").write_text("\t".join(keys) + "\n" + "".join(
            "\t".join("" if r[k] is None else str(r[k]) for k in keys) + "\n" for r in rows))
        written.append(out / "membrane_validation_lipids.tsv")
    written += render_figures(out, record, stages)
    return written


MATRIX_ROWS = [("Ring penetration / clashes", "pipeline (ring_piercing.json, lipid scan, audit)", ("slice", "backmap", "minimize")),
               ("Voronoi APL", "measure", ("cg_frame", "embed", "slice", "backmap", "minimize", "equilibrated", "final")),
               ("Lateral RDF", "rdfs", ("cg_frame", "embed", "slice")),
               ("Lipids around the protein", "protein_density", ("embed", "slice", "backmap", "minimize", "equilibrated", "final")),
               ("Bilayer thickness", "thickness", ("cg_frame", "embed", "slice", "backmap", "minimize", "equilibrated", "final")),
               ("Composition / asymmetry", "measure", ("cg_frame", "embed", "slice", "backmap", "minimize", "equilibrated", "final")),
               ("Z-density / core hydration", "core_hydration", ("minimize", "equilibrated", "final")),
               ("Tail interdigitation", "interdigitation", ("backmap", "minimize", "equilibrated", "final")),
               ("Protein orientation", "orientation", ("cg_frame", "embed", "slice", "backmap", "minimize", "equilibrated", "final")),
               ("Tail order S_CD", None, ("equilibrated", "final")),
               ("Area / volume convergence", None, ("equilibrated", "final")),
               ("Lateral diffusion", None, ("final",)), ("K_A", None, ("final",)), ("Leaflet tension", None, ("final",))]
MATRIX_COLUMNS = ["cg_frame", "embed", "slice", "backmap", "minimize", "equilibrated", "final"]


def stage_matrix(stages: dict) -> list[dict]:
    """Which validation ran at which stage in this build ('done'), is planned there ('planned') or does not apply ('-')."""
    rows = []
    for label, key, applies in MATRIX_ROWS:
        cells = {}
        for column in MATRIX_COLUMNS:
            if column not in applies:
                cells[column] = "-"
            elif key is None or column in ("equilibrated", "final"):
                cells[column] = "planned" if column not in stages else ("done" if key and stages[column].get(key) is not None else "planned")
            elif column in stages and (key.startswith("pipeline") or stages[column].get(key) is not None):
                cells[column] = "done"
            else:
                cells[column] = "planned"
        rows.append({"validation": label, **cells})
    return rows


def matrix_markdown(rows: list[dict]) -> str:
    """The stage matrix as a Markdown table (checkmark = ran in this build, o = runs at that stage, blank = not applicable)."""
    mark = {"done": "✓", "planned": "o", "-": ""}
    head = "| Validation | " + " | ".join(STAGE_LABELS.get(c, c) for c in MATRIX_COLUMNS) + " |\n|---|" + "---|" * len(MATRIX_COLUMNS) + "\n"
    return head + "\n".join(f"| {r['validation']} | " + " | ".join(mark[r[c]] for c in MATRIX_COLUMNS) + " |" for r in rows) + "\n"


def record_stage_order(stages: dict) -> list[str]:
    """Stages in build order."""
    order = ["cg_frame", "embed", "slice", "backmap", "minimize"]
    return [n for n in order if n in stages] + [n for n in stages if n not in order]


def render_figures(out: Path, record: dict, stages: dict) -> list[Path]:
    """Four PNG figures (global APL, species APL, RDF, composition) when matplotlib is available; else nothing."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    written = []
    names = record_stage_order(stages)
    summaries = {s["name"]: s for s in record["stages"]}
    colours = [STAGE_COLOURS.get(n, "#4a3aa7") for n in names]
    # A. global APL by leaflet and stage
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    width = 0.8 / max(len(names), 1)
    for i, name in enumerate(names):
        for j, leaflet in enumerate(("upper", "lower")):
            value = summaries[name]["leaflets"][leaflet]["apl_A2"]
            ax.bar(j + (i - (len(names) - 1) / 2) * width, value, width * 0.92, color=colours[i], hatch=LEAFLET_HATCH[leaflet],
                   edgecolor="white", linewidth=0.5, label=summaries[name]["label"] if j == 0 else None)
            ax.text(j + (i - (len(names) - 1) / 2) * width, value + 0.4, f"{value:.1f}\nN={summaries[name]['leaflets'][leaflet]['lipids']}",
                    ha="center", va="bottom", fontsize=6.5, color="#0b0b0b")
    old = record.get("old_rule_estimate")
    if old and old["leaflets"]["upper"]["apl_A2"]:
        for j, leaflet in enumerate(("upper", "lower")):
            ax.axhline(old["leaflets"][leaflet]["apl_A2"], xmin=0.08 + 0.5 * j, xmax=0.42 + 0.5 * j, color=STAGE_COLOURS["old_rule"],
                       linestyle="--", linewidth=1.2, label="old all-beads-inside slice (estimate)" if j == 0 else None)
    ax.set_xticks([0, 1], ["upper leaflet", "lower leaflet"])
    ax.set_ylabel("area per lipid (Å²)")
    ax.set_title("Global leaflet APL by stage (periodic Voronoi, protein cells discarded)", fontsize=10)
    ax.legend(fontsize=7, frameon=False, loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.12))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_ylim(0, max(s["leaflets"][l]["apl_A2"] for s in summaries.values() for l in LEAFLETS) * 1.3)
    fig.tight_layout()
    fig.savefig(out / "membrane_validation_apl.png", dpi=150)
    plt.close(fig)
    written.append(out / "membrane_validation_apl.png")
    # B. species APL by leaflet
    species = sorted({r["lipid"] for s in record["stages"] for r in s["species"]})
    fig, axes = plt.subplots(2, 1, figsize=(max(7.5, 1.1 * len(species) + 2), 6.2), sharex=True)
    for ax, leaflet in zip(axes, ("upper", "lower")):
        for i, name in enumerate(names):
            rows = {r["lipid"]: r for r in summaries[name]["species"] if r["leaflet"] == leaflet}
            x = np.arange(len(species)) + (i - (len(names) - 1) / 2) * width
            medians = np.array([rows[sp]["apl_median_A2"] if sp in rows else np.nan for sp in species])
            half = np.array([0.5 * rows[sp]["apl_iqr_A2"] if sp in rows else 0.0 for sp in species])   # the interquartile range
            ax.bar(x, medians, width * 0.92, yerr=half, color=colours[i], edgecolor="white", linewidth=0.5, error_kw={"elinewidth": 0.8},
                   label=summaries[name]["label"])
            for xi, sp in zip(x, species):
                if sp in rows:
                    ax.text(xi, medians[list(species).index(sp)] + half[list(species).index(sp)] + 1.0, f"N={rows[sp]['n']}", ha="center",
                            va="bottom", fontsize=5.5, rotation=90, color="#0b0b0b")
        ax.set_ylabel(f"{leaflet} leaflet\nmedian local area (Å²)")
        ax.set_ylim(0, None)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[1].set_xticks(np.arange(len(species)), [display_name(sp) for sp in species], fontsize=8)
    axes[0].set_title("Species-specific local APL: median of the per-lipid Voronoi cells, bars = interquartile range, N per bar", fontsize=10)
    axes[1].legend(fontsize=7, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    fig.tight_layout()
    fig.savefig(out / "membrane_validation_species.png", dpi=150)
    plt.close(fig)
    written.append(out / "membrane_validation_species.png")
    # C. RDF: embed vs slice (primary), all-atom stages dashed
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
    for ax, leaflet in zip(axes, ("upper", "lower")):
        for i, name in enumerate(names):
            curve = stages[name]["rdfs"].get(leaflet, {}).get("all")
            if curve is None or name == "cg_frame":
                continue
            ax.plot(10.0 * curve["r_nm"], curve["g"], color=colours[i], linewidth=2.0 if stages[name]["stage"]["resolution"] == "cg" else 1.4,
                    linestyle="-" if stages[name]["stage"]["resolution"] == "cg" else "--", label=summaries[name]["label"])
        ax.axhline(1.0, color="#c3c2b7", linewidth=0.8)
        ax.set_title(f"{leaflet} leaflet: headgroup anchor g(r)", fontsize=10)
        ax.set_xlabel("r (Å)")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("g(r)")
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle("Lateral RDF: embedded vs sliced CG (solid, the slicing test); all-atom stages dashed (cross-resolution)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "membrane_validation_rdf.png", dpi=150)
    plt.close(fig)
    written.append(out / "membrane_validation_rdf.png")
    density = density_figure(out, record, plt)
    if density is not None:
        written.append(density)
    # D. composition: counts by species and leaflet, embed vs slice
    pair = [n for n in ("embed", "slice") if n in summaries]
    if len(pair) == 2:
        fig, axes = plt.subplots(1, 2, figsize=(max(8, 1.2 * len(species) + 2), 3.6))
        for ax, leaflet in zip(axes, ("upper", "lower")):
            for i, name in enumerate(pair):
                counts = Counter(summaries[name]["leaflets"][leaflet]["composition"])
                total = sum(counts.values()) or 1
                x = np.arange(len(species)) + (i - 0.5) * 0.4
                ax.bar(x, [counts.get(sp, 0) for sp in species], 0.38, color=STAGE_COLOURS[name], edgecolor="white", label=summaries[name]["label"])
                for xi, sp in zip(x, species):
                    ax.text(xi, counts.get(sp, 0) + 0.5, f"{100 * counts.get(sp, 0) / total:.0f}%", ha="center", va="bottom", fontsize=6,
                            color="#0b0b0b")
            ax.set_xticks(np.arange(len(species)), [display_name(sp) for sp in species], fontsize=7, rotation=30)
            ax.set_title(f"{leaflet} leaflet lipid counts (labels: mole %)", fontsize=10)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
        axes[0].set_ylabel("lipids")
        axes[0].legend(fontsize=7, frameon=False)
        fig.tight_layout()
        fig.savefig(out / "membrane_validation_composition.png", dpi=150)
        plt.close(fig)
        written.append(out / "membrane_validation_composition.png")
    return written


def protein_of_cg_frame(cg_protein: list[list[dict]]) -> np.ndarray:
    """Bead coordinates (nm) of the frame's own coarse-grained protein, the footprint of the uncut frame."""
    return np.array([[b["x"], b["y"], b["z"]] for res in cg_protein for b in res], dtype=float)


__all__ += ['protein_of_cg_frame', 'residues_in_order_by_chain']
_ = (AMINO, composition_table, residues_in_order)  # imported for re-use by callers of this module
