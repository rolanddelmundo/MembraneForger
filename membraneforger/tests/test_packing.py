"""Packing, composition, RDF and QC: Voronoi areas per lipid, global and species APL, composition, lateral RDF under PBC.

Synthetic membranes with known answers, plus the real 6WHC frame cut and compared with itself.
"""
import json
import tempfile
import unittest
from collections import Counter
from pathlib import Path

import numpy as np

from .common import CG_GRO, MAPPING, REPO, mf

CELL = [8.0, 8.0, 12.0]
MIDPLANE = 6.0
SETTINGS = mf.Settings()
# the configured RDF criteria: first-shell shift PASS within 0.5 A, WARNING within 1.0 A, FAIL beyond or above 3 x counting noise
RDF_CRITERIA = (SETTINGS.rdf_max_peak_shift_a, SETTINGS.rdf_max_noise_units, SETTINGS.rdf_peak_warning_a)


def lattice(nx, ny, spacing, z, name="POPC", jitter=0.0, seed=0, spread=0.2):
    """A square lattice of two-bead lipids (PO4 head at z, tail bead 1 nm towards the midplane and `spread` along x)."""
    rng = np.random.default_rng(seed)
    out = []
    for i in range(nx):
        for j in range(ny):
            x, y = (i + 0.5) * spacing + jitter * rng.normal(), (j + 0.5) * spacing + jitter * rng.normal()
            tail = z - 1.0 if z > MIDPLANE else z + 1.0
            beads = [{"atom": "PO4", "x": x, "y": y, "z": z}, {"atom": "C1A", "x": x + spread, "y": y, "z": tail}]
            out.append({"cg": name, "aa": name, "cls": "phospholipid", "beads": beads, "xyz": np.array([[b["x"], b["y"], b["z"]] for b in beads])})
    return out


def bilayer(nx=10, ny=10, spacing=0.8, jitter=0.0, species=("POPC",), spread=0.2):
    """Two lattice leaflets; species are assigned round-robin so every species has a known count."""
    upper = lattice(nx, ny, spacing, MIDPLANE + 2.0, jitter=jitter, seed=1, spread=spread)
    lower = lattice(nx, ny, spacing, MIDPLANE - 2.0, jitter=jitter, seed=2, spread=spread)
    for k, mol in enumerate(upper + lower):
        mol["cg"] = mol["aa"] = species[k % len(species)]
    return upper + lower


class Voronoi(unittest.TestCase):
    def test_a_square_lattice_gives_one_cell_of_spacing_squared_per_lipid(self):
        box = np.array([8.0, 8.0])
        points = np.array([[(i + 0.5) * 0.8, (j + 0.5) * 0.8] for i in range(10) for j in range(10)])
        areas = mf.periodic_voronoi_areas(points, box)
        self.assertEqual(len(areas), 100)                                     # one value per lipid
        self.assertTrue(np.allclose(areas, 0.64, atol=1e-6))
        self.assertAlmostEqual(float(areas.sum()), 64.0, places=6)            # the cells tile the periodic cell exactly

    def test_cells_tile_the_cell_with_jittered_points_and_extra_generators(self):
        rng = np.random.default_rng(3)
        points = rng.uniform(0, 8, size=(60, 2))
        extra = rng.uniform(3, 5, size=(200, 2))                               # a dense "protein" in the middle
        areas = mf.periodic_voronoi_areas(points, np.array([8.0, 8.0]), extra)
        protein = mf.periodic_voronoi_areas(extra, np.array([8.0, 8.0]), points)
        self.assertAlmostEqual(float(areas.sum() + protein.sum()), 64.0, places=5)  # lipid cells + discarded protein cells = A_xy
        self.assertTrue((areas > 0).all())

    def test_a_lattice_bilayer_has_the_lattice_apl_in_both_leaflets(self):
        stage = mf.cg_stage("t", bilayer(), [], CELL, MIDPLANE)
        measure = mf.measure_packing(stage)
        for leaflet in ("upper", "lower"):
            self.assertEqual(measure[leaflet]["lipids"], 100)
            self.assertAlmostEqual(measure[leaflet]["apl_A2"], 64.0, places=4)
            self.assertAlmostEqual(measure[leaflet]["protein_area_nm2"], 0.0, places=6)
        self.assertEqual(sum("voronoi_area_nm2" in r for r in stage["lipids"]), 200)   # one area per lipid

    def test_the_protein_footprint_is_removed_from_the_lipid_area(self):
        membrane = [m for m in bilayer() if not (3.0 < m["beads"][0]["x"] < 5.0 and 3.0 < m["beads"][0]["y"] < 5.0)]  # a 2 x 2 nm hole, 4 lipids
        rng = np.random.default_rng(0)
        protein = [{"atom": "CA", "x": 10 * x, "y": 10 * y, "z": 10 * z}
                   for x, y, z in np.column_stack([rng.uniform(3.1, 4.9, 400), rng.uniform(3.1, 4.9, 400), rng.uniform(3.0, 9.0, 400)])]
        stage = mf.cg_stage("t", membrane, protein, CELL, MIDPLANE)
        measure = mf.measure_packing(stage)
        for leaflet in ("upper", "lower"):
            self.assertEqual(measure[leaflet]["lipids"], 96)
            self.assertGreater(measure[leaflet]["protein_area_nm2"], 2.0)      # about the 4 nm^2 hole, less half a cell at its rim
            self.assertLess(measure[leaflet]["protein_area_nm2"], 5.5)
            self.assertLess(abs(measure[leaflet]["apl_A2"] - 64.0), 3.0)      # lipids next to the protein are not inflated by the hole
        without = mf.measure_packing(mf.cg_stage("t", membrane, [], CELL, MIDPLANE))
        self.assertAlmostEqual(without["upper"]["apl_A2"], 6400.0 / 96, places=3)   # without the protein the hole is assigned to lipids


class Species(unittest.TestCase):
    def test_species_apl_is_the_mean_cell_of_that_species_and_n_is_reported(self):
        stage = mf.cg_stage("t", bilayer(species=("POPC", "CHOL", "SAP6", "POPC"), jitter=0.01), [], CELL, MIDPLANE)
        measure = mf.measure_packing(stage)
        rows = mf.species_table(stage)
        by = {(r["leaflet"], r["lipid"]): r for r in rows}
        self.assertEqual(by[("upper", "POPC")]["n"], 50)
        self.assertEqual(by[("upper", "CHOL")]["n"], 25)
        self.assertAlmostEqual(by[("upper", "CHOL")]["apl_mean_A2"], 64.0, delta=1.0)      # the local cell, not 6400 / 25
        self.assertAlmostEqual(by[("lower", "SAP6")]["mole_fraction"], 0.25, places=6)
        self.assertIsNone(by[("upper", "CHOL")]["head_apl_mean_A2"])                        # sterols have no head-plane cell
        self.assertAlmostEqual(measure["upper"]["head_apl_A2"], 6400.0 / 75, places=3)         # the 75 heads alone share the plane
        self.assertIsNotNone(by[("upper", "POPC")]["head_apl_mean_A2"])
        self.assertTrue(all(r["well_sampled"] for r in rows))

    def test_composition_distance_and_finite_crop_noise(self):
        reference = Counter({"POPC": 60, "CHOL": 30, "SAP6": 10})
        self.assertAlmostEqual(mf.composition_distance(reference, Counter({"POPC": 6, "CHOL": 3, "SAP6": 1})), 0.0)
        self.assertAlmostEqual(mf.composition_distance(reference, Counter({"POPC": 10})), 0.4)
        self.assertGreater(mf.expected_composition_distance(reference, 20), mf.expected_composition_distance(reference, 200))
        table = mf.composition_table({"upper": {"composition": dict(reference)}, "lower": {"composition": dict(reference)}},
                                     {"upper": {"composition": {"POPC": 30, "CHOL": 15, "SAP6": 5}}, "lower": {"composition": {"POPC": 50}}})
        self.assertAlmostEqual(table["upper"]["distance"], 0.0)
        self.assertTrue(table["upper"]["within_finite_crop_variability"])
        self.assertFalse(table["lower"]["within_finite_crop_variability"])
        self.assertIn("POPC", table["lower"]["flagged_species"])


class LateralRDF(unittest.TestCase):
    def test_rdf_of_a_periodic_lattice_peaks_at_the_lattice_spacing_and_is_the_same_in_every_image(self):
        rng = np.random.default_rng(7)
        points = np.array([[(i + 0.5) * 0.8, (j + 0.5) * 0.8] for i in range(10) for j in range(10)]) + 0.05 * rng.normal(size=(100, 2))
        curve = mf.lateral_rdf(points, np.array([8.0, 8.0]))
        features = mf.rdf_features(curve)
        self.assertAlmostEqual(features["first_peak_A"], 8.0, delta=1.5)                    # the smoothed first shell sits at the spacing
        self.assertGreater(features["first_peak_height"], 1.5)
        shifted = mf.lateral_rdf(np.mod(points + 3.3, 8.0), np.array([8.0, 8.0]))      # the same arrangement, moved across the boundary
        self.assertTrue(np.allclose(curve["g"], shifted["g"], atol=1e-9))

    def test_rdf_normalisation_of_random_points_is_one(self):
        rng = np.random.default_rng(5)
        points = rng.uniform(0, 20, size=(4000, 2))
        curve = mf.lateral_rdf(points, np.array([20.0, 20.0]))
        self.assertAlmostEqual(float(curve["g"][5:].mean()), 1.0, delta=0.03)

    def test_an_unchanged_periodic_arrangement_passes_the_comparison(self):
        stage = mf.cg_stage("a", bilayer(jitter=0.08), [], CELL, MIDPLANE)
        same = mf.cg_stage("b", bilayer(jitter=0.08), [], CELL, MIDPLANE)
        comparison = mf.rdf_comparison(mf.leaflet_rdfs(stage), mf.leaflet_rdfs(same), *RDF_CRITERIA)
        self.assertEqual(comparison["status"], "PASS")
        for leaflet in ("upper", "lower"):
            self.assertEqual(comparison["leaflets"][leaflet]["all"]["first_peak_shift_A"], 0.0)
            self.assertAlmostEqual(comparison["leaflets"][leaflet]["all"]["rms_difference"], 0.0)

    def test_a_dilated_arrangement_fails_the_comparison(self):
        stage = mf.cg_stage("a", bilayer(jitter=0.05), [], CELL, MIDPLANE)
        wide = mf.cg_stage("b", bilayer(nx=8, ny=8, spacing=1.0, jitter=0.05), [], CELL, MIDPLANE)   # the same cell, 1.0 nm lattice
        comparison = mf.rdf_comparison(mf.leaflet_rdfs(stage), mf.leaflet_rdfs(wide), *RDF_CRITERIA)
        self.assertEqual(comparison["status"], "FAIL")
        self.assertGreater(abs(comparison["leaflets"]["upper"]["all"]["first_peak_shift_A"]), 1.0)


class QC(unittest.TestCase):
    def test_classification_tiers(self):
        self.assertEqual(mf.classify(2.0, 3.0, 5.0), "PASS")
        self.assertEqual(mf.classify(-4.0, 3.0, 5.0), "WARNING")
        self.assertEqual(mf.classify(8.0, 3.0, 5.0), "FAIL")
        self.assertEqual(mf.classify(None, 3.0, 5.0), "INSUFFICIENT SAMPLING")
        control = np.linspace(-1, 1, 200)
        self.assertEqual(mf.classify_percentile(0.0, control)[0], "PASS")
        self.assertEqual(mf.classify_percentile(0.97, control)[0], "WARNING")
        self.assertEqual(mf.classify_percentile(5.0, control)[0], "FAIL")
        self.assertEqual(mf.classify_percentile(0.0, control[:5])[0], "INSUFFICIENT SAMPLING")

    def test_overall_status(self):
        rec = lambda status: mf.metric("m", "s", 1.0, "u", status=status)
        self.assertEqual(mf.overall_status([rec("PASS"), rec("PASS")])["status"], "PASS")
        self.assertEqual(mf.overall_status([rec("PASS"), rec("WARNING")])["status"], "PASS WITH WARNINGS")
        self.assertEqual(mf.overall_status([rec("PASS"), rec("INSUFFICIENT SAMPLING")])["status"], "INSUFFICIENT SAMPLING")
        self.assertEqual(mf.overall_status([rec("WARNING"), rec("FAIL")])["status"], "FAIL")
        self.assertEqual(mf.overall_status([rec("PASS")], required=("other",))["status"], "INSUFFICIENT SAMPLING")


class SliceGate(unittest.TestCase):
    """A lattice bilayer cut to half its size: the anchor rule keeps the lattice density, the old rule would not."""

    def test_anchor_crop_keeps_the_density_and_the_old_rule_loses_it(self):
        membrane = bilayer(nx=20, ny=20, spacing=0.8, jitter=0.03, species=("POPC", "DOPC", "CHOL", "POPE"), spread=0.5)  # tails cross edges
        cell = [16.0, 16.0, 12.0]
        solute = [{"atom": "CA", "resname": "ALA", "chain": "A", "resid": 1, "x": 60.0, "y": 60.0, "z": 60.0},
                  {"atom": "CA", "resname": "ALA", "chain": "A", "resid": 2, "x": 100.0, "y": 100.0, "z": 60.0}]
        embed = mf.cg_stage("embed", membrane, solute, cell, MIDPLANE)
        reference = mf.reference_packing(embed)
        centred = mf.Settings(box_xy_buffer_nm=2.0, slice_optimize_offset=False)
        cut = mf.slice_membrane_cg(membrane, solute, cell, MIDPLANE, centred, reference=reference)
        self.assertEqual(cut["box"][:2], [8.0, 8.0])
        shifted = mf.slice_membrane_cg(membrane, solute, cell, MIDPLANE, mf.Settings(box_xy_buffer_nm=2.0), reference=reference)
        self.assertLess(abs(shifted["report"]["lipids_after"] - 200), 12)                        # an offset window keeps the density too
        stage = mf.cg_stage("slice", cut["membrane"], cut["placed"], cut["box"], MIDPLANE)
        comparison = mf.compare_packing(embed, stage, 5.0)
        for leaflet in ("upper", "lower"):
            self.assertEqual(comparison["leaflets"][leaflet]["lipids"], 100)                    # 8 x 8 nm of a 0.8 nm lattice
            self.assertLess(abs(comparison["leaflets"][leaflet]["delta_apl_percent"]), 1.0)
            self.assertTrue(comparison["composition"][leaflet]["within_finite_crop_variability"])
        selection = cut["report"]["selection"]
        self.assertLess(selection["kept_by_all_beads_inside_rule"], selection["kept_by_anchor_rule"])   # the old rule drops edge lipids
        self.assertEqual(len(cut["kept_indices"]), len(set(cut["kept_indices"])))                        # no duplicates
        region = mf.region_reference(embed, cut["kept_indices"])
        self.assertEqual(region["upper"]["lipids"], 100)
        rdf = mf.rdf_comparison(mf.leaflet_rdfs(embed), mf.leaflet_rdfs(stage), *RDF_CRITERIA)
        self.assertEqual(rdf["status"], "PASS")

    def test_window_distribution_of_a_lattice_is_narrow(self):
        membrane = bilayer(nx=20, ny=20, spacing=0.8, jitter=0.03)
        embed = mf.cg_stage("embed", membrane, [], [16.0, 16.0, 12.0], MIDPLANE)
        mf.measure_packing(embed)
        windows = mf.window_distribution(embed, np.array([8.0, 8.0]), [True, True], samples=6)
        self.assertEqual(windows["upper"]["windows"], 36)
        self.assertLess(float(np.ptp(windows["upper"]["window_apl_A2"])), 2.0)


class RealFrameGate(unittest.TestCase):
    """The real 6WHC frame: embedded membrane -> slice -> gate, with the report written to a temporary directory."""

    def test_the_6whc_slice_passes_the_gate_and_writes_the_report(self):
        aa, _ = mf.read_all_atom(REPO / "examples/preeq_cg_cellmem/6WHC_MTZP_prot-lig.pdb", REPO / "forcefield")
        raw, box = mf.read_cg(CG_GRO)
        cg = mf.classify_cg(raw, MAPPING)
        slab = mf.make_membrane_whole(cg["membrane"], box)
        midplane = mf.bilayer_midplane(cg["membrane"], slab)[0]
        settings = mf.Settings()
        fit = mf.map_all_atom_to_cg(aa, cg["protein"], box, slab, settings)
        placed = mf.apply_transform(aa, fit["R"], fit["t"])
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / mf.LOG_NAME).write_text("")
            analysis = mf.MembraneValidation(out, settings)
            analysis.add_cg("cg_frame", cg["membrane"], [], box, midplane, mf.protein_of_cg_frame(cg["protein"]))
            embedded = analysis.add_cg("embed", cg["membrane"], placed, box, midplane)
            for leaflet in ("upper", "lower"):
                self.assertTrue(44.0 < embedded["leaflets"][leaflet]["apl_A2"] < 58.0)        # about 50 A^2 per lipid, protein removed
            cut = mf.slice_membrane_cg(cg["membrane"], placed, box, midplane, settings, None, analysis.slice_reference())
            analysis.add_cg("slice", cut["membrane"], cut["placed"], cut["box"], midplane)
            check = analysis.validate_slice(cut, cg["membrane"], placed, box)
            self.assertEqual(check["status"], "PASS")                                         # the configured gate as a whole
            self.assertTrue(all(check["verdict"].values()))
            for leaflet in ("upper", "lower"):
                v = check["leaflets"][leaflet]
                self.assertEqual(v["construction_status"], "PASS")                         # <= apl_slice_warning_percent (3 %)
                self.assertLessEqual(abs(v["construction_delta_percent"]), settings.apl_slice_warning_percent)
                self.assertIn(v["composition_status"], ("PASS",))
            records = {(r["metric"], r["leaflet"]): r for r in check["records"]}
            self.assertEqual(records[("slice integrity", None)]["status"], "PASS")
            self.assertTrue(all(records[("headgroup RDF first-shell peak shift", l)]["status"] == "PASS" for l in ("upper", "lower")))
            files = analysis.write(cut)
            names = {f.name for f in files}
            self.assertTrue({"membrane_validation.json", "membrane_validation.md", "membrane_validation_lipids.tsv"} <= names)
            record = json.loads((out / "membrane_validation.json").read_text())
            self.assertEqual([s["name"] for s in record["stages"]], ["cg_frame", "embed", "slice"])
            self.assertIn("Old versus new slicing algorithm", (out / "membrane_validation.md").read_text())
            rows = (out / "membrane_validation_lipids.tsv").read_text().splitlines()
            self.assertEqual(len(rows) - 1, sum(s["leaflets"]["upper"]["lipids"] + s["leaflets"]["lower"]["lipids"] for s in record["stages"]))


if __name__ == "__main__":
    unittest.main()


class ParentWindowRDF(unittest.TestCase):
    def test_window_rdf_distribution_of_a_lattice_is_narrow_and_a_window_of_it_passes(self):
        membrane = bilayer(nx=20, ny=20, spacing=0.8, jitter=0.04)
        embed = mf.cg_stage("embed", membrane, [], [16.0, 16.0, 12.0], MIDPLANE)
        windows = mf.window_rdf_distribution(embed, np.array([8.0, 8.0]), [True, True], 2.5, samples=6)
        self.assertEqual(windows["upper"]["windows"], 36)
        self.assertLess(float(windows["upper"]["rms"].max()), 0.6)
        whole = mf.leaflet_rdfs(embed)["upper"]["all"]
        inside = [r for r in embed["lipids"] if r["leaflet"] == "upper" and r["x_nm"] < 8.0 and r["y_nm"] < 8.0]
        one = mf.lateral_rdf(np.array([[r["x_nm"], r["y_nm"]] for r in inside]), np.array([16.0, 16.0]), r_max=2.5, area=64.0)
        rms = float(np.sqrt(np.mean((one["g"][3:20] - whole["g"][3:20]) ** 2)))
        self.assertNotEqual(mf.classify_percentile(rms, windows["upper"]["rms"])[0], "FAIL")   # a parent window is within its own distribution
