"""Accessible-shell-area normalization of g(r): validated against random anchors in a geometry-matched open area."""
import unittest

import numpy as np

from .common import mf

BOX = np.array([6.0, 6.0])


def disk(radius, centre=(3.0, 3.0)):
    """Heavy atoms filling a disk of the given radius (a stand-in for the protein cross-section)."""
    atoms = [centre]
    for r in np.arange(0.15, radius + 1e-9, 0.15):
        n = max(6, int(2 * np.pi * r / 0.15))
        atoms += [(centre[0] + r * np.cos(a), centre[1] + r * np.sin(a)) for a in np.linspace(0, 2 * np.pi, n, endpoint=False)]
    return np.array(atoms)


class OpenArea(unittest.TestCase):
    def test_no_protein_leaves_the_whole_cell_open(self):
        grid = mf.open_area_grid(np.zeros((0, 2)), BOX)
        self.assertAlmostEqual(grid["accessible_nm2"], 36.0, places=6)

    def test_the_protein_closes_its_footprint_plus_the_accessible_margin(self):
        grid = mf.open_area_grid(disk(1.0), BOX)
        closed = 36.0 - grid["accessible_nm2"]
        self.assertAlmostEqual(closed, np.pi * 1.2 ** 2, delta=0.2)        # outermost atoms at 0.9 nm + 0.3 nm accessible margin

    def test_shell_areas_of_a_free_anchor_are_full_annuli(self):
        grid = mf.open_area_grid(np.zeros((0, 2)), BOX)
        edges = np.arange(0.0, 1.01, 0.1)
        area = mf.shell_areas(np.array([[3.0, 3.0]]), grid, edges)[0]
        expected = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
        np.testing.assert_allclose(np.cumsum(area)[4:], np.cumsum(expected)[4:], rtol=0.03)   # the 0.05 nm grid: cumulative area within r

    def test_an_anchor_beside_the_protein_has_less_room(self):
        grid = mf.open_area_grid(disk(1.0), BOX)
        edges = np.arange(0.0, 2.01, 0.1)
        near = mf.shell_areas(np.array([[3.0 + 1.5, 3.0]]), grid, edges)[0]
        far = mf.shell_areas(np.array([[0.2, 0.2]]), grid, edges)[0]
        self.assertLess(near[10:].sum(), 0.9 * far[10:].sum())


class AgainstRandomAnchors(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grid = mf.open_area_grid(disk(1.2), BOX)
        cls.ideal = mf.random_reference_rdfs(cls.grid, 120, replicates=40, seed=1, r_max=2.5)

    def test_the_shell_curve_is_one_for_an_ideal_gas_and_the_existing_one_is_not(self):
        r = self.ideal["r_nm"]
        band = (r >= 0.3) & (r <= 2.5)
        self.assertAlmostEqual(float(self.ideal["shell_mean"][band].mean()), 1.0, delta=0.015)
        self.assertLess(float(self.ideal["raw_mean"][band].mean()), 0.96)                    # the offset the shell area removes

    def test_the_deficit_of_the_existing_curve_grows_with_the_protein_share(self):
        offsets = []
        for radius in (0.6, 1.2, 1.8):
            grid = mf.open_area_grid(disk(radius), BOX)
            ref = mf.random_reference_rdfs(grid, 100, replicates=25, seed=2, r_max=2.5)
            band = (ref["r_nm"] >= 0.5) & (ref["r_nm"] <= 2.5)
            offsets.append(float(ref["raw_mean"][band].mean()) - 1.0)
            self.assertAlmostEqual(float(ref["shell_mean"][band].mean()), 1.0, delta=0.02)
        self.assertGreater(offsets[0], offsets[1])
        self.assertGreater(offsets[1], offsets[2])

    def test_hard_core_anchors_are_one_beyond_the_core_with_shell_normalization(self):
        ref = mf.random_reference_rdfs(self.grid, 70, replicates=30, hard_core_nm=0.45, seed=3, r_max=2.5)
        r = ref["r_nm"]
        self.assertLess(float(ref["shell_mean"][r < 0.4].max()), 0.05)                          # nothing inside the core
        beyond = (r >= 0.9) & (r <= 2.5)
        self.assertAlmostEqual(float(ref["shell_mean"][beyond].mean()), 1.0, delta=0.05)
        self.assertLess(float(ref["raw_mean"][beyond].mean()), float(ref["shell_mean"][beyond].mean()) - 0.03)

    def test_an_impossible_packing_is_refused_not_looped_on(self):
        grid = mf.open_area_grid(np.zeros((0, 2)), np.array([2.0, 2.0]))
        with self.assertRaises(SystemExit):
            mf.random_reference_rdfs(grid, 60, replicates=1, hard_core_nm=0.45)

    def test_without_a_protein_both_normalizations_agree(self):
        grid = mf.open_area_grid(np.zeros((0, 2)), BOX)
        rng = np.random.default_rng(5)
        points = rng.uniform(0, 6, size=(150, 2))
        raw = mf.lateral_rdf(points, BOX, area=36.0, r_max=2.5)
        shell = mf.lateral_rdf_shell(points, grid, r_max=2.5)
        band = raw["r_nm"] >= 0.5
        np.testing.assert_allclose(shell["g"][band], raw["g"][band], atol=0.15)
        self.assertAlmostEqual(float(shell["g"][band].mean()), float(raw["g"][band].mean()), delta=0.02)

    def test_the_partner_form_is_one_for_independent_sets(self):
        rng = np.random.default_rng(6)
        pool = self.grid["open_xy"]
        pick = lambda n: pool[rng.integers(len(pool), size=n)] + rng.uniform(-0.025, 0.025, size=(n, 2))
        curves = [mf.lateral_rdf_shell(pick(80), self.grid, r_max=2.5, partners=pick(80))["g"] for _ in range(30)]
        r = mf.lateral_rdf_shell(pick(5), self.grid, r_max=2.5, partners=pick(5))["r_nm"]
        self.assertAlmostEqual(float(np.mean(curves, axis=0)[r >= 0.3].mean()), 1.0, delta=0.02)

    def test_compare_normalizations_reports_the_offsets(self):
        raw = {"r_nm": self.ideal["r_nm"], "g": self.ideal["raw_mean"]}
        shell = {"r_nm": self.ideal["r_nm"], "g": self.ideal["shell_mean"]}
        result = mf.compare_normalizations(raw, shell)
        self.assertLess(result["offset_existing_percent"], -3.0)
        self.assertLess(abs(result["offset_shell_percent"]), 1.5)


class LeafletCurves(unittest.TestCase):
    def test_a_stage_gives_a_curve_per_leaflet_without_touching_the_existing_ones(self):
        rng = np.random.default_rng(7)
        pool = mf.open_area_grid(disk(1.2), BOX)["open_xy"]
        lipids = [{"leaflet": leaflet, "resname": "POPC", "x_nm": float(p[0]), "y_nm": float(p[1]), "z_nm": 1.0 if leaflet == "upper" else -1.0,
                   "index": i + (0 if leaflet == "upper" else 1000)}
                  for leaflet in ("upper", "lower") for i, p in enumerate(pool[rng.integers(len(pool), size=80)])]
        protein = np.array([[x, y, z] for x, y in disk(1.2) for z in (-1.0, 1.0)])
        stage = {"box_nm": [6.0, 6.0, 8.0], "lipids": lipids, "protein_nm": protein}
        shell = mf.leaflet_rdfs_shell(stage, r_max_nm=2.5)
        self.assertEqual(set(shell), {"lower", "upper"})
        for leaflet in shell:
            curve = shell[leaflet]["all"]
            self.assertEqual(curve["normalization"], "accessible shell area")
            self.assertIn("features", curve)
            self.assertAlmostEqual(curve["accessible_nm2"], mf.open_area_grid(disk(1.2), BOX)["accessible_nm2"], places=6)


if __name__ == "__main__":
    unittest.main()
