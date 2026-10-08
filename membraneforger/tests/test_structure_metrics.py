"""Structural metrics on synthetic membranes: thickness, protein tilt and depth, z profiles, core hydration, interdigitation.

Every system here is built from a few lines of numpy so that the expected answer is known exactly up to binning.
"""
import unittest

import numpy as np

from .common import mf

BOX = [8.0, 8.0, 12.0]
MIDPLANE = 4.0


def plane(z, n=8, resname="POPC", leaflet="upper", box=BOX):
    """n x n anchors on a square lattice at height z."""
    step = box[0] / n
    rows = []
    for i in range(n):
        for j in range(n):
            rows.append({"index": len(rows), "resid": len(rows) + 1, "resname": resname, "leaflet": leaflet,
                         "x_nm": (i + 0.5) * step, "y_nm": (j + 0.5) * step, "z_nm": z, "z_range_nm": (z - 2.0, z)})
    return rows


def stage(upper_z=6.0, lower_z=2.0, extra=()):
    """A flat two-plane stage with the planes at upper_z and lower_z (midplane 4 nm)."""
    lipids = plane(upper_z, leaflet="upper") + plane(lower_z, leaflet="lower") + list(extra)
    for k, r in enumerate(lipids):
        r["index"], r["resid"] = k, k + 1
    return {"name": "synthetic", "resolution": "cg", "box_nm": list(BOX), "midplane_nm": MIDPLANE, "lipids": lipids, "protein_nm": np.zeros((0, 3))}


def rod(tilt_deg=0.0, length=6.0, n=61, centre_z=MIDPLANE, x=4.0, y=4.0):
    """A straight line of points along z (or tilted by tilt_deg about y) through (x, y, centre_z), with a little width."""
    t = np.linspace(-0.5 * length, 0.5 * length, n)
    axis = np.array([np.sin(np.radians(tilt_deg)), 0.0, np.cos(np.radians(tilt_deg))])
    rng = np.random.default_rng(1)
    pts = np.outer(t, axis) + rng.normal(scale=0.05, size=(n, 3))
    return pts + np.array([x, y, centre_z])


class Thickness(unittest.TestCase):
    def test_two_flat_planes_give_their_separation_and_a_flat_map(self):
        out = mf.bilayer_thickness(stage(6.0, 2.0), bins_nm=1.0)
        self.assertAlmostEqual(out["thickness_nm"], 4.0, places=6)
        self.assertAlmostEqual(out["thickness_A"], 40.0, places=5)
        self.assertEqual(out["local_map"].shape, (8, 8))
        self.assertEqual(out["cells_filled"], 64)
        self.assertTrue(np.allclose(out["local_map"], 4.0))
        self.assertLess(out["local_sd_nm"], 1e-9)
        self.assertEqual(out["n_per_leaflet"], {"lower": 64, "upper": 64})

    def test_sterol_anchors_below_the_phosphates_do_not_thin_the_bilayer(self):
        chol = [{"index": 0, "resid": 0, "resname": "CHOL", "leaflet": "upper", "x_nm": 0.3 + 0.1 * k, "y_nm": 0.3, "z_nm": 5.5,
                 "z_range_nm": (4.0, 5.5)} for k in range(30)]
        out = mf.bilayer_thickness(stage(6.0, 2.0, extra=chol))
        self.assertAlmostEqual(out["thickness_nm"], 4.0, places=6)
        self.assertEqual(out["n_per_leaflet"]["upper"], 64)

    def test_an_empty_cell_is_nan_not_zero(self):
        s = stage(6.0, 2.0)
        s["lipids"] = [r for r in s["lipids"] if not (r["leaflet"] == "upper" and r["x_nm"] < 1.0 and r["y_nm"] < 1.0)]
        out = mf.bilayer_thickness(s, bins_nm=1.0)
        self.assertTrue(np.isnan(out["local_map"][0, 0]))
        self.assertEqual(out["cells_filled"], 63)

    def test_thickness_change_grades_with_the_supplied_ceilings(self):
        ref = mf.bilayer_thickness(stage(6.0, 2.0))
        for factor, expected in ((1.02, "PASS"), (0.96, "WARNING"), (1.08, "FAIL")):
            new = mf.bilayer_thickness(stage(4.0 + 2.0 * factor, 4.0 - 2.0 * factor))
            change = mf.thickness_change(ref, new, 3.0, 5.0)
            self.assertEqual(change["status"], expected, factor)
            self.assertAlmostEqual(change["percent"], 100.0 * (factor - 1.0), places=6)
            self.assertEqual(change["records"][0]["status"], expected)
            self.assertEqual(change["records"][0]["units"], "A")


class Orientation(unittest.TestCase):
    def test_a_rod_along_z_has_no_tilt(self):
        out = mf.protein_orientation(rod(0.0), MIDPLANE)
        self.assertLess(out["tilt_deg"], 2.0)
        self.assertGreater(out["n_slab"], 10)

    def test_a_rod_tilted_by_20_degrees_reads_about_20(self):
        out = mf.protein_orientation(rod(20.0), MIDPLANE)
        self.assertAlmostEqual(out["tilt_deg"], 20.0, delta=2.0)

    def test_depth_follows_the_transmembrane_atoms(self):
        # A 2 nm rod centred 0.3 nm above the midplane (inside the +-1.5 nm slab) with a big extramembrane blob far
        # above: the depth is that of the slab atoms, the centre-of-mass z is dragged up by the blob, polarity up.
        blob = np.random.default_rng(2).normal(scale=0.3, size=(400, 3)) + np.array([4.0, 4.0, MIDPLANE + 6.0])
        pts = np.vstack([rod(0.0, length=2.0, centre_z=MIDPLANE + 0.3), blob])
        out = mf.protein_orientation(pts, MIDPLANE)
        self.assertAlmostEqual(out["depth_nm"], 0.3, delta=0.1)
        self.assertGreater(out["com_z_nm"], 3.0)
        self.assertEqual(out["polarity"], 1)

    def test_ca_mask_restricts_the_atoms_used(self):
        pts = np.vstack([rod(0.0), rod(40.0, x=6.0)])
        mask = np.zeros(len(pts), dtype=bool)
        mask[: len(pts) // 2] = True
        self.assertLess(mf.protein_orientation(pts, MIDPLANE, ca_mask=mask)["tilt_deg"], 2.0)

    def test_inversion_is_seen_against_a_reference(self):
        blob = np.random.default_rng(3).normal(scale=0.3, size=(400, 3)) + np.array([4.0, 4.0, MIDPLANE + 6.0])
        up = mf.protein_orientation(np.vstack([rod(0.0), blob]), MIDPLANE)
        flipped = np.vstack([rod(0.0), blob]).copy()
        flipped[:, 2] = 2.0 * MIDPLANE - flipped[:, 2]
        down = mf.protein_orientation(flipped, MIDPLANE, reference=up)
        self.assertTrue(down["inverted"])
        same = mf.protein_orientation(np.vstack([rod(0.0), blob]), MIDPLANE, reference=up)
        self.assertFalse(same["inverted"])

    def test_hydrophobic_mismatch_counts_slab_atoms_outside_the_planes(self):
        out = mf.protein_orientation(rod(0.0), MIDPLANE, headgroup_planes_nm=(MIDPLANE - 1.0, MIDPLANE + 1.0))
        self.assertAlmostEqual(out["hydrophobic_mismatch_fraction"], 1.0 / 3.0, delta=0.1)   # slab +-1.5 nm, planes +-1 nm
        self.assertIsNone(mf.protein_orientation(rod(0.0), MIDPLANE)["hydrophobic_mismatch_fraction"])

    def test_orientation_change_grades_tilt_and_depth(self):
        ref = mf.protein_orientation(rod(0.0), MIDPLANE)
        tilted = mf.protein_orientation(rod(4.0), MIDPLANE)
        change = mf.orientation_change(ref, tilted, 1.0, 3.0, 0.5, 1.5)
        by_name = {r["metric"]: r for r in change["records"]}
        self.assertEqual(by_name["protein tilt"]["status"], "FAIL")
        self.assertAlmostEqual(change["tilt_change_deg"], 4.0, delta=1.5)
        short = mf.protein_orientation(rod(0.0, length=2.0), MIDPLANE)
        shifted = mf.protein_orientation(rod(0.0, length=2.0, centre_z=MIDPLANE + 0.1), MIDPLANE)
        change = mf.orientation_change(short, shifted, 1.0, 3.0, 0.5, 1.5)
        by_name = {r["metric"]: r for r in change["records"]}
        self.assertAlmostEqual(change["depth_change_A"], 1.0, delta=0.3)
        self.assertEqual(by_name["protein insertion depth"]["status"], "WARNING")
        self.assertEqual(by_name["protein tilt"]["status"], "PASS")


class Profiles(unittest.TestCase):
    def test_densities_conserve_the_particle_count_and_wrap_z(self):
        rng = np.random.default_rng(4)
        water = rng.uniform(0, BOX[2], size=(500, 3)) * np.array([BOX[0] / BOX[2], BOX[1] / BOX[2], 1.0])
        water[:100, 2] += BOX[2]      # images above the cell
        water[100:150, 2] -= BOX[2]   # and below
        heads = np.column_stack([rng.uniform(0, 8, 64), rng.uniform(0, 8, 64), np.full(64, 6.05)])
        out = mf.z_density_profiles({"water": water, "headgroups": heads}, BOX, bin_nm=0.1)
        area_bin = out["area_nm2"] * out["bin_nm"]
        self.assertAlmostEqual(float(out["density_nm3"]["water"].sum() * area_bin), 500.0, places=6)
        self.assertAlmostEqual(float(out["density_nm3"]["headgroups"].sum() * area_bin), 64.0, places=6)
        self.assertAlmostEqual(float(out["z_nm"][np.argmax(out["density_nm3"]["headgroups"])]), 6.05, delta=0.1)
        self.assertTrue(np.all(out["z_nm"] >= 0) and np.all(out["z_nm"] <= BOX[2]))
        self.assertEqual(out["n"]["water"], 500)

    def test_aa_groups_are_picked_by_residue_and_atom_name(self):
        atom = lambda rn, name, z: {"resname": rn, "atom": name, "x": 1.0, "y": 1.0, "z": z}
        atoms = [atom("POPC", "P", 6.0), atom("POPC", "N", 6.3), atom("POPC", "C22", 5.0), atom("POPC", "C316", 4.2),
                 atom("POPC", "H2R", 5.0), atom("POPC", "C1", 5.9), atom("POPC", "C13", 6.2),
                 atom("CHL1", "O3", 5.5), atom("CHL1", "C25", 4.5), atom("CHL1", "H3", 5.5),
                 atom("PSM", "P", 6.0), atom("PSM", "C5S", 5.0), atom("PSM", "C10F", 4.5),
                 atom("GLPA", "NF", 6.5), atom("GLPA", "C8S", 5.0),
                 atom("TIP3", "OH2", 9.0), atom("TIP3", "H1", 9.0), atom("SOD", "SOD", 9.5), atom("CLA", "CLA", 9.5),
                 atom("ALA", "CA", 4.0), atom("ALA", "HA", 4.0), atom("ALA", "N", 4.1)]
        groups = mf.aa_groups_from_atoms(atoms, {"POPC", "CHL1", "PSM", "GLPA"})
        self.assertEqual(len(groups["headgroups"]), 3)   # POPC P, PSM P, GLPA NF; the choline N and C13 are not anchors
        self.assertEqual(len(groups["tails"]), 5)        # C22, C316, C5S, C10F, C8S; C1 and C13 are glycerol/choline
        self.assertEqual(len(groups["cholesterol"]), 2)  # O3 and C25 (its C25 is not mistaken for a chain carbon)
        self.assertEqual(len(groups["water"]), 1)
        self.assertEqual(len(groups["ions"]), 2)
        self.assertEqual(len(groups["protein"]), 2)
        self.assertEqual(groups["water"].shape, (1, 3))


class Hydration(unittest.TestCase):
    def bulk(self, rng, n=2000):
        """Water filling the two bulk slabs |z - midplane| > 3 nm (planes at +-2 nm, one nm margin)."""
        xyz = rng.uniform(0, 1, size=(n, 3)) * np.array([BOX[0], BOX[1], 2.0])
        xyz[:, 2] += np.where(rng.uniform(size=n) < 0.5, MIDPLANE + 3.0, MIDPLANE - 5.0)
        return xyz

    def test_bulk_water_alone_gives_zero_core_water_and_passes(self):
        out = mf.core_hydration(self.bulk(np.random.default_rng(5)), np.zeros((0, 3)), MIDPLANE, 0.8, BOX, thickness_nm=4.0)
        self.assertEqual(out["n_core_free"], 0)
        self.assertEqual(out["ratio_percent"], 0.0)
        self.assertEqual(out["status"], "PASS")
        self.assertGreater(out["bulk_density_nm3"], 4.0)
        self.assertFalse(out["transmembrane_water_path"])

    def test_water_in_the_core_is_a_high_percentage_and_fails(self):
        rng = np.random.default_rng(6)
        core = rng.uniform(0, 1, size=(400, 3)) * np.array([BOX[0], BOX[1], 1.6]) + np.array([0, 0, MIDPLANE - 0.8])
        out = mf.core_hydration(np.vstack([self.bulk(rng), core]), np.zeros((0, 3)), MIDPLANE, 0.8, BOX, thickness_nm=4.0)
        self.assertEqual(out["n_core_free"], 400)
        self.assertGreater(out["ratio_percent"], 20.0)
        self.assertEqual(out["status"], "FAIL")
        self.assertEqual(out["records"][0]["status"], "FAIL")

    def test_waters_near_the_protein_are_excluded(self):
        rng = np.random.default_rng(7)
        protein = np.array([[4.0, 4.0, MIDPLANE]])
        near = protein + rng.normal(scale=0.15, size=(20, 3))        # within 0.5 nm of the protein atom
        far = np.array([[1.0, 1.0, MIDPLANE], [7.0, 7.0, MIDPLANE + 0.3]])
        out = mf.core_hydration(np.vstack([self.bulk(rng), near, far]), protein, MIDPLANE, 0.8, BOX, thickness_nm=4.0, protein_exclusion_nm=0.5)
        self.assertEqual(out["n_protein_associated"], 20)
        self.assertEqual(out["n_core_free"], 2)

    def test_a_water_column_spanning_the_core_is_a_path(self):
        rng = np.random.default_rng(8)
        column = np.column_stack([np.full(14, 2.0), np.full(14, 2.0), MIDPLANE + np.linspace(-1.1, 1.1, 14)])  # 0.17 nm apart
        out = mf.core_hydration(np.vstack([self.bulk(rng), column]), np.zeros((0, 3)), MIDPLANE, 0.8, BOX, thickness_nm=4.0)
        self.assertTrue(out["transmembrane_water_path"])
        self.assertEqual(out["records"][1]["status"], "FAIL")
        self.assertGreaterEqual(max(out["component_sizes"]), 12)

    def test_scattered_isolated_core_waters_are_not_a_path(self):
        rng = np.random.default_rng(9)
        scattered = np.array([[1.0, 1.0, MIDPLANE - 0.7], [3.0, 5.0, MIDPLANE - 0.3], [6.0, 2.0, MIDPLANE + 0.1], [2.0, 7.0, MIDPLANE + 0.5],
                              [5.0, 5.0, MIDPLANE + 0.75]])
        out = mf.core_hydration(np.vstack([self.bulk(rng), scattered]), np.zeros((0, 3)), MIDPLANE, 0.8, BOX, thickness_nm=4.0)
        self.assertFalse(out["transmembrane_water_path"])
        self.assertEqual(out["records"][1]["status"], "PASS")
        self.assertEqual(out["n_core_free"], 5)
        self.assertTrue(all(s == 1 for s in out["component_sizes"]))

    def test_the_bulk_region_can_be_given_explicitly(self):
        rng = np.random.default_rng(10)
        out = mf.core_hydration(self.bulk(rng), np.zeros((0, 3)), MIDPLANE, 0.8, BOX, bulk_region=(MIDPLANE + 3.0, MIDPLANE + 5.0))
        self.assertGreater(out["n_bulk"], 800)
        self.assertEqual(out["status"], "PASS")


class Interdigitation(unittest.TestCase):
    def test_separated_blocks_do_not_overlap_and_identical_blocks_overlap_fully(self):
        rng = np.random.default_rng(11)
        upper = rng.uniform(0, 1, size=(2000, 3)) * np.array([8.0, 8.0, 1.0]) + np.array([0, 0, MIDPLANE + 0.2])
        lower = rng.uniform(0, 1, size=(2000, 3)) * np.array([8.0, 8.0, 1.0]) + np.array([0, 0, MIDPLANE - 1.2])
        apart = mf.tail_interdigitation(upper, lower, MIDPLANE)
        self.assertLess(apart["overlap"], 0.02)
        same = mf.tail_interdigitation(upper, upper.copy(), MIDPLANE)
        self.assertAlmostEqual(same["overlap"], 1.0, places=6)
        self.assertEqual(len(same["z_nm"]), len(same["rho_upper"]))
        partial = mf.tail_interdigitation(upper, upper - np.array([0, 0, 0.5]), MIDPLANE)
        self.assertTrue(0.3 < partial["overlap"] < 0.7)

    def test_interdigitation_change_grades_or_asks_for_a_reference(self):
        self.assertEqual(mf.interdigitation_change(None, 0.3)["status"], "REFERENCE NEEDED")
        self.assertEqual(mf.interdigitation_change(0.30, 0.315)["status"], "PASS")
        self.assertEqual(mf.interdigitation_change(0.30, 0.345)["status"], "WARNING")
        self.assertEqual(mf.interdigitation_change(0.30, 0.40)["status"], "FAIL")
        self.assertAlmostEqual(mf.interdigitation_change(0.30, 0.33)["relative_percent"], 10.0, places=6)

    def test_aa_tail_split_assigns_by_the_anchor_not_the_tail_atom(self):
        def residue(resname, p_z, tail_z):
            return [{"resname": resname, "atom": "P", "x": 1.0, "y": 1.0, "z": p_z},
                    {"resname": resname, "atom": "C22", "x": 1.0, "y": 1.0, "z": tail_z},
                    {"resname": resname, "atom": "C316", "x": 1.0, "y": 1.0, "z": tail_z - 0.1},
                    {"resname": resname, "atom": "H2R", "x": 1.0, "y": 1.0, "z": tail_z}]
        chol = [{"resname": "CHL1", "atom": "O3", "x": 1, "y": 1, "z": 5.5}, {"resname": "CHL1", "atom": "C25", "x": 1, "y": 1, "z": 4.3}]
        upper, lower = mf.aa_leaflet_tail_split([residue("POPC", 6.0, MIDPLANE - 0.3), residue("DOPC", 2.0, MIDPLANE + 0.2), chol], MIDPLANE)
        self.assertEqual(upper.shape, (2, 3))   # the POPC tails dipping below the midplane stay with the upper leaflet
        self.assertEqual(lower.shape, (2, 3))
        self.assertLess(upper[:, 2].max(), MIDPLANE)


class Asymmetry(unittest.TestCase):
    def test_mole_fractions_sum_to_one_per_leaflet(self):
        chol = [{"index": 0, "resid": 0, "resname": "CHOL", "leaflet": "lower", "x_nm": 0.3 + 0.1 * k, "y_nm": 0.3, "z_nm": 2.5,
                 "z_range_nm": (2.5, 4.0)} for k in range(16)]
        out = mf.leaflet_asymmetry(stage(6.0, 2.0, extra=chol))
        self.assertEqual(out["n"], {"lower": 80, "upper": 64})
        self.assertAlmostEqual(out["ratio_upper_over_lower"], 0.8, places=9)
        self.assertAlmostEqual(sum(s["x_lower"] for s in out["species"].values()), 1.0, places=9)
        self.assertAlmostEqual(sum(s["x_upper"] for s in out["species"].values()), 1.0, places=9)
        self.assertEqual(out["species"]["CHOL"]["n_upper"], 0)
        self.assertAlmostEqual(out["species"]["CHOL"]["x_lower"], 0.2, places=9)



class LayeredOrder(unittest.TestCase):
    def profiles(self, head_z, tail_z, water_z):
        """Groups of points at given |z| offsets on both sides of a midplane at z = 5 in a 10 nm box."""
        rng = np.random.default_rng(0)
        make = lambda offsets: np.array([[rng.uniform(0, 6), rng.uniform(0, 6), 5.0 + sign * o]
                                         for o in offsets for sign in (-1, 1) for _ in range(50)])
        groups = {"headgroups": make([head_z]), "tails": make([tail_z]), "water": make([water_z]), "cholesterol": make([0.5 * (head_z + tail_z)])}
        return mf.z_density_profiles(groups, [6.0, 6.0, 10.0])

    def test_a_bilayer_architecture_passes_and_an_inverted_one_fails(self):
        self.assertEqual(mf.layered_order(self.profiles(2.0, 0.7, 3.5), 5.0)["status"], "PASS")
        self.assertEqual(mf.layered_order(self.profiles(0.7, 2.0, 3.5), 5.0)["status"], "FAIL")       # heads inside the tails

    def test_misplaced_cholesterol_warns(self):
        rng = np.random.default_rng(1)
        make = lambda o: np.array([[rng.uniform(0, 6), rng.uniform(0, 6), 5.0 + sign * o] for sign in (-1, 1) for _ in range(50)])
        groups = {"headgroups": make(2.0), "tails": make(0.7), "water": make(3.5), "cholesterol": make(3.2)}
        self.assertEqual(mf.layered_order(mf.z_density_profiles(groups, [6.0, 6.0, 10.0]), 5.0)["status"], "WARNING")


if __name__ == "__main__":
    unittest.main()
