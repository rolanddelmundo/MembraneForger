"""Ion numbers: the added NaCl concentration refers to the water, and neutralizing counterions are counted apart."""
import math
import unittest

from .common import failure, mf

# The 38527384 build: 65217 waters kept by solvate in a 2755.6 nm^3 box, system charge -217 e.
WATERS, BOX_NM3, CHARGE = 65217, 2755.6, -217


class SaltVolume(unittest.TestCase):
    def test_water_count_converts_to_the_spc216_volume(self):
        self.assertAlmostEqual(mf.WATERS_PER_NM3, 216 / 1.86206 ** 3, places=9)
        self.assertAlmostEqual(mf.WATERS_PER_NM3, 33.46, places=2)
        self.assertAlmostEqual(mf.salt_volume("water", [10.0, 10.0, 10.0], WATERS), 1949.3, places=1)

    def test_box_basis_is_the_cell_volume_and_a_bad_basis_is_refused(self):
        self.assertAlmostEqual(mf.salt_volume("box", [10.0, 10.0, 27.556], WATERS), BOX_NM3, places=6)
        self.assertIn("water", failure(mf.salt_volume, "total", [1.0, 1.0, 1.0], 10))

    def test_the_default_measures_the_water(self):
        self.assertEqual(mf.Settings().salt_volume, "water")


class IonCounts(unittest.TestCase):
    def additive(self, volume, charge):
        return mf.ion_counts(volume, charge, mode="additive")[:3]

    def test_box_and_additive_reproduce_the_earlier_rule(self):
        box = mf.salt_volume("box", [10.0, 10.0, 27.556], WATERS)
        self.assertEqual(self.additive(box, CHARGE), (249, 466, 249))      # as logged by the Oct 9 and Oct 17 builds
        self.assertEqual(self.additive(2737.8, -183), (247, 430, 247))

    def test_water_and_additive_on_the_cluster_build(self):
        water = mf.salt_volume("water", [10.0, 10.0, 27.556], WATERS)
        self.assertEqual(self.additive(water, CHARGE), (176, 393, 176))    # 0.15 M of the 1949 nm^3 of water
        self.assertEqual(self.additive(64784 / mf.WATERS_PER_NM3, -183), (175, 358, 175))

    def test_additive_counterions_follow_the_charge_sign_and_never_the_salt(self):
        for charge in (-217, -1, 0, 1, 5):
            pairs, pos, neg = self.additive(1949.3, charge)
            self.assertEqual(pairs, 176)
            self.assertEqual(pos - neg, -charge)
            self.assertEqual(min(pos, neg), pairs)

    def test_inclusive_negative_charge_the_counterion_is_the_target(self):
        target, pos, neg, short = mf.ion_counts(1949.3, -20, mode="inclusive")
        self.assertEqual((target, pos, neg, short), (176, 176, 156, 0))

    def test_inclusive_positive_charge_the_counterion_is_the_target(self):
        target, pos, neg, short = mf.ion_counts(1949.3, 20, mode="inclusive")
        self.assertEqual((target, pos, neg, short), (176, 156, 176, 0))

    def test_inclusive_neutral_system_is_the_plain_salt(self):
        self.assertEqual(mf.ion_counts(1949.3, 0, mode="inclusive"), (176, 176, 176, 0))

    def test_inclusive_always_neutralizes(self):
        for charge in (-217, -176, -175, -1, 0, 1, 175, 176, 217):
            _, pos, neg, _ = mf.ion_counts(1949.3, charge, mode="inclusive")
            self.assertEqual(pos - neg, -charge)
            self.assertGreaterEqual(min(pos, neg), 0)

    def test_inclusive_clamps_at_zero_when_the_charge_exceeds_the_target(self):
        target, pos, neg, short = mf.ion_counts(1949.3, CHARGE, mode="inclusive")   # the GPR9 build
        self.assertEqual((target, pos, neg, short), (176, 217, 0, 41))
        self.assertAlmostEqual(mf.salt_summary(WATERS, pos, neg)["na_M"], 0.185, places=3)
        _, pos, neg, short = mf.ion_counts(1949.3, 176, mode="inclusive")           # exactly the target: no shortfall
        self.assertEqual((pos, neg, short), (0, 176, 0))

    def test_the_water_needed_for_the_target_and_the_box_that_gives_it(self):
        needed = mf.water_volume_for_salt(CHARGE)
        self.assertAlmostEqual(needed, 217 / (0.15 * mf.MOLAR_PER_NM3), places=6)
        self.assertEqual(mf.ion_counts(needed, CHARGE, mode="inclusive")[3], 0)    # at that volume the shortfall is gone
        self.assertEqual(mf.ion_counts(needed, CHARGE, mode="inclusive")[1:3], (217, 0))

    def test_the_warning_states_the_shortfall_and_the_taller_box(self):
        box = [15.8, 9.6, 18.17]
        target, pos, neg, short = mf.ion_counts(1949.3, CHARGE, mode="inclusive")
        text = mf.shortfall_warning("water", 1949.3, box, CHARGE, target, short)
        self.assertTrue(text.startswith("WARNING"))
        self.assertIn("41 ions short of the 176 requested", text)
        dz = (mf.water_volume_for_salt(CHARGE) - 1949.3) / (box[0] * box[1])
        self.assertIn(f"{dz:.2f} nm taller", text)
        self.assertGreater(dz, 2.0)

    def test_a_bad_mode_is_refused_and_the_defaults(self):
        self.assertIn("counterion_mode", failure(mf.ion_counts, 1000.0, -1, mode="both"))
        self.assertEqual(mf.Settings().counterion_mode, "inclusive")

    def test_rounding_is_within_half_a_pair_of_the_target(self):
        for nwater in range(1000, 90000, 777):
            volume = nwater / mf.WATERS_PER_NM3
            target = mf.ion_counts(volume, 0)[0]
            self.assertLessEqual(abs(target - 0.15 * mf.MOLAR_PER_NM3 * volume), 0.5 + 1e-9)

    def test_no_volume_no_salt(self):
        self.assertEqual(mf.ion_counts(0.0, -3, mode="additive"), (0, 3, 0, 0))
        self.assertEqual(mf.ion_counts(0.0, -3, mode="inclusive"), (0, 3, 0, 3))


class SaltSummary(unittest.TestCase):
    def test_added_salt_counterions_and_ionic_strength_are_reported_apart(self):
        summary = mf.salt_summary(WATERS, 393, 176)
        self.assertAlmostEqual(summary["nacl_added_M"], 0.150, places=3)
        self.assertAlmostEqual(summary["cl_M"], 0.150, places=3)
        self.assertAlmostEqual(summary["na_M"], 0.335, places=3)           # Cl- plus 217 neutralizing Na+
        self.assertAlmostEqual(summary["ionic_strength_M"], 0.242, places=3)
        self.assertGreater(summary["ionic_strength_M"], summary["nacl_added_M"])

    def test_the_earlier_rule_put_about_0_21_m_in_the_water(self):
        summary = mf.salt_summary(WATERS, 466, 249)
        self.assertAlmostEqual(summary["nacl_added_M"], 0.212, places=3)

    def test_a_neutral_system_has_only_the_salt(self):
        pairs, pos, neg, _ = mf.ion_counts(WATERS / mf.WATERS_PER_NM3, 0)
        summary = mf.salt_summary(WATERS, pos, neg)
        self.assertAlmostEqual(summary["ionic_strength_M"], summary["nacl_added_M"], places=9)
        self.assertLess(abs(summary["nacl_added_M"] - 0.15), 0.5 / (WATERS / mf.WATERS_PER_NM3 * mf.MOLAR_PER_NM3))

    def test_ions_replace_waters_one_for_one(self):
        pairs, pos, neg, _ = mf.ion_counts(WATERS / mf.WATERS_PER_NM3, CHARGE, mode="additive")
        self.assertEqual(WATERS - pos - neg, 64648)                        # TIP3 left for the topology, never negative
        self.assertTrue(math.isfinite(mf.salt_summary(WATERS, pos, neg)["na_M"]))


class Census(unittest.TestCase):
    names = ["N", "CA"] + ["OH2", "H1", "H2"] * 4 + ["SOD"] * 3 + ["CLA"] * 2

    def test_counts_come_from_the_atom_names(self):
        self.assertEqual(mf.ion_census(self.names), {"water": 4, "na": 3, "cl": 2})

    def test_gro_and_topology_counts_agree(self):
        rows = [("PROA", 1), ("TIP3", 4), ("SOD", 3), ("CLA", 2)]
        self.assertEqual(mf.check_ion_counts(mf.ion_census(self.names), rows), [])

    def test_a_count_that_differs_is_named(self):
        rows = [("PROA", 1), ("TIP3", 4), ("SOD", 4), ("CLA", 2)]
        self.assertEqual(mf.check_ion_counts(mf.ion_census(self.names), rows), ["SOD: 3 in the structure, 4 in the topology"])
        rows = [("PROA", 1), ("TIP3", 5), ("SOD", 3)]
        problems = mf.check_ion_counts(mf.ion_census(self.names), rows)
        self.assertEqual(len(problems), 2)

    def test_the_salt_recomputed_from_the_final_structure_matches_the_request(self):
        pairs, pos, neg, _ = mf.ion_counts(WATERS / mf.WATERS_PER_NM3, -20, mode="inclusive")
        names = ["OH2", "H1", "H2"] * (WATERS - pos - neg) + ["SOD"] * pos + ["CLA"] * neg
        census = mf.ion_census(names)
        salt = mf.salt_summary(census["water"] + census["na"] + census["cl"], census["na"], census["cl"])
        self.assertLess(abs(salt["na_M"] - 0.15), 0.5 / (WATERS / mf.WATERS_PER_NM3 * mf.MOLAR_PER_NM3))   # within one ion
        self.assertAlmostEqual(salt["na_M"] - salt["cl_M"], 20 / (WATERS / mf.WATERS_PER_NM3 * mf.MOLAR_PER_NM3), places=9)


if __name__ == "__main__":
    unittest.main()
