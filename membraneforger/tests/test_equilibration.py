"""Gates 3 and 4 from files: box series, stress profile, report writing, and the honest statuses when data are missing."""
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np

from .common import mf


def write_xvg(path: Path, time, x, y, z):
    """A gmx energy box series."""
    lines = ["# gmx energy", '@    title "Energy"', '@ s0 legend "Box-X"', '@ s1 legend "Box-Y"', '@ s2 legend "Box-Z"']
    lines += [f"{t:10.3f} {a:10.5f} {b:10.5f} {c:10.5f}" for t, a, b, c in zip(time, x, y, z)]
    path.write_text("\n".join(lines) + "\n")


class EquilibrationFromFiles(unittest.TestCase):
    def setUp(self):
        self.out = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.out, True)
        rng = np.random.default_rng(0)
        self.time = np.arange(0, 20000.0, 10.0)                                   # 20 ns, 10 ps frames
        n = len(self.time)
        self.flat = 10.0 + 0.02 * rng.normal(size=n)                                # a converged box edge (nm)
        self.ramp = np.linspace(10.0, 9.0, n) + 0.01 * rng.normal(size=n)          # a shrinking box edge
        (self.out / "membrane_validation.json").write_text(json.dumps({"stages": [{"name": "slice"}], "records": [], "gates": {}}))
        (self.out / "membrane_validation.md").write_text("# x\n## 3. Equilibration convergence (gate 3)\n"
                                                         "## 4. Equilibrium physical validation (gate 4)\n")

    def test_a_converged_box_series_passes_gate_3_and_ka_needs_a_reference(self):
        write_xvg(self.out / "ener.xvg", self.time, self.flat, self.flat, 12.0 + 0.0 * self.flat)
        code = mf.equilibration.main(["--out", str(self.out), "--energy", str(self.out / "ener.xvg")]) if hasattr(mf, "equilibration") else None
        if code is None:
            from membraneforger import equilibration
            code = equilibration.main(["--out", str(self.out), "--energy", str(self.out / "ener.xvg")])
        self.assertEqual(code, 0)
        result = json.loads((self.out / "equilibration_validation.json").read_text())
        by = {(r["metric"], r["stage"]): r for r in result["records"]}
        self.assertEqual(by[("box area", "equilibrated")]["status"], "PASS")
        self.assertEqual(by[("box volume", "equilibrated")]["status"], "PASS")
        self.assertEqual(by[("area compressibility K_A", "final")]["status"], "REFERENCE NEEDED")   # no --reference-ka
        self.assertGreater(by[("area compressibility K_A", "final")]["value"], 0.0)
        self.assertEqual(by[("APL", "equilibrated")]["status"], "NOT RUN")                           # no frames given
        self.assertEqual(by[("leaflet differential tension", "final")]["status"], "NOT RUN")
        merged = json.loads((self.out / "membrane_validation.json").read_text())
        self.assertIn("3 equilibration convergence", merged["gates"])
        self.assertIn("| box area |", (self.out / "membrane_validation.md").read_text())

    def test_a_shrinking_box_fails_gate_3(self):
        from membraneforger import equilibration
        write_xvg(self.out / "ener.xvg", self.time, self.ramp, self.ramp, 12.0 + 0.0 * self.ramp)
        equilibration.main(["--out", str(self.out), "--energy", str(self.out / "ener.xvg")])
        result = json.loads((self.out / "equilibration_validation.json").read_text())
        area = next(r for r in result["records"] if r["metric"] == "box area")
        self.assertEqual(area["status"], "FAIL")
        self.assertIn("monotonic", area["reason"].lower())
        self.assertEqual(result["gates"]["3 equilibration convergence"]["status"], "FAIL")

    def test_a_stress_profile_gives_leaflet_tensions(self):
        from membraneforger import equilibration
        z = np.linspace(-3.0, 3.0, 61)
        pi = 150.0 * np.sign(z) * np.exp(-((np.abs(z) - 1.5) ** 2) / 0.2)           # antisymmetric: equal and opposite leaflets
        (self.out / "stress.dat").write_text("# z pi\n" + "\n".join(f"{a:.3f} {b:.4f}" for a, b in zip(z, pi)) + "\n")
        equilibration.main(["--out", str(self.out), "--stress", str(self.out / "stress.dat")])
        result = json.loads((self.out / "equilibration_validation.json").read_text())
        delta = next(r for r in result["records"] if r["metric"] == "leaflet differential tension")
        upper = next(r for r in result["records"] if r["metric"] == "leaflet tension" and r["leaflet"] == "upper")
        lower = next(r for r in result["records"] if r["metric"] == "leaflet tension" and r["leaflet"] == "lower")
        self.assertAlmostEqual(upper["value"], -lower["value"], places=3)
        self.assertAlmostEqual(delta["value"], upper["value"] - lower["value"], places=3)
        self.assertIn(delta["status"], ("FAIL", "WARNING", "INSUFFICIENT SAMPLING"))    # a large imbalance never passes


class OrientationDriftAndWaterPath(unittest.TestCase):
    def test_tilt_drift_is_graded_in_degrees(self):
        from membraneforger import equilibration
        time = np.arange(0, 10000.0, 10.0)
        steady = equilibration.absolute_drift_record("protein tilt", "equilibrated", mf.convergence(time, 70.0 + 0.2 * np.sin(time / 300.0)),
                                                     "deg", equilibration.ORIENTATION_DRIFT["tilt_deg"])
        self.assertEqual(steady["status"], "PASS")
        drifting = equilibration.absolute_drift_record("protein tilt", "equilibrated", mf.convergence(time, 70.0 + 0.002 * time), "deg",
                                                       equilibration.ORIENTATION_DRIFT["tilt_deg"])
        self.assertEqual(drifting["status"], "FAIL")                                            # 10 degrees over the final window, monotonic
        self.assertGreater(abs(drifting["deviation"]), 7.0)

    def test_water_path_persistence(self):
        from membraneforger import equilibration
        self.assertEqual(equilibration.water_path_persistence([False] * 10)["status"], "PASS")
        self.assertEqual(equilibration.water_path_persistence([True] + [False] * 9)["status"], "WARNING")
        self.assertEqual(equilibration.water_path_persistence([True] * 6 + [False] * 4)["status"], "FAIL")
        self.assertEqual(equilibration.water_path_persistence([None, None])["status"], "NOT RUN")


class CoreWaterRecord(unittest.TestCase):
    """Core hydration over frames: occupancy first, a trend only with enough frames and waters."""
    def record(self, counts, bulk=31.5, area=130.0):
        from membraneforger.equilibration import core_water_record
        counts = np.asarray(counts, dtype=float)
        return core_water_record(100.0 * counts / (area * 1.6) / bulk, counts)

    def test_a_few_frames_with_a_handful_of_waters_is_inconclusive_not_a_fail(self):
        record = self.record([5, 9, 5, 6, 4, 5, 3, 0, 4, 6, 5])                                   # the KOR config A series
        self.assertEqual(record["status"], "INSUFFICIENT SAMPLING")
        self.assertIn("water in 10 of 11 frames", record["reason"])

    def test_a_steady_well_sampled_core_passes_and_a_rising_one_warns(self):
        steady = self.record(np.full(40, 12.0) + np.tile([-1.0, 1.0], 20))
        self.assertEqual(steady["status"], "PASS")
        rising = self.record(np.linspace(8.0, 20.0, 40))
        self.assertEqual(rising["status"], "WARNING")

    def test_a_high_level_fails_whatever_the_sampling(self):
        self.assertEqual(self.record([400, 410, 390])["status"], "FAIL")

    def test_no_measurable_water_is_not_run(self):
        from membraneforger.equilibration import core_water_record
        self.assertEqual(core_water_record([np.nan, np.nan])["status"], "NOT RUN")

    def test_the_report_keeps_the_free_water_counts_the_trend_verdict_used(self):
        import json
        import tempfile
        from pathlib import Path

        from membraneforger.equilibration import write_equilibration_report
        keys = ("frame", "apl_upper_A2", "apl_lower_A2", "thickness_A", "core_water_percent", "core_water_free", "tilt_deg", "depth_A", "area_nm2")
        series = {k: [1.0, 2.0] for k in keys}
        series["core_water_free"] = [5.0, float("nan")]
        with tempfile.TemporaryDirectory() as tmp:
            write_equilibration_report(Path(tmp), [], series)
            saved = json.loads((Path(tmp) / "equilibration_validation.json").read_text())["series"]
        self.assertEqual(saved["core_water_free"], [5.0, None])


if __name__ == "__main__":
    unittest.main()


class PeriodicImageClearance(unittest.TestCase):
    def test_image_gap_of_a_slab_and_of_the_protein(self):
        from membraneforger.equilibration import image_gap
        protein = np.array([[1.0, 1.0, 5.0], [2.0, 1.0, 5.0], [0.5, 1.0, 9.0]])
        box = [6.0, 6.0, 12.0]
        full = image_gap(protein, box, None, mf.Settings())
        slab = image_gap(protein, box, 5.0, mf.Settings())
        self.assertAlmostEqual(full, 4.5, places=6)   # the atom at x = 0.5 is 4.5 nm from the image of the one at x = 2
        self.assertAlmostEqual(slab, 5.0, places=6)   # only the two atoms at the midplane count: 2 to 1 + 6


class ProfileChange(unittest.TestCase):
    def profile(self, length, shift):
        n = int(round(length / 0.1))
        z = (np.arange(n) + 0.5) * (length / n)
        return {"z_nm": z, "bin_nm": length / n, "density_nm3": {"water": 1.0 / (1.0 + np.exp(-(z - length / 2.0 - shift) * 4.0))}}

    def test_profiles_of_different_box_height_and_position_compare_after_shifting_to_the_midplane(self):
        from membraneforger.equilibration import profile_change
        first, last = self.profile(20.9, 0.0), self.profile(23.7, 1.0)  # taller box, membrane moved up by 1 nm
        change = profile_change(first, last, 20.9 / 2.0, 23.7 / 2.0 + 1.0)
        self.assertLess(change["water"], 0.05)  # the same shape relative to the midplane
        self.assertGreater(profile_change(first, last, 20.9 / 2.0, 23.7 / 2.0)["water"], 0.1)  # misplaced midplane: a real change


class ScdRecord(unittest.TestCase):
    def test_per_carbon_counts_become_one_count_per_chain(self):
        from membraneforger.equilibration import scd_record
        values = {"carbon": ["C22", "C23", "C24"], "index": [2, 3, 4], "scd": [0.20, 0.21, 0.19], "sem": [0.01] * 3, "n": [120, 118, 120]}
        record = scd_record("POPC", "sn1", values, None)
        self.assertEqual((record["n"], record["status"]), (120, "REFERENCE NEEDED"))
        self.assertAlmostEqual(record["value"], 0.20, places=6)
