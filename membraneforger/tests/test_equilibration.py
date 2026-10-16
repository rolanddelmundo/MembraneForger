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


if __name__ == "__main__":
    unittest.main()
