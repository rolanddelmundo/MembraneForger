"""Trajectory analysis on synthetic data: readers, convergence, order parameters, lateral diffusion, K_A, leaflet tension."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from .common import mf

XVG = """# gmx energy
@    title "GROMACS Energies"
@    xaxis  label "Time (ps)"
@    yaxis  label "(nm)"
@ s0 legend "Box-X"
@ s1 legend "Box-Y"
@ s2 legend "Box-Z"
    0.000000   10.000000   10.000000    8.000000
   10.000000   10.100000    9.900000    8.100000
   20.000000   10.200000    9.800000    8.200000
"""


class Readers(unittest.TestCase):
    def test_read_xvg_parses_legends_and_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "box.xvg"
            path.write_text(XVG)
            xvg = mf.read_xvg(path)
            self.assertEqual(xvg["legends"], ["Box-X", "Box-Y", "Box-Z"])
            self.assertEqual(xvg["data"].shape, (3, 4))
            self.assertEqual(xvg["labels"]["xaxis"], "Time (ps)")
            box = mf.box_series_from_xvg(path)
            self.assertTrue(np.allclose(box["time_ps"], [0.0, 10.0, 20.0]))
            self.assertTrue(np.allclose(box["area_nm2"], [100.0, 99.99, 99.96]))
            self.assertAlmostEqual(box["volume_nm3"][0], 800.0)

    def test_read_stress_profile_parses_two_and_three_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "stress.dat"
            path.write_text("# z pi err\n-1.0 5.0 0.5\n0.0 0.0 0.5  # midplane\n1.0 -5.0 0.5\n")
            profile = mf.read_stress_profile(path)
            self.assertTrue(np.allclose(profile["z_nm"], [-1.0, 0.0, 1.0]) and np.allclose(profile["error_bar"], 0.5))
            path.write_text("-1.0 5.0\n1.0 -5.0\n")
            self.assertIsNone(mf.read_stress_profile(path)["error_bar"])

    def test_trajectory_frames_without_mdanalysis_is_a_clear_refusal(self):
        try:
            import MDAnalysis  # noqa: F401
        except ImportError:
            with self.assertRaises(SystemExit) as ctx:
                mf.trajectory_frames(Path("top.gro"), Path("traj.xtc"))
            self.assertIn("pip install MDAnalysis", str(ctx.exception))


class Convergence(unittest.TestCase):
    def test_a_flat_noisy_series_passes(self):
        rng = np.random.default_rng(1)
        t = np.arange(1000.0)
        result = mf.convergence(t, 64.0 + 0.3 * rng.standard_normal(1000))
        self.assertEqual(result["status"], "PASS", result["reason"])
        self.assertFalse(result["monotonic"])
        self.assertTrue(np.isfinite(result["sem"]) and result["sem"] < 0.1)

    def test_a_ten_percent_ramp_fails_with_the_monotonic_flag(self):
        t = np.arange(1000.0)
        result = mf.convergence(t, 60.0 + 6.0 * t / t[-1])
        self.assertEqual(result["status"], "FAIL")
        self.assertTrue(result["monotonic"])
        self.assertGreater(result["drift_percent"], 1.0)

    def test_a_short_series_is_insufficient_sampling(self):
        self.assertEqual(mf.convergence(np.arange(6.0), np.ones(6))["status"], "INSUFFICIENT SAMPLING")

    def test_record_uses_the_qc_metric_shape(self):
        result = mf.convergence(np.arange(100.0), np.full(100, 3.9))
        record = mf.convergence_record("thickness", "equilibration", result, "nm")
        self.assertEqual((record["metric"], record["value"], record["status"]), ("thickness", 3.9, "PASS"))

    def test_block_sem_of_white_noise(self):
        rng = np.random.default_rng(2)
        sem = mf.block_sem(rng.standard_normal(10000), blocks=10)
        self.assertTrue(0.002 < sem < 0.03)


class OrderParameters(unittest.TestCase):
    def test_along_z_is_one_and_perpendicular_is_minus_a_half(self):
        along = np.zeros((50, 2, 2, 3))
        along[..., 2] = 1.0
        self.assertTrue(np.allclose(mf.order_parameters(along), 1.0))
        across = np.zeros((50, 3))
        across[:, 0] = 0.109
        self.assertAlmostEqual(mf.order_parameters([across])[0], -0.5)

    def test_isotropic_vectors_average_to_zero(self):
        rng = np.random.default_rng(3)
        v = rng.standard_normal((20000, 3))
        self.assertLess(abs(mf.order_parameters([v])[0]), 0.05)

    def test_profile_from_residues_uses_the_hydrogens_that_exist(self):
        def residue(scale):
            atoms = [{"atom": "C22", "x": 0.0, "y": 0.0, "z": 0.0}, {"atom": "H2R", "x": 0.0, "y": 0.0, "z": 0.109},
                     {"atom": "H2S", "x": 0.0, "y": 0.0, "z": -0.109},
                     {"atom": "C29", "x": 0.0, "y": 0.0, "z": 1.0}, {"atom": "H91", "x": 0.109 * scale, "y": 0.0, "z": 1.0}]
            return atoms
        profile = mf.scd_profile_from_residues([residue(1.0), residue(2.0)], mf.CHARMM36_CHAINS)
        sn2 = profile["sn2"]
        self.assertEqual(sn2["carbon"], ["C22", "C29"])
        self.assertEqual(sn2["index"], [2, 9])
        self.assertTrue(np.allclose(sn2["scd"], [1.0, -0.5]) and sn2["n"] == [2, 2])
        self.assertEqual(profile["sn1"]["carbon"], [])

    def test_profile_rmse_grading(self):
        reference = np.full(14, 0.20)
        self.assertEqual(mf.scd_profile_rmse(reference, reference + 0.02)["status"], "PASS")
        self.assertEqual(mf.scd_profile_rmse(reference, reference - 0.04)["status"], "WARNING")
        self.assertEqual(mf.scd_profile_rmse(reference, reference + 0.10)["status"], "FAIL")
        self.assertEqual(mf.scd_profile_rmse(None, reference)["status"], "REFERENCE NEEDED")


class Diffusion(unittest.TestCase):
    def test_unwrap_undoes_periodic_jumps(self):
        rng = np.random.default_rng(4)
        box = np.array([6.0, 6.0])
        true = np.cumsum(0.3 * rng.standard_normal((400, 20, 2)), axis=0) + 3.0
        unwrapped = mf.unwrap_xy(np.mod(true, box), box)
        self.assertTrue(np.allclose(unwrapped - unwrapped[:1], true - true[:1]))
        self.assertTrue(np.allclose(mf.unwrap_xy(np.mod(true, box), np.tile(box, (400, 1))), unwrapped))

    def test_random_walk_recovers_d_and_alpha_one(self):
        rng = np.random.default_rng(5)
        d_nm2_per_ps, dt = 0.5e-5, 10.0         # 5e-6 nm^2/ps = 5 um^2/s, a fluid-phase lipid
        steps = np.sqrt(2.0 * d_nm2_per_ps * dt) * rng.standard_normal((2000, 200, 2))
        positions = np.cumsum(steps, axis=0)
        t = dt * np.arange(2000)
        msd = mf.msd_xy(positions, t)
        self.assertEqual(len(msd["lag_ps"]), 1000)
        fit = mf.diffusion_fit(msd["lag_ps"], msd["msd_nm2"])
        self.assertLess(abs(fit["alpha"] - 1.0), 0.1, fit)
        self.assertLess(abs(fit["d_um2_per_s"] - 5.0) / 5.0, 0.2, fit)
        self.assertAlmostEqual(fit["d_um2_per_s"], fit["d_1e8_cm2_per_s"])
        self.assertEqual(fit["status"], "PASS", fit["reason"])

    def test_direct_and_fft_msd_agree(self):
        rng = np.random.default_rng(6)
        p = np.cumsum(rng.standard_normal((60, 7, 2)), axis=0)
        t = np.arange(60.0)
        msd = mf.msd_xy(p, t, remove_drift=False)
        direct = [np.mean(((p[lag:] - p[:-lag]) ** 2).sum(axis=2)) for lag in range(1, 31)]
        self.assertTrue(np.allclose(msd["msd_nm2"], direct))

    def test_ballistic_motion_is_not_diffusive(self):
        t = np.arange(500.0)
        velocity = np.array([[0.01, -0.02], [0.03, 0.01], [-0.01, 0.0]])
        positions = t[:, None, None] * velocity[None]
        msd = mf.msd_xy(positions, t, remove_drift=False)
        fit = mf.diffusion_fit(msd["lag_ps"], msd["msd_nm2"])
        self.assertLess(abs(fit["alpha"] - 2.0), 0.05)
        self.assertEqual(fit["status"], "INSUFFICIENT SAMPLING")
        self.assertIn("not in a diffusive regime", fit["reason"])

    def test_diffusion_against_a_reference(self):
        self.assertEqual(mf.diffusion_vs_reference(6.0, 5.0)["status"], "PASS")
        self.assertEqual(mf.diffusion_vs_reference(1.5, 5.0)["status"], "WARNING")
        self.assertEqual(mf.diffusion_vs_reference(40.0, 5.0)["status"], "FAIL")
        self.assertEqual(mf.diffusion_vs_reference(6.0, None)["status"], "REFERENCE NEEDED")


class Compressibility(unittest.TestCase):
    def test_gaussian_area_series_recovers_ka(self):
        rng = np.random.default_rng(7)
        mean, sd, temperature = 40.0, 0.25, 310.0
        areas = mean + sd * rng.standard_normal(5000)
        expected = mf.KB_J_PER_K * temperature * mean * 1e-18 / (sd ** 2 * 1e-36) * 1000.0  # mN/m
        result = mf.area_compressibility(areas, temperature)
        self.assertLess(abs(result["ka_mn_per_m"] - expected) / expected, 0.15, result)
        self.assertTrue(np.isfinite(result["uncertainty_mn_per_m"]) and result["uncertainty_mn_per_m"] > 0)
        self.assertIn(result["status"], ("PASS", "WARNING"))
        self.assertEqual(mf.area_compressibility(areas[:6], temperature)["status"], "INSUFFICIENT SAMPLING")

    def test_ka_against_a_reference(self):
        self.assertEqual(mf.ka_vs_reference(275.0, 250.0)["status"], "PASS")
        self.assertEqual(mf.ka_vs_reference(325.0, 250.0)["status"], "WARNING")
        self.assertEqual(mf.ka_vs_reference(375.0, 250.0)["status"], "FAIL")
        self.assertEqual(mf.ka_vs_reference(375.0, None)["status"], "REFERENCE NEEDED")


class Tension(unittest.TestCase):
    Z = np.arange(-2.95, 3.0, 0.1)

    def test_antisymmetric_profile_gives_opposite_leaflet_tensions(self):
        pi = 100.0 * np.sin(np.pi * self.Z / 3.0)
        result = mf.leaflet_tension(self.Z, pi, 0.0)
        self.assertAlmostEqual(result["gamma_upper_mn_per_m"], -result["gamma_lower_mn_per_m"], places=6)
        self.assertAlmostEqual(result["delta_mn_per_m"], 2.0 * result["gamma_upper_mn_per_m"], places=6)
        self.assertEqual(result["status"], "FAIL")

    def test_symmetric_profile_has_no_difference(self):
        pi = 100.0 * np.cos(np.pi * self.Z / 3.0) - 100.0 * np.cos(2.0 * np.pi * self.Z / 3.0)
        result = mf.leaflet_tension(self.Z, pi, 0.0)
        self.assertLess(abs(result["delta_mn_per_m"]), 1e-6)
        self.assertEqual(result["status"], "PASS" if abs(result["gamma_upper_mn_per_m"]) <= 10.0 else "WARNING")

    def test_ten_bar_over_one_nanometre_is_one_millinewton_per_metre(self):
        z = np.arange(0.05, 1.0, 0.1)
        result = mf.leaflet_tension(np.concatenate([-z, z]), np.concatenate([np.zeros(10), np.full(10, 10.0)]), 0.0)
        self.assertAlmostEqual(result["gamma_upper_mn_per_m"], 1.0)
        self.assertAlmostEqual(result["gamma_lower_mn_per_m"], 0.0)

    def test_a_large_but_uncertain_difference_is_insufficient_sampling(self):
        z = np.arange(0.05, 1.0, 0.1)
        pi = np.concatenate([np.zeros(10), np.full(10, 150.0)])        # delta = 15 mN/m
        noisy = mf.leaflet_tension(np.concatenate([-z, z]), pi, 0.0, uncertainty_bar=np.full(20, 600.0))
        self.assertEqual(noisy["status"], "INSUFFICIENT SAMPLING")
        self.assertEqual(mf.leaflet_tension(np.concatenate([-z, z]), pi, 0.0, uncertainty_bar=np.full(20, 1.0))["status"], "FAIL")


if __name__ == "__main__":
    unittest.main()
