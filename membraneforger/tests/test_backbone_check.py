"""Backbone check after EM: peptide planes, centre handedness, phi/psi drift and strict segments, on synthetic chains of known geometry."""
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

from .common import failure, mf

BOX = np.array([6.0, 6.0, 6.0])   # nm


def place(a, b, c, bond, angle, torsion):
    """The point d with |cd| = bond, angle b-c-d and torsion a-b-c-d (degrees)."""
    angle, torsion = math.radians(angle), math.radians(torsion)
    bc = c - b
    bc = bc / np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n = n / np.linalg.norm(n)
    m = np.cross(n, bc)
    return c - bond * math.cos(angle) * bc + bond * math.sin(angle) * (math.cos(torsion) * m + math.sin(torsion) * n)


def volume(centre, a, b, c):
    return float((a - centre) @ np.cross(b - centre, c - centre))


def chain(names, angles=None, hand=None, thr_inverted=()):
    """A backbone of ideal bond lengths and angles (nm). angles[i] = (phi, psi, omega); omega joins residue i to i+1.
    hand[i] = +1 (L) or -1 (D). Residues named THR/ILE get OG1|CG1 and CG2 with the natural C-beta handedness unless listed in thr_inverted."""
    n = len(names)
    angles = angles or [(-63.0, -41.0, 180.0)] * n
    hand = hand or [1] * n
    N = np.zeros(3)
    CA = N + np.array([1.458, 0, 0])
    C = place(np.array([0, 1.0, 0]), N, CA, 1.525, 111.2, 180.0)
    res = []
    for i in range(n):
        phi, psi, omega = angles[i]
        atoms = {"N": N, "CA": CA, "C": C}
        cb = [place(C, N, CA, 1.53, 110.5, t) for t in (122.6, -122.6)]
        atoms["CB"] = next(x for x in cb if (volume(CA, N, C, x) > 0) == (hand[i] > 0))
        if names[i] in ("THR", "ILE"):
            og = place(N, CA, atoms["CB"], 1.43, 109.5, 60.0)
            cg = place(N, CA, atoms["CB"], 1.53, 109.5, -60.0)
            if (volume(atoms["CB"], CA, og, cg) > 0) == (i in thr_inverted):
                og, cg = cg, og
            atoms["OG1" if names[i] == "THR" else "CG1"], atoms["CG2"] = og, cg
        res.append({"name": names[i], "resid": i + 1, "atoms": atoms})
        if i + 1 < n:
            N2 = place(N, CA, C, 1.329, 116.2, psi)
            CA2 = place(CA, C, N2, 1.458, 121.7, omega)
            C2 = place(C, N2, CA2, 1.525, 111.2, angles[i + 1][0])
            N, CA, C = N2, CA2, C2
    for r in res:                                  # Angstrom -> nm
        r["atoms"] = {k: v / 10.0 for k, v in r["atoms"].items()}
    return {"box": BOX, "residues": res}


NAMES = ["ALA", "LEU", "ASP", "THR", "ALA", "PRO", "GLY", "ILE", "ALA", "SER"]


def check(before, after, **kw):
    return mf.check_backbone(before, after, kw.pop("labels", None), kw.pop("restrained", True), kw.pop("margin_deg", 30.0),
                             kw.pop("strict", ()), kw.pop("mode", "yes"))


class Omega(unittest.TestCase):
    def angles(self, **changes):
        a = [(-63.0, -41.0, 180.0)] * len(NAMES)
        a[5] = (-65.0, 150.0, 180.0)            # the proline: phi in its window
        for i, v in changes.items():
            a[int(i[1:])] = v
        return a

    def test_identical_structures_pass(self):
        s = chain(NAMES, self.angles())
        r = check(s, chain(NAMES, self.angles()))
        self.assertEqual(r["status"], "PASS")
        self.assertEqual(r["failures"], [])
        self.assertAlmostEqual(r["max_drift_deg"], 0.0, places=0)

    def test_a_trans_bond_that_em_twisted_fails(self):
        r = check(chain(NAMES, self.angles()), chain(NAMES, self.angles(i1=(-63.0, -41.0, 100.0))))
        self.assertEqual(r["status"], "FAIL")
        self.assertTrue(any("LEU2-ASP3" in f and "trans -> twisted" in f for f in r["failures"]), r["failures"])

    def test_report_mode_lists_the_same_as_warnings(self):
        r = check(chain(NAMES, self.angles()), chain(NAMES, self.angles(i1=(-63.0, -41.0, 100.0))), mode="report")
        self.assertEqual(r["status"], "PASS")
        self.assertTrue(any(w.startswith("(would fail)") for w in r["warnings"]))

    def test_a_bond_twisted_before_and_after_is_inherited_not_a_failure(self):
        twisted = self.angles(i1=(-63.0, -41.0, 100.0))
        r = check(chain(NAMES, twisted), chain(NAMES, twisted))
        self.assertEqual(r["status"], "PASS")
        self.assertEqual(r["omega_twisted_before_and_after"], 1)

    def test_a_twisted_bond_that_relaxes_is_no_failure(self):
        r = check(chain(NAMES, self.angles(i1=(-63.0, -41.0, 100.0))), chain(NAMES, self.angles()))
        self.assertEqual(r["status"], "PASS")

    def test_a_small_loss_of_planarity_only_warns(self):
        r = check(chain(NAMES, self.angles()), chain(NAMES, self.angles(i1=(-63.0, -41.0, 140.0))))
        self.assertEqual(r["status"], "PASS")
        self.assertTrue(any("below 150" in w for w in r["warnings"]))

    def test_xaa_pro_keeps_its_isomer(self):
        cis = self.angles(i4=(-63.0, -41.0, 0.0))                # ALA5-PRO6 cis
        self.assertEqual(check(chain(NAMES, cis), chain(NAMES, cis))["status"], "PASS")
        flipped = check(chain(NAMES, self.angles()), chain(NAMES, cis))
        self.assertEqual(flipped["status"], "FAIL")
        self.assertIn("trans -> cis", flipped["failures"][0])
        back = check(chain(NAMES, cis), chain(NAMES, self.angles()))
        self.assertEqual(back["status"], "FAIL")
        self.assertIn("cis -> trans", back["failures"][0])


class Chirality(unittest.TestCase):
    def test_an_inverted_calpha_fails(self):
        hand = [1] * len(NAMES)
        hand[1] = -1
        r = check(chain(NAMES), chain(NAMES, hand=hand))
        self.assertEqual(r["status"], "FAIL")
        self.assertTrue(any("LEU2 C-alpha L -> D" in f for f in r["failures"]), r["failures"])

    def test_a_d_residue_in_the_input_is_not_an_error(self):
        hand = [1] * len(NAMES)
        hand[1] = -1
        self.assertEqual(check(chain(NAMES, hand=hand), chain(NAMES, hand=hand))["status"], "PASS")

    def test_an_inverted_thr_or_ile_cbeta_fails(self):
        r = check(chain(NAMES), chain(NAMES, thr_inverted=(3,)))
        self.assertEqual(r["status"], "FAIL")
        self.assertTrue(any("THR4 C-beta inverted" in f for f in r["failures"]), r["failures"])
        self.assertEqual(check(chain(NAMES, thr_inverted=(3, 7)), chain(NAMES, thr_inverted=(3, 7)))["status"], "PASS")

    def test_glycine_and_unknown_residues_have_no_centre_to_invert(self):
        names = ["ALA", "GLY", "AIB", "ALA"]
        before, after = chain(names), chain(names, hand=[1, -1, 1, 1])      # glycine "flipped": nothing to flip
        for structure in (before, after):
            structure["residues"][2]["atoms"].pop("CB")                     # AIB has CB1 and CB2, no atom called CB
        self.assertEqual(check(before, after)["status"], "PASS")


class Drift(unittest.TestCase):
    def moved(self, delta):
        a = [(-63.0, -41.0, 180.0)] * len(NAMES)
        a[5] = (-65.0, 150.0, 180.0)
        b = list(a)
        b[2] = (a[2][0] - delta, a[2][1], 180.0)
        return chain(NAMES, a), chain(NAMES, b)

    def test_with_the_restraints_on_a_drift_beyond_the_restraint_width_plus_margin_fails(self):
        r = check(*self.moved(45.0))
        self.assertEqual(r["status"], "PASS")
        self.assertAlmostEqual(r["max_drift_deg"], 45.0, places=0)
        r = check(*self.moved(75.0))
        self.assertEqual(r["status"], "FAIL")
        self.assertIn("moved more than 60", r["failures"][0])
        self.assertEqual(check(*self.moved(75.0), margin_deg=50.0)["status"], "PASS")

    def test_without_the_restraints_a_drift_only_warns(self):
        r = check(*self.moved(75.0), restrained=False)
        self.assertEqual(r["status"], "PASS")
        self.assertTrue(any("no backbone restraints" in w for w in r["warnings"]))


class StrictSegments(unittest.TestCase):
    NAMES = ["ALA", "ASP", "ALA", "ALA", "ALA", "ALA"]
    LABELS = [("A", 9 + i, n) for i, n in enumerate(NAMES)]

    def pair(self):
        a = [(-63.0, -41.0, 180.0)] * len(self.NAMES)
        b = list(a)
        b[1] = (-11.0, -83.0, 180.0)             # Asp10 leaves the helix for the phi > -20 band, 32 degrees of phi
        return chain(self.NAMES, a), chain(self.NAMES, b)

    def test_a_modelled_residue_in_an_outlier_region_fails_although_the_drift_is_small(self):
        r = check(*self.pair(), labels=self.LABELS, strict=mf.parse_segments("A:9-13"))
        self.assertEqual(r["status"], "FAIL")
        self.assertIn("A:ASP10", r["failures"][0])

    def test_the_same_residue_outside_the_segment_only_warns(self):
        r = check(*self.pair(), labels=self.LABELS, strict=mf.parse_segments("B:9-13,A:20-30"))
        self.assertEqual(r["status"], "PASS")
        self.assertTrue(any("became Ramachandran outliers" in w for w in r["warnings"]))

    def test_segment_syntax(self):
        self.assertEqual(mf.parse_segments("R:338-345, A:5-9,B:12"), (("R", 338, 345), ("A", 5, 9), ("B", 12, 12)))
        self.assertEqual(mf.parse_segments(""), ())
        for bad in ("R338-345", "R:x-5", "R:9-5", ":1-2"):
            with self.assertRaises(ValueError, msg=bad):
                mf.parse_segments(bad)


class Pairing(unittest.TestCase):
    def test_residue_lists_that_do_not_correspond_are_skipped_not_failed(self):
        r = check(chain(NAMES), chain(NAMES[:-1]))
        self.assertEqual(r["status"], "skipped")
        r = check(chain(NAMES), chain(["ALA"] * len(NAMES)))
        self.assertEqual(r["status"], "skipped")

    def test_labels_of_the_wrong_length_fall_back_to_the_gro_numbers(self):
        r = check(chain(NAMES), chain(NAMES, hand=[1, -1] + [1] * 8), labels=[("A", 1, "ALA")])
        self.assertEqual(r["status"], "FAIL")
        self.assertIn("LEU2", r["failures"][0])
        self.assertEqual(r["labels"], "residue numbers of the .gro")


class GroReading(unittest.TestCase):
    def write(self, path, structure, wrap):
        atoms = []
        for r in structure["residues"]:
            for name, x in r["atoms"].items():
                x = np.mod(x + 2.9, BOX) if wrap else x
                atoms.append({"resid": r["resid"], "resname": r["name"], "atom": name, "x": float(x[0]), "y": float(x[1]), "z": float(x[2])})
        atoms.append({"resid": 99, "resname": "TIP3", "atom": "OH2", "x": 1.0, "y": 1.0, "z": 1.0})
        mf.write_gro(atoms, list(BOX), path, "synthetic")

    def test_a_chain_wrapped_across_the_periodic_box_gives_the_same_angles(self):
        s = chain(NAMES)
        with tempfile.TemporaryDirectory() as tmp:
            self.write(Path(tmp) / "whole.gro", s, False)
            self.write(Path(tmp) / "wrapped.gro", s, True)
            whole, wrapped = (mf.backbone_residues(Path(tmp) / f) for f in ("whole.gro", "wrapped.gro"))
        self.assertEqual(len(whole["residues"]), len(NAMES))      # water is not a protein residue
        r = check(whole, wrapped)
        self.assertEqual(r["status"], "PASS", r["failures"])
        self.assertLess(r["max_drift_deg"], 1.0)


class Config(unittest.TestCase):
    def test_defaults_and_generated_names(self):
        s = mf.Settings()
        self.assertEqual(s.backbone_check, "report")
        self.assertEqual(s.backbone_drift_margin_deg, 30.0)
        self.assertIn("backbone_check.json", mf.GENERATED)

    def test_the_stage_has_a_hint_that_names_the_report(self):
        self.assertIn("backbone", mf.pipeline.STAGE_HINTS)
        self.assertIn("backbone_check.json", mf.pipeline.STAGE_HINTS["backbone"])


class StageIntegration(unittest.TestCase):
    """review_backbone on files in an output directory, as finish_system calls it after validate_em."""

    def session(self, out, **settings):
        return mf.Session(out=out, name="t", gmx="gmx", forcefield=Path("."), data=None, python="python", ntomp=1, nsteps=1,
                          settings=mf.Settings(**settings))

    def files(self, out, after_hand=None):
        for name, hand in (("solv_ions.gro", None), ("em.unverified.gro", after_hand)):
            s = chain(NAMES, hand=hand)
            atoms = [{"resid": r["resid"], "resname": r["name"], "atom": k, "x": float(v[0]), "y": float(v[1]), "z": float(v[2])}
                     for r in s["residues"] for k, v in r["atoms"].items()]
            mf.write_gro(atoms, list(BOX), out / name, "synthetic")
        (out / "toppar").mkdir()
        (out / "toppar" / "PROA.itp").write_text("[ dihedral_restraints ]\n#ifdef BBDIHRES\n#endif\n")

    system = {"chains": {"A": [{"resid": i + 1, "resname": n, "atom": "CA"} for i, n in enumerate(NAMES)]}}

    def test_pass_writes_the_report_and_returns_a_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            self.files(out)
            summary = mf.review_backbone(self.session(out), self.system)
            self.assertEqual(summary["status"], "PASS")
            self.assertTrue(summary["restraints_active"])
            self.assertTrue((out / "backbone_check.json").is_file())
            self.assertIn("backbone after EM", (out / "membranebuilder.log").read_text())

    def test_a_failure_stops_the_build_with_the_residue_named(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            self.files(out, after_hand=[1, -1] + [1] * 8)
            session = self.session(out, backbone_check="yes")
            message = failure(mf.review_backbone, session, self.system)
            self.assertIn("backbone defect", message)
            self.assertIn("A:LEU2", message)
            stage = None
            try:
                mf.run_stage(session, "backbone", mf.review_backbone, session, self.system)
            except mf.StageFailure as exc:
                stage = exc.stage
            self.assertEqual(stage, "backbone")

    def test_off_and_report_modes(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            self.files(out, after_hand=[1, -1] + [1] * 8)
            self.assertEqual(mf.review_backbone(self.session(out, backbone_check="no"), self.system), {"status": "off"})
            self.assertFalse((out / "backbone_check.json").exists())
            self.assertEqual(mf.review_backbone(self.session(out, backbone_check="report"), self.system)["status"], "PASS")


if __name__ == "__main__":
    unittest.main()
