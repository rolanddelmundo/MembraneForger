"""Membrane orientation: anchor selection, providers, the rigid transform, its validation and the CG registration.

The PPM 3.0 provider is exercised through a stand-in executable that reproduces the real program's I/O contract
(input on stdin, res.lib in the working directory, anchorout.pdb + datapar1 + datasub1 as outputs, dummy atoms,
hydrogens dropped). Tests that need the real program run only when MEMBRANEFORGER_PPM points at it, or `immers`
is on PATH; tests that need the network (OPM download) run only when MEMBRANEFORGER_TEST_NETWORK is set.
"""
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

import numpy as np

from .common import AA_PDB, CG_GRO, FORCEFIELD, MAPPING, REPO, failure, mf

OPM_6WHC = Path(__file__).parent / "data" / "6whc_opm_chainR.pdb"  # chain R and the dummy atoms of the OPM entry 6WHC


def rotation(axis: np.ndarray, degrees: float) -> np.ndarray:
    """Rotation matrix about an axis."""
    axis = np.asarray(axis, dtype=float) / np.linalg.norm(axis)
    a = np.radians(degrees)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * K @ K


FAKE_PPM = textwrap.dedent('''\
    #!/usr/bin/env python3
    """Stand-in for PPM 3.0 (immers): same I/O contract, a fixed rigid transform instead of the energy optimisation."""
    import math, os, sys
    from pathlib import Path
    mode = os.environ.get("FAKE_PPM_MODE", "ok")
    text = sys.stdin.read().splitlines()
    assert text[0] == "2", text
    pdb = text[2]
    topo = text[6].strip()
    assert Path("res.lib").is_file(), "res.lib must be in the working directory"
    out = Path(pdb).stem + "out.pdb"
    if mode == "crash":
        sys.stderr.write("segmentation fault (fake)\\n"); sys.exit(139)
    if mode == "silent":
        sys.exit(0)  # exit 0, no output at all
    lines = [l for l in Path(pdb).read_text().splitlines() if l[:6] in ("ATOM  ", "HETATM")]
    th = math.radians(float(os.environ.get("FAKE_PPM_ANGLE", "33.0")))
    ax = [0.3, -0.5, 0.81]; n = math.sqrt(sum(v * v for v in ax)); ax = [v / n for v in ax]
    c, s, C = math.cos(th), math.sin(th), 1 - math.cos(th)
    R = [[c + ax[0]*ax[0]*C, ax[0]*ax[1]*C - ax[2]*s, ax[0]*ax[2]*C + ax[1]*s],
         [ax[1]*ax[0]*C + ax[2]*s, c + ax[1]*ax[1]*C, ax[1]*ax[2]*C - ax[0]*s],
         [ax[2]*ax[0]*C - ax[1]*s, ax[2]*ax[1]*C + ax[0]*s, c + ax[2]*ax[2]*C]]
    t = [12.5, -7.25, 3.0]
    rows = []
    for l in lines:
        if l[13] == "H" or l[12] == "H":
            continue  # PPM drops hydrogens from its output
        x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
        p = [R[i][0]*x + R[i][1]*y + R[i][2]*z + t[i] for i in range(3)]
        if mode == "reflect":
            p[0] = -p[0]
        if mode == "shuffle" and l[17:20] == "LEU" and l[12:16].strip() == "CA":
            p[2] += 2.5  # a non-rigid change
        rows.append(l[:30] + "%8.3f%8.3f%8.3f" % tuple(p) + l[54:])
    half = 15.7
    if mode == "malformed":
        body = "REMARK junk\\n" + "\\n".join(rows[:5]) + "\\n"
    elif mode == "curved":
        body = ""  # datapar1 stays empty: PPM found a curved membrane
    else:
        dummies = []
        k = 90000
        for i in range(-3, 4):
            for j in range(-3, 4):
                for sign, name in ((-1, "N"), (1, "O")):
                    k += 1
                    dummies.append("HETATM%5d  %s   DUM %5d    %8.3f%8.3f%8.3f" % (k % 100000, name, k % 10000, 10.0*i, 10.0*j, sign*half))
        body = "REMARK      1/2 of bilayer thickness:%7.1f\\n" % half + "\\n".join(rows) + "\\n" + "\\n".join(dummies) + "\\nEND\\n"
    Path(out).write_text(body)
    Path("datapar1").write_text("" if mode == "curved" else "%-14s;%4.1f;%4.1f;%4d;%4d;%6.1f;%6.1f;\\n" % (pdb, 2*half, 1.1, 9, 2, -45.6, 0.0))
    Path("datapar2").write_text("%-14s;%4.1f;%5.0f;%4d;%6.1f;%6.1f;\\n" % (pdb, 2*half, 90.0, 9, -50.0, -45.6) if mode == "curved" else "")
    Path("datasub1").write_text("%-14s;A; 9; 1(  10-  30)\\n" % pdb)
    print("anchor.pdb emin=   -45.6 thickn= %4.1f+- 1.1 tilt=     9.+-    2." % (2*half))
    sys.exit(0)
''')


class Fixtures:
    """Shared, lazily built fixtures: the example complex, a fake PPM, and OPM-like references."""
    aa = None

    @classmethod
    def complex(cls):
        if cls.aa is None:
            cls.aa, _ = mf.read_all_atom(AA_PDB, FORCEFIELD)
        return cls.aa

    @staticmethod
    def fake_ppm(root: Path) -> Path:
        exe_dir = root / "ppm"
        exe_dir.mkdir(exist_ok=True)
        exe = exe_dir / "immers"
        exe.write_text(FAKE_PPM)
        exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
        (exe_dir / "res.lib").write_text("fake residue library\n")
        return exe

    @staticmethod
    def session_dirs(root: Path):
        out, work = root / "out", root / "out" / "work" / "orientation"
        work.mkdir(parents=True, exist_ok=True)
        (out / mf.LOG_NAME).write_text("")
        return out, work


class AnchorSelection(unittest.TestCase):
    def setUp(self):
        self.chains = mf.protein_chains(Fixtures.complex())

    def test_protein_chains_of_the_example(self):
        self.assertEqual(list(self.chains), ["R", "P", "A", "B", "C"])
        self.assertEqual(len(self.chains["R"]), 401)

    def test_explicit_chain_and_chains(self):
        self.assertEqual(mf.select_anchor_chains(self.chains, ("R",)), (["R"], "command line"))
        self.assertEqual(mf.select_anchor_chains(self.chains, ("A", "B"))[0], ["A", "B"])
        self.assertIn("not protein chains", failure(mf.select_anchor_chains, self.chains, ("Z",)))
        self.assertIn("not protein chains", failure(mf.select_anchor_chains, self.chains, ("L",)))  # the ligand chain

    def test_single_chain_is_selected_automatically(self):
        single = {"R": self.chains["R"]}
        self.assertEqual(mf.select_anchor_chains(single, ()), (["R"], "the only protein chain"))

    def test_multiple_chains_without_a_choice_fail_with_instructions(self):
        message = failure(mf.select_anchor_chains, self.chains, ())
        self.assertIn("Multiple protein chains were detected: R, P, A, B, C", message)
        self.assertIn("--orient-chain R", message)
        self.assertIn("--orient-chains R,P", message)

    def test_the_file_name_is_never_a_pdb_id(self):
        self.assertEqual(mf.resolve_pdb_id(None, AA_PDB), (None, "none"))  # prot-lig.pdb carries no HEADER
        self.assertEqual(mf.resolve_pdb_id("6whc", AA_PDB), ("6WHC", "command line"))
        self.assertIn("not a PDB ID", failure(mf.resolve_pdb_id, "GLP1R", AA_PDB))

    def test_header_record_gives_the_pdb_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.pdb"
            path.write_text("HEADER    MEMBRANE PROTEIN                        07-APR-20   6WHC              \n"
                            "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C\n")
            self.assertEqual(mf.pdb_id_from_header(path), "6WHC")
            path.write_text("HEADER    MEMBRANE PROTEIN                        07-APR-20   XXXX              \n")
            self.assertIsNone(mf.pdb_id_from_header(path))

    def test_cli_chain_parsing(self):
        self.assertEqual(mf.parse_chains("R"), ("R",))
        self.assertEqual(mf.parse_chains("A,B", "C"), ("A", "B", "C"))
        self.assertEqual(mf.parse_chains(None, ""), ())


class OPMReferenceMode(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = mf.parse_opm_file(OPM_6WHC, "6WHC", str(OPM_6WHC))
        cls.chains = mf.protein_chains(Fixtures.complex())

    def test_reference_is_parsed_without_dummy_atoms(self):
        self.assertEqual(self.reference.half_thickness_a, 14.3)
        self.assertEqual(list(self.reference.chains), ["R"])
        self.assertEqual(mf.spanning_chains(self.reference), {"R"})
        self.assertTrue(all(name != "DUM" for _, name, _ in self.reference.chains["R"]))

    def test_mismatched_dummy_planes_or_missing_remark_are_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            text = OPM_6WHC.read_text()
            bad = Path(tmp) / "bad.pdb"
            bad.write_text(text.replace("1/2 of bilayer thickness:   14.3", "1/2 of bilayer thickness:   20.0"))
            self.assertIn("do not match the half thickness", failure(mf.parse_opm_file, bad, "6WHC", "x"))
            bad.write_text("\n".join(l for l in text.splitlines() if "1/2 of" not in l))
            self.assertIn("not an OPM/OPRLM", failure(mf.parse_opm_file, bad, "6WHC", "x"))

    def test_exact_hit_identifies_the_spanning_chain_and_the_nterm_side(self):
        anchors, how = mf.select_anchor_chains(self.chains, (), self.reference)
        self.assertEqual(anchors, ["R"])
        self.assertIn("membrane-spanning", how)
        side, why = mf.nterm_side_from_reference(self.reference, "R")
        self.assertEqual(side, "out")
        self.assertIn("GLN27", why)

    def test_orientation_from_the_reference_fits_the_embedded_core(self):
        derived = mf.orientation_from_opm({"R": self.chains["R"]}, self.reference, mf.Settings())
        fit = derived["fit"]
        self.assertLessEqual(fit["core_rmsd_A"], 1.0)
        self.assertGreaterEqual(fit["core_fraction"], 0.7)
        self.assertGreater(fit["matched_atoms"], 400)
        self.assertAlmostEqual(fit["determinant"], 1.0, places=6)
        self.assertEqual(derived["matches"][0]["identity"], 1.0)
        # the whole chain deviates more than the core: the example model has moved outside the membrane
        self.assertGreater(fit["whole_chain_rmsd_A"], fit["core_rmsd_A"])

    def test_pdb_id_hit_but_incompatible_chain_is_rejected(self):
        message = failure(mf.orientation_from_opm, {"A": self.chains["A"]}, self.reference, mf.Settings())
        self.assertIn("no chain of OPM 6WHC", message)  # G alpha is not in the single-chain reference
        bent = [dict(a, z=a["z"] * 0.6) if a["chain"] == "R" else a for a in Fixtures.complex()]
        message = failure(mf.orientation_from_opm, {"R": mf.protein_chains(bent)["R"]}, self.reference, mf.Settings())
        self.assertIn("does not superpose", message)

    def test_full_complex_from_reference_file_in_opm_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, work = Fixtures.session_dirs(Path(tmp))
            request = mf.OrientationRequest(mode="opm", chains=("R",), pdb_id="6WHC", opm_file=OPM_6WHC)
            result = mf.orient_complex(Fixtures.complex(), AA_PDB, request, mf.Settings(), out, work)
            report = result["report"]
            self.assertEqual((report["provider"], report["nterm_side"], report["status"]), ("opm", "out", "PASS"))
            self.assertEqual(report["validation"]["nterm"]["status"], "PASS")
            self.assertLess(report["validation"]["pairwise_distance_max_change_A"], 1e-6)
            self.assertTrue((out / "oriented.pdb").is_file())
            # the complex sits in the membrane frame: the receptor straddles z = 0, the G protein lies below (IN)
            z = {label: np.array([a["z"] for a in result["oriented"] if a["chain"] == label and a["atom"] == "CA"]) for label in "RA"}
            self.assertTrue((z["R"] > 14.3).any() and (z["R"] < -14.3).any())
            self.assertLess(z["A"].mean(), -14.3)

    def test_opm_mode_without_a_reference_fails_before_ppm(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, work = Fixtures.session_dirs(Path(tmp))
            request = mf.OrientationRequest(mode="opm", chains=("R",))
            self.assertIn("needs an exact OPM/OPRLM entry", failure(mf.orient_complex, Fixtures.complex(), AA_PDB, request, mf.Settings(), out, work))

    def test_contradicting_nterm_side_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, work = Fixtures.session_dirs(Path(tmp))
            request = mf.OrientationRequest(mode="opm", chains=("R",), pdb_id="6WHC", opm_file=OPM_6WHC, nterm_side="in")
            message = failure(mf.orient_complex, Fixtures.complex(), AA_PDB, request, mf.Settings(), out, work)
            self.assertIn("contradicts the exact reference", message)


class LocalPPMProvider(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.exe = Fixtures.fake_ppm(self.tmp)
        self.out, self.work = Fixtures.session_dirs(self.tmp)
        self.aa = Fixtures.complex()
        self.chains = mf.protein_chains(self.aa)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)
        os.environ.pop("FAKE_PPM_MODE", None)
        os.environ.pop("FAKE_PPM_ANGLE", None)

    def run_ppm(self, anchors=("R",), nterm="out", **kwargs):
        ppm = mf.LocalPPM(self.exe, **kwargs)
        return ppm.run(self.aa, list(anchors), nterm, self.out, self.work)

    def test_executable_discovery(self):
        self.assertEqual(mf.find_ppm_executable(str(self.exe)), self.exe.resolve())
        if not shutil.which("immers") and not os.environ.get("MEMBRANEFORGER_PPM"):
            self.assertIn("not found", failure(mf.find_ppm_executable, None))
        (self.exe.parent / "res.lib").unlink()
        self.assertIn("res.lib is missing", failure(mf.find_ppm_executable, str(self.exe)))
        self.assertIn("not an executable", failure(mf.find_ppm_executable, str(self.tmp / "absent")))

    def test_input_file_is_the_documented_planar_single_membrane_form(self):
        ppm = mf.LocalPPM(self.exe, membrane="PMm")
        self.assertEqual(ppm.input_text(["A", "B"], "in"), "2\nno\nanchor.pdb\n1\nPMm\nplanar\nin \nA,B\n")
        self.assertEqual(mf.LocalPPM(self.exe).input_text(["A"], "out").splitlines()[4], "   ")  # undefined membrane
        self.assertIn("unknown PPM membrane code", failure(mf.LocalPPM, self.exe, membrane="XYZ"))

    def test_only_the_anchor_chain_is_submitted_and_outputs_are_recorded(self):
        result = self.run_ppm()
        submitted, _ = mf.read_pdb(Path(result["workdir"]) / "anchor.pdb")
        self.assertEqual({a["chain"] for a in submitted}, {"A"})  # relabelled; the relabel map is recorded
        self.assertEqual(result["chain_relabel"], {"R": "A"})
        self.assertEqual(len(submitted), sum(a["chain"] == "R" and a["resname"] in mf.AMINO for a in self.aa))
        self.assertEqual(result["residues_submitted"], 401)
        self.assertEqual(result["half_thickness_a"], 15.7)
        self.assertEqual((result["hydrophobic_thickness_a"], result["tilt_deg"], result["transfer_energy_kcal_mol"]), (31.4, 9.0, -45.6))
        self.assertGreater(result["dummy_atoms_ignored"], 0)
        self.assertEqual(result["input"].splitlines()[1], "no")
        self.assertTrue(result["executable_sha256"] and result["res_lib_sha256"] and result["output_sha256"])

    def test_heteroatoms_are_submitted_only_on_request(self):
        result = self.run_ppm(anchors=("R", "L")) if "L" in self.chains else None  # L is not a protein chain
        self.assertIsNone(result)
        with_hetero = self.run_ppm(anchors=("R",), heteroatoms=True)
        self.assertEqual(with_hetero["input"].splitlines()[1], "yes")

    def test_transform_is_derived_despite_dummy_atoms_and_missing_hydrogens(self):
        result = self.run_ppm()
        derived = mf.orientation_from_ppm({"R": self.chains["R"]}, result, mf.Settings())
        self.assertLess(derived["fit"]["rmsd_A"], 0.002)
        self.assertGreaterEqual(derived["fit"]["matched_fraction"], 0.99)
        expected = rotation([0.3, -0.5, 0.81], 33.0)
        self.assertTrue(np.allclose(derived["R"], expected, atol=1e-4))
        self.assertTrue(np.allclose(derived["t"], [12.5, -7.25, 3.0], atol=2e-3))

    def test_missing_executable(self):
        self.assertIn("not an executable", failure(mf.find_ppm_executable, str(self.tmp / "nope")))

    def test_nonzero_exit(self):
        os.environ["FAKE_PPM_MODE"] = "crash"
        self.assertIn("exit 139", failure(self.run_ppm))

    def test_zero_exit_without_output(self):
        os.environ["FAKE_PPM_MODE"] = "silent"
        self.assertIn("did not write", failure(self.run_ppm))

    def test_zero_exit_with_malformed_output(self):
        os.environ["FAKE_PPM_MODE"] = "malformed"
        self.assertIn("no thickness REMARK", failure(self.run_ppm))

    def test_curved_result_is_refused(self):
        os.environ["FAKE_PPM_MODE"] = "curved"
        self.assertIn("datapar1", failure(self.run_ppm))

    def test_stale_output_cannot_be_reused(self):
        self.run_ppm()
        stale = Path(self.work / "ppm" / "anchorout.pdb")
        self.assertTrue(stale.is_file())
        os.environ["FAKE_PPM_MODE"] = "silent"
        self.assertIn("did not write", failure(self.run_ppm))
        self.assertFalse((self.work / "ppm" / "anchorout.pdb").exists())  # the directory was rebuilt, nothing stale survives

    def test_reflection_and_non_rigid_output_are_refused(self):
        os.environ["FAKE_PPM_MODE"] = "reflect"
        result = self.run_ppm()
        self.assertIn("not a rigid copy", failure(mf.orientation_from_ppm, {"R": self.chains["R"]}, result, mf.Settings()))
        os.environ["FAKE_PPM_MODE"] = "shuffle"
        result = self.run_ppm()
        self.assertIn("not a rigid copy", failure(mf.orientation_from_ppm, {"R": self.chains["R"]}, result, mf.Settings()))

    def test_nterm_side_is_required_without_a_trustworthy_source(self):
        request = mf.OrientationRequest(mode="ppm", chains=("R",), ppm_exe=str(self.exe))
        message = failure(mf.orient_complex, self.aa, AA_PDB, request, mf.Settings(), self.out, self.work)
        self.assertIn("sidedness cannot be determined safely", message)
        self.assertIn("--nterm-side in", message)
        self.assertFalse((self.work / "ppm").exists())  # nothing expensive ran

    def test_nterm_side_in_and_out_are_passed_through(self):
        for side in ("in", "out"):
            request = mf.OrientationRequest(mode="ppm", chains=("R",), ppm_exe=str(self.exe), nterm_side=side)
            try:
                result = mf.orient_complex(self.aa, AA_PDB, request, mf.Settings(), self.out, self.work)
            except SystemExit as exc:  # the fake transform does not know sides: an "in" request fails the sidedness check
                self.assertEqual(side, "in")
                self.assertIn("N terminus out", str(exc))
                continue
            self.assertEqual(result["report"]["ppm"]["nterm_side"], side)
            self.assertEqual(result["report"]["ppm"]["input"].splitlines()[6].strip(), side)

    def test_nterm_auto_with_trustworthy_reference_then_ppm(self):
        request = mf.OrientationRequest(mode="ppm", chains=("R",), ppm_exe=str(self.exe), pdb_id="6WHC", opm_file=OPM_6WHC)
        result = mf.orient_complex(self.aa, AA_PDB, request, mf.Settings(), self.out, self.work)
        report = result["report"]
        self.assertEqual((report["provider"], report["nterm_side"]), ("ppm", "out"))
        self.assertIn("OPM 6WHC chain R", report["nterm_side_source"])

    def test_auto_mode_falls_back_to_ppm_when_the_reference_chain_is_incompatible(self):
        bent = [dict(a, z=a["z"] * 0.6) if a["chain"] == "R" else a for a in self.aa]
        request = mf.OrientationRequest(mode="auto", chains=("R",), ppm_exe=str(self.exe), pdb_id="6WHC", opm_file=OPM_6WHC)
        result = mf.orient_complex(bent, AA_PDB, request, mf.Settings(), self.out, self.work)
        self.assertEqual(result["report"]["provider"], "ppm")
        self.assertIn("does not superpose", result["report"]["reference_rejected"])

    def test_auto_mode_without_pdb_id_uses_ppm(self):
        request = mf.OrientationRequest(mode="auto", chains=("R",), ppm_exe=str(self.exe), nterm_side="out")
        result = mf.orient_complex(self.aa, AA_PDB, request, mf.Settings(), self.out, self.work)
        self.assertEqual(result["report"]["provider"], "ppm")
        self.assertEqual(result["report"]["pdb_id_source"], "none")

    def test_two_anchor_chains(self):
        request = mf.OrientationRequest(mode="ppm", chains=("R", "P"), ppm_exe=str(self.exe), nterm_side="out")
        result = mf.orient_complex(self.aa, AA_PDB, request, mf.Settings(), self.out, self.work)
        report = result["report"]
        self.assertEqual(report["anchor_chains"], ["R", "P"])
        self.assertEqual(report["ppm"]["chain_relabel"], {"R": "A", "P": "B"})
        self.assertEqual(report["ppm"]["input"].splitlines()[7], "A,B")
        self.assertEqual(len(report["ppm"]["chain_matches"]), 2)


class RigidTransformOfTheWholeComplex(unittest.TestCase):
    """The whole input receives one identical proper rigid transform; nothing is independently moved or rebuilt."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp())
        cls.exe = Fixtures.fake_ppm(cls.tmp)
        cls.out, cls.work = Fixtures.session_dirs(cls.tmp)
        cls.aa = Fixtures.complex()
        request = mf.OrientationRequest(mode="ppm", chains=("R",), ppm_exe=str(cls.exe), nterm_side="out")
        cls.result = mf.orient_complex(cls.aa, AA_PDB, request, mf.Settings(), cls.out, cls.work)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_atoms_names_residues_chains_and_order_are_preserved(self):
        keys = lambda atoms: [(a["atom"], a["resname"], a["resid"], a["chain"], a["segid"], a["elem"]) for a in atoms]
        self.assertEqual(keys(self.result["oriented"]), keys(self.aa))
        self.assertEqual(len(self.result["oriented"]), 10987)
        self.assertTrue(any(a["resname"] == "GDP" for a in self.result["oriented"]))  # ligands rode along

    def test_pairwise_distances_are_invariant(self):
        X, Y = mf.xyz_nm(self.aa), mf.xyz_nm(self.result["oriented"])
        rng = np.random.default_rng(7)
        i, j = rng.integers(0, len(X), size=(2, 20000))
        self.assertLess(np.abs(np.linalg.norm(X[i] - X[j], axis=1) - np.linalg.norm(Y[i] - Y[j], axis=1)).max(), 1e-9)

    def test_every_chain_received_the_same_transform(self):
        R, t = self.result["R"], self.result["t"]
        self.assertAlmostEqual(np.linalg.det(R), 1.0, places=9)
        for label in "RPABCL":
            idx = [k for k, a in enumerate(self.aa) if a["chain"] == label]
            P, Q = mf.xyz_nm([self.aa[k] for k in idx]), mf.xyz_nm([self.result["oriented"][k] for k in idx])
            self.assertLess(np.abs(P @ R.T + t - Q).max(), 1e-9, label)
            own = mf.rigid_fit(P, Q)  # the chain's own best fit is the global transform
            self.assertTrue(np.allclose(own["R"], R, atol=1e-6) and np.allclose(own["t"], t, atol=1e-6), label)
        chains = self.result["report"]["validation"]["chains"]
        self.assertEqual(set(chains), set("RPABCL"))
        self.assertTrue(all(c["rmsd_to_global_transform_A"] < 1e-9 for c in chains.values()))

    def test_anchor_matches_its_ppm_target(self):
        ppm, _ = mf.read_pdb(Path(self.result["report"]["ppm"]["workdir"]) / "anchorout.pdb")
        target = {(a["resid"], a["resname"], a["atom"]): (a["x"], a["y"], a["z"]) for a in ppm if a["resname"] != "DUM"}
        moved = [a for a in self.result["oriented"] if a["chain"] == "R" and (a["resid"], a["resname"], a["atom"]) in target]
        self.assertGreater(len(moved), 3000)
        error = np.array([np.linalg.norm(np.array([a["x"], a["y"], a["z"]]) - target[(a["resid"], a["resname"], a["atom"])]) for a in moved])
        self.assertLess(error.max(), 0.002)

    def test_report_is_complete_and_reproducible(self):
        report = json.loads(json.dumps(self.result["report"], default=str))
        for key in ("provider", "mode", "anchor_chains", "pdb_id", "nterm_side", "ppm", "fit", "rotation_matrix",
                    "translation_A", "membrane_normal", "membrane_center_A", "validation", "input_sha256", "oriented_sha256"):
            self.assertIn(key, report)
        for key in ("version", "executable", "executable_sha256", "membrane_code", "curvature", "input", "atoms_submitted"):
            self.assertIn(key, report["ppm"])
        R, t = np.array(report["rotation_matrix"]), np.array(report["translation_A"])
        self.assertLess(np.abs(mf.xyz_nm(self.aa) @ R.T + t - mf.xyz_nm(self.result["oriented"])).max(), 1e-9)

    def test_relabelled_chains_and_renumbered_residues_give_the_same_orientation(self):
        swap = {"R": "X", "A": "Q", "P": "Z"}
        other = [dict(a, chain=swap.get(a["chain"], a["chain"]), resid=a["resid"] + 1000 if a["resname"] in mf.AMINO else a["resid"])
                 for a in self.aa]
        request = mf.OrientationRequest(mode="ppm", chains=("X",), ppm_exe=str(self.exe), nterm_side="out")
        result = mf.orient_complex(other, AA_PDB, request, mf.Settings(), self.out, self.work)
        self.assertTrue(np.allclose(result["R"], self.result["R"], atol=1e-6) and np.allclose(result["t"], self.result["t"], atol=1e-3))

    def test_missing_terminal_residues_still_orient(self):
        trimmed = [a for a in self.aa if not (a["chain"] == "R" and (a["resid"] < 40 or a["resid"] > 410))]
        request = mf.OrientationRequest(mode="ppm", chains=("R",), ppm_exe=str(self.exe), nterm_side="out")
        result = mf.orient_complex(trimmed, AA_PDB, request, mf.Settings(), self.out, self.work)
        self.assertEqual(result["report"]["status"], "PASS")
        self.assertTrue(np.allclose(result["R"], self.result["R"], atol=1e-6))

    def test_custom_structure_without_pdb_id_orients_through_ppm(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "model.pdb"
            lines = [l for l in AA_PDB.read_text().splitlines() if l[:6] in ("ATOM  ", "HETATM") and l[21] == "R"]
            path.write_text("\n".join(lines) + "\nEND\n")
            atoms, _ = mf.read_all_atom(path, FORCEFIELD)
            request = mf.OrientationRequest(mode="auto", ppm_exe=str(self.exe), nterm_side="out")  # single chain: no --orient-chain
            result = mf.orient_complex(atoms, path, request, mf.Settings(), self.out, self.work)
            self.assertEqual((result["report"]["provider"], result["report"]["anchor_selection"]), ("ppm", "the only protein chain"))

    def test_orientation_none_leaves_the_input_untouched(self):
        result = mf.orient_complex(self.aa, AA_PDB, mf.OrientationRequest(mode="none"), mf.Settings(), self.out, self.work)
        self.assertIs(result["oriented"], self.aa)
        self.assertFalse(result["report"]["enabled"])
        self.assertTrue(np.allclose(result["R"], np.eye(3)))


class RegistrationKeepsTheOrientation(unittest.TestCase):
    """After CG registration the membrane normal is still z and the complex is still at its membrane depth."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp())
        cls.out, cls.work = Fixtures.session_dirs(cls.tmp)
        raw, cls.box = mf.read_cg(CG_GRO)
        cls.cg = mf.classify_cg(raw, MAPPING)
        cls.slab = mf.make_membrane_whole(cls.cg["membrane"], cls.box)
        cls.midplane, cls.how, cls.leaflets = mf.bilayer_midplane(cls.cg["membrane"], cls.slab)
        request = mf.OrientationRequest(mode="opm", chains=("R",), pdb_id="6WHC", opm_file=OPM_6WHC)
        cls.oriented = mf.orient_complex(Fixtures.complex(), AA_PDB, request, mf.Settings(), cls.out, cls.work)["oriented"]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_cg_midplane_from_phosphate_planes(self):
        self.assertIn("PO4 leaflet planes", self.how)
        self.assertLess(abs(self.midplane - 0.5 * sum(self.leaflets["po4_planes_nm"])), 1e-9)
        self.assertLess(self.leaflets["po4_normal_tilt_from_z_deg"], 5.0)
        self.assertTrue(self.slab[0] < self.midplane < self.slab[1])

    def test_constrained_registration(self):
        fit = mf.map_all_atom_to_cg(self.oriented, self.cg["protein"], self.box, self.slab, mf.Settings(), self.midplane, ("R",), 14.3)
        R, t = fit["R"], fit["t"]
        self.assertTrue(np.allclose(R[2], [0, 0, 1]) and np.allclose(R[:, 2], [0, 0, 1]))  # a rotation about z only
        self.assertAlmostEqual(np.linalg.det(R), 1.0, places=9)
        registration = fit["metrics"]["registration"]
        self.assertIn("inside +-14.3 A", registration["fitted_on"])
        self.assertLess(registration["tilt_between_cg_anchor_pose_and_orientation_deg"], mf.Settings().register_max_tilt_deg)
        self.assertLessEqual(registration["anchor_constrained_core_rmsd_A"], mf.Settings().fit_max_core_rmsd_a)
        placed = mf.apply_transform(self.oriented, R, t)
        check = mf.check_orientation_preserved(self.oriented, placed, ["R"], self.midplane * 10.0)
        self.assertEqual(check["status"], "PASS")
        self.assertLess(check["normal_tilt_deg"], 1e-6)
        self.assertAlmostEqual(check["anchor_ca_depth_A"], check["anchor_ca_depth_after_registration_A"], places=6)
        # the oriented midplane z = 0 lands on the CG midplane (modulo the periodic z image chosen for the complex)
        cell_z = self.box[2] * 10.0
        offset = (t[2] - self.midplane * 10.0) / cell_z
        self.assertLess(abs(offset - round(offset)), 1e-9)
        # correspondence rows still describe every matched chain
        self.assertEqual(fit["assigned"], [("A", 1), ("P", 2), ("R", 0)])

    def test_a_tilted_cg_pose_beyond_the_limit_is_refused(self):
        centre = np.array([[b["x"], b["y"], b["z"]] for r in self.cg["protein"] for b in r]).mean(axis=0)
        Rt = rotation([1.0, 0.0, 0.0], 30.0)
        rotated = []
        for residue in self.cg["protein"]:
            new = []
            for b in residue:
                p = Rt @ (np.array([b["x"], b["y"], b["z"]]) - centre) + centre
                new.append(dict(b, x=float(p[0]), y=float(p[1]), z=float(p[2])))
            rotated.append(new)
        message = failure(mf.map_all_atom_to_cg, self.oriented, rotated, self.box, self.slab, mf.Settings(), self.midplane, ("R",), 14.3)
        self.assertIn("tilted", message)
        self.assertIn("no scientifically defensible registration", message)

    def test_registration_cannot_tilt_the_complex_even_when_the_free_fit_would(self):
        settings = mf.Settings(register_max_tilt_deg=90.0)
        fit = mf.map_all_atom_to_cg(self.oriented, self.cg["protein"], self.box, self.slab, settings, self.midplane, ("R",), 14.3)
        self.assertTrue(np.allclose(fit["R"][2], [0, 0, 1]))
        free = mf.map_all_atom_to_cg(self.oriented, self.cg["protein"], self.box, self.slab, settings)
        self.assertGreater(np.degrees(np.arccos(free["R"][2, 2])), 1.0)  # the plain workflow would have tilted it

    def test_check_detects_a_tilt(self):
        placed = mf.apply_transform(self.oriented, rotation([0.0, 1.0, 0.0], 0.5), np.zeros(3))
        self.assertIn("changed the orientation", failure(mf.check_orientation_preserved, self.oriented, placed, ["R"], 0.0))


class BoxSizing(unittest.TestCase):
    def make_system(self):
        coords = []
        lipid = lambda atom, k, z: {"atom": atom, "group": "MEMB", "resid": 100 + k, "resname": "POPC", "x": 10.0 * k, "y": 5.0, "z": z}
        for k in range(40):  # two phosphate planes 4 nm apart around z = 60 A, plus lipid tails between them
            coords.append(lipid("P", k, 40.0 + (k % 2) * 40.0))
            coords.append(lipid("C2", k, 60.0 + (k % 3) * 2.0))
        coords.append(lipid("C1", 99, 95.0))  # one stray tall head group
        for n, z in enumerate(np.linspace(20.0, 130.0, 50)):  # the solute spans 2-13 nm
            for atom, dz in (("CA", 0.0), ("HA", 1.0)):
                coords.append({"atom": atom, "group": "Protein_LIG", "resid": n + 1, "resname": "GLY", "x": 50.0, "y": 50.0, "z": float(z) + dz})
        return coords

    def rebox(self, requested_z=None, pad=1.5):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / mf.LOG_NAME).write_text("")
            system = {"cryst1": "CRYST1  100.000  100.000  150.000  90.00  90.00  90.00 P 1           1", "out": out, "name": "t"}
            box = mf.rebox_system(system, {"coords": self.make_system()}, requested_z, pad)
            atoms, gro_box = mf.read_gro(out / "boxed.gro")
            return box, system["box_report"], atoms, gro_box

    def test_auto_box_is_centred_on_the_phosphate_midplane(self):
        box, report, atoms, gro_box = self.rebox()
        self.assertEqual(report["mode"], "auto")
        self.assertEqual(box[:2], [10.0, 10.0])
        self.assertAlmostEqual(report["bilayer_midplane_nm"], 6.0)  # not shifted by the stray head group at 9.5 nm
        # the solute reaches 13.1 nm (hydrogen) / 13.0 heavy: 7.0 nm above the midplane + 1.5 pad -> half 8.5 -> 17.0 nm
        self.assertAlmostEqual(box[2], 17.0, places=3)
        self.assertEqual(gro_box[:3], box)
        z = np.array([a["z"] for a in atoms])
        self.assertAlmostEqual(z[[a["atom"] == "P" for a in atoms]].mean(), box[2] / 2.0, places=6)

    def test_user_box_is_opt_in_and_validated(self):
        box, report, _, _ = self.rebox(requested_z=20.0)
        self.assertEqual((report["mode"], box[2]), ("user", 20.0))
        self.assertEqual(box[:2], [10.0, 10.0])                    # x and y are the cell this stage receives
        self.assertIn("is below the 170 A", failure(self.rebox, requested_z=15.0))

    def test_cli_box_parsing(self):
        self.assertIsNone(mf.parse_box("auto"))
        self.assertIsNone(mf.parse_box(None))
        self.assertEqual(mf.parse_box("11.18,11.18,20.3"), (11.18, 11.18, 20.3))
        self.assertIn("three positive edge lengths in A", failure(mf.parse_box, "11,20"))
        self.assertIn("three positive edge lengths in A", failure(mf.parse_box, "a,b,c"))


class CommandLine(unittest.TestCase):
    def run_cli(self, *args):
        result = subprocess.run([sys.executable, "-m", "membraneforger", *map(str, args)], cwd=REPO, text=True, capture_output=True)
        return result.returncode, result.stdout + result.stderr

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.gmx = os.environ.get("MEMBRANEFORGER_GMX") or shutil.which("gmx") or "gmx"

    def test_option_validation(self):
        self.assertIn("invalid choice", self.run_cli("--all-atom", AA_PDB, "--coarse-grain", CG_GRO, "--orientation", "maybe")[1])
        self.assertIn("expected 3 arguments", self.run_cli("--all-atom", AA_PDB, "--coarse-grain", CG_GRO, "--box", "100", "100")[1])
        self.assertIn("three positive edge lengths", self.run_cli("--all-atom", AA_PDB, "--coarse-grain", CG_GRO, "--box", "100", "-5", "100")[1])
        self.assertIn("--bilayer-z applies with --orientation none",
                      self.run_cli("--all-atom", AA_PDB, "--coarse-grain", CG_GRO, "--orient-chain", "R", "--bilayer-z", "5")[1])
        self.assertIn("do not apply", self.run_cli("--membrane", AA_PDB, "--orient-chain", "R")[1])

    def test_multiple_chains_fail_early_with_instructions(self):
        if not shutil.which(self.gmx):
            self.skipTest("GROMACS is needed to reach the orientation stage")
        code, output = self.run_cli("--all-atom", AA_PDB, "--coarse-grain", CG_GRO, "--out", self.tmp / "o", "--gmx", self.gmx)
        self.assertEqual(code, 1)
        self.assertIn("ERROR: orient: Multiple protein chains were detected", output)
        self.assertIn("--orient-chain R", output)
        self.assertFalse((self.tmp / "o" / "em.gro").exists())
        manifest = json.loads((self.tmp / "o" / "run_manifest.json").read_text())
        self.assertEqual(manifest["error"]["stage"], "orient")

    def test_nterm_side_is_requested_before_anything_expensive(self):
        if not shutil.which(self.gmx):
            self.skipTest("GROMACS is needed to reach the orientation stage")
        exe = Fixtures.fake_ppm(self.tmp)
        code, output = self.run_cli("--all-atom", AA_PDB, "--coarse-grain", CG_GRO, "--out", self.tmp / "o", "--gmx", self.gmx,
                                    "--orientation", "ppm", "--orient-chain", "R", "--ppm-exe", exe)
        self.assertEqual(code, 1)
        self.assertIn("--nterm-side out", output)
        self.assertNotIn("stage mstool", output)


@unittest.skipUnless(os.environ.get("MEMBRANEFORGER_PPM") or shutil.which("immers"), "real PPM 3.0 not installed")
class RealPPM(unittest.TestCase):
    """With the real program: the output is an exact rigid copy and the receptor straddles the slab with its N terminus out."""

    def test_real_ppm_orients_the_example_receptor(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, work = Fixtures.session_dirs(Path(tmp))
            request = mf.OrientationRequest(mode="ppm", chains=("R",), nterm_side="out")
            result = mf.orient_complex(Fixtures.complex(), AA_PDB, request, mf.Settings(), out, work)
            report = result["report"]
            self.assertLess(report["fit"]["rmsd_A"], 0.01)
            self.assertEqual(report["validation"]["nterm"]["status"], "PASS")
            frame = report["validation"]["frame"]
            self.assertGreater(frame["anchor_ca_inside_slab"], 100)
            self.assertTrue(frame["anchor_ca_above_slab"] > 0 and frame["anchor_ca_below_slab"] > 0)
            self.assertTrue(28.0 <= report["ppm"]["hydrophobic_thickness_a"] <= 36.0)


@unittest.skipUnless(os.environ.get("MEMBRANEFORGER_TEST_NETWORK"), "set MEMBRANEFORGER_TEST_NETWORK to download from OPM")
class OPMDownload(unittest.TestCase):
    def test_download_and_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            reference = mf.fetch_opm_reference("6WHC", Path(tmp))
            self.assertEqual(reference.half_thickness_a, 14.3)
            self.assertTrue((Path(tmp) / "6whc.pdb").is_file())
            self.assertIsNone(mf.fetch_opm_reference("9ZZZ", Path(tmp)))


if __name__ == "__main__":
    unittest.main()
