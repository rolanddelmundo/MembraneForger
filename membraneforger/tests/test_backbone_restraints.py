"""EM-only backbone phi/psi restraints: the clipped windows stay inside the backbone check's allowed region, contain the input
angle, wrap through +-180, and leave residues the check does not judge (or that sit too close to a boundary) at +-30 degrees."""
import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from .common import AA_PDB, mf
from .test_backbone_check import chain

HALF = mf.BB_DIHRES_DPHI


def residue(name, phi, psi, resid=1):
    """A residue for restraint_windows: angles in degrees (None for a chain end), dummy atom numbers."""
    return {"resid": resid, "resname": name, "phi": None if phi is None else ((1, 2, 3, 4), phi),
            "psi": None if psi is None else ((2, 3, 4, 5), psi)}


def rows_of(name, phi, psi, **kw):
    return {r["angle"]: r for r in mf.restraint_windows([residue(name, phi, psi)], **kw)}


def wrap(x):
    return (np.asarray(x, float) + 180.0) % 360.0 - 180.0


class AllowedRegion(unittest.TestCase):
    def test_vectorized_region_is_the_backbone_check_region(self):
        grid = np.arange(-180.0, 180.0, 1.0)
        phi, psi = np.meshgrid(grid, grid, indexing="ij")
        for parent in ("ALA", "PRO"):
            for inset in (0.0, 5.0, 12.5):
                mask = mf.rama_allowed(parent, phi, psi, inset)
                expected = np.array([[not mf.rama_outlier(parent, f, s, margin=-inset) for s in grid] for f in grid])
                self.assertTrue(np.array_equal(mask, expected), f"{parent} inset {inset}")


class Windows(unittest.TestCase):
    def test_window_near_the_eclipsed_band_is_clipped_and_keeps_the_input(self):
        rows = rows_of("PHE", -28.83, 138.37)            # 8.8 deg from the phi = -20 boundary
        phi, psi = rows["phi"], rows["psi"]
        self.assertEqual(phi["status"], "clipped")
        self.assertAlmostEqual(phi["low"], -58.83)
        self.assertAlmostEqual(phi["high"], -25.33)        # the last 0.5 deg step at or below -25 (= -20 + inset)
        self.assertAlmostEqual(phi["centre"], -42.08)
        self.assertAlmostEqual(phi["half_width"], 16.75)
        self.assertEqual(psi["status"], "unchanged")
        self.assertAlmostEqual(psi["half_width"], HALF)

    def test_window_through_180_wraps_its_centre(self):
        rows = rows_of("ASP", 156.98, 104.5)              # beta region across phi = +-180
        phi = rows["phi"]
        self.assertEqual(phi["status"], "clipped")
        self.assertGreater(phi["high"], 180.0)             # unwrapped around the input
        self.assertTrue(-180.0 <= phi["centre"] < 180.0)
        self.assertAlmostEqual(float(wrap(phi["centre"] - (phi["low"] + phi["high"]) / 2)), 0.0)
        self.assertLessEqual(phi["low"], 156.98)
        self.assertGreaterEqual(phi["low"], 155.0)         # 150 + inset

    def test_proline_only_phi_is_limited(self):
        rows = rows_of("PRO", -45.0, 150.0)
        self.assertEqual(rows["psi"]["status"], "unchanged")
        self.assertEqual(rows["phi"]["status"], "clipped")
        self.assertLessEqual(rows["phi"]["high"], -40.0)

    def test_residues_that_keep_thirty_degrees(self):
        cases = {"input_outlier": ("ALA", 10.0, 0.0), "flag_margin": ("ALA", -22.0, 140.0),
                 "not_checked": ("GLY", -28.0, 138.0)}
        for status, (name, phi, psi) in cases.items():
            for row in rows_of(name, phi, psi).values():
                self.assertEqual(row["status"], status)
                self.assertAlmostEqual(row["half_width"], HALF)
                self.assertAlmostEqual(row["centre"], row["input"])
        self.assertEqual(rows_of("XYZ", -28.0, 138.0)["phi"]["status"], "not_checked")   # no standard parent residue
        self.assertEqual(rows_of("ALA", -28.0, None)["phi"]["status"], "not_checked")    # chain end: no psi
        self.assertEqual(rows_of("PHE", -28.83, 138.37, min_half=20.0)["phi"]["status"], "flag_narrow")

    def test_fixed_mode_and_bad_arguments(self):
        for row in rows_of("PHE", -28.83, 138.37, mode="fixed").values():
            self.assertEqual((row["status"], row["half_width"]), ("fixed", HALF))
        with self.assertRaises(ValueError):
            rows_of("ALA", -60.0, -40.0, mode="off")
        with self.assertRaises(ValueError):
            rows_of("ALA", -60.0, -40.0, inset=HALF)

    def test_every_window_is_inside_the_region_and_contains_its_input(self):
        rng = np.random.default_rng(7)
        inset, fine = 5.0, np.linspace(0.0, 1.0, 81)
        for name in ("ALA", "PRO", "HSD", "KSM"):          # HSD and KSM map to their parents HIS and LYS
            parent = mf.PARENT_RESIDUE[name]
            for phi0, psi0 in rng.uniform(-180.0, 180.0, size=(250, 2)):
                rows = rows_of(name, float(phi0), float(psi0), inset=inset)
                p, s = rows["phi"], rows["psi"]
                for r in (p, s):
                    self.assertTrue(r["low"] - 1e-9 <= r["input"] <= r["high"] + 1e-9)
                    self.assertGreaterEqual(r["half_width"], mf.BB_DIHRES_MIN_HALF_DEG)
                if p["status"] in ("clipped", "unchanged"):
                    f = wrap(p["low"] + fine * (p["high"] - p["low"]))
                    g = wrap(s["low"] + fine * (s["high"] - s["low"]))
                    self.assertTrue(mf.rama_allowed(parent, f[:, None], g[None, :], inset).all(), f"{name} {phi0:.1f}/{psi0:.1f}")
                else:
                    self.assertIn(p["status"], ("input_outlier", "flag_margin", "flag_narrow"))
                    self.assertAlmostEqual(p["half_width"], HALF)


def write_chain(directory: Path, names, angles, shift_from=None):
    """A synthetic chain as an .itp (atoms only) and the matching PDB in Angstrom; residues from shift_from on move 10 A away."""
    s = chain(names, angles)
    itp, atoms, nr = ["[ moleculetype ]\nPROA 3\n\n[ atoms ]\n"], [], 0
    for i, r in enumerate(s["residues"]):
        for k in ("N", "CA", "C", "CB"):
            nr += 1
            xyz = r["atoms"][k] * 10.0 + (np.array([10.0, 0, 0]) if shift_from is not None and i >= shift_from else 0)
            itp.append(f"{nr:6d} XX {r['resid']:5d} {r['name']:>4s} {k:>4s} {nr:5d} 0.0 12.0\n")
            atoms.append({"atom": k, "resname": r["name"], "resid": r["resid"], "chain": "A", "x": xyz[0], "y": xyz[1], "z": xyz[2]})
    (directory / "PROA.itp").write_text("".join(itp))
    mf.write_pdb(atoms, directory / "PROA.pdb")
    return directory / "PROA.itp", directory / "PROA.pdb"


class Topology(unittest.TestCase):
    NAMES = ["ALA", "PHE", "ASP", "PRO", "GLY", "ALA"]
    ANGLES = [(-63.0, -41.0, 180.0), (-28.83, 138.37, 180.0), (156.98, 104.5, 180.0), (-65.0, 150.0, 180.0),
              (80.0, 10.0, 180.0), (-63.0, -41.0, 180.0)]

    def test_angles_rows_and_itp_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            itp, pdb = write_chain(Path(tmp), self.NAMES, self.ANGLES)
            residues = mf.backbone_angles(mf.read_itp(itp)["atoms"], mf.read_pdb(pdb)[0])
            self.assertIsNone(residues[0]["phi"])
            self.assertIsNone(residues[-1]["psi"])
            for r, (phi, psi, _) in zip(residues[1:-1], self.ANGLES[1:-1]):
                self.assertAlmostEqual(r["phi"][1], phi, delta=0.1)   # PDB coordinates keep 3 decimals
                self.assertAlmostEqual(r["psi"][1], psi, delta=0.1)
            rows = mf.write_backbone_dihres(itp, pdb, "clipped", 5.0)
            self.assertEqual(len(rows), 2 * len(self.NAMES) - 2)
            self.assertEqual({r["molecule"] for r in rows}, {"PROA"})
            text = itp.read_text()
            block = text[text.index("#ifdef BBDIHRES"):text.index("#endif")]
            restraint_rows = [line.split() for line in block.splitlines() if line.strip()[:1].isdigit()]
            self.assertEqual(len(restraint_rows), len(rows))
            for line, row in zip(restraint_rows, rows):
                self.assertEqual(tuple(int(x) for x in line[:4]), row["atoms"])
                self.assertAlmostEqual(float(line[5]), row["centre"], places=2)
                self.assertAlmostEqual(float(line[6]), row["half_width"], places=2)
                self.assertEqual(line[7], "DIHRES_FC")
            status = {(r["resname"], r["angle"]): r["status"] for r in rows}
            self.assertEqual(status["PHE", "phi"], "clipped")
            self.assertEqual(status["GLY", "phi"], "not_checked")
            path = mf.write_restraint_record(Path(tmp), rows)
            self.assertEqual(len(path.read_text().splitlines()), len(rows) + 1)
            self.assertEqual(path.parent, Path(tmp))          # never inside toppar/

    def test_chain_break_and_repeated_residue_numbers(self):
        with tempfile.TemporaryDirectory() as tmp:
            itp, pdb = write_chain(Path(tmp), self.NAMES, self.ANGLES, shift_from=3)
            residues = mf.backbone_angles(mf.read_itp(itp)["atoms"], mf.read_pdb(pdb)[0])
            self.assertIsNone(residues[2]["psi"])           # no peptide bond across the break
            self.assertIsNone(residues[3]["phi"])
            atoms = mf.read_itp(itp)["atoms"]
            for a in atoms:                                  # every residue numbered 7: still six residues
                a["resnr"], a["residue"] = 7, "ALA"
            self.assertEqual(len(mf.backbone_angles(atoms, mf.read_pdb(pdb)[0])), len(self.NAMES))

    def test_audit_ignores_the_backbone_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            itp = Path(tmp) / "MOL.itp"
            itp.write_text("[ moleculetype ]\nMOL 3\n\n[ dihedral_restraints ]\n1 2 3 4 1 60.0 0.0 DIHRES_FC\n"
                           "\n#ifdef BBDIHRES\n[ dihedral_restraints ]\n5 6 7 8 1 -42.0 16.75 DIHRES_FC\n#endif\n")
            self.assertEqual(len(mf.itp_dihedral_restraints([itp])["MOL"]), 1)

    def test_audit_adds_posres_to_an_existing_define_line(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            (run / "em.mdp").write_text("define          = -DDIHRES_FC=1000 -DBBDIHRES\nintegrator      = steep\n")
            for name in ("em.gro", "topol.top", "index_ini.ndx"):
                (run / name).write_text("")
            with mock.patch("membraneforger.audit.gromacs", return_value=(0, "")):
                mf.check_grompp(run, "gmx", "em.gro")
            posres = (mf.audit.audit_scratch(run) / "posres.mdp").read_text()
            defines = [line for line in posres.splitlines() if line.startswith("define")]
            self.assertEqual(defines, ["define          = -DDIHRES_FC=1000 -DBBDIHRES -DPOSRES"])   # grompp rejects a repeated key

    def test_minimization_defines_bbdihres_in_both_stages(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / "toppar").mkdir()
            (out / "toppar" / "PROA.itp").write_text("\n#ifdef BBDIHRES\n[ dihedral_restraints ]\n1 2 3 4 1 -60 30 DIHRES_FC\n#endif\n")
            (out / "toppar" / "LIP.itp").write_text("[ dihedral_restraints ]\n1 2 3 4 1 60 0 DIHRES_FC\n")
            with mock.patch("membraneforger.minimization._minimize") as run:
                mf.run_em(out, {"charge": 0.0}, "gmx", 1, 10)
            defines = {c.args[2]: c.args[7] for c in run.call_args_list}
            self.assertIn("-DBBDIHRES", defines["emres"])
            self.assertIn("-DBBDIHRES", defines["em"])
            self.assertIn("-DDIHRES_FC=", defines["em"])
            self.assertNotIn("-DDIHRES ", defines["em"] + " ")


class CommandLine(unittest.TestCase):
    def test_default_is_off_and_the_inset_is_checked(self):
        self.assertEqual((mf.Settings().backbone_restraints, mf.Settings().backbone_restraint_inset_deg), ("off", 5.0))
        args = mf.make_parser().parse_args(["--aa", "x.pdb", "--backbone-restraints", "clipped"])
        self.assertEqual((args.backbone_restraints, args.backbone_restraint_inset), ("clipped", 5.0))
        for bad in ("-1", "25"):
            err = io.StringIO()
            with contextlib.redirect_stderr(err), self.assertRaises(SystemExit):
                mf.main(["--aa", str(AA_PDB), "--gmx", sys.executable, "--orientation", "none", "--backbone-restraints", "clipped",
                         "--backbone-restraint-inset", bad])
            self.assertIn("--backbone-restraint-inset", err.getvalue())


if __name__ == "__main__":
    unittest.main()
