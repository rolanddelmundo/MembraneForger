"""Command-line defaults: lengths in nm, output naming from the input structure, and refusing to overwrite earlier work."""
import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from .common import AA_PDB, FORCEFIELD, failure, mf


def parse(*argv):
    parser = mf.make_parser()
    args = parser.parse_args(["--aa", "x.pdb", *argv])
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        mf.resolve_lengths(args, parser)
    return args, err.getvalue()


class Lengths(unittest.TestCase):
    def test_conversion_is_ten_angstrom_per_nm(self):
        self.assertAlmostEqual(mf.angstrom_to_nm(125.0), 12.5)
        self.assertAlmostEqual(mf.nm_to_angstrom(12.5), 125.0)
        self.assertAlmostEqual(mf.nm_to_angstrom(mf.angstrom_to_nm(37.3)), 37.3)

    def test_nm_options_become_the_internal_angstrom_values(self):
        args, note = parse("--box-nm", "12", "12", "20", "--bilayer-z-nm", "6.5", "--fit-max-core-rmsd-nm", "0.4")
        self.assertEqual(args.box, [120.0, 120.0, 200.0])
        self.assertAlmostEqual(args.bilayer_z, 65.0)
        self.assertAlmostEqual(args.fit_max_core_rmsd, 4.0)
        self.assertEqual(note, "")                                                    # nothing to warn about

    def test_legacy_angstrom_spellings_still_work_and_say_what_they_were_read_as(self):
        args, note = parse("--box", "120", "120", "200", "--bilayer-z", "65", "--fit-max-core-rmsd", "4")
        self.assertEqual((args.box, args.bilayer_z, args.fit_max_core_rmsd), ([120.0, 120.0, 200.0], 65.0, 4.0))
        self.assertIn("--box is in angstrom", note)
        self.assertIn("--box-nm 12 12 20 nm", note)
        self.assertIn("--bilayer-z-nm 6.5 nm", note)
        self.assertIn("--fit-max-core-rmsd-nm 0.4 nm", note)

    def test_both_spellings_of_one_length_are_refused(self):
        parser = mf.make_parser()
        args = parser.parse_args(["--aa", "x.pdb", "--box-nm", "12", "12", "20", "--box", "120", "120", "200"])
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            mf.resolve_lengths(args, parser)

    def test_defaults_are_unchanged_in_angstrom(self):
        args, _ = parse()
        self.assertIsNone(args.box)
        self.assertIsNone(args.bilayer_z)
        self.assertEqual(args.fit_max_core_rmsd, mf.Settings().fit_max_core_rmsd_a)
        self.assertEqual(args.xy_buffer, mf.Settings().box_xy_buffer_nm)                # nm all along
        self.assertEqual(mf.make_parser().parse_args(["--aa", "x", "--xy-buffer-nm", "2"]).xy_buffer, 2.0)

    def test_the_box_header_is_nm(self):
        self.assertIsNone(mf.parse_box("auto"))
        self.assertEqual(mf.parse_box("14,14,22"), (140.0, 140.0, 220.0))
        self.assertIn("nm", failure(mf.parse_box, "14,14"))
        self.assertIn("nm", failure(mf.parse_box, "14,-1,22"))


class OutputNaming(unittest.TestCase):
    def build_with(self, *argv):
        sessions = []
        with mock.patch.object(mf.cli, "build", side_effect=lambda session, *a: sessions.append(session) or 0), \
                mock.patch.object(mf.cli, "run_sessions", side_effect=lambda session, *a: sessions.append(session) or 0), \
                tempfile.TemporaryDirectory() as tmp, mock.patch("pathlib.Path.cwd", return_value=Path(tmp).resolve()):
            code = mf.main(["--aa", str(AA_PDB), "--gmx", sys.executable, "--orientation", "none", *argv])
            return code, sessions[0], Path(tmp).resolve()

    def test_the_default_is_named_after_the_structure_whatever_the_membrane(self):
        outs = set()
        for cg in ("1", "2", "KOR_Membrane1", "GPR139_Membrane5"):
            code, session, cwd = self.build_with("--cg", cg)
            outs.add(session.out.name)
            self.assertEqual(session.out.parent, cwd)
        self.assertEqual(outs, {f"{AA_PDB.stem}_membraneforger"})

    def test_runs_use_the_same_name_and_a_folder_per_frame(self):
        code, session, cwd = self.build_with("--cg", "2", "--runs", "2")
        self.assertEqual(session.out, cwd / f"{AA_PDB.stem}_membraneforger")

    def test_an_explicit_path_wins(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, session, _ = self.build_with("--cg", "KOR_Membrane1", "--out", str(Path(tmp) / "mine"))
            self.assertEqual(session.out, (Path(tmp) / "mine").resolve())

    def test_default_out_helper(self):
        self.assertEqual(mf.default_out(Path("/data/6WHC.pdb"), Path("/work")), Path("/work/6WHC_membraneforger"))

    def test_runs_refuse_two_frames_with_one_folder(self):
        frames = [Path("/x/KOR_Membrane1_cg_cellmem.gro"), Path("/y/KOR_Membrane1_cg_cellmem.gro")]
        session = mf.Session(out=Path("/nonexistent/out"), name="t", gmx="gmx", forcefield=FORCEFIELD, data=None, python=sys.executable,
                             ntomp=1, nsteps=1)
        self.assertIn("share one output folder", failure(mf.run_sessions, session, AA_PDB, frames, []))


class RestartSafety(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "out"
        self.out.mkdir()

    def touch(self, *names):
        for name in names:
            path = self.out / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("x")

    def test_an_empty_or_new_directory_is_fine(self):
        mf.refuse_overwrite(self.out)
        mf.refuse_overwrite(self.out / "not_yet")
        self.touch(mf.LOG_NAME)                                                        # a log alone is not a build
        mf.refuse_overwrite(self.out)

    def test_an_earlier_build_needs_overwrite(self):
        self.touch("em.gro", "topol.top")
        message = failure(mf.refuse_overwrite, self.out)
        self.assertIn("already holds the output of an earlier build", message)
        self.assertIn("em.gro", message)
        self.assertIn("--overwrite", message)
        mf.refuse_overwrite(self.out, overwrite=True)

    def test_a_failed_build_leaves_outputs_too(self):
        self.touch("int/em.unverified.gro", "run_manifest.json")
        self.assertIn("earlier build", failure(mf.refuse_overwrite, self.out))

    def test_checkpoints_and_trajectories_need_the_restart_flag_as_well(self):
        self.touch("em.gro", "md/state.cpt", "md/md.xtc")
        self.assertEqual(mf.existing_outputs(self.out)["simulation"], ["md/md.xtc", "md/state.cpt"])
        self.assertIn("simulation restart files", failure(mf.refuse_overwrite, self.out))
        self.assertIn("simulation restart files", failure(mf.refuse_overwrite, self.out, overwrite=True))
        self.assertIn("simulation restart files", failure(mf.refuse_overwrite, self.out, overwrite_restart=True))
        mf.refuse_overwrite(self.out, overwrite=True, overwrite_restart=True)

    def test_restart_files_alone_are_refused_too(self):
        self.touch("state.cpt")
        self.assertIn("simulation restart files", failure(mf.refuse_overwrite, self.out))

    def test_the_minimizations_own_files_are_not_simulation_restart_files(self):
        self.touch("em.gro", "em.cpt", "em_prev.cpt", "em.trr", "emres.cpt", "emres_prev.cpt", "int/em.trr")
        self.assertEqual(mf.existing_outputs(self.out)["simulation"], [])
        mf.refuse_overwrite(self.out, overwrite=True)

    def test_clearing_a_build_never_deletes_simulation_files(self):
        self.touch("em.gro", "topol.top", "state.cpt", "state_prev.cpt", "md.xtc", "md.log")
        mf.clear_stale_outputs(self.out)
        self.assertFalse((self.out / "em.gro").exists())
        for name in ("state.cpt", "state_prev.cpt", "md.xtc", "md.log"):
            self.assertTrue((self.out / name).exists(), name)

    def test_the_membrane_source_does_not_count_as_an_earlier_build(self):
        self.touch("int/membrane.pdb")
        keep = self.out / "int" / "membrane.pdb"
        mf.refuse_overwrite(self.out, keep=keep)
        self.touch("int/boxed.gro")
        self.assertIn("earlier build", failure(mf.refuse_overwrite, self.out, keep=keep))

    def test_build_refuses_before_touching_anything(self):
        self.touch("em.gro", "topol.top", "state.cpt")
        session = mf.Session(out=self.out, name="t", gmx="gmx", forcefield=FORCEFIELD, data=None, python=sys.executable, ntomp=1, nsteps=1)
        with contextlib.redirect_stdout(io.StringIO()) as shown:
            code = mf.build(session, AA_PDB, None, None, [])
        self.assertEqual(code, 1)
        self.assertIn("simulation restart files", shown.getvalue())
        self.assertEqual(sorted(p.name for p in self.out.iterdir()), ["em.gro", "state.cpt", "topol.top"])   # no log, no deletion
        self.assertEqual((self.out / "em.gro").read_text(), "x")

    def test_runs_check_every_folder_before_the_first_build(self):
        (self.out / "KOR_Membrane2").mkdir()
        (self.out / "KOR_Membrane2" / "em.gro").write_text("x")
        session = mf.Session(out=self.out, name="t", gmx="gmx", forcefield=FORCEFIELD, data=None, python=sys.executable, ntomp=1, nsteps=1)
        frames = [Path("/x/KOR_Membrane1_cg_cellmem.gro"), Path("/x/KOR_Membrane2_cg_cellmem.gro")]
        with mock.patch.object(mf.cli, "build", return_value=0) as built:
            self.assertIn("earlier build", failure(mf.run_sessions, session, AA_PDB, frames, []))
            built.assert_not_called()

    def test_the_cli_flags_reach_the_session(self):
        sessions = []
        with mock.patch.object(mf.cli, "build", side_effect=lambda session, *a: sessions.append(session) or 0):
            mf.main(["--aa", str(AA_PDB), "--cg", "KOR_Membrane1", "--gmx", sys.executable, "--orientation", "none",
                     "--out", str(self.out), "--overwrite", "--overwrite-restart"])
            mf.main(["--aa", str(AA_PDB), "--cg", "KOR_Membrane1", "--gmx", sys.executable, "--orientation", "none", "--out", str(self.out)])
        self.assertEqual([(s.overwrite, s.overwrite_restart) for s in sessions], [(True, True), (False, False)])

    def test_a_gromacs_continuation_in_the_directory_is_not_blocked_by_this_code(self):
        # The build runs steepest descent only and never passes -cpi, -append or -noappend to mdrun, so GROMACS continuation
        # (mdrun -cpi state.cpt) in a build directory is a user command this code neither starts nor interferes with.
        source = Path(mf.minimization.__file__).read_text()
        for flag in ("-cpi", "-append", "-noappend", "-cpo"):
            self.assertNotIn(f'"{flag}"', source)
        self.assertNotIn("*.cpt", " ".join(mf.GENERATED))
        self.assertIn("*.cpt", mf.SIMULATION_PATTERNS)


if __name__ == "__main__":
    unittest.main()
