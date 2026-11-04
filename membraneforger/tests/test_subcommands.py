"""The command word in front of the options: build (the default), check, slice, backmap, solvate and equilibration."""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from .common import AA_PDB, CG_GRO, FORCEFIELD, POPC, mf, residue


def quiet(function, *args):
    with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
        return function(*args)


class CommandWord(unittest.TestCase):
    def test_without_a_command_the_old_command_line_is_a_build(self):
        self.assertEqual(mf.split_command(["--aa", "x.pdb", "--out", "d"]), ("build", ["--aa", "x.pdb", "--out", "d"]))
        self.assertEqual(mf.split_command([]), ("build", []))

    def test_a_command_word_is_taken_off_the_front(self):
        for name in mf.COMMANDS:
            self.assertEqual(mf.split_command([name, "--aa", "x.pdb"]), (name, ["--aa", "x.pdb"]))

    def test_the_help_lists_every_command(self):
        text = mf.make_parser().format_help()
        for name in mf.COMMANDS:
            self.assertIn(f"\n  {name}", text)

    def session_for(self, *argv):
        sessions = []
        with mock.patch.object(mf.cli, "build", side_effect=lambda session, *a: sessions.append(session) or 0), \
                mock.patch.object(mf.cli, "run_sessions", side_effect=lambda session, *a: sessions.append(session) or 0), \
                tempfile.TemporaryDirectory() as tmp:
            code = mf.main([*argv, "--gmx", sys.executable, "--orientation", "none", "--out", str(Path(tmp) / "b")])
        self.assertEqual(code, 0)
        return sessions[0]

    def test_slice_and_backmap_stop_the_pipeline_after_their_stage(self):
        self.assertIsNone(self.session_for("--aa", str(AA_PDB), "--cg", str(CG_GRO)).stop_after)
        self.assertIsNone(self.session_for("build", "--aa", str(AA_PDB), "--cg", str(CG_GRO)).stop_after)
        self.assertEqual(self.session_for("slice", "--aa", str(AA_PDB), "--cg", str(CG_GRO)).stop_after, "slice")
        self.assertEqual(self.session_for("backmap", "--aa", str(AA_PDB), "--cg", "KOR_Membrane1").stop_after, "backmap")
        self.assertEqual(self.session_for("slice", "--aa", str(AA_PDB), "--cg", "2", "--runs", "2").stop_after, "slice")

    def test_check_is_the_check_flag(self):
        with mock.patch.object(mf.cli, "run_check", return_value=0) as checked, tempfile.TemporaryDirectory() as tmp:
            code = mf.main(["check", "--aa", str(AA_PDB), "--gmx", sys.executable, "--orientation", "none", "--out", str(Path(tmp) / "b")])
        self.assertEqual(code, 0)
        checked.assert_called_once()

    def test_solvate_needs_an_assembled_membrane_and_slice_refuses_one(self):
        with self.assertRaises(SystemExit):
            quiet(mf.main, ["solvate", "--aa", str(AA_PDB), "--gmx", sys.executable])
        with self.assertRaises(SystemExit):
            quiet(mf.main, ["slice", "--membrane", str(AA_PDB), "--gmx", sys.executable])
        with self.assertRaises(SystemExit):
            quiet(mf.main, ["backmap", "--membrane", str(AA_PDB), "--gmx", sys.executable])

    def test_solvate_reaches_build_with_the_membrane_and_no_stop(self):
        calls = []
        with mock.patch.object(mf.cli, "build", side_effect=lambda session, aa, cg, membrane, cmd: calls.append((session, aa, cg, membrane)) or 0), \
                tempfile.TemporaryDirectory() as tmp:
            code = mf.main(["solvate", "--membrane", str(AA_PDB), "--gmx", sys.executable, "--out", str(Path(tmp) / "b")])
        self.assertEqual(code, 0)
        session, aa, cg, membrane = calls[0]
        self.assertEqual((aa, cg, membrane, session.stop_after), (None, None, AA_PDB.resolve(), None))

    def test_equilibration_is_handed_to_its_own_parser(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out), self.assertRaises(SystemExit) as stop:
            mf.main(["equilibration", "--help"])
        self.assertEqual(stop.exception.code, 0)
        self.assertIn("--energy", out.getvalue())


def fake_prepared(n_lipids=3):
    """What prepare_inputs returns, for a three-lipid membrane around a two-atom complex."""
    membrane = []
    for i in range(n_lipids):
        beads = [dict(b, resname="POPC") for b in residue(i + 1, "POPC", POPC)]
        membrane.append({"cg": "POPC", "cls": "lipid", "aa": "POPC", "beads": beads, "dropped": []})
    placed = [{"atom": "CA", "altloc": " ", "resname": "ALA", "chain": "R", "resid": 1, "x": 12.0, "y": 13.0, "z": 14.0, "segid": "", "elem": "C"},
              {"atom": "N", "altloc": " ", "resname": "ALA", "chain": "R", "resid": 1, "x": 11.0, "y": 13.0, "z": 14.0, "segid": "", "elem": "N"}]
    return {"placed": placed, "membrane": membrane, "box": [5.0, 6.0, 7.0], "mapping": {}, "composition": {"POPC": n_lipids}, "ligands": set()}


class SliceOutputs(unittest.TestCase):
    def test_the_cut_membrane_and_the_placed_complex_are_written_in_the_cut_box(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            written = mf.write_slice_outputs(out, fake_prepared())
            beads, box = mf.read_gro(out / "slice_cg.gro")
            atoms, cryst1 = mf.read_pdb(out / "placed_complex.pdb")
        self.assertEqual((written["membrane_molecules"], written["beads"], written["complex_atoms"]), (3, 3 * len(POPC), 2))
        self.assertEqual(len(beads), 3 * len(POPC))
        self.assertEqual({b["resname"] for b in beads}, {"POPC"})
        self.assertEqual([round(v, 3) for v in box], [5.0, 6.0, 7.0])
        self.assertEqual(len(atoms), 2)
        self.assertIn("50.000   60.000   70.000", cryst1)

    def test_slice_cg_gro_and_placed_complex_pdb_are_generated_files(self):
        self.assertIn("slice_cg.gro", mf.GENERATED)
        self.assertIn("placed_complex.pdb", mf.GENERATED)


class StoppedBuild(unittest.TestCase):
    """build() with stop_after: the pipeline ends after the stage, PASS, with the stage outputs at the top level."""

    def run_build(self, stop_after, prepared):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        out = Path(tmp.name) / "b"
        session = mf.Session(out=out, name="b", gmx=sys.executable, forcefield=FORCEFIELD, data=None, python=sys.executable,
                             ntomp=1, nsteps=1, stop_after=stop_after)
        with mock.patch.object(mf.pipeline, "gromacs_version", return_value="none"), \
                mock.patch.object(mf.pipeline, "prepare_inputs", return_value=prepared), \
                mock.patch.object(mf.pipeline, "finish_system", side_effect=AssertionError("must not run")) as finished, \
                contextlib.redirect_stdout(io.StringIO()):
            code = mf.build(session, AA_PDB, CG_GRO, None, ["x"])
        finished.assert_not_called()
        return code, out

    def test_slice_stops_after_writing_the_cut(self):
        code, out = self.run_build("slice", fake_prepared())
        self.assertEqual(code, 0)
        self.assertTrue((out / "slice_cg.gro").is_file() and (out / "placed_complex.pdb").is_file())
        manifest = json.loads((out / "run_manifest.json").read_text())
        self.assertEqual(manifest["status"], "PASS")
        self.assertEqual((manifest["stop_after"], manifest["slice_outputs"]["membrane_molecules"]), ("slice", 3))
        log = (out / mf.LOG_NAME).read_text().splitlines()
        self.assertTrue(log[-1].startswith("PASS:") and "slice_cg.gro" in log[-1], log[-1])

    def test_backmap_stops_with_membrane_pdb_at_the_top_level(self):
        prepared = fake_prepared()

        def assemble(session, prepared, seed):
            (session.out / "membrane.pdb").write_text("CRYST1\nEND\n")
            session.record["backmap_isomer_review"] = {}
            return session.out / "membrane.pdb"
        staged = {"solute_atoms": 2, "rings": {"pierced": []}, "contact": {"distance_nm": 1.0, "atoms": ["a", "b"]}, "overlaps": None}
        with mock.patch.object(mf.pipeline, "backmap_and_assemble", side_effect=assemble), \
                mock.patch.object(mf.pipeline, "build_topology_and_box", return_value=staged), \
                mock.patch.object(mf.pipeline, "backmap_verdict", return_value={"accept": True, "problem": None, "stage": "rings", "counts": {}}):
            code, out = self.run_build("backmap", prepared)
        self.assertEqual(code, 0)
        self.assertTrue((out / "membrane.pdb").is_file(), "membrane.pdb is the result and stays at the top level")
        self.assertFalse((out / "int" / "membrane.pdb").exists())
        self.assertEqual(json.loads((out / "run_manifest.json").read_text())["status"], "PASS")
