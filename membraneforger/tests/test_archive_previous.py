"""--archive-previous: an earlier build's reports (and structures) are copied aside before --overwrite replaces it."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from .common import mf


def earlier_build(out, status="FAIL", started="2026-09-03T08:15:54-0700", passed=False):
    """A directory shaped like the output of a finished build (files only; the contents are marked with their own name)."""
    out.mkdir(exist_ok=True)
    (out / "run_manifest.json").write_text(json.dumps({"status": status, "started": started}))
    for name in ("membranebuilder.log", "audit.json", "backbone_check.json", "topol.top", "index_ini.ndx"):
        (out / name).write_text(f"{name} of the earlier build\n")
    (out / "int").mkdir(exist_ok=True)
    (out / "em.gro" if passed else out / "int" / "em.unverified.gro").write_text("structure of the earlier build\n")
    return out


class Archive(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.out = earlier_build(Path(self._tmp.name) / "run")

    def tearDown(self):
        self._tmp.cleanup()

    def test_reports_mode_copies_the_explanation_not_the_structures(self):
        folder = mf.archive_previous_outputs(self.out, "reports")
        self.assertEqual(folder.parent, self.out / mf.ARCHIVE_DIR)
        self.assertEqual(folder.name, "20260903T081554_FAIL")          # month 09: no "-0" splitting accident
        names = {p.name for p in folder.iterdir()}
        self.assertTrue({"membranebuilder.log", "run_manifest.json", "audit.json", "backbone_check.json", "ARCHIVED.json"} <= names)
        self.assertFalse(names & {"em.gro", "em.unverified.gro", "topol.top", "index_ini.ndx"})
        self.assertEqual((folder / "audit.json").read_text(), "audit.json of the earlier build\n")
        self.assertEqual(json.loads((folder / "ARCHIVED.json").read_text())["mode"], "reports")

    def test_all_mode_adds_the_structures_from_int_when_the_build_failed(self):
        folder = mf.archive_previous_outputs(self.out, "all")
        self.assertEqual((folder / "em.unverified.gro").read_text(), "structure of the earlier build\n")
        self.assertTrue((folder / "topol.top").is_file() and (folder / "index_ini.ndx").is_file())

    def test_all_mode_takes_em_gro_of_a_build_that_passed(self):
        out = earlier_build(Path(self._tmp.name) / "passed", status="PASS", passed=True)
        folder = mf.archive_previous_outputs(out, "all")
        self.assertTrue((folder / "em.gro").is_file())
        self.assertTrue(folder.name.endswith("_PASS"))

    def test_nothing_to_archive_means_no_folder(self):
        self.assertIsNone(mf.archive_previous_outputs(self.out, "no"))
        fresh = Path(self._tmp.name) / "fresh"
        fresh.mkdir()
        self.assertIsNone(mf.archive_previous_outputs(fresh, "reports"))
        self.assertFalse((fresh / mf.ARCHIVE_DIR).exists())

    def test_the_archive_survives_the_rerun_and_a_second_one_does_not_overwrite_the_first(self):
        first = mf.archive_previous_outputs(self.out, "reports")
        removed = mf.clear_stale_outputs(self.out)
        self.assertIn("run_manifest.json", removed)
        self.assertTrue((first / "run_manifest.json").is_file())
        earlier_build(self.out)                                    # the same start time again
        second = mf.archive_previous_outputs(self.out, "reports")
        self.assertNotEqual(first, second)
        self.assertTrue(second.name.endswith("_FAIL_2"))
        self.assertTrue((first / "membranebuilder.log").is_file() and (second / "membranebuilder.log").is_file())

    def test_a_manifest_that_cannot_be_read_still_archives_under_the_file_time(self):
        (self.out / "run_manifest.json").write_text("{not json")
        folder = mf.archive_previous_outputs(self.out, "reports")
        self.assertTrue(folder.name.endswith("_unknown"))
        self.assertTrue((folder / "membranebuilder.log").is_file())

    def test_a_copy_that_fails_raises_before_anything_is_replaced(self):
        with mock.patch("membraneforger.validation.shutil.copy2", side_effect=OSError("disk full")):
            with self.assertRaises(OSError):
                mf.archive_previous_outputs(self.out, "reports")
        self.assertTrue((self.out / "run_manifest.json").is_file() and (self.out / "audit.json").is_file())   # nothing was cleared

    def test_the_option_needs_overwrite(self):
        parser = mf.make_parser()
        args = parser.parse_args(["--aa", "x.pdb", "--archive-previous", "all"])
        self.assertEqual((args.archive_previous, args.overwrite), ("all", False))
        self.assertEqual(parser.parse_args(["--aa", "x.pdb"]).archive_previous, "no")


if __name__ == "__main__":
    unittest.main()
