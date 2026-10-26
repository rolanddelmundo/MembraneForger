"""add_lipidated_residue appends to the private force field of the build: what setup_forcefield or an earlier residue put there survives."""
import tempfile
import unittest
from pathlib import Path

from .common import AA_PDB, FORCEFIELD, mf


class PrivateForceField(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)
        self.body = mf.read_pdb(AA_PDB)[0]                                   # the bundled GCGR complex carries a KTZ at the peptide
        mf.setup_forcefield(self.work, {"CYSG", "KTZ"}, FORCEFIELD)
        self.ff = self.work / "charmm36.ff"

    def earlier_additions(self):
        """What a residue added before this one leaves in the working copy (setup_forcefield's CYSM block, an earlier lipidated residue)."""
        for name, marker in (("merged.rtp", "\n[ EARLIER ]\n  [ atoms ]\n     N       NH1   -0.470  0\n"), ("ffbonded.itp", "\n; EARLIER 1 2\n")):
            text = (self.ff / name).read_text()
            (self.ff / name).unlink()                                        # the entry was a symlink to the pristine file
            (self.ff / name).write_text(text + marker)

    def test_earlier_additions_survive_and_the_residue_is_added_once(self):
        self.earlier_additions()
        rtp0, bonded0 = (self.ff / "merged.rtp").read_text(), (self.ff / "ffbonded.itp").read_text()
        note = mf.add_lipidated_residue(self.work, "KTZ", self.body, FORCEFIELD)
        self.assertIn("KTZ polymer RTP", note)
        rtp, bonded = (self.ff / "merged.rtp").read_text(), (self.ff / "ffbonded.itp").read_text()
        self.assertTrue(rtp.startswith(rtp0) and bonded.startswith(bonded0))   # appended, not rebuilt from the pristine files
        self.assertEqual((rtp.count("[ EARLIER ]"), rtp.count("[ KTZ ]")), (1, 1))
        self.assertEqual(bonded.count("; EARLIER 1 2"), 1)
        self.assertGreater(len(bonded), len(bonded0))                           # the supplemental CGenFF angle and dihedrals
        mf.add_lipidated_residue(self.work, "KTZ", self.body, FORCEFIELD)       # a second chain with the same residue changes nothing
        self.assertEqual((self.ff / "merged.rtp").read_text(), rtp)
        self.assertEqual((self.ff / "ffbonded.itp").read_text(), bonded)

    def test_without_earlier_additions_the_pristine_files_are_the_base(self):
        pristine = FORCEFIELD / "charmm36.ff"
        mf.add_lipidated_residue(self.work, "KTZ", self.body, FORCEFIELD)
        self.assertTrue((self.ff / "merged.rtp").read_text().startswith((pristine / "merged.rtp").read_text()))
        self.assertTrue((self.ff / "ffbonded.itp").read_text().startswith((pristine / "ffbonded.itp").read_text()))
        self.assertEqual((self.ff / "merged.rtp").read_text().count("[ KTZ ]"), 1)
        self.assertNotIn("[ KTZ ]", (pristine / "merged.rtp").read_text())


if __name__ == "__main__":
    unittest.main()
