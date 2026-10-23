"""The ligand is frozen in both minimizations: index group, mdp lines, the check after EM, and the audit's index rules."""
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from .common import mf

ATOM = lambda name, charge=0.0: {"atom": name, "type": "C", "charge": charge}


def topology(ligands=("LIGX",)):
    """Protein PROA (3 atoms), ligand LIGX (C1, H1, C2) x2, POPC (2 atoms) x2, TIP3 (3 atoms) x1."""
    return {"molecules": [("PROA", 1, [ATOM("N"), ATOM("CA"), ATOM("C")], "Protein_LIG"),
                          ("LIGX", 2, [ATOM("C1"), ATOM("H1"), ATOM("C2")], "Protein_LIG"),
                          ("POPC", 2, [ATOM("C1"), ATOM("C2")], "MEMB"),
                          ("TIP3", 1, [ATOM("OH2"), ATOM("H1"), ATOM("H2")], "SOL_ION")],
            "ligands": list(ligands), "charge": 0.0}


class IndexGroup(unittest.TestCase):
    def setUp(self):
        self.out = Path(tempfile.mkdtemp())

    def groups(self):
        import re
        text = (self.out / "index_ini.ndx").read_text()
        return {n: [int(v) for v in b.split()] for n, b in zip(re.findall(r"\[ (\S+) \]", text), re.split(r"\[ \S+ \]", text)[1:])}

    def test_ligand_atoms_get_their_own_trailing_group(self):
        sizes = mf.make_index(self.out, topology())
        groups = self.groups()
        self.assertEqual(list(groups), ["System", "Protein_LIG", "MEMB", "SOL_ION", "LIG"])
        self.assertEqual(groups["Protein_LIG"], [1, 2, 3, 4, 5, 6, 7, 8, 9])
        self.assertEqual(groups["LIG"], [4, 6, 7, 9])  # both copies, heavy atoms only: the hydrogens (5, 8) stay free
        self.assertEqual(sizes["LIG"], 4)
        self.assertEqual(sizes["Protein_LIG"] + sizes["MEMB"] + sizes["SOL_ION"], sizes["System"])

    def test_no_ligand_keeps_the_four_groups(self):
        sizes = mf.make_index(self.out, topology(ligands=()))
        self.assertEqual(list(self.groups()), ["System", "Protein_LIG", "MEMB", "SOL_ION"])
        self.assertNotIn("LIG", sizes)


class FrozenInMinimization(unittest.TestCase):
    def run_em(self, top, restraints):
        out = Path(tempfile.mkdtemp())
        with mock.patch.object(mf.minimization, "run_command"), \
                mock.patch.object(mf.minimization, "has_dihedral_restraints", return_value=restraints):
            mf.run_em(out, top, "gmx", 1, 10)
        return out

    def test_both_stages_freeze_the_ligand_group(self):
        out = self.run_em(topology(), True)
        for name in ("emres", "em"):
            mdp = (out / f"{name}.mdp").read_text()
            self.assertIn("freezegrps      = LIG\n", mdp)
            self.assertIn("freezedim       = Y Y Y\n", mdp)
        self.assertIn("define          = -DDIHRES", (out / "emres.mdp").read_text())

    def test_no_ligand_no_freeze(self):
        out = self.run_em(topology(ligands=()), False)
        self.assertNotIn("freeze", (out / "em.mdp").read_text())

    def test_the_ions_run_is_unchanged(self):
        self.assertNotIn("freeze", mf.EM_MDP.format(emtol=500, nsteps=0, coulombtype="PME", define="", freeze=""))


class LigandStaysPut(unittest.TestCase):
    def setUp(self):
        self.out = Path(tempfile.mkdtemp())
        atoms = [{"resid": i + 1, "resname": "X", "atom": f"A{i}", "x": 0.1 * i, "y": 1.0, "z": 1.0} for i in range(6)]
        self.before = (atoms, [3.0, 3.0, 3.0])
        mf.make_index(self.out, {"molecules": [("PROA", 1, [ATOM("a")] * 4, "Protein_LIG"), ("LIGX", 1, [ATOM("a")] * 2, "Protein_LIG")],
                                 "ligands": ["LIGX"]})
        self.index = {"LIG": 2}

    def xyz(self, shift=None):
        xyz = np.array([[0.1 * i, 1.0, 1.0] for i in range(6)])
        if shift is not None:
            xyz[4:] += shift
        return xyz

    def test_unmoved_ligand_passes_even_in_another_periodic_image(self):
        mf.minimization.check_ligand_frozen(self.out, self.index, self.xyz(), self.before)
        mf.minimization.check_ligand_frozen(self.out, self.index, self.xyz([3.0, 0.0, -3.0]), self.before)

    def test_receptor_atoms_may_move(self):
        xyz = self.xyz()
        xyz[:4] += 0.05
        mf.minimization.check_ligand_frozen(self.out, self.index, xyz, self.before)

    def test_moved_ligand_is_refused(self):
        with self.assertRaises(SystemExit) as caught:
            mf.minimization.check_ligand_frozen(self.out, self.index, self.xyz([0.0, 0.05, 0.0]), self.before)
        self.assertIn("frozen ligand moved during EM", str(caught.exception))

    def test_no_ligand_group_is_not_checked(self):
        mf.minimization.check_ligand_frozen(self.out, {}, self.xyz([1.0, 0.0, 0.0]), self.before)


if __name__ == "__main__":
    unittest.main()
