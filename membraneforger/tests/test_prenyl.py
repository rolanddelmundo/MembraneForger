"""E/Z configuration of the S-geranylgeranyl cysteine: detection, the rigid conversion, the input policy and the audit check.

The bundled GCGR-Gs input carries its Ggamma CYSG with all three CH=C(CH3) bonds in Z (dihedral 0 to 5 degrees), so it is the real test case."""
import copy
import tempfile
import unittest
from pathlib import Path

import numpy as np

from .common import AA_PDB, failure, mf
from membraneforger.audit import check_prenyl_configuration
from membraneforger.config import Settings


def cysg(atoms):
    return [a for a in atoms if a["resname"] == "CYSG"]


def coords(res):
    return {a["atom"]: np.array([a["x"], a["y"], a["z"]]) for a in res}


def bond_lengths(c):
    """Every bond of the prenyl chain, and each hydrogen to its nearest heavy atom."""
    pairs = [(a, b) for a, b in mf.prenyl.PRENYL_BONDS if a in c and b in c]
    heavy = [n for n in c if n in mf.prenyl.HEAVY_NAMES]
    out = [float(np.linalg.norm(c[a] - c[b])) for a, b in pairs]
    for n in c:
        if n not in mf.prenyl.HEAVY_NAMES:
            out.append(min(float(np.linalg.norm(c[n] - c[h])) for h in heavy))
    return np.array(out)


class Prenyl(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.atoms = mf.read_pdb(AA_PDB)[0]

    def setUp(self):
        self.atoms = copy.deepcopy(self.__class__.atoms)

    def test_the_bundled_input_is_all_z(self):
        (record,) = mf.review_prenyl(self.atoms)
        self.assertEqual(record["z_bonds"], ["C2=C3", "C7=C8", "C12=C13"])
        self.assertEqual(record["missing"], [])
        self.assertTrue(all(abs(d) < 10 for d in record["dihedral_deg"]), record["dihedral_deg"])

    def test_conversion_gives_e_and_keeps_every_bond_exactly(self):
        c = coords(cysg(self.atoms))
        new, converted = mf.convert_to_e(c)
        self.assertEqual(converted, ["C2=C3", "C7=C8", "C12=C13"])
        self.assertTrue(all(abs(d) > 170 for d in mf.prenyl_dihedrals(new)), mf.prenyl_dihedrals(new))
        self.assertLess(float(np.abs(bond_lengths(c) - bond_lengths(new)).max()), 1e-6)
        for name in ("N", "CA", "C", "O", "CB", "SG", "C1", "C2", "C3"):             # the proximal part and the two axis carbons stay
            if name in c:
                self.assertTrue(np.allclose(c[name], new[name]), name)
        self.assertFalse(np.allclose(c["C4"], new["C4"]))                          # both substituents of C3 swap sides
        again, converted_again = mf.convert_to_e(new)
        self.assertEqual(converted_again, [])
        self.assertTrue(all(np.allclose(new[n], again[n]) for n in new))

    def test_hydrogens_follow_their_carbon(self):
        c = coords(cysg(self.atoms))
        new, _ = mf.convert_to_e(c)
        hydrogens = [n for n in c if n not in mf.prenyl.HEAVY_NAMES and n[0] == "H"]
        self.assertGreater(len(hydrogens), 20)
        self.assertLess(float(np.abs(bond_lengths(c) - bond_lengths(new)).max()), 1e-6)

    def test_report_warns_and_leaves_the_atoms(self):
        before = copy.deepcopy(self.atoms)
        atoms, notes, record = mf.apply_prenyl_configuration(self.atoms, "report", "in.pdb")
        self.assertEqual(atoms, before)
        self.assertTrue(any(level == "WARN" and "Z configuration" in text for level, text in notes))
        self.assertEqual(record["z_residues"], ["C:CYSG68"] if record["z_residues"] else [])

    def test_yes_stops_the_build_and_names_the_fix(self):
        message = failure(mf.apply_prenyl_configuration, self.atoms, "yes", "in.pdb")
        self.assertIn("Z configuration", message)
        self.assertIn("--prenyl-configuration fix", message)

    def test_fix_converts_and_a_second_pass_finds_nothing(self):
        atoms, notes, record = mf.apply_prenyl_configuration(self.atoms, "fix", "in.pdb")
        self.assertTrue(record["converted"])
        self.assertTrue(any("converted from Z to E" in text for _, text in notes))
        self.assertFalse(any(r["z_bonds"] for r in mf.review_prenyl(atoms)))
        _, notes2, record2 = mf.apply_prenyl_configuration(atoms, "yes", "in.pdb")          # an E group passes the strict mode
        self.assertEqual(record2["z_residues"], [])

    def test_no_does_nothing(self):
        before = copy.deepcopy(self.atoms)
        atoms, notes, record = mf.apply_prenyl_configuration(self.atoms, "no", "in.pdb")
        self.assertEqual((atoms, notes, record), (before, [], {"mode": "no"}))

    def test_fix_refuses_a_conversion_that_lands_on_another_residue(self):
        new, _ = mf.convert_to_e(coords(cysg(self.atoms)))
        target = new["C10"]
        blocker = {"atom": "CA", "altloc": " ", "resname": "ALA", "chain": "Z", "resid": 1, "x": float(target[0] + 1.0), "y": float(target[1]),
                   "z": float(target[2]), "segid": "", "elem": "C"}
        message = failure(mf.apply_prenyl_configuration, self.atoms + [blocker], "fix", "in.pdb")
        self.assertIn("leaves a prenyl atom", message)

    def test_a_group_without_its_tail_is_reported_not_crashed(self):
        atoms = [a for a in self.atoms if not (a["resname"] == "CYSG" and a["atom"] in ("C16", "C17", "C18", "C19", "C20"))]
        (record,) = mf.review_prenyl(atoms)
        self.assertEqual(record["missing"], ["C16", "C17", "C18", "C19", "C20"])
        _, notes, _ = mf.apply_prenyl_configuration(atoms, "report", "in.pdb")
        self.assertTrue(any("lacks prenyl atom" in text for _, text in notes))


class AuditPrenyl(unittest.TestCase):
    def gro(self, directory, atoms):
        res = cysg(atoms)
        with (Path(directory) / "em.gro").open("w") as fh:
            fh.write(f"test\n{len(res)}\n")
            for i, a in enumerate(res, 1):
                fh.write(f"{68:5d}{'CYSG':<5s}{a['atom']:>5s}{i:5d}{a['x'] / 10:8.3f}{a['y'] / 10:8.3f}{a['z'] / 10:8.3f}\n")
            fh.write("   10.00000   10.00000   10.00000\n")

    def test_enforced_modes_fail_on_z_and_pass_on_e(self):
        atoms = mf.read_pdb(AA_PDB)[0]
        with tempfile.TemporaryDirectory() as tmp:
            self.gro(tmp, atoms)
            report, fails = check_prenyl_configuration(Path(tmp), "em.gro", Settings(prenyl_configuration="yes"))
            self.assertEqual((report["residues"], report["z_residues"], bool(fails)), (1, 1, True))
            report, fails = check_prenyl_configuration(Path(tmp), "em.gro", Settings(prenyl_configuration="report"))
            self.assertEqual((report["z_residues"], fails), (1, []))                                  # listed, not failed
            fixed, _, _ = mf.apply_prenyl_configuration(copy.deepcopy(atoms), "fix", "in.pdb")
            self.gro(tmp, fixed)
            report, fails = check_prenyl_configuration(Path(tmp), "em.gro", Settings(prenyl_configuration="yes"))
            self.assertEqual((report["z_residues"], fails), (0, []))


if __name__ == "__main__":
    unittest.main()
