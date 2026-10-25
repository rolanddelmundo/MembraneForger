"""Position restraints written for molecules without their own: a ligand follows the backbone schedule, a lipid the lipid one."""
import tempfile
import unittest
from pathlib import Path

from .common import FORCEFIELD, mf


class LigandRestraints(unittest.TestCase):
    def written(self, name: str, ligands: set) -> tuple[str, str]:
        """The posre file and the itp that molecule_itp writes into a fresh toppar/."""
        with tempfile.TemporaryDirectory() as tmp:
            toppar = Path(tmp)
            itp = mf.molecule_itp(name, toppar, FORCEFIELD, ligands)
            return (toppar / f"posre_{name}.itp").read_text(), itp.read_text()

    def test_a_ligand_without_restraints_gets_backbone_force_constants_on_its_heavy_atoms(self):
        for name in ("GDP", "MG"):
            posre, itp = self.written(name, {name})
            rows = [line.split() for line in posre.splitlines() if line and line[0] != "["]
            self.assertTrue(rows, name)
            self.assertTrue(all(r[2:] == ["POSRES_FC_BB"] * 3 for r in rows), name)
            self.assertNotIn("POSRES_FC_LIG", posre)
            self.assertIn('#ifdef POSRES\n#include "posre_%s.itp"\n#endif' % name, itp)
        data = mf.read_itp(FORCEFIELD / "toppar" / "GDP.itp")
        heavy = [a["nr"] for a in data["atoms"] if mf.element(a["atom"]) != "H"]
        self.assertEqual([int(line.split()[0]) for line in self.written("GDP", {"GDP"})[0].splitlines()[1:]], heavy)

    def test_the_same_molecule_not_given_as_a_ligand_keeps_the_lipid_schedule(self):
        posre, _ = self.written("GDP", set())
        self.assertIn("POSRES_FC_LIPID", posre)
        self.assertNotIn("POSRES_FC_BB", posre)


if __name__ == "__main__":
    unittest.main()
