"""Metamorphic tests: equivalent inputs must give the same placement within numerical tolerance."""
import unittest

import numpy as np

from .common import AA_PDB, CG_GRO, FORCEFIELD, MAPPING, mf


class PlacementInvariance(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.aa, _ = mf.read_all_atom(AA_PDB, FORCEFIELD)
        cls.raw, cls.box = mf.read_cg(CG_GRO)
        cls.cg = mf.classify_cg(cls.raw, MAPPING)
        cls.slab = mf.make_membrane_whole(cls.cg["membrane"], cls.box)
        cls.fit = mf.map_all_atom_to_cg(cls.aa, cls.cg["protein"], cls.box, cls.slab)
        cls.placed = mf.xyz_nm(cls.aa) @ cls.fit["R"].T + cls.fit["t"]

    def test_reference_mapping(self):
        self.assertEqual(self.fit["assigned"], [("A", 1), ("P", 2), ("R", 0)])
        self.assertEqual(self.fit["metrics"]["unmatched_aa_chains"], ["B", "C"])
        self.assertAlmostEqual(np.linalg.det(self.fit["R"]), 1.0, places=6)
        self.assertEqual(self.fit["metrics"]["backbone_pairs"], 600)
        self.assertTrue(self.fit["metrics"]["ambiguity"].startswith("none"))

    def test_chain_relabeling_and_residue_renumbering(self):
        swap = {"R": "X", "A": "Q", "P": "Z"}
        other = [dict(a, chain=swap.get(a["chain"], a["chain"]), resid=a["resid"] + 1000 if a["resname"] in mf.AMINO else a["resid"])
                 for a in self.aa]
        fit = mf.map_all_atom_to_cg(other, self.cg["protein"], self.box, self.slab)
        self.assertEqual(fit["assigned"], [("Q", 1), ("X", 0), ("Z", 2)])
        self.assertTrue(np.allclose(fit["R"], self.fit["R"]) and np.allclose(fit["t"], self.fit["t"]))

    def test_rigid_rotation_and_translation_of_the_all_atom_input(self):
        angle = 1.1
        Rz = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        M = Rz @ np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
        X = mf.xyz_nm(self.aa) @ M.T + np.array([50.0, -30.0, 12.0])
        moved = [dict(a, x=float(x), y=float(y), z=float(z)) for a, (x, y, z) in zip(self.aa, X)]
        fit = mf.map_all_atom_to_cg(moved, self.cg["protein"], self.box, self.slab)
        self.assertLess(np.abs(X @ fit["R"].T + fit["t"] - self.placed).max(), 1e-6)

    def test_periodic_wrapping_of_the_cg_input(self):
        shift = np.array([self.box[0] / 2, 0.0, self.box[2] * 0.45])
        wrapped = [dict(a, x=(a["x"] + shift[0]) % self.box[0], z=(a["z"] + shift[2]) % self.box[2]) for a in self.raw]
        cg = mf.classify_cg(wrapped, MAPPING)
        slab = mf.make_membrane_whole(cg["membrane"], self.box)
        fit = mf.map_all_atom_to_cg(self.aa, cg["protein"], self.box, slab)
        self.assertAlmostEqual(fit["metrics"]["core_rmsd_A"], self.fit["metrics"]["core_rmsd_A"], places=2)
        self.assertAlmostEqual(slab[1] - slab[0], self.slab[1] - self.slab[0], places=3)
        # the complex sits at the same place relative to the membrane, modulo the lattice
        placed = mf.xyz_nm(self.aa) @ fit["R"].T + fit["t"]
        delta = (placed - self.placed).mean(axis=0) / 10.0 - shift
        cell = np.array(self.box)
        self.assertLess(np.abs(delta - cell * np.round(delta / cell)).max(), 0.01)
        self.assertEqual(fit["metrics"]["bb_beads_in_bilayer"], self.fit["metrics"]["bb_beads_in_bilayer"])

    def test_cg_bead_order_within_the_file_is_irrelevant_for_membrane_molecules(self):
        protein = [a for r in self.cg["protein"] for a in r]
        rest = [a for a in self.raw if a["resname"] not in mf.MARTINI3_PROTEIN]
        residues = mf.cg_residues(rest)
        lipids = [r for r in residues if r[0]["resname"] not in ("W", "ION")]
        other = [r for r in residues if r[0]["resname"] in ("W", "ION")]
        reordered = protein + [a for r in other for a in r] + [a for r in lipids for a in r]
        cg = mf.classify_cg(reordered, MAPPING)
        self.assertEqual(cg["composition"], self.cg["composition"])
        self.assertEqual(cg["counts"], self.cg["counts"])

    def test_classification_is_independent_of_residue_numbering(self):
        renumbered = [dict(a, resid=(a["resid"] + 7) % 100000) for a in self.raw if a["resname"] not in mf.MARTINI3_PROTEIN]
        protein = [a for r in self.cg["protein"] for a in r]
        self.assertEqual(mf.classify_cg(protein + renumbered, MAPPING)["composition"], self.cg["composition"])


if __name__ == "__main__":
    unittest.main()
