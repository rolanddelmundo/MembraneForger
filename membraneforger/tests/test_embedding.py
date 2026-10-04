"""Embedding tests: placing a complex into a membrane simulated with another protein, trimming the patch, lipid edits."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .common import AA_PDB, CG_GRO, FORCEFIELD, MAPPING, REPO, failure, mf, residue, small_system

KOR1 = REPO / "examples/preeq_cg_cellmem/KOR1_cg_cellmem.gro"


def load(cg_path):
    """The example complex plus one classified coarse-grained frame, its box and the bilayer z range."""
    aa = mf.read_all_atom(AA_PDB, FORCEFIELD)[0]
    beads, box = mf.read_cg(cg_path)
    cg = mf.classify_cg(beads, MAPPING)
    slab = mf.make_membrane_whole(cg["membrane"], box)
    return aa, cg, box, slab


class HydrophobicBelt(unittest.TestCase):
    def test_belt_of_the_fitted_example_sits_on_the_bilayer_midplane(self):
        aa, cg, box, slab = load(CG_GRO)
        fit = mf.map_all_atom_to_cg(aa, cg["protein"], box, slab)
        moved = mf.xyz_nm(aa) @ fit["R"].T + fit["t"]
        placed = [dict(a, x=float(x), y=float(y), z=float(z)) for a, (x, y, z) in zip(aa, moved)]
        belt = mf.hydrophobic_belt(placed)
        self.assertLess(abs(belt["z_a"] / 10.0 - 0.5 * (slab[0] + slab[1])), 0.6)
        self.assertGreater(belt["hydrophobic"], belt["charged"])

    def test_belt_moves_with_the_structure(self):
        aa = mf.read_all_atom(AA_PDB, FORCEFIELD)[0]
        base = mf.hydrophobic_belt(aa)["z_a"]
        shifted = mf.hydrophobic_belt([dict(a, z=a["z"] + 25.0) for a in aa])["z_a"]
        self.assertAlmostEqual(shifted - base, 25.0, delta=1.0)

    def test_too_few_residues_are_refused(self):
        one = [dict(atom="CA", resname="ALA", chain="A", resid=1, x=0.0, y=0.0, z=1.0)]
        self.assertIn("fewer than 10 protein residues", failure(mf.hydrophobic_belt, one))


class Embedding(unittest.TestCase):
    def test_complex_is_placed_in_the_frames_protein_hole_and_overlaps_are_removed(self):
        aa, cg, box, slab = load(KOR1)
        before = len(cg["membrane"])
        result = mf.embed_complex(aa, cg["protein"], cg["membrane"], box, slab)
        placed, kept = result["placed"], result["membrane"]
        removed = sum(result["removed"].values())
        self.assertEqual(len(kept) + removed, before)
        self.assertGreater(removed, 0)
        self.assertLess(removed, 0.25 * before)
        heavy = [a for a in placed if mf.element(a["atom"]) != "H"]
        tree = cKDTree(mf.wrap(mf.xyz_nm(heavy) / 10.0, box), boxsize=box)
        beads = np.array([[b["x"], b["y"], b["z"]] for mol in kept for b in mol["beads"]])
        self.assertFalse(any(tree.query_ball_point(mf.wrap(beads, box), mf.OVERLAP_NM)))
        centre_z = mf.hydrophobic_belt(placed)["z_a"] / 10.0
        self.assertLess(abs(centre_z - result["metrics"]["midplane_nm"]), 0.3)
        self.assertEqual(result["metrics"]["mode"], "embed")
        self.assertTrue(np.allclose(result["R"], np.eye(3)))

    def test_given_bilayer_centre_is_used(self):
        aa, cg, box, slab = load(KOR1)
        auto = mf.embed_complex(aa, cg["protein"], cg["membrane"], box, slab)
        given = mf.embed_complex(aa, cg["protein"], cg["membrane"], box, slab, bilayer_z_a=auto["metrics"]["bilayer_centre_input_A"] + 10.0)
        self.assertAlmostEqual(given["t"][2] - auto["t"][2], -10.0, places=3)
        self.assertTrue(given["metrics"]["hydrophobic_belt"].get("given"))

    def test_given_midplane_puts_the_bilayer_centre_on_it(self):
        aa, cg, box, slab = load(KOR1)
        default = mf.embed_complex(aa, cg["protein"], cg["membrane"], box, slab, bilayer_z_a=0.0)
        moved = mf.embed_complex(aa, cg["protein"], cg["membrane"], box, slab, bilayer_z_a=0.0, midplane_nm=default["metrics"]["midplane_nm"] + 0.5)
        self.assertAlmostEqual(moved["t"][2] - default["t"][2], 5.0, places=3)
        self.assertAlmostEqual(moved["metrics"]["midplane_nm"], default["metrics"]["midplane_nm"] + 0.5, places=3)

    def test_periodic_mean_crosses_the_boundary(self):
        self.assertAlmostEqual(mf.periodic_mean(np.array([0.1, 9.9]), 10.0), 0.0, places=6)
        self.assertAlmostEqual(mf.periodic_mean(np.array([4.0, 6.0]), 10.0), 5.0, places=6)


def grid_membrane(box, spacing=1.0):
    """A synthetic membrane: one POPC per grid point, each a column of beads 0.5 nm long in z."""
    membrane = []
    for i in range(int(box[0] / spacing)):
        for j in range(int(box[1] / spacing)):
            x, y = (i + 0.5) * spacing, (j + 0.5) * spacing
            beads = [{"resid": 1, "resname": "POPC", "atom": n, "chain": "", "x": x, "y": y, "z": 5.0 + 0.05 * k}
                     for k, n in enumerate(MAPPING["POPC"]["beads"])]
            membrane.append({"cg": "POPC", "cls": "phospholipid", "aa": "POPC", "beads": beads, "dropped": []})
    return membrane


class Trimming(unittest.TestCase):
    def test_trim_keeps_the_molecules_around_the_complex_and_centres_it(self):
        box = [10.0, 10.0, 12.0]
        placed = [{"atom": "CA", "resname": "ALA", "x": 70.0 + dx, "y": 70.0, "z": 50.0} for dx in (-1.0, 1.0)]
        result = mf.trim_membrane(grid_membrane(box), placed, box, (4.0, 4.0))
        self.assertEqual(result["box"], [4.0, 4.0, 12.0])
        self.assertEqual(len(result["membrane"]), 16)
        self.assertEqual(sum(result["removed"].values()), 84)
        for a in result["placed"]:
            self.assertAlmostEqual(a["y"], 20.0, places=6)
        self.assertAlmostEqual(np.mean([a["x"] for a in result["placed"]]), 20.0, places=6)
        for mol in result["membrane"]:
            for b in mol["beads"]:
                self.assertTrue(0.0 <= b["x"] < 4.0 and 0.0 <= b["y"] < 4.0)

    def test_trim_wraps_across_the_periodic_boundary(self):
        box = [10.0, 10.0, 12.0]
        placed = [{"atom": "CA", "resname": "ALA", "x": 5.0, "y": 5.0, "z": 50.0}]  # complex at the corner (0.5, 0.5) nm
        result = mf.trim_membrane(grid_membrane(box), placed, box, (4.0, 4.0))
        self.assertEqual(len(result["membrane"]), 16)

    def test_box_larger_than_the_patch_is_refused(self):
        box = [10.0, 10.0, 12.0]
        placed = [{"atom": "CA", "resname": "ALA", "x": 50.0, "y": 50.0, "z": 50.0}]
        self.assertIn("exceeds the membrane patch", failure(mf.trim_membrane, grid_membrane(box), placed, box, (12.0, 4.0)))


class LipidEdits(unittest.TestCase):
    def setUp(self):
        self.membrane = mf.classify_cg(small_system() + residue(8, "POPC", MAPPING["POPC"]["beads"]), MAPPING)["membrane"]

    def test_names_and_aliases(self):
        self.assertEqual(mf.lipid_name("DPG3"), "GM3")
        self.assertEqual(mf.lipid_name("pip2"), "SAP6")
        self.assertEqual(mf.lipid_name("chol"), "CHOL")
        self.assertIn("unknown lipid DPPC", failure(mf.lipid_name, "DPPC"))

    def test_delete_removes_every_molecule_of_the_species(self):
        kept, edits = mf.edit_lipids(self.membrane, ["CHOL"], None, MAPPING)
        self.assertEqual([m["cg"] for m in kept], ["POPC", "POPC"])
        self.assertEqual(edits["deleted"], {"CHOL": 1})

    def test_convert_renames_beads_in_place(self):
        kept, edits = mf.edit_lipids(self.membrane, ["POPC"], "DOPE", MAPPING)
        self.assertEqual(edits["converted_to_DOPE"], {"POPC": 2})
        dope = [m for m in kept if m["cg"] == "DOPE"]
        self.assertEqual(len(dope), 2)
        self.assertEqual([b["atom"] for b in dope[0]["beads"]], MAPPING["DOPE"]["beads"])
        self.assertEqual({b["resname"] for b in dope[0]["beads"]}, {"DOPE"})
        self.assertEqual(dope[0]["aa"], "DOPE")
        self.assertEqual(len(kept), 3)

    def test_refusals(self):
        self.assertIn("say which lipid to replace", failure(mf.edit_lipids, self.membrane, [], "POPC", MAPPING))
        self.assertIn("cannot turn CHOL into POPC", failure(mf.edit_lipids, self.membrane, ["CHOL"], "POPC", MAPPING))
        self.assertIn("no DOPS to remove", failure(mf.edit_lipids, self.membrane, ["DOPS"], None, MAPPING))
        self.assertIn("also named by --dellipid", failure(mf.edit_lipids, self.membrane, ["POPC"], "POPC", MAPPING))
        self.assertIn("every membrane molecule was removed", failure(mf.edit_lipids, self.membrane, ["POPC", "CHOL"], None, MAPPING))


class RequestedBoxZ(unittest.TestCase):
    def system(self):
        coords = [{"resid": 1, "resname": "ALA", "atom": "CA", "group": "Protein_LIG", "x": 10.0, "y": 10.0, "z": 50.0 + dz}
                  for dz in (-30.0, 30.0)]
        coords += [{"resid": 2, "resname": "POPC", "atom": "P", "group": "MEMB", "x": 10.0, "y": 10.0, "z": z} for z in (30.0, 70.0)]
        out = Path(tempfile.mkdtemp())
        system = {"cryst1": "CRYST1   40.000   40.000  100.000  90.00  90.00  90.00 P 1           1", "out": out, "name": "t"}
        return system, {"coords": coords}

    def test_requested_z_is_used_when_large_enough_and_refused_otherwise(self):
        system, topology = self.system()
        auto = mf.rebox_system(system, topology)
        self.assertAlmostEqual(auto[2], 2.0 * (3.0 + mf.SLAB_Z_PAD_NM))
        self.assertIn("is below the", failure(mf.rebox_system, system, topology, auto[2] - 1.0))
        self.assertAlmostEqual(mf.rebox_system(system, topology, 20.0)[2], 20.0)
