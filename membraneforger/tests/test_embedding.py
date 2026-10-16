"""Embedding tests: placing a complex into a membrane simulated with another protein, trimming the patch, lipid edits."""
import tempfile
import unittest
from collections import Counter
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


GPR1 = REPO / "examples/preeq_cg_cellmem/GPR1_cg_cellmem.gro"
FAST = {"relax_steps": 0, "max_void_nm3": float("inf")}  # placement only: no push, no pocket gate


def heavy_tree(placed, box):
    heavy = [a for a in placed if mf.element(a["atom"]) != "H"]
    return heavy, cKDTree(mf.wrap(mf.xyz_nm(heavy) / 10.0, box), boxsize=box)


def lipid_anchor_segments(atom):
    """True for the G protein's membrane anchors in the example: the geranylgeranylated Ggamma C-terminus and the Galpha N-terminus."""
    return (atom["chain"] == "C" and atom["resid"] >= 55) or (atom["chain"] == "A" and atom["resid"] <= 31)


class Embedding(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.aa, cls.cg, cls.box, cls.slab = load(KOR1)
        cls.result = mf.embed_complex(cls.aa, cls.cg["protein"], cls.cg["membrane"], cls.box, cls.slab)

    def test_complex_is_placed_in_the_frames_protein_hole_and_room_is_made(self):
        before = len(self.cg["membrane"])
        result, box = self.result, self.box
        placed, kept = result["placed"], result["membrane"]
        removed = sum(result["removed"].values())
        self.assertEqual(len(kept) + removed, before)
        self.assertEqual(sorted(result["removed_indices"]), result["removed_indices"])
        self.assertEqual(len(result["removed_indices"]), removed)
        self.assertGreater(removed, 0)
        self.assertLess(removed, 0.1 * before)
        _, tree = heavy_tree(placed, box)
        beads = np.array([[b["x"], b["y"], b["z"]] for mol in kept for b in mol["beads"]])
        self.assertFalse(any(tree.query_ball_point(mf.wrap(beads, box), mf.HARD_CORE_NM)))
        self.assertGreaterEqual(result["metrics"]["closest_bead_to_complex_nm"], mf.HARD_CORE_NM)
        gone = set(result["removed_indices"])
        originals = [mol for k, mol in enumerate(self.cg["membrane"]) if k not in gone]
        worst_shape = worst_move = 0.0
        for old, mol in zip(originals, kept):  # same lipids in the same order; beads and xyz agree; shape is kept
            self.assertEqual([b["atom"] for b in old["beads"]], [b["atom"] for b in mol["beads"]])
            self.assertTrue(np.allclose(mol["xyz"], [[b["x"], b["y"], b["z"]] for b in mol["beads"]]))
            a, b = np.asarray(old["xyz"]), mol["xyz"]
            i, j = np.triu_indices(len(a), 1)
            worst_shape = max(worst_shape, np.abs(np.linalg.norm(a[i] - a[j], axis=1) - np.linalg.norm(b[i] - b[j], axis=1)).max())
            worst_move = max(worst_move, np.linalg.norm(a - b, axis=1).max())
        self.assertLess(worst_shape, 0.1)
        self.assertLess(worst_move, 1.0)
        self.assertAlmostEqual(result["metrics"]["relaxation"]["max_intramolecular_distance_change_nm"], worst_shape, places=3)
        self.assertAlmostEqual(result["metrics"]["relaxation"]["max_bead_displacement_nm"], worst_move, places=3)
        centre_z = mf.hydrophobic_belt(placed)["z_a"] / 10.0
        self.assertLess(abs(centre_z - result["metrics"]["midplane_nm"]), 0.3)
        self.assertEqual(result["metrics"]["mode"], "embed")
        self.assertTrue(np.allclose(result["R"], np.eye(3)))
        self.assertLessEqual(result["metrics"]["voids"]["largest_pocket_nm3"], mf.MAX_VOID_NM3)

    def test_input_membrane_is_not_modified(self):
        fresh = load(KOR1)[1]["membrane"]
        self.assertEqual(len(fresh), len(self.cg["membrane"]))
        for a, b in zip(fresh, self.cg["membrane"]):
            self.assertTrue(np.allclose(a["xyz"], b["xyz"]))

    def test_given_bilayer_centre_is_used(self):
        aa, cg, box, slab = self.aa, self.cg, self.box, self.slab
        auto = mf.embed_complex(aa, cg["protein"], cg["membrane"], box, slab, **FAST)
        given = mf.embed_complex(aa, cg["protein"], cg["membrane"], box, slab, bilayer_z_a=auto["metrics"]["bilayer_centre_input_A"] + 10.0, **FAST)
        self.assertAlmostEqual(given["t"][2] - auto["t"][2], -10.0, places=3)
        self.assertTrue(given["metrics"]["hydrophobic_belt"].get("given"))

    def test_given_midplane_puts_the_bilayer_centre_on_it(self):
        aa, cg, box, slab = self.aa, self.cg, self.box, self.slab
        default = mf.embed_complex(aa, cg["protein"], cg["membrane"], box, slab, bilayer_z_a=0.0, **FAST)
        moved = mf.embed_complex(aa, cg["protein"], cg["membrane"], box, slab, bilayer_z_a=0.0,
                                 midplane_nm=default["metrics"]["midplane_nm"] + 0.5, **FAST)
        self.assertAlmostEqual(moved["t"][2] - default["t"][2], 5.0, places=3)
        self.assertAlmostEqual(moved["metrics"]["midplane_nm"], default["metrics"]["midplane_nm"] + 0.5, places=3)

class LipidAnchoredComplex(unittest.TestCase):
    """The receptor-Gs complex in a GPR139 frame: the Ggamma geranylgeranyl chain and the Galpha N-terminus reach into
    the inner leaflet. Removing every lipid within 0.40 nm of them emptied a column of that leaflet (the old rule)."""

    @classmethod
    def setUpClass(cls):
        cls.aa, cls.cg, cls.box, cls.slab = load(GPR1)
        cls.result = mf.embed_complex(cls.aa, cls.cg["protein"], cls.cg["membrane"], cls.box, cls.slab)
        cls.midplane = cls.result["metrics"]["midplane_nm"]
        cls.heavy, cls.tree = heavy_tree(cls.result["placed"], cls.box)
        touching = lambda mol: any(cls.tree.query_ball_point(mf.wrap(np.asarray(mol["xyz"]), cls.box), mf.OVERLAP_NM))
        cls.old_rule = [mol for mol in cls.cg["membrane"] if not touching(mol)]

    def test_the_old_contact_rule_leaves_an_inner_leaflet_pocket_that_the_gate_would_stop(self):
        protein = mf.xyz_nm(self.heavy) / 10.0
        old = mf.membrane_voids(self.old_rule, protein, self.box, self.midplane)
        self.assertGreater(old["lower"]["largest_pocket_nm3"], 2 * mf.MAX_VOID_NM3)
        self.assertGreater(len(self.cg["membrane"]) - len(self.old_rule), 80)

    def test_no_inner_leaflet_pocket_after_embedding(self):
        voids = self.result["metrics"]["voids"]
        self.assertLessEqual(voids["lower"]["largest_pocket_nm3"], mf.MAX_VOID_NM3)
        self.assertLessEqual(voids["upper"]["largest_pocket_nm3"], mf.MAX_VOID_NM3)
        again = mf.membrane_voids(self.result["membrane"], mf.xyz_nm(self.heavy) / 10.0, self.box, self.midplane)
        self.assertEqual(again["largest_pocket_nm3"], voids["largest_pocket_nm3"])

    def test_lipid_anchors_displace_lipids_instead_of_deleting_them(self):
        removed = set(self.result["removed_indices"])
        self.assertLess(len(removed), 0.6 * (len(self.cg["membrane"]) - len(self.old_rule)))
        anchors = np.array([lipid_anchor_segments(a) for a in self.heavy])
        protein = mf.xyz_nm(self.heavy) / 10.0
        rest = cKDTree(mf.wrap(protein[~anchors], self.box), boxsize=self.box)
        only_anchors = [k for k, mol in enumerate(self.cg["membrane"])
                        if any(self.tree.query_ball_point(mf.wrap(np.asarray(mol["xyz"]), self.box), mf.OVERLAP_NM))
                        and not any(rest.query_ball_point(mf.wrap(np.asarray(mol["xyz"]), self.box), mf.OVERLAP_NM))]
        self.assertGreater(len(only_anchors), 25)  # the old rule deleted all of these
        self.assertLess(len(removed & set(only_anchors)), 0.25 * len(only_anchors))
        moved = {k for k in only_anchors if k not in removed}
        self.assertTrue(moved)

    def test_bookkeeping_by_species_and_leaflet(self):
        result = self.result
        self.assertEqual(sum(result["metrics"]["removed_by_leaflet"].values()), sum(result["removed"].values()))
        before = Counter(m["cg"] for m in self.cg["membrane"])
        after = Counter(m["cg"] for m in result["membrane"])
        self.assertEqual(before - after, result["removed"])
        self.assertEqual(Counter(self.cg["membrane"][k]["cg"] for k in result["removed_indices"]), result["removed"])

    def test_a_pocket_above_the_limit_stops_the_build(self):
        # Without the push the hard-core rule alone deletes lipids around the anchors and leaves a pocket.
        message = failure(mf.embed_complex, self.aa, self.cg["protein"], self.cg["membrane"], self.box, self.slab, relax_steps=0)
        self.assertIn("empty pocket", message)
        self.assertIn("lower leaflet", message)


TM6 = range(349, 377)  # chain R of the example receptor: one helix crossing the bilayer, CA z from -18 to +15 A


def single_helix(aa):
    """TM6 of the example receptor alone: a stand-in for a single-pass membrane protein. Returns it and the belt z (A)."""
    return [a for a in aa if a["chain"] == "R" and a["resid"] in TM6], mf.hydrophobic_belt(aa)["z_a"]


class SinglePassProteinAwayFromTheHole(unittest.TestCase):
    """One transmembrane helix in a frame equilibrated around a GPCR: the receptor's hole cannot be filled by it, so the
    'free' site puts it on the unbroken bilayer farthest from the frame's protein."""

    @classmethod
    def setUpClass(cls):
        aa, cls.cg, cls.box, cls.slab = load(GPR1)
        cls.helix, cls.centre_a = single_helix(aa)
        cls.result = mf.embed_complex(cls.helix, cls.cg["protein"], cls.cg["membrane"], cls.box, cls.slab,
                                      bilayer_z_a=cls.centre_a, site="free")

    def test_the_helix_is_one_transmembrane_segment(self):
        ca = [a["z"] - self.centre_a for a in self.helix if a["atom"] == "CA"]
        self.assertEqual(len(ca), len(TM6))
        self.assertLess(ca[0], -15.0)
        self.assertGreater(ca[-1], 10.0)
        self.assertGreater(np.corrcoef(np.arange(len(ca)), ca)[0, 1], 0.95)  # climbs through the bilayer once

    def test_free_site_is_far_from_the_frames_protein(self):
        site, distance = mf.free_site(self.cg["protein"], self.slab, self.box)
        self.assertGreater(distance, 6.0)  # an 18.3 nm cell around a ~4 nm receptor
        beads = mf.frame_protein_beads(self.cg["protein"], self.slab, self.box)
        gap = np.abs((beads[:, :2] - site + 0.5 * np.array(self.box[:2])) % np.array(self.box[:2]) - 0.5 * np.array(self.box[:2]))
        self.assertAlmostEqual(float(np.linalg.norm(gap, axis=1).min()), distance, places=6)

    def test_the_site_keeps_the_slice_window_clear_of_the_receptor(self):
        metrics = self.result["metrics"]
        self.assertEqual(metrics["site"], "free")
        self.assertGreaterEqual(metrics["slice_window_clearance_nm"], mf.FRAME_PROTEIN_CLEARANCE_NM)
        self.assertGreater(metrics["distance_to_frame_protein_nm"], metrics["slice_window_clearance_nm"])
        xyz = mf.xyz_nm(self.helix) / 10.0
        tm = np.abs(xyz[:, 2] - self.centre_a / 10.0) <= mf.BELT_NM / 2.0
        window = mf.slice_window(xyz, xyz[tm, :2].mean(axis=0), np.array(self.box[:2]), 1.0, 0.5)
        site, clear = mf.free_site(self.cg["protein"], self.slab, self.box, window)
        self.assertTrue(np.allclose(metrics["target_xy_nm"], np.round(site, 3)))
        self.assertAlmostEqual(metrics["slice_window_clearance_nm"], round(clear, 3))

    def test_a_cut_too_large_to_clear_the_receptor_is_refused_before_the_push(self):
        message = failure(mf.embed_complex, self.helix, self.cg["protein"], self.cg["membrane"], self.box, self.slab,
                          bilayer_z_a=self.centre_a, site="free", requested_xy_nm=(self.box[0] - 1.0, self.box[1] - 1.0))
        self.assertIn("no place in this", message)
        self.assertIn("--embed-site hole", message)

    def test_the_push_makes_room_with_few_removals_and_no_pocket(self):
        result = self.result
        self.assertLess(sum(result["removed"].values()), 15)
        self.assertGreater(result["metrics"]["lipids_touching_complex"], sum(result["removed"].values()))
        self.assertLessEqual(result["metrics"]["voids"]["largest_pocket_nm3"], mf.MAX_VOID_NM3)
        self.assertGreaterEqual(result["metrics"]["closest_bead_to_complex_nm"], mf.HARD_CORE_NM)

    def test_the_receptor_hole_stops_the_build_and_names_the_free_site(self):
        message = failure(mf.embed_complex, self.helix, self.cg["protein"], self.cg["membrane"], self.box, self.slab,
                          bilayer_z_a=self.centre_a, relax_steps=300)
        self.assertIn("empty pocket", message)
        self.assertIn("--embed-site free", message)

    def test_unknown_site_is_refused(self):
        self.assertIn("unknown embedding site", failure(mf.embed_complex, self.helix, self.cg["protein"], self.cg["membrane"],
                                                        self.box, self.slab, bilayer_z_a=self.centre_a, site="middle"))


class FreeSiteHelpers(unittest.TestCase):
    box = [10.0, 10.0, 12.0]
    slab = (4.0, 8.0)

    def protein(self, *xyz):
        return [[{"atom": "BB", "x": x, "y": y, "z": z} for x, y, z in xyz]]

    def test_without_protein_in_the_bilayer_the_patch_centre_is_used(self):
        site, distance = mf.free_site(self.protein((1.0, 1.0, 11.0)), self.slab, self.box)  # above the bilayer only
        self.assertTrue(np.allclose(site, [5.0, 5.0]))
        self.assertIsNone(distance)

    def test_the_site_is_opposite_the_protein_under_periodicity(self):
        site, distance = mf.free_site(self.protein((0.5, 0.5, 6.0)), self.slab, self.box, grid_nm=0.5)
        self.assertTrue(np.allclose(site, [5.25, 5.25]) or np.allclose(site, [5.75, 5.75]))
        self.assertAlmostEqual(distance, float(np.hypot(4.75, 4.75)), places=6)

    def test_protein_beads_inside_the_slice_window_are_counted(self):
        protein = self.protein((1.0, 1.0, 6.0), (6.0, 6.0, 6.0), (9.5, 1.0, 6.0))
        placed = [{"x": 50.0, "y": 50.0}]
        cut = lambda lower, size: {"box": [size[0], size[1], 12.0], "placed": [{"x": 50.0 - 10 * lower[0], "y": 50.0 - 10 * lower[1]}]}
        count = lambda lower, size, clearance: mf.frame_protein_in_slice(protein, self.slab, self.box, placed, cut(lower, size), clearance)
        self.assertEqual(count((4.0, 4.0), (3.0, 3.0), 0.0), 1)
        self.assertEqual(count((3.0, 3.0), (2.0, 2.0), 0.0), 0)
        self.assertEqual(count((9.0, 0.5), (3.0, 1.0), 0.0), 2)  # across the periodic edge
        self.assertEqual(count((3.0, 0.0), (2.0, 10.0), 0.0), 0)
        self.assertEqual(count((0.0, 0.0), (10.0, 10.0), 0.0), 3)  # uncut: the whole cell
        # the receptor's hole is wider than its bead centres: a bead just outside the window counts with a clearance
        self.assertEqual(count((3.0, 3.0), (2.5, 2.5), 0.0), 0)                  # (6, 6) is 0.5 nm off the edge
        self.assertEqual(count((3.0, 3.0), (2.5, 2.5), 1.0), 1)
        self.assertEqual(mf.frame_protein_in_slice(protein, self.slab, self.box, placed, cut((3.0, 3.0), (2.5, 2.5))), 1)  # default 1 nm
        self.assertEqual(count((3.0, 3.0), (1.5, 1.5), 1.0), 0)                  # 1.5 nm off
        self.assertIn("comes within 1.0 nm of 1 bead(s) of the frame's own protein",
                      failure(mf.refuse_frame_protein_in_slice, protein, self.slab, self.box, placed, cut((3.0, 3.0), (2.5, 2.5))))
        self.assertEqual(mf.refuse_frame_protein_in_slice(protein, self.slab, self.box, placed, cut((3.0, 3.0), (1.5, 1.5))), 0)

    def test_clearance_is_euclidean_at_a_corner_as_in_the_placement_check(self):
        corner = self.protein((5.8, 5.8, 6.0))  # 0.8 nm past both edges of the window's corner: 1.13 nm away
        edge = self.protein((5.8, 4.0, 6.0))    # 0.8 nm past one edge: 0.8 nm away
        placed = [{"x": 50.0, "y": 50.0}]
        cut = {"box": [2.0, 2.0, 12.0], "placed": [{"x": 50.0 - 30.0, "y": 50.0 - 30.0}]}  # window [3, 5) x [3, 5)
        self.assertEqual(mf.frame_protein_in_slice(corner, self.slab, self.box, placed, cut), 0)
        self.assertEqual(mf.frame_protein_in_slice(edge, self.slab, self.box, placed, cut), 1)
        distance = mf.rectangle_distance(np.array([[5.8, 5.8], [5.8, 4.0]]), np.array([3.0, 3.0]), np.array([2.0, 2.0]),
                                         np.array(self.box[:2]))[0]
        self.assertTrue(np.allclose(distance, [np.hypot(0.8, 0.8), 0.8]))
        # the placement check measures the same distance: one candidate point (5, 5) with the window [3, 5) x [3, 5) there
        self.assertAlmostEqual(mf.free_site(corner, self.slab, self.box, window=((-2.0, -2.0), (0.0, 0.0)), grid_nm=10.0)[1],
                               float(np.hypot(0.8, 0.8)), places=9)

    def test_protein_beads_are_imaged_onto_the_bilayer_in_z(self):
        # a bilayer made whole across the z boundary (slab 9-13 nm in a 12 nm cell) and a protein bead left at z = 0.5
        beads = mf.frame_protein_beads(self.protein((2.0, 2.0, 0.5), (2.0, 2.0, 6.0)), (9.0, 13.0), self.box)
        self.assertTrue(np.allclose(beads, [[2.0, 2.0, 12.5]]))
        self.assertTrue(np.allclose(mf.image_onto_slab(np.array([[0.0, 0.0, -11.0]]), (4.0, 8.0), 12.0), [[0.0, 0.0, 1.0]]))

    def test_the_window_decides_the_site(self):
        protein = self.protein((5.0, 5.0, 6.0))
        point, _ = mf.free_site(protein, self.slab, self.box, grid_nm=0.5)
        # a window 8 nm wide in x and 1 nm in y leaves at most 1 nm of room in x but 4.5 nm in y: the site moves away
        # from the receptor mostly along y (on the 0.5 nm grid: 0.75 nm gap in x, 4.25 nm in y)
        site, clear = mf.free_site(protein, self.slab, self.box, window=((-4.0, -0.5), (4.0, 0.5)), grid_nm=0.5)
        self.assertAlmostEqual(abs((site[1] - 5.0 + 5.0) % 10.0 - 5.0), 4.25 + 0.5, places=6)
        self.assertAlmostEqual(abs((site[0] - 5.0 + 5.0) % 10.0 - 5.0), 0.75 + 4.0, places=6)
        self.assertAlmostEqual(clear, float(np.hypot(0.75, 4.25)), places=6)
        self.assertGreater(np.hypot(*(point - 5.0)), 6.0)

    def test_slice_window_mirrors_the_slicer(self):
        xyz = np.array([[1.0, 2.0, 0.0], [3.0, 2.5, 0.0]])
        lower, upper = mf.slice_window(xyz, np.array([2.0, 2.0]), np.array([10.0, 10.0]), 1.0, 0.5)
        # extent 2 x 0.5 nm centred at (2, 2.25); + 1 nm buffer and 0.5 nm slack on each side = 5 x 3.5 nm
        self.assertTrue(np.allclose(lower, [-2.5, -1.5]) and np.allclose(upper, [2.5, 2.0]))
        lower, upper = mf.slice_window(xyz, np.array([2.0, 2.0]), np.array([10.0, 10.0]), 1.0, 0.5, requested_xy_nm=(6.0, 10.0))
        self.assertTrue(np.allclose(upper - lower, [7.0, 10.0]))  # requested x plus the slack; y is the whole cell


def lattice_membrane(box, spacing=0.4, half_thickness=2.0, midplane=5.0):
    """A dense synthetic bilayer: per leaflet, one 'lipid' per lattice column, a PO4 bead on top and tail beads down to the midplane."""
    names = ("PO4", "GL1", "C1A", "C2A", "C3A", "C4A")
    membrane = []
    for sign in (-1.0, 1.0):
        for i in range(int(round(box[0] / spacing))):
            for j in range(int(round(box[1] / spacing))):
                x, y = (i + 0.5) * spacing, (j + 0.5) * spacing
                zs = [midplane + sign * (half_thickness - spacing * k) for k in range(int(half_thickness / spacing))]
                beads = [{"resid": 1, "resname": "POPC", "atom": names[k], "chain": "", "x": x, "y": y, "z": z} for k, z in enumerate(zs)]
                membrane.append({"cg": "POPC", "cls": "phospholipid", "aa": "POPC", "beads": beads, "dropped": [],
                                 "xyz": np.array([[x, y, z] for z in zs])})
    return membrane


class Voids(unittest.TestCase):
    box = [8.0, 8.0, 10.0]

    def test_an_intact_bilayer_has_no_pocket(self):
        voids = mf.membrane_voids(lattice_membrane(self.box), np.zeros((0, 3)), self.box, 5.0)
        self.assertEqual(voids["largest_pocket_nm3"], 0.0)
        self.assertEqual(voids["lower"]["empty_points"] + voids["upper"]["empty_points"], 0)

    def test_missing_lipids_make_a_pocket_in_their_leaflet_only(self):
        membrane = lattice_membrane(self.box)
        centre = np.array([0.2, 4.2])  # across the periodic edge in x
        gone = lambda m: m["xyz"][0, 2] < 5.0 and np.all(np.abs((m["xyz"][0, :2] - centre + 4.0) % 8.0 - 4.0) < 1.0)
        holed = [m for m in membrane if not gone(m)]
        self.assertEqual(len(membrane) - len(holed), 25)
        voids = mf.membrane_voids(holed, np.zeros((0, 3)), self.box, 5.0)
        self.assertGreater(voids["lower"]["largest_pocket_nm3"], 0.5)
        self.assertEqual(voids["upper"]["largest_pocket_nm3"], 0.0)
        self.assertEqual(voids["lower"]["pockets"], 1)
        x, y, z = voids["lower"]["largest_pocket_centre_nm"]
        self.assertLess(min(abs(x - 0.2), abs(x - 8.2)), 0.15)
        self.assertAlmostEqual(y, 4.2, delta=0.15)
        self.assertLess(z, 5.0)

    def test_solute_in_the_hole_fills_it(self):
        membrane = lattice_membrane(self.box)
        gone = lambda m: m["xyz"][0, 2] < 5.0 and np.all(np.abs(m["xyz"][0, :2] - 4.2) < 1.0)
        holed = [m for m in membrane if not gone(m)]
        g = np.arange(3.2, 5.21, 0.3)
        solute = np.array([[x, y, z] for x in g for y in g for z in np.arange(3.0, 5.0, 0.3)])
        self.assertEqual(mf.membrane_voids(holed, solute, self.box, 5.0)["largest_pocket_nm3"], 0.0)


class ProteinPush(unittest.TestCase):
    def test_relaxation_pushes_a_lipid_off_the_complex_only_when_asked(self):
        lipid = np.array([[2.0, 2.0, 5.0 - 0.4 * k] for k in range(5)])
        protein = np.array([[2.1, 2.0, 4.2]])
        seam = mf.relax_seam([lipid], np.array([6.0, 6.0]), 10.0, [False, False], protein, 0.30)
        self.assertEqual(seam["steps"], 1)  # seam mode: nothing was created by a cut, so it stops at once
        self.assertTrue(np.allclose(seam["xyz"][0], lipid))
        pushed = mf.relax_seam([lipid], np.array([6.0, 6.0]), 10.0, [False, False], protein, 0.30, 400,
                               protein_repulsion_nm=0.45, protein_clearance_nm=0.40)
        self.assertTrue(pushed["converged"])
        self.assertGreaterEqual(pushed["closest_bead_to_complex_nm"], 0.40)
        self.assertLess(pushed["max_intramolecular_distance_change_nm"], 0.05)


class PeriodicMean(unittest.TestCase):
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
    def test_too_small_z_is_refused_before_backmapping(self):
        placed = [{"atom": "CA", "resname": "ALA", "x": 0.0, "y": 0.0, "z": z} for z in (20.0, 80.0)]
        slab = (3.0, 7.0)  # centre 5 nm; complex 2-8 nm -> 2 x (3 + 1.5) = 9 nm
        self.assertAlmostEqual(mf.check_box_z(placed, slab, 9.0), 9.0)
        self.assertIn("below the 90 A", failure(mf.check_box_z, placed, slab, 8.0))

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
