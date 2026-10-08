"""BOX = auto slicing of the coarse-grained membrane: window, whole molecules, periodic images, seam, user box, pipeline.

None of these tests backmaps anything: they cut the real coarse-grained frames and stop before mstool runs.
"""
import copy
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from .common import AA_PDB, CG_GRO, FORCEFIELD, MAPPING, REPO, failure, mf

OPM_6WHC = Path(__file__).parent / "data" / "6whc_opm_chainR.pdb"
KOR_FRAME = REPO / "examples" / "preeq_cg_cellmem" / "KOR1_cg_cellmem.gro"
GPR139_FRAME = REPO / "examples" / "preeq_cg_cellmem" / "GPR3_cg_cellmem.gro"
CELL = [20.0, 20.0, 12.0]
MIDPLANE = 6.0


def lipid(name, x, y, z=MIDPLANE, spread=0.0):
    """A two-bead 'lipid' (PO4 head and a tail bead) with its head at x, y, z; spread makes it wider in x."""
    beads = [{"atom": "PO4", "x": x, "y": y, "z": z + 1.0}, {"atom": "C1A", "x": x + spread, "y": y, "z": z - 1.0}]
    return {"cg": name, "aa": name, "cls": "phospholipid", "beads": beads,
            "xyz": np.array([[b["x"], b["y"], b["z"]] for b in beads])}


def solute(x_nm, y_nm):
    """Two solute atoms spanning the given extents (nm), in Angstrom as the placed complex is."""
    return [{"atom": "CA", "resname": "ALA", "chain": "A", "resid": 1, "x": 10.0 * x_nm[0], "y": 10.0 * y_nm[0], "z": 60.0},
            {"atom": "CA", "resname": "ALA", "chain": "A", "resid": 2, "x": 10.0 * x_nm[1], "y": 10.0 * y_nm[1], "z": 60.0}]


def fill(n_per_leaflet=40):
    """Enough lipids, in both leaflets, to pass the minimum-lipids check; positions far from the test windows."""
    return [lipid("POPC", 0.5 + 0.4 * k, 0.5, MIDPLANE + (3.0 if k % 2 else -3.0)) for k in range(2 * n_per_leaflet)]


SETTINGS = mf.Settings(slice_min_lipids=4)
SMALL = mf.Settings(slice_min_lipids=2)


class WindowAndMolecules(unittest.TestCase):
    def test_window_is_the_solute_extent_plus_the_buffer(self):
        low, high = np.array([8.0, 9.0]), np.array([12.0, 11.0])
        windows = mf.crop_windows(low, high, np.array([20.0, 20.0]), 1.0)
        self.assertTrue(np.allclose(windows["lower"], [7.0, 8.0]) and np.allclose(windows["size"], [6.0, 4.0]))
        self.assertTrue(all(w["cropped"] for w in windows["axes"]))

    def test_an_axis_as_wide_as_the_cell_is_not_cut(self):
        windows = mf.crop_windows(np.array([1.0, 9.0]), np.array([19.0, 11.0]), np.array([20.0, 20.0]), 1.0)
        self.assertEqual([w["cropped"] for w in windows["axes"]], [False, True])
        self.assertEqual(windows["size"].tolist(), [20.0, 4.0])

    def test_a_lipid_is_kept_by_its_headgroup_anchor_even_when_a_tail_bead_crosses_the_window(self):
        inside = [lipid("POPC", 9.0, 10.0, MIDPLANE + 3), lipid("POPC", 10.0, 10.0, MIDPLANE - 3)]
        straddling = lipid("DOPC", 12.9, 10.0, spread=0.3)       # PO4 anchor at x = 12.9 inside, the tail bead at 13.2 outside
        far = lipid("POPE", 2.0, 10.0)
        above = lipid("POPS", 9.0, 12.5)                          # anchor outside the y window 8..12
        cut = mf.slice_membrane_cg(inside + [straddling, far, above] + fill(), solute((8, 12), (9, 11)), CELL, MIDPLANE, SMALL)
        self.assertEqual(cut["report"]["composition_after"], {"DOPC": 1, "POPC": 2})  # the fill lipids sit at y = 0.5, outside
        self.assertEqual(cut["report"]["lipids_outside_window"], len(fill()) + 2)
        dopc = next(m for m in cut["membrane"] if m["cg"] == "DOPC")
        self.assertEqual(len(dopc["beads"]), 2)                                       # the complete residue, tail bead included
        self.assertAlmostEqual(dopc["beads"][1]["x"], 13.2 - 7.0)                     # the tail still crosses the new cell edge
        selection = cut["report"]["selection"]
        self.assertEqual((selection["kept_by_anchor_rule"], selection["kept_by_all_beads_inside_rule"]), (3, 2))
        self.assertEqual(selection["rejected_by_old_rule_with_anchor_inside_by_species"], {"DOPC": 1})

    def test_a_lipid_whose_anchor_is_outside_is_not_kept_even_when_a_tail_bead_is_inside(self):
        outside_head = lipid("DOPE", 13.1, 10.0, spread=-0.4)    # anchor at 13.1 beyond x = 13, tail bead at 12.7 inside
        cut = mf.slice_membrane_cg([outside_head, lipid("POPC", 9.0, 10.0, MIDPLANE + 3), lipid("POPE", 10.0, 10.0, MIDPLANE - 3)],
                                   solute((8, 12), (9, 11)), CELL, MIDPLANE, SMALL)
        self.assertEqual(cut["report"]["composition_after"], {"POPC": 1, "POPE": 1})

    def test_half_open_window_counts_a_boundary_lipid_once(self):
        on_edge = lipid("DOPC", 7.0, 10.0, MIDPLANE + 3)          # anchor exactly on the lower edge: inside (lower <= a)
        on_far_edge = lipid("DOPE", 13.0, 10.0, MIDPLANE - 3)     # anchor exactly on the upper edge: outside (a < lower + width)
        cut = mf.slice_membrane_cg([on_edge, on_far_edge, lipid("POPC", 9.0, 10.0, MIDPLANE + 3), lipid("POPE", 10.0, 10.0, MIDPLANE - 3)],
                                   solute((8, 12), (9, 11)), CELL, MIDPLANE, SMALL)
        self.assertEqual(cut["report"]["composition_after"], {"DOPC": 1, "POPC": 1, "POPE": 1})
        self.assertEqual(len(cut["membrane"]), 3)                                     # no duplicate from the periodic image

    def test_the_crop_is_deterministic(self):
        members = [lipid("POPC", 9.0 + 0.3 * k, 10.0 + 0.1 * (k % 3), MIDPLANE + (3 if k % 2 else -3)) for k in range(12)]
        first = mf.slice_membrane_cg(members, solute((8, 12), (9, 11)), CELL, MIDPLANE, SMALL)
        second = mf.slice_membrane_cg(copy.deepcopy(members), solute((8, 12), (9, 11)), CELL, MIDPLANE, SMALL)
        self.assertEqual(first["report"]["composition_after"], second["report"]["composition_after"])
        self.assertTrue(all(np.allclose(a["xyz"], b["xyz"]) for a, b in zip(first["membrane"], second["membrane"])))

    def test_window_origin_new_cell_and_complex_shift(self):
        placed = solute((8, 12), (9, 11))
        members = [lipid("POPC", 9.0, 10.0, MIDPLANE + 3), lipid("POPE", 11.0, 9.5, MIDPLANE - 3)]
        cut = mf.slice_membrane_cg(members, placed, CELL, MIDPLANE, SMALL)
        self.assertEqual(cut["box"], [6.0, 4.0, 12.0])
        first = cut["membrane"][0]["beads"][0]
        self.assertAlmostEqual(first["x"], 9.0 - 7.0)
        self.assertAlmostEqual(first["y"], 10.0 - 8.0)
        self.assertAlmostEqual(first["z"], MIDPLANE + 4.0)  # z is untouched
        moved = np.array([[a["x"], a["y"], a["z"]] for a in cut["placed"]])
        original = np.array([[a["x"], a["y"], a["z"]] for a in placed])
        self.assertTrue(np.allclose(moved - original, [-70.0, -80.0, 0.0]))     # one rigid translation, in Angstrom
        self.assertAlmostEqual(np.linalg.norm(moved[0] - moved[1]), np.linalg.norm(original[0] - original[1]))
        self.assertEqual(cut["report"]["membrane_margin_around_complex_nm"], {"low_x_nm": 1.0, "high_x_nm": 1.0, "low_y_nm": 1.0, "high_y_nm": 1.0})

    def test_the_inputs_are_not_modified(self):
        members = [lipid("POPC", 9.0, 10.0, MIDPLANE + 3), lipid("POPE", 11.0, 9.5, MIDPLANE - 3)]
        before, placed = copy.deepcopy(members), solute((8, 12), (9, 11))
        reference = copy.deepcopy(placed)
        mf.slice_membrane_cg(members, placed, CELL, MIDPLANE, SMALL)
        self.assertEqual(members[0]["beads"], before[0]["beads"])
        self.assertEqual(placed, reference)

    def test_periodic_image_nearest_the_complex_is_used(self):
        # the complex sits at the cell edge (x 18.2..19.6): the window is 17.2..20.6 and a lipid at x = 0.3 belongs to it
        wrapped = lipid("POPC", 0.3, 10.0, MIDPLANE + 3)
        below = lipid("POPE", 18.5, 10.0, MIDPLANE - 3)
        cut = mf.slice_membrane_cg([wrapped, below], solute((18.2, 19.6), (9, 11)), CELL, MIDPLANE, SMALL)
        self.assertEqual(cut["report"]["composition_after"], {"POPC": 1, "POPE": 1})
        by_name = {m["cg"]: m for m in cut["membrane"]}
        self.assertAlmostEqual(by_name["POPC"]["beads"][0]["x"], 20.3 - 17.2)      # imaged to x = 20.3, then shifted
        self.assertAlmostEqual(by_name["POPE"]["beads"][0]["x"], 18.5 - 17.2)

    def test_nothing_is_cut_when_the_window_reaches_the_cell(self):
        members = [lipid("POPC", 3.0, 3.0, MIDPLANE + 3), lipid("POPE", 15.0, 15.0, MIDPLANE - 3), lipid("DOPC", 1.0, 18.0, MIDPLANE - 3),
                   lipid("DOPE", 17.0, 2.0, MIDPLANE + 3)]
        cut = mf.slice_membrane_cg(members, solute((1, 19), (1, 19)), CELL, MIDPLANE, SMALL)
        self.assertFalse(cut["report"]["cropped"])
        self.assertEqual(cut["box"], CELL)
        self.assertEqual(len(cut["membrane"]), 4)
        self.assertEqual(cut["report"]["seam"]["lipids_removed"], 0)

    def test_leaflets_and_composition_are_reported(self):
        members = [lipid("POPC", 9.0, 10.0, MIDPLANE + 3), lipid("POPC", 10.0, 10.0, MIDPLANE + 3), lipid("POPE", 11.0, 10.0, MIDPLANE - 3),
                   lipid("DOPS", 12.0, 9.5, MIDPLANE - 3)]
        cut = mf.slice_membrane_cg(members, solute((8, 12), (9, 11)), CELL, MIDPLANE, SETTINGS)
        self.assertEqual(cut["report"]["leaflets_after"], {"lower": {"DOPS": 1, "POPE": 1}, "upper": {"POPC": 2}})
        self.assertEqual(cut["report"]["lipids_before"], 4)

    def test_a_membrane_that_is_too_small_or_one_sided_is_refused(self):
        one_sided = [lipid("POPC", 9.0 + 0.3 * k, 10.0, MIDPLANE + 3) for k in range(6)]
        self.assertIn("leaves 6 lipids", failure(mf.slice_membrane_cg, one_sided, solute((8, 12), (9, 11)), CELL, MIDPLANE, SETTINGS))
        few = [lipid("POPC", 9.0, 10.0, MIDPLANE + 3), lipid("POPE", 10.0, 10.0, MIDPLANE - 3)]
        message = failure(mf.slice_membrane_cg, few, solute((8, 12), (9, 11)), CELL, MIDPLANE, mf.Settings())
        self.assertIn("at least 30 are needed", message)
        self.assertIn("--xy-buffer", message)


class Seam(unittest.TestCase):
    def test_overlaps_created_only_by_the_new_periodicity_are_relaxed_not_deleted(self):
        # window x 7..13: heads at x = 7.05 and 12.95 are 0.1 nm apart across the NEW boundary, 5.9 nm apart in the source
        members = [lipid("POPC", 7.05, 10.0, MIDPLANE + 3), lipid("POPE", 12.95, 10.0, MIDPLANE + 3),
                   lipid("DOPC", 9.0, 9.0, MIDPLANE - 3), lipid("DOPE", 10.0, 11.0, MIDPLANE - 3)]
        cut = mf.slice_membrane_cg(members, solute((8, 12), (9, 11)), CELL, MIDPLANE, SMALL)
        seam = cut["report"]["seam"]
        self.assertEqual(seam["lipids_removed"], 0)                                   # nothing is deleted
        self.assertEqual(len(cut["membrane"]), 4)
        self.assertTrue(seam["relaxation"]["converged"])
        self.assertEqual(seam["pairs_created_by_new_periodicity"], 0)                 # no pair closer than the threshold is left
        self.assertEqual(seam["relaxation"]["lipids_moved"], 2)                        # only the two seam lipids moved
        self.assertLess(seam["relaxation"]["max_bead_displacement_nm"], 0.3)
        points = np.vstack([m["xyz"] for m in cut["membrane"]])
        owners = np.concatenate([np.full(len(m["xyz"]), i) for i, m in enumerate(cut["membrane"])])
        from scipy.spatial import cKDTree
        box = np.array(cut["box"])
        pairs = cKDTree(np.mod(points - points.min(axis=0), box), boxsize=box).query_pairs(SETTINGS.seam_min_bead_nm, output_type="ndarray")
        self.assertFalse((owners[pairs[:, 0]] != owners[pairs[:, 1]]).any())          # no seam clash remains
        for mol, before in zip(cut["membrane"], members):                               # every lipid kept its shape
            length = lambda m: float(np.linalg.norm(m["xyz"][0] - m["xyz"][1]))
            self.assertAlmostEqual(length(mol), length(before), places=1)

    def test_a_hard_core_overlap_that_cannot_be_relaxed_costs_one_lipid(self):
        members = [lipid("POPC", 7.05, 10.0, MIDPLANE + 3), lipid("POPE", 12.95, 10.0, MIDPLANE + 3),
                   lipid("DOPC", 9.0, 9.0, MIDPLANE - 3), lipid("DOPE", 10.0, 11.0, MIDPLANE - 3)]
        cut = mf.slice_membrane_cg(members, solute((8, 12), (9, 11)), CELL, MIDPLANE, mf.Settings(slice_min_lipids=2, seam_relax_steps=0))
        seam = cut["report"]["seam"]
        self.assertEqual((seam["pairs_created_by_new_periodicity"], seam["lipids_removed"]), (2, 1))
        self.assertEqual(len(cut["membrane"]), 3)

    def test_contacts_that_were_already_close_in_the_source_are_left_alone(self):
        members = [lipid("POPC", 9.0, 10.0, MIDPLANE + 3), lipid("POPE", 9.1, 10.0, MIDPLANE + 3),   # 0.1 nm apart in the source
                   lipid("DOPC", 11.0, 9.5, MIDPLANE - 3), lipid("DOPE", 10.0, 11.0, MIDPLANE - 3)]
        cut = mf.slice_membrane_cg(members, solute((8, 12), (9, 11)), CELL, MIDPLANE, SMALL)
        seam = cut["report"]["seam"]
        self.assertEqual((seam["pairs_created_by_new_periodicity"], seam["lipids_removed"]), (0, 0))
        self.assertGreater(seam["pairs_already_close_in_source_frame"], 0)
        self.assertEqual(len(cut["membrane"]), 4)

    def test_pairs_across_an_uncut_axis_keep_their_source_periodicity(self):
        # x is not cut (window = cell), y is: beads at x = 0.05 and 19.95 are neighbours in the source cell, not a seam
        members = [lipid("POPC", 0.05, 10.0, MIDPLANE + 3), lipid("POPE", 19.95, 10.0, MIDPLANE + 3),
                   lipid("DOPC", 9.0, 9.0, MIDPLANE - 3), lipid("DOPE", 10.0, 11.0, MIDPLANE - 3)]
        cut = mf.slice_membrane_cg(members, solute((1, 19), (9, 11)), CELL, MIDPLANE, SMALL)
        self.assertEqual(cut["report"]["seam"]["lipids_removed"], 0)
        self.assertEqual(len(cut["membrane"]), 4)


class UserBox(unittest.TestCase):
    def test_requested_xy_cuts_the_membrane_around_the_complex(self):
        members = [lipid("POPC", 9.0, 10.0, MIDPLANE + 3), lipid("POPE", 11.0, 9.5, MIDPLANE - 3), lipid("DOPC", 14.2, 10.0, MIDPLANE - 3)]
        cut = mf.slice_membrane_cg(members, solute((9, 11), (9.5, 10.5)), CELL, MIDPLANE, SMALL, requested=(8.0, 6.0, 12.0))
        self.assertEqual(cut["box"][:2], [8.0, 6.0])
        self.assertEqual(cut["report"]["mode"], "user")
        self.assertEqual(cut["report"]["requested_xy_nm"], [8.0, 6.0])
        self.assertEqual(cut["report"]["composition_after"], {"POPC": 1, "POPE": 1})   # DOPC at x = 14.2 is outside the 6..14 window

    def test_requested_box_limits(self):
        args = ([lipid("POPC", 9.0, 10.0, MIDPLANE + 3), lipid("POPE", 11.0, 9.5, MIDPLANE - 3)], solute((9, 11), (9.5, 10.5)), CELL, MIDPLANE, SMALL)
        self.assertIn("larger than the 20.000 nm", failure(mf.slice_membrane_cg, *args, requested=(25.0, 6.0, 12.0)))
        self.assertIn("leaves less than 0.5 nm", failure(mf.slice_membrane_cg, *args, requested=(2.5, 6.0, 12.0)))
        whole = mf.slice_membrane_cg(*args, requested=(20.0, 20.0, 12.0))
        self.assertFalse(whole["report"]["cropped"])
        self.assertEqual(whole["box"], CELL)


class RealFrames(unittest.TestCase):
    """The real coarse-grained frames, cut but never backmapped."""

    @classmethod
    def setUpClass(cls):
        cls.settings = mf.Settings()
        aa, _ = mf.read_all_atom(AA_PDB, FORCEFIELD)
        raw, cls.box = mf.read_cg(CG_GRO)
        cg = mf.classify_cg(raw, MAPPING)
        slab = mf.make_membrane_whole(cg["membrane"], cls.box)
        cls.midplane = mf.bilayer_midplane(cg["membrane"], slab)[0]
        fit = mf.map_all_atom_to_cg(aa, cg["protein"], cls.box, slab, cls.settings)
        cls.membrane, cls.placed = cg["membrane"], mf.apply_transform(aa, fit["R"], fit["t"])

    def check_whole(self, cut):
        beads = {m["cg"]: len(m["beads"]) for m in self.membrane}
        for mol in cut["membrane"]:
            self.assertEqual(len(mol["beads"]), beads[mol["cg"]])                 # no lipid lost a bead
        for axis, window in enumerate(cut["report"]["axes"]):
            if window["cropped"]:
                anchors = np.array([mf.anchor_xyz(m["xyz"], [b["atom"] for b in m["beads"]], mf.anchor_bead(m["cg"])) for m in cut["membrane"]])
                self.assertGreaterEqual(anchors[:, axis].min(), 0.0)                 # every anchor inside the half-open cell
                self.assertLess(anchors[:, axis].max(), cut["box"][axis])
                low = min(m["xyz"][:, axis].min() for m in cut["membrane"])
                high = max(m["xyz"][:, axis].max() for m in cut["membrane"])
                self.assertGreater(low, -2.0)                                          # tails cross the edge by less than a lipid length
                self.assertLess(high, cut["box"][axis] + 2.0)

    def test_6whc_example(self):
        cut = mf.slice_membrane_cg(self.membrane, self.placed, self.box, self.midplane, self.settings)
        report = cut["report"]
        self.assertTrue(report["axes"][0]["cropped"] and not report["axes"][1]["cropped"])  # the complex nearly spans y
        self.assertLess(cut["box"][0], self.box[0])
        self.assertEqual(cut["box"][1:], self.box[1:])
        self.assertEqual(report["lipids_before"], len(self.membrane))
        self.assertEqual(report["lipids_before"] - report["lipids_after"],
                         report["lipids_outside_window"] + report["seam"]["lipids_removed"])
        self.assertGreater(report["selection"]["rejected_by_old_rule_with_anchor_inside"], 0)   # the old rule would have lost lipids
        self.assertEqual(sum(report["composition_after"].values()), report["lipids_after"])
        self.assertTrue(set(report["composition_after"]) <= set(report["composition_before"]))
        self.assertTrue(all(report["composition_after"][k] <= v for k, v in report["composition_before"].items() if k in report["composition_after"]))
        self.assertEqual(report["seam"]["lipids_removed"], 0)
        self.assertLess(report["area_fraction_kept"], 1.0)
        self.check_whole(cut)
        # the complex moved rigidly and stays inside the new cell with the buffer on each cropped side
        moved = mf.xyz_nm(cut["placed"]) / 10.0
        self.assertGreaterEqual(moved[:, 0].min(), self.settings.box_xy_buffer_nm - 1e-6)
        self.assertLessEqual(moved[:, 0].max(), cut["box"][0] - self.settings.box_xy_buffer_nm + 1e-6)
        before = mf.xyz_nm(self.placed)
        after = mf.xyz_nm(cut["placed"])
        self.assertAlmostEqual(float(np.linalg.norm(before[0] - before[-1])), float(np.linalg.norm(after[0] - after[-1])), places=6)

    def test_a_larger_buffer_keeps_more_membrane(self):
        small = mf.slice_membrane_cg(self.membrane, self.placed, self.box, self.midplane, mf.Settings(box_xy_buffer_nm=0.5))
        large = mf.slice_membrane_cg(self.membrane, self.placed, self.box, self.midplane, mf.Settings(box_xy_buffer_nm=1.0))
        self.assertLess(small["report"]["lipids_after"], large["report"]["lipids_after"])
        full = mf.slice_membrane_cg(self.membrane, self.placed, self.box, self.midplane, mf.Settings(box_xy_buffer_nm=5.0))
        self.assertFalse(full["report"]["cropped"])
        self.assertEqual(full["report"]["lipids_after"], len(self.membrane))

    def test_18_nm_gpcr_frames_are_cut_to_a_fraction_of_their_area(self):
        for frame in (KOR_FRAME, GPR139_FRAME):
            raw, box = mf.read_cg(frame)
            raw, alignment = mf.align_frame_to_box(raw, box)                 # the bundled frames are periodic in their box as shipped
            self.assertEqual(alignment["rotation_about_z_deg"], 0.0)
            self.assertEqual(alignment["overlapping_pairs_as_read"], 0)
            cg = mf.classify_cg(raw, MAPPING)
            slab = mf.make_membrane_whole(cg["membrane"], box)
            midplane = mf.bilayer_midplane(cg["membrane"], slab)[0]
            cell = np.array(box[:2])
            xy = np.array([[b["x"], b["y"]] for res in cg["protein"] for b in res])
            xy = xy - cell * np.round((xy - xy[0]) / cell)                      # unwrap the protein across the boundary
            stand_in = [{"atom": "BB", "resname": "X", "chain": "A", "resid": i, "x": 10 * x, "y": 10 * y, "z": 100.0} for i, (x, y) in enumerate(xy)]
            cut = mf.slice_membrane_cg(cg["membrane"], stand_in, box, midplane, self.settings)
            report = cut["report"]
            with self.subTest(frame=frame.name):
                self.assertTrue(report["cropped"])
                self.assertLess(report["area_fraction_kept"], 0.2)
                self.assertGreater(report["lipids_after"], 40)
                self.assertLess(report["lipids_after"], report["lipids_before"] // 5)
                self.assertGreater(min(report["membrane_margin_around_complex_nm"].values()), 0.99)
                self.assertTrue(report["leaflets_after"]["lower"] and report["leaflets_after"]["upper"])
                self.assertEqual(cut["box"][2], box[2])
                self.check_whole_for(cg["membrane"], cut)
                # no hard-core seam overlap survives in the new periodic cell, and few contacts at all are left to backmapping
                from scipy.spatial import cKDTree
                pts = np.vstack([m["xyz"] for m in cut["membrane"]])
                own = np.concatenate([np.full(len(m["xyz"]), i) for i, m in enumerate(cut["membrane"])])
                new_box = np.array(cut["box"])
                tree = cKDTree(np.mod(pts - pts.min(axis=0), new_box), boxsize=new_box)
                for limit, allowed in ((self.settings.seam_hard_core_nm, 0), (self.settings.seam_min_bead_nm, 5)):
                    pairs = tree.query_pairs(limit, output_type="ndarray")
                    flat = cKDTree(pts).query_pairs(limit, output_type="ndarray")
                    created = {tuple(p) for p in pairs.tolist() if own[p[0]] != own[p[1]]} - {tuple(p) for p in flat.tolist()}
                    self.assertLessEqual(len(created), allowed)

    def check_whole_for(self, membrane, cut):
        beads = {m["cg"]: len(m["beads"]) for m in membrane}
        self.assertTrue(all(len(m["beads"]) == beads[m["cg"]] for m in cut["membrane"]))


class OffsetSlack(unittest.TestCase):
    def test_slack_is_the_search_range_clamped_to_the_minimum_margin(self):
        self.assertAlmostEqual(mf.offset_slack(mf.Settings()), 0.5)                              # 1.0 buffer, 0.5 margin
        self.assertAlmostEqual(mf.offset_slack(mf.Settings(box_xy_buffer_nm=0.7)), 0.2)
        self.assertEqual(mf.offset_slack(mf.Settings(box_xy_buffer_nm=0.5)), 0.0)                # cannot shift at all
        self.assertEqual(mf.offset_slack(mf.Settings(box_xy_buffer_nm=0.3)), 0.0)
        self.assertEqual(mf.offset_slack(mf.Settings(slice_optimize_offset=False)), 0.0)


class Integration(unittest.TestCase):
    """prepare_inputs end to end up to (not including) backmapping: orientation, registration, slicing, box."""

    @classmethod
    def setUpClass(cls):
        if subprocess.run([sys.executable, "-c", "import mstool"], capture_output=True).returncode:
            raise unittest.SkipTest("mstool is not importable in this interpreter")

    def prepare(self, cg=CG_GRO, aa=AA_PDB, **kwargs):
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, True)
        out = tmp / "out"
        (out / "work").mkdir(parents=True)
        (out / mf.LOG_NAME).write_text("")
        request = kwargs.pop("orientation", mf.OrientationRequest(mode="none"))
        session = mf.Session(out=out, name="t", gmx="gmx", forcefield=FORCEFIELD, data=mf.locate_data(None), python=sys.executable,
                             ntomp=1, nsteps=1, orientation=request, work=out / "work", **kwargs)
        return session, mf.prepare_inputs(session, aa, cg)

    def test_auto_box_slices_before_backmapping_and_records_it(self):
        session, prepared = self.prepare()
        report = session.record["slice"]
        self.assertTrue(report["cropped"])
        self.assertEqual(prepared["box"], report["new_cell_nm"])
        self.assertLess(prepared["box"][0], 11.18)
        self.assertEqual(sum(prepared["composition"].values()), report["lipids_after"])
        self.assertEqual(len(prepared["membrane"]), report["lipids_after"])
        moved = mf.xyz_nm(prepared["placed"]) / 10.0
        # the window may be shifted by the density/composition search, never below the minimum margin; both margins add up to 2 x buffer
        margins = (float(moved[:, 0].min()), float(prepared["box"][0] - moved[:, 0].max()))
        self.assertGreaterEqual(min(margins), session.settings.box_xy_min_buffer_nm - 1e-6)
        self.assertAlmostEqual(sum(margins), 2.0 * session.settings.box_xy_buffer_nm, places=5)
        self.assertEqual(report["crop_offset"]["offset_nm"][0], round(1.0 - margins[0], 3))
        self.assertTrue(session.record["slice_check"]["pass"])
        self.assertTrue((Path(session.out) / "membrane_validation.md").is_file())

    def test_oriented_complex_slices_and_keeps_its_orientation(self):
        request = mf.OrientationRequest(mode="opm", chains=("R",), pdb_id="6WHC", opm_file=OPM_6WHC)
        session, prepared = self.prepare(orientation=request)
        registration = session.record["orientation"]["registration"]
        self.assertEqual(registration["check"]["status"], "PASS")
        self.assertLess(registration["check"]["normal_tilt_deg"], 1e-4)
        self.assertTrue(session.record["slice"]["cropped"])
        # the complex is still the oriented one up to a rotation about z and a translation: its z coordinates only shift
        oriented, _ = mf.read_pdb(Path(session.out) / "oriented.pdb")
        z_before, z_after = mf.xyz_nm(oriented)[:, 2], mf.xyz_nm(prepared["placed"])[:, 2]
        self.assertLess(float(np.ptp(z_after - z_before)), 1e-2)                       # PDB rounding only: a pure z translation
        self.assertTrue((Path(session.out) / "orientation_report.json").is_file())

    def test_oriented_complex_embeds_into_a_bundled_membrane_and_is_sliced(self):
        request = mf.OrientationRequest(mode="opm", chains=("R",), pdb_id="6WHC", opm_file=OPM_6WHC)
        session, prepared = self.prepare(cg=KOR_FRAME, orientation=request, embed=True)
        report = session.record["slice"]
        self.assertTrue(report["cropped"])
        self.assertLess(prepared["box"][0] * prepared["box"][1], 0.5 * 18.28 * 18.28)       # a fraction of the 18 nm patch
        self.assertLess(report["lipids_after"], 0.5 * report["lipids_before"])
        registration = session.record["orientation"]["registration"]
        self.assertEqual(registration["mode"], "embed")
        self.assertEqual(registration["check"]["status"], "PASS")
        self.assertLess(registration["check"]["normal_tilt_deg"], 1e-4)
        self.assertEqual(session.record["embedding"]["bilayer_centre_input_A"], 0.0)          # the oriented frame: bilayer centre at z = 0
        self.assertAlmostEqual(session.record["embedding"]["midplane_nm"], session.record["coarse_grain"]["bilayer_midplane_nm"], places=3)
        oriented, _ = mf.read_pdb(Path(session.out) / "oriented.pdb")
        z_shift = mf.xyz_nm(prepared["placed"])[:, 2] - mf.xyz_nm(oriented)[:, 2]
        self.assertLess(float(np.ptp(z_shift)), 1e-2)                                         # embedding only translated the complex
        self.assertAlmostEqual(float(z_shift.mean()) / 10.0, session.record["coarse_grain"]["bilayer_midplane_nm"], places=2)
        self.assertEqual(sum(prepared["composition"].values()), report["lipids_after"])

    def test_single_pass_protein_embeds_away_from_the_frames_receptor_hole(self):
        # TM6 of the example receptor alone stands in for a single-pass protein; the GPR139 frame's receptor hole is
        # far larger than one helix, so it goes on the unbroken bilayer away from that hole ("free" site)
        aa = mf.read_all_atom(AA_PDB, FORCEFIELD)[0]
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, True)
        helix = tmp / "tm6.pdb"
        mf.write_pdb([a for a in aa if a["chain"] == "R" and 349 <= a["resid"] <= 376], helix)
        # one helix plus the default 1 nm buffer is a ~4 nm cell, where the cut edges alone move the APL past the gate
        settings = mf.Settings(box_xy_buffer_nm=3.0)
        session, prepared = self.prepare(cg=GPR139_FRAME, aa=helix, embed=True, embed_site="free", settings=settings)
        embedding, report = session.record["embedding"], session.record["slice"]
        self.assertEqual(embedding["site"], "free")
        self.assertGreater(embedding["distance_to_frame_protein_nm"], 6.0)
        self.assertLess(sum(embedding["removed_lipids"].values()), 15)
        self.assertLessEqual(embedding["voids"]["largest_pocket_nm3"], mf.MAX_VOID_NM3)
        self.assertTrue(report["cropped"])
        self.assertEqual(report["frame_protein_beads_in_window"], 0)
        self.assertTrue(session.record["slice_check"]["pass"])
        self.assertEqual(sum(prepared["composition"].values()), report["lipids_after"])
        self.assertGreater(min(prepared["box"][:2]), 7.0)
        # the segment named with --orientation none sets the bilayer centre the embedding uses
        z = [a["z"] for a in aa if a["chain"] == "R" and 349 <= a["resid"] <= 376 and a["atom"] == "CA"]
        session, _ = self.prepare(cg=GPR139_FRAME, aa=helix, embed=True, embed_site="free", settings=settings,
                                  orientation=mf.OrientationRequest(mode="none", residues=((349, 376),)))
        self.assertAlmostEqual(session.record["embedding"]["bilayer_centre_input_A"], round(0.5 * (min(z) + max(z)), 3), places=3)
        self.assertTrue(session.record["embedding"]["hydrophobic_belt"]["given"])
        with self.assertRaises(mf.StageFailure) as caught:  # the same helix in the receptor's hole leaves a pocket
            self.prepare(cg=GPR139_FRAME, aa=helix, embed=True, settings=settings)
        self.assertEqual(caught.exception.stage, "mapping")
        self.assertIn("--embed-site free", caught.exception.what)

    def test_user_box_is_validated_against_the_cell_and_cuts_the_membrane(self):
        session, prepared = self.prepare(box_a=(96.0, 111.8, 250.0))
        self.assertEqual(session.record["slice"]["mode"], "user")
        self.assertAlmostEqual(prepared["box"][0], 9.6)
        with self.assertRaises(mf.StageFailure) as caught:
            self.prepare(box_a=(300.0, 111.8, 250.0))
        self.assertEqual(caught.exception.stage, "slice")
        self.assertIn("larger than the 11.180 nm", caught.exception.what)

    def test_xy_buffer_option_and_help(self):
        result = subprocess.run([sys.executable, "-m", "membraneforger", "--help"], cwd=REPO, text=True, capture_output=True)
        self.assertIn("--xy-buffer", result.stdout)
        result = subprocess.run([sys.executable, "-m", "membraneforger", "--all-atom", str(AA_PDB), "--coarse-grain", str(CG_GRO),
                                 "--xy-buffer", "0"], cwd=REPO, text=True, capture_output=True)
        self.assertIn("--xy-buffer must be positive", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
