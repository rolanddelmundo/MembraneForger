"""The coarse-grained frame gate (MembraneValidation.validate_cg_frame): the bundled frames pass, damaged copies of one fail.

None of these tests backmaps anything: they read and classify the real frames and stop before embedding.
"""
import copy
import tempfile
import unittest
from pathlib import Path

import numpy as np

from membraneforger.membrane_report import MembraneValidation, protein_of_cg_frame

from .common import MAPPING, failure, mf


def check(membrane, protein, box, midplane, detail, cg_validate=True):
    """Run the gate on a copy of a membrane; returns (gate result or None, the SystemExit message or None)."""
    validation = MembraneValidation(Path(tempfile.mkdtemp()), mf.Settings(cg_validate=cg_validate))
    validation.add_cg("cg_frame", membrane, [], box, midplane, protein_of_cg_frame(protein))
    result = []
    message = failure(lambda: result.append(validation.validate_cg_frame(membrane, protein, box, midplane, detail)))
    return (result[0] if result else None), message


def statuses(result):
    """metric name -> list of statuses of the gate's records."""
    out = {}
    for r in result["records"]:
        out.setdefault(r["metric"], []).append(r["status"])
    return out


def move(mol, bead_index, delta):
    """Move one bead of a whole molecule by delta (nm), in both the bead records and the whole-molecule coordinates."""
    mol["xyz"][bead_index] += delta
    b = mol["beads"][bead_index]
    b["x"], b["y"], b["z"] = (float(v) for v in mol["xyz"][bead_index])


class BundledFrames(unittest.TestCase):
    def test_every_shipped_frame_passes_the_gate(self):
        for frame in mf.bundled_frames().values():
            atoms, box = mf.read_cg(frame)
            cg = mf.classify_cg(atoms, MAPPING)
            slab = mf.make_membrane_whole(cg["membrane"], box)
            midplane, _, detail = mf.bilayer_midplane(cg["membrane"], slab)
            result, message = check(cg["membrane"], cg["protein"], box, midplane, detail)
            with self.subTest(frame=frame.name):
                self.assertIsNone(message)
                self.assertEqual(result["status"], "PASS", result["overall"]["warnings"])


class DamagedFrame(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        atoms, cls.box = mf.read_cg(mf.bundled_frames()["GPR5"])
        cls.cg = mf.classify_cg(atoms, MAPPING)
        slab = mf.make_membrane_whole(cls.cg["membrane"], cls.box)
        cls.midplane, _, cls.detail = mf.bilayer_midplane(cls.cg["membrane"], slab)

    def run_gate(self, membrane, protein=None, detail=None, **kw):
        return check(membrane, protein or self.cg["protein"], self.box, self.midplane, self.detail if detail is None else detail, **kw)

    def test_a_tilted_bilayer_fails(self):
        membrane, protein = copy.deepcopy(self.cg["membrane"]), copy.deepcopy(self.cg["protein"])
        angle = np.radians(15.0)
        rotate = np.array([[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]])
        centre = np.array([0.0, 0.5 * self.box[1], self.midplane])
        for mol in membrane:
            mol["xyz"] = (mol["xyz"] - centre) @ rotate.T + centre
            for b, (x, y, z) in zip(mol["beads"], mol["xyz"]):
                b["x"], b["y"], b["z"] = float(x), float(y), float(z)
        tilted = mf.bilayer_midplane(membrane, (self.midplane - 3.0, self.midplane + 3.0))[2]
        _, message = self.run_gate(membrane, protein, tilted)
        self.assertIn("bilayer normal tilt", message)
        _, message = self.run_gate(membrane, protein, tilted, cg_validate=False)
        self.assertIsNone(message)

    def test_a_frame_with_one_leaflet_fails(self):
        membrane = copy.deepcopy(self.cg["membrane"])
        lower = [m for m in membrane if float(m["xyz"][:, 2].mean()) < self.midplane]
        _, message = self.run_gate(lower)
        self.assertIn("not a usable bilayer", message)

    def test_a_broken_molecule_fails(self):
        membrane = copy.deepcopy(self.cg["membrane"])
        move(membrane[0], 0, np.array([3.0, 0.0, 0.0]))
        _, message = self.run_gate(membrane)
        self.assertIn("longest bead link inside a lipid", message)

    def test_lipids_with_their_head_on_the_wrong_side_are_flagged(self):
        membrane = copy.deepcopy(self.cg["membrane"])
        popc = [m for m in membrane if m["cg"] == "POPC"][:10]
        for mol in popc:  # lying flat across the midplane: the head on one side, most of the beads on the other
            side = 1.0 if float(mol["xyz"][:, 2].max()) - self.midplane > self.midplane - float(mol["xyz"][:, 2].min()) else -1.0
            for i, b in enumerate(mol["beads"]):
                head = b["atom"] in ("NC3", "PO4", "GL1", "GL2")
                move(mol, i, np.array([0.0, 0.0, self.midplane + side * (0.2 if head else -0.2) - mol["xyz"][i, 2]]))
        result, message = self.run_gate(membrane)
        self.assertIsNone(message)  # reported, never a stop
        self.assertEqual(statuses(result)["lipids on the wrong leaflet"], ["WARNING"])

    def test_lipid_beads_on_top_of_the_protein_fail(self):
        membrane = copy.deepcopy(self.cg["membrane"])
        beads = [b for res in self.cg["protein"] for b in res]
        rng = np.random.default_rng(0)
        for mol in membrane[:150]:
            target = beads[int(rng.integers(len(beads)))]
            move(mol, 0, np.array([target["x"], target["y"], target["z"]]) + 0.05 - mol["xyz"][0])
        _, message = self.run_gate(membrane)
        self.assertIn("protein-lipid bead overlaps", message)

    def test_the_gate_is_in_the_record_and_the_stage_matrix_stays_whole(self):
        membrane = copy.deepcopy(self.cg["membrane"])
        validation = MembraneValidation(Path(tempfile.mkdtemp()), mf.Settings())
        validation.add_cg("cg_frame", membrane, [], self.box, self.midplane, protein_of_cg_frame(self.cg["protein"]))
        validation.validate_cg_frame(membrane, self.cg["protein"], self.box, self.midplane, self.detail)
        record = validation.record()
        self.assertEqual(record["gates"]["0 CG frame quality"]["status"], "PASS")
        self.assertTrue(any(r["stage"] == "cg_frame" and r["status"] == "PASS" for r in record["records"]))


if __name__ == "__main__":
    unittest.main()
