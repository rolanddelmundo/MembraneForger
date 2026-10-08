"""A frame rotated about z relative to its box is detected and rotated back before anything periodic is done with it."""
import math
import unittest

import numpy as np

from .common import CG_GRO, REPO, failure, mf

KOR1 = REPO / "examples" / "preeq_cg_cellmem" / "KOR1_cg_cellmem.gro"


def liquid_frame(box=18.0, spacing=0.42, rotate_deg=0.0, seed=0):
    """A periodic 2D liquid of two-bead 'lipids' (random positions at least `spacing` apart under PBC), optionally rotated about the centre.

    A lattice would not do: a square lattice folds onto itself without overlaps at every Pythagorean angle (36.87 degrees, ...),
    which a liquid membrane never does.
    """
    rng = np.random.default_rng(seed)
    points = []
    for _ in range(40000):
        p = rng.uniform(0.0, box, 2)
        if not points or np.all(np.linalg.norm(np.mod(np.array(points) - p + box / 2, box) - box / 2, axis=1) >= spacing):
            points.append(p)
        if len(points) >= int(0.6 * (box / spacing) ** 2):
            break
    t = math.radians(rotate_deg)
    atoms = []
    for k, (x, y) in enumerate(points):
        for dx, dz, name in ((0.0, 1.0, "PO4"), (0.15, 0.0, "C1A")):
            px, py = x + dx - box / 2, y - box / 2
            rx, ry = px * math.cos(t) - py * math.sin(t) + box / 2, px * math.sin(t) + py * math.cos(t) + box / 2
            atoms.append({"resid": k + 1, "resname": "POPC", "atom": name, "chain": "", "x": rx, "y": ry, "z": 5.0 + dz})
    return atoms, [box, box, 10.0]


class FrameAlignment(unittest.TestCase):
    def test_a_periodic_frame_is_left_alone(self):
        atoms, box = liquid_frame()
        fixed, report = mf.align_frame_to_box(atoms, box)
        self.assertEqual(report["rotation_about_z_deg"], 0.0)
        self.assertEqual(report["overlapping_pairs_as_read"], 0)
        self.assertIs(fixed, atoms)

    def test_a_rotated_frame_is_rotated_back(self):
        atoms, box = liquid_frame(rotate_deg=30.0)
        self.assertGreater(mf.frame_overlaps(np.array([[a["x"], a["y"], a["z"]] for a in atoms]), np.arange(len(atoms)) // 2, np.array(box)), 100)
        fixed, report = mf.align_frame_to_box(atoms, box)
        angle = abs(report["rotation_about_z_deg"]) % 90.0
        self.assertLess(min(abs(angle - 30.0), abs(angle - 60.0)), 0.1)   # 30 degrees back, or 60 forward: the same square
        self.assertEqual(report["overlapping_pairs_after"], 0)
        xy = np.mod(np.array([[a["x"], a["y"]] for a in fixed if a["atom"] == "PO4"]), box[:2])
        nearest = np.sort(np.linalg.norm(xy[:, None, :] - xy[None, :, :], axis=2), axis=1)[:, 1]
        self.assertGreater(float(nearest.min()), 0.3)                    # no folded corner: the liquid is intact everywhere, edges included
        self.assertAlmostEqual(fixed[0]["z"], atoms[0]["z"])              # z is untouched

    def test_a_frame_that_no_rotation_makes_periodic_is_refused(self):
        atoms, box = liquid_frame()
        squeezed = [dict(a, x=a["x"] * 1.3) for a in atoms]              # content wider than the box
        self.assertIn("not periodic in its", failure(mf.align_frame_to_box, squeezed, box))

    def test_the_6whc_frame_is_periodic_and_a_bundled_frame_is_rotated(self):
        atoms, box = mf.read_cg(CG_GRO)
        _, report = mf.align_frame_to_box(atoms, box)
        self.assertEqual(report["rotation_about_z_deg"], 0.0)
        atoms, box = mf.read_cg(KOR1)
        fixed, report = mf.align_frame_to_box(atoms, box)
        self.assertGreater(abs(report["rotation_about_z_deg"]), 1.0)
        self.assertLess(report["overlapping_pairs_after"], 0.02 * report["overlapping_pairs_as_read"])
        self.assertEqual(len(fixed), len(atoms))


if __name__ == "__main__":
    unittest.main()
