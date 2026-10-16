"""GM3 stereochemistry: the 16 sugar centres and 2 ceramide trans bonds are judged geometrically, so an inverted sugar carbon fails."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from .common import DATA, REPO, mf

REFERENCE = json.loads((Path(__file__).parent / "data" / "gm3_sugar_reference.json").read_text())
SUGAR_CENTRES = {"glucose (beta-D-Glc)": ("C1", "C2", "C3", "C4", "C5"), "galactose (beta-D-Gal)": ("C7", "C8", "C10", "C11", "C12"),
                 "sialic acid (alpha-Neu5Ac)": ("C15", "C17", "C20", "C18", "C19", "C21")}


def reference_xyz() -> dict:
    """GM3 atom name -> xyz (A) of the three reference sugar models, joined at their shared linking oxygens.

    Each sugar is modelled in its own frame and carries the linking oxygen of the next bond under the GM3 name it shares
    with its neighbour (O4 between Glc and Gal, O10 between Gal and Neu5Ac). Translating Gal so that its O4 lands on
    Glc's, and Neu5Ac so that its O10 lands on Gal's, keeps every chirality (a translation) and makes the shared names
    consistent; the sugars may overlap sterically, which the chirality test does not care about.
    """
    glc, gal, neu = (dict((n, np.array(x)) for n, x in sugar.items()) for sugar in REFERENCE["sugars"].values())
    shift = glc["O4"] - gal["O4"]
    gal = {n: x + shift for n, x in gal.items()}
    shift = gal["O10"] - neu["O10"]
    neu = {n: x + shift for n, x in neu.items()}
    return {**neu, **gal, **glc}


class Definitions(unittest.TestCase):
    def test_map_dat_defines_all_sixteen_sugar_centres_and_the_trans_bonds(self):
        chirals = mf.gm3_chiral_definitions(DATA)
        self.assertEqual(sorted(d[1] for d in chirals), sorted(mf.GM3_POSITIONS))
        self.assertEqual(len(mf.gm3_trans_definitions(DATA)), 2)
        self.assertEqual(sum(len(v) for v in SUGAR_CENTRES.values()), 16)

    def test_the_reference_sugars_satisfy_every_definition_and_their_mirror_images_none(self):
        named = reference_xyz()
        chirals = [d for d in mf.gm3_chiral_definitions(DATA) if all(n in named for n in d)]
        self.assertEqual(len(chirals), 16)                                            # every centre is testable on the models
        inverted, _ = mf.wrong_gm3_centres(named, chirals, [])
        self.assertEqual(inverted, [])
        mirrored = {k: v * np.array([-1.0, 1.0, 1.0]) for k, v in named.items()}
        inverted, _ = mf.wrong_gm3_centres(mirrored, chirals, [])
        self.assertEqual(sorted(inverted), sorted(mf.GM3_POSITIONS))                   # a mirror image inverts all 16

    def test_one_inverted_carbon_in_a_correct_chair_is_caught(self):
        # Swap the two substituents of Glc C4 (H4 and O4) through the carbon: the ring stays as it is, the centre inverts.
        named = reference_xyz()
        c4 = named["C4"]
        named["H4"], named["O4"] = c4 + (named["O4"] - c4) * np.linalg.norm(named["H4"] - c4) / np.linalg.norm(named["O4"] - c4), \
            c4 + (named["H4"] - c4) * np.linalg.norm(named["O4"] - c4) / np.linalg.norm(named["H4"] - c4)
        chirals = [d for d in mf.gm3_chiral_definitions(DATA) if all(n in named for n in d)]
        inverted, _ = mf.wrong_gm3_centres(named, chirals, [])
        self.assertEqual(inverted, ["C4"])

    def test_trans_and_cis_bonds(self):
        trans = [("A", "B", "C", "D")]
        planar = {"A": np.array([1.0, 0.0, 0.0]), "B": np.array([0.0, 0.0, 0.0]), "C": np.array([0.0, 1.0, 0.0])}
        self.assertEqual(mf.wrong_gm3_centres({**planar, "D": np.array([-1.0, 1.0, 0.0])}, [], trans)[1], [])   # 180: trans
        self.assertEqual(mf.wrong_gm3_centres({**planar, "D": np.array([1.0, 1.0, 0.0])}, [], trans)[1], trans)  # 0: cis


class Structures(unittest.TestCase):
    def test_a_structure_without_gm3_passes_and_a_mirrored_gm3_fails(self):
        xml = (DATA / "GM3.xml").read_text() if (DATA / "GM3.xml").is_file() else ""
        self.assertTrue(xml)
        with tempfile.TemporaryDirectory() as tmp:
            empty = Path(tmp) / "empty.pdb"
            empty.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n")
            report, fails = mf.check_gm3_stereo(empty, DATA)
            self.assertEqual((report["molecules"], fails), (0, []))
            # a GM3 residue written by name from the reference sugars plus a straight ceramide: correct, then mirrored
            named = reference_xyz()
            names = [name for name, _ in mf.GM3_XML_TO_GLPA]
            ceramide = {n: np.array([0.5 * i, 0.0, -5.0 - 0.3 * (i % 2)]) for i, n in enumerate(n for n in names if n not in named)}
            full = {**named, **ceramide}
            # a genuinely trans ceramide: C3S-C4S-C5S-C6S and C2S-NF-C1F-C2F at 180 degrees
            for a, b, c, d in mf.gm3_trans_definitions(DATA):
                full[a], full[b], full[c], full[d] = (np.array(p) for p in ((1.0, 0.0, -9.0), (0.0, 0.0, -9.0), (0.0, 1.0, -9.0), (-1.0, 1.0, -9.0)))
            for label, sign in (("right", 1.0), ("mirror", -1.0)):
                path = Path(tmp) / f"{label}.pdb"
                rows = [f"HETATM{i + 1:5d} {n:<4s} GM3 M   1    {sign * full[n][0]:8.3f}{full[n][1]:8.3f}{full[n][2]:8.3f}  1.00  0.00\n"
                        for i, n in enumerate(names)]
                path.write_text("".join(rows) + "END\n")
                report, fails = mf.check_gm3_stereo(path, DATA)
                self.assertEqual(report["molecules"], 1)
                if label == "right":
                    self.assertEqual(fails, [])
                else:
                    self.assertEqual(len(fails), 1)
                    self.assertEqual(sum(report["wrong_by_centre"].values()), 16)
                    self.assertIn("inverted sugar stereocentre", fails[0])

    def test_the_example_script_uses_the_same_check(self):
        text = (REPO / "examples" / "check_gm3_stereo.py").read_text()
        self.assertIn("from membraneforger.stereo import check_gm3_stereo", text)


if __name__ == "__main__":
    unittest.main()
