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


def correct_gm3() -> tuple[list, dict]:
    """GM3 atom names in topology order and a correct whole molecule: the reference sugars plus a straight trans ceramide."""
    named = reference_xyz()
    names = [name for name, _ in mf.GM3_XML_TO_GLPA]
    ceramide = {n: np.array([0.5 * i, 0.0, -5.0 - 0.3 * (i % 2)]) for i, n in enumerate(n for n in names if n not in named)}
    full = {**named, **ceramide}
    # a genuinely trans ceramide: C3S-C4S-C5S-C6S and C2S-NF-C1F-C2F at 180 degrees
    for a, b, c, d in mf.gm3_trans_definitions(DATA):
        full[a], full[b], full[c], full[d] = (np.array(p) for p in ((1.0, 0.0, -9.0), (0.0, 0.0, -9.0), (0.0, 1.0, -9.0), (-1.0, 1.0, -9.0)))
    return names, full


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


class MapDatBlocks(unittest.TestCase):
    def test_a_residue_block_is_found_with_trailing_blanks_and_windows_line_ends(self):
        text = (DATA / "map.dat").read_text().replace("RESI GM3\n", "RESI GM3  \n").replace("RESI SAPI25\n", "RESI SAPI25 \n")
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "map.dat").write_text(text.replace("\n", "\r\n"), newline="")
            self.assertEqual(sorted(d[1] for d in mf.gm3_chiral_definitions(Path(tmp))), sorted(mf.GM3_POSITIONS))
            self.assertEqual(sorted(d[1] for d in mf.pi_chiral_definitions(Path(tmp))), sorted(mf.PI_POSITIONS))
            self.assertEqual([d for d in mf.gm3_trans_definitions(Path(tmp))], [tuple(d) for d in mf.gm3_trans_definitions(DATA)])

    def test_a_missing_block_is_named(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "map.dat").write_text("RESI PSM\n\n[ chiral ]\n")
            with self.assertRaisesRegex(SystemExit, "no RESI SAPI25 block"):
                mf.pi_chiral_definitions(Path(tmp))


class AmideRestraint(unittest.TestCase):
    """The CHARMM36 term C2S-NF-C1F-C2F (CTL1 NHL C CTL2) has k = 0, so emres holds the amide trans only through its DIHRES row."""

    def test_the_amide_row_is_the_map_dat_trans_definition_in_both_gm3_topologies(self):
        order = [gm3 for gm3, _ in mf.GM3_XML_TO_GLPA]
        amide = next(d for d in mf.gm3_trans_definitions(DATA) if mf.GM3_TRANS_NAMES.get(d) == "ceramide amide")
        quad = tuple(order.index(name) for name in amide)  # zero-based: the topology lists the atoms in GM3_XML_TO_GLPA order
        for itp, molecule in (("GLPA.itp", "GLPA"), ("GM3.itp", "GM3")):
            rows = mf.itp_dihedral_restraints([REPO / "membraneforger" / "forcefield" / "toppar" / itp])[molecule]
            self.assertEqual([t for q, t in rows if q == quad], [180.0], itp)

    def test_the_psm_amide_row_is_the_map_dat_trans_definition(self):
        itp = REPO / "membraneforger" / "forcefield" / "toppar" / "PSM.itp"
        names = [name for name, _ in mf.itp_atoms([itp])["PSM"]]
        amide = next(tuple(d) for d in mf.map_section(DATA, "trans", "PSM") if d[1] == "NF")
        quad = tuple(names.index(name) for name in amide)
        self.assertEqual(amide, ("C2S", "NF", "C1F", "C2F"))
        rows = mf.itp_dihedral_restraints([itp])["PSM"]
        self.assertEqual([t for q, t in rows if q == quad], [180.0])
        self.assertEqual(len(rows), 4)  # the two centres, the C4=C5 bond and the amide

    def test_three_bands(self):
        self.assertEqual([mf.trans_band(p) for p in (180.0, -150.0, 149.9, -90.0, 30.1, 30.0, 0.0, -12.0)],
                         ["trans", "trans", "twisted", "twisted", "twisted", "cis", "cis", "cis"])

    def test_twisted_between_90_and_150_is_reported_but_not_refused_and_under_90_is_refused(self):
        _, full = correct_gm3()
        amide = ("C2S", "NF", "C1F", "C2F")
        a, b, c = (full[n] for n in amide[:3])
        axis = (c - b) / np.linalg.norm(c - b)
        start = full["C2F"] - c
        for angle, band, refused in ((120.0, "twisted", False), (60.0, "twisted", True), (0.0, "cis", True), (175.0, "trans", False)):
            theta = np.radians(mf.dihedral(a, b, c, full["C2F"]) - angle)  # rotate C2F about the NF-C1F axis to the wanted dihedral
            named = {**full, "C2F": c + start * np.cos(theta) + np.cross(axis, start) * np.sin(theta) + axis * (axis @ start) * (1 - np.cos(theta))}
            phi = mf.dihedral(a, b, c, named["C2F"])
            self.assertAlmostEqual(abs(phi), angle, delta=0.5)
            self.assertEqual(mf.trans_band(phi), band)
            self.assertEqual(bool(mf.wrong_gm3_centres(named, [], [amide])[1]), refused, angle)


def correct_inositol() -> dict:
    """A chair of six carbons (C11..C16) with a hydrogen and an oxygen on each, every one in the configuration of map.dat.

    Axial and equatorial positions are assigned to H and O alternately, then H and O are exchanged at each carbon the
    SAPI25 definitions read as inverted: that exchange flips exactly that carbon, which no ring-chair restraint sees."""
    ring = ["C11", "C12", "C13", "C14", "C15", "C16"]
    named = {}
    for k, c in enumerate(ring):
        angle = np.radians(60.0 * k)
        up = 1.0 if k % 2 == 0 else -1.0
        radial = np.array([np.cos(angle), np.sin(angle), 0.0])
        named[c] = np.array([1.45 * radial[0], 1.45 * radial[1], 0.25 * up])
        axial, equatorial = named[c] + np.array([0.0, 0.0, 1.0 * up]), named[c] + 1.0 * radial + np.array([0.0, 0.0, -0.3 * up])
        h, o = (f"H{k + 1}", "O12" if k == 0 else f"O{k + 1}")
        named[h], named[o] = axial, equatorial
    inverted, _ = mf.wrong_gm3_centres(named, mf.pi_chiral_definitions(DATA), [])
    for centre in inverted:
        h, o = f"H{ring.index(centre) + 1}", "O12" if centre == "C11" else f"O{ring.index(centre) + 1}"
        named[h], named[o] = named[o], named[h]
    return named


def write_inositol(path: Path, named: dict, resname: str = "SAPI", box: float | None = None) -> None:
    """The inositol atoms as one SAPI25 residue in a PDB (or a .gro when the path ends in .gro), wrapped into a cell if asked."""
    xyz = {n: (np.mod(x, box) if box else x) for n, x in named.items()}
    if path.suffix == ".gro":
        lines = [f"{1:5d}{resname:<5s}{n:>5s}{i + 1:5d}{x[0] / 10:8.3f}{x[1] / 10:8.3f}{x[2] / 10:8.3f}\n" for i, (n, x) in enumerate(xyz.items())]
        path.write_text(f"inositol\n{len(lines)}\n" + "".join(lines) + f"{(box or 100.0) / 10:10.5f}" * 3 + "\n")
        return
    rows = [f"HETATM{i + 1:5d} {n:<4s} {resname:<4s}M   1    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00\n" for i, (n, x) in enumerate(xyz.items())]
    head = f"CRYST1{box:9.3f}{box:9.3f}{box:9.3f}  90.00  90.00  90.00 P 1           1\n" if box else ""
    path.write_text(head + "".join(rows) + "END\n")


class Inositol(unittest.TestCase):
    """SAPI25: the six ring-chair rows are satisfied when a hydrogen and a hydroxyl swap places, so the centres are tested geometrically."""

    def test_map_dat_defines_the_six_inositol_centres(self):
        self.assertEqual(sorted(d[1] for d in mf.pi_chiral_definitions(DATA)), sorted(mf.PI_POSITIONS))

    def test_a_correct_inositol_passes_and_each_swapped_carbon_is_the_one_reported(self):
        good = correct_inositol()
        definitions = mf.pi_chiral_definitions(DATA)
        self.assertEqual(mf.wrong_gm3_centres(good, definitions, [])[0], [])
        with tempfile.TemporaryDirectory() as tmp:
            for suffix, resname in ((".pdb", "SAPI"), (".gro", "SAPI2")):
                path = Path(tmp) / f"good{suffix}"
                write_inositol(path, good, resname)
                report, fails = mf.check_pi_stereo(path, DATA)
                self.assertEqual((report["molecules"], fails), (1, []))
                for k, centre in enumerate(("C11", "C12", "C13", "C14", "C15", "C16")):
                    bad = dict(good)
                    h, o = f"H{k + 1}", "O12" if k == 0 else f"O{k + 1}"
                    bad[h], bad[o] = good[o], good[h]
                    broken = Path(tmp) / f"bad{k}{suffix}"
                    write_inositol(broken, bad, resname)
                    report, fails = mf.check_pi_stereo(broken, DATA)
                    self.assertEqual([c for c, n in report["wrong_by_centre"].items() if n], [mf.PI_POSITIONS[centre]])
                    self.assertEqual(len(fails), 1)
                    self.assertIn("inverted inositol stereocentre", fails[0])

    def test_a_structure_without_sapi25_has_nothing_to_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            empty = Path(tmp) / "empty.pdb"
            empty.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n")
            report, fails = mf.check_pi_stereo(empty, DATA)
            self.assertEqual((report["molecules"], fails), (0, []))

    def test_a_correct_inositol_split_by_the_periodic_boundary_is_made_whole(self):
        good = {n: x + 3.0 for n, x in correct_inositol().items()}
        box = 6.0  # coordinates run from 0.7 to 5.2 A and are shifted by -3: negative ones jump to the far side of the cell
        shifted = {n: x - 3.0 for n, x in good.items()}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wrapped.pdb"
            write_inositol(path, shifted, box=box)
            self.assertGreater(float(np.ptp([np.mod(x, box) for x in shifted.values()], axis=0).max()), box / 2)  # genuinely split
            report, fails = mf.check_pi_stereo(path, DATA)
            self.assertEqual((report["molecules"], fails), (1, []))


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
            names, full = correct_gm3()
            for label, sign in (("right", 1.0), ("mirror", -1.0)):
                path = Path(tmp) / f"{label}.pdb"
                rows = [f"HETATM{i + 1:5d} {n:<4s} GM3 M   1    {sign * full[n][0]:8.3f}{full[n][1]:8.3f}{full[n][2]:8.3f}  1.00  0.00\n"
                        for i, n in enumerate(names)]
                path.write_text("".join(rows) + "END\n")
                report, fails = mf.check_gm3_stereo(path, DATA)
                self.assertEqual(report["molecules"], 1)
                if label == "right":
                    self.assertEqual(fails, [])
                    self.assertEqual(report["bands_by_bond"], {"ceramide C4=C5": {"trans": 1, "twisted": 0, "cis": 0},
                                                               "ceramide amide": {"trans": 1, "twisted": 0, "cis": 0}})
                else:
                    self.assertEqual(len(fails), 1)
                    self.assertEqual(sum(report["wrong_by_centre"].values()), 16)
                    self.assertIn("inverted sugar stereocentre", fails[0])

    def test_a_correct_gm3_wrapped_across_the_periodic_boundary_is_made_whole(self):
        # mdrun may wrap each atom into the cell independently (GPU and CPU builds wrap differently), splitting a
        # molecule across the boundary; on raw coordinates every such molecule reads as inverted - a false positive.
        names, full = correct_gm3()
        box = 200.0
        wrapped = {n: np.mod(x, box) for n, x in full.items()}  # negative coordinates jump to the far side of the cell
        self.assertGreater(float(np.ptp([wrapped[n] for n in names], axis=0).max()), box / 2)  # genuinely split
        with tempfile.TemporaryDirectory() as tmp:
            pdb = Path(tmp) / "wrapped.pdb"
            rows = [f"HETATM{i + 1:5d} {n:<4s} GM3 M   1    {wrapped[n][0]:8.3f}{wrapped[n][1]:8.3f}{wrapped[n][2]:8.3f}  1.00  0.00\n"
                    for i, n in enumerate(names)]
            pdb.write_text(f"CRYST1{box:9.3f}{box:9.3f}{box:9.3f}  90.00  90.00  90.00 P 1           1\n" + "".join(rows) + "END\n")
            report, fails = mf.check_gm3_stereo(pdb, DATA)
            self.assertEqual((report["molecules"], fails), (1, []))
            gro = Path(tmp) / "wrapped.gro"
            lines = [f"{1:5d}{'GM3':<5s}{n:>5s}{(i + 1) % 100000:5d}{wrapped[n][0] / 10:8.3f}{wrapped[n][1] / 10:8.3f}{wrapped[n][2] / 10:8.3f}\n"
                     for i, n in enumerate(names)]
            gro.write_text(f"wrapped GM3\n{len(names)}\n" + "".join(lines) + f"{box / 10:10.5f}{box / 10:10.5f}{box / 10:10.5f}\n")
            report, fails = mf.check_gm3_stereo(gro, DATA)
            self.assertEqual((report["molecules"], fails), (1, []))
            # without a box (no CRYST1) the same split coordinates must still be reported, not silently unwrapped
            raw = Path(tmp) / "raw.pdb"
            raw.write_text("".join(rows) + "END\n")
            report, _ = mf.check_gm3_stereo(raw, DATA)
            self.assertGreater(sum(report["wrong_by_centre"].values()) + sum(report["cis_by_bond"].values()), 0)

    def test_the_example_script_uses_the_same_check(self):
        text = (REPO / "membraneforger" / "examples" / "check_gm3_stereo.py").read_text()
        self.assertIn("from membraneforger.stereo import check_gm3_stereo", text)


if __name__ == "__main__":
    unittest.main()
