"""Unit tests for pure logic: classification, aliases, I/O, alignment, disulfides, topology parsing, subprocesses."""
import copy
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from .common import AA_PDB, CG_GRO, DATA, FORCEFIELD, MAPPING, POPC, failure, mf, residue, small_system


class MartiniClassification(unittest.TestCase):
    def test_classes_of_a_small_system(self):
        cg = mf.classify_cg(small_system(), MAPPING)
        self.assertEqual(dict(cg["counts"]), {"protein": 3, "phospholipid": 1, "sterol": 1, "solvent": 1, "ion": 1})

    def test_single_residue_classifier(self):
        self.assertEqual(mf.classify_cg_residue(residue(1, "W", ["W"]), MAPPING), ("solvent", None))
        self.assertEqual(mf.classify_cg_residue(residue(1, "NA", ["NA"]), MAPPING), ("ion", None))
        self.assertEqual(mf.classify_cg_residue(residue(1, "HIS", ["BB", "SC1", "SC2", "SC3"]), MAPPING)[0], "protein")
        self.assertEqual(mf.classify_cg_residue(residue(1, "POPC", POPC), MAPPING)[0], "phospholipid")
        self.assertEqual(mf.classify_cg_residue(residue(1, "QQQ", ["X1"]), MAPPING)[0], "unsupported")

    def test_cholesterol_maps_to_chl1_and_reports_r6(self):
        cg = mf.classify_cg(small_system(), MAPPING)
        sterol = next(m for m in cg["membrane"] if m["cls"] == "sterol")
        self.assertEqual(sterol["aa"], "CHL1")
        self.assertNotIn("R6", [b["atom"] for b in sterol["beads"]])
        self.assertEqual(cg["dropped"], {"CHOL": ["R6"]})

    def test_real_system_composition_comes_from_the_file(self):
        atoms, _ = mf.read_cg(CG_GRO)
        cg = mf.classify_cg(atoms, MAPPING)
        self.assertEqual(dict(cg["composition"]), {"CHOL": 99, "POPC": 49, "DOPC": 49, "POPE": 49, "DOPE": 49, "PSM": 30,
                                                   "GM3": 20, "POPS": 15, "DOPS": 13, "SAP6": 19})
        self.assertEqual(cg["counts"]["glycolipid"], 20)
        self.assertTrue(all(len(m["beads"]) == 20 for m in cg["membrane"] if m["cg"] == "GM3"))
        self.assertEqual(cg["dropped"]["SAP6"], ["C4"])
        self.assertEqual(cg["dropped"]["GM3"], ["GLC:V", "GAL:V", "NMC:V"])
        # every bead is accounted for: classified residues plus the 3 extra residues merged into each GM3
        self.assertEqual(sum(cg["counts"].values()) + 3 * 20, len(mf.cg_residues(atoms)))

    def test_pip2_maps_to_sapi25(self):
        atoms, _ = mf.read_cg(CG_GRO)
        pip2 = next(m for m in mf.classify_cg(atoms, MAPPING)["membrane"] if m["cg"] == "SAP6")
        self.assertEqual((pip2["aa"], pip2["cls"], len(pip2["beads"])), ("SAPI25", "phospholipid", 17))


class Aliases(unittest.TestCase):
    def test_gro_alias_round_trip(self):
        names = ["CHL1", "POPC", "SAPI25", "GM3", "LONGNAME7"]
        alias = mf.gro_safe_aliases(names)
        self.assertEqual(set(alias), {"SAPI25", "LONGNAME7"})
        self.assertTrue(all(len(short) <= 5 for short in alias.values()))
        back = {short: name for name, short in alias.items()}
        self.assertEqual([back.get(alias.get(n, n), alias.get(n, n)) for n in names], names)

    def test_gro_field_truncates_without_alias(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "t.gro"
            mf.write_gro([dict(residue(1, "SAPI25", ["PO4"])[0])], [5, 5, 5], path, "t")
            self.assertEqual(mf.read_gro(path)[0][0]["resname"], "SAPI2")  # why the alias layer exists

    def test_alias_collision_is_refused(self):
        self.assertIn("collide", failure(mf.gro_safe_aliases, ["SAPI25", "M000"]))

    def test_pdb_resname_round_trip(self):
        self.assertEqual(mf.structure_resname("SAPI25"), "SAPI")
        self.assertEqual(mf.RENAME_MOLECULE[mf.structure_resname("SAPI25")], "SAPI25")
        self.assertEqual(mf.structure_resname("GM3"), "GM3")
        self.assertIn("does not survive", failure(mf.structure_resname, "SAPI24"))


class MstoolCompatibility(unittest.TestCase):
    """Regression: an all-atom input without hydrogens aborted mstool's isomer review (wave-1 integration failure)."""

    def test_rock_never_carries_a_mapped_residue_name(self):
        atoms, _ = mf.read_pdb(AA_PDB)
        rock = mf.rock_atoms(atoms)[:len(atoms)]
        self.assertEqual({a["resname"] for a in mf.rock_atoms(atoms)}, {mf.ROCK_RESNAME})
        self.assertNotIn(mf.ROCK_RESNAME, MAPPING)
        self.assertTrue({a["resname"] for a in atoms} & set(MAPPING))  # the original names would have been reviewed
        self.assertEqual([(a["atom"], a["chain"], a["resid"], a["x"]) for a in rock],
                         [(a["atom"], a["chain"], a["resid"], a["x"]) for a in atoms])
        self.assertEqual(atoms[0]["resname"] in mf.AMINO, True)  # the input records themselves are untouched

    def test_flipped_centres_are_counted_only_for_well_formed_definitions(self):
        import pandas

        from membraneforger.mstool_worker import flipped_chirals
        rows = []
        for hand in (1.0, -1.0):  # two residues, the second is the mirror image
            for name, xyz in (("H", (0, 0, hand)), ("C", (0, 0, 0)), ("A", (1, 0, -0.3)), ("B", (-0.5, 0.87, -0.3)), ("D", (-0.5, -0.87, -0.3))):
                rows.append({"name": name, "x": xyz[0], "y": xyz[1], "z": xyz[2]})
        atoms = pandas.DataFrame(rows)
        definition = [["H", "C", "A", "B", "D"]]
        self.assertEqual(sum(flipped_chirals(atoms, definition, []).values()), 1)
        self.assertEqual(flipped_chirals(atoms, definition, ["C"]), {})

    def test_ring_centres_block_every_protein_ring(self):
        atoms, _ = mf.read_pdb(AA_PDB)
        rings = {"PHE": 1, "TYR": 1, "TRP": 2, "HSD": 1, "HSE": 1, "HSP": 1, "HIS": 1, "PRO": 1}
        residues = {(a["chain"], a["resid"], a["resname"]) for a in atoms if a["resname"] in rings}
        expected = sum(rings[name] for _, _, name in residues)
        centres = mf.ring_centres([a for a in atoms if a["resname"] in rings])
        self.assertEqual(len(centres), expected)
        phe = [a for a in atoms if a["resname"] == "PHE"]
        first = [a for a in phe if (a["chain"], a["resid"]) == (phe[0]["chain"], phe[0]["resid"])]
        ring = mf.xyz_nm([a for a in first if a["atom"] in ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")]).mean(axis=0)
        self.assertLess(np.linalg.norm(np.array(mf.ring_centres(first)[0]) - ring), 1e-6)
        extra = mf.rock_atoms(atoms)[len(atoms):]
        self.assertEqual(len(extra), len(mf.ring_centres(atoms)))
        self.assertTrue(all(a["atom"].startswith("C") for a in extra))  # mstool types rock atoms by their first letter

    def test_isomer_review_is_parsed(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "log.txt"
            log.write_text("x\nIn summary, the number of residues with the flipped isomers:\nchiral    :         26\ndihedral  :          6\n####\n")
            self.assertEqual(mf.read_isomer_review(log), {"chiral": 26, "dihedral": 6})
            self.assertEqual(mf.read_isomer_review(Path(tmp) / "absent.txt"), {})
            chirality = Path(tmp) / "c.json"
            chirality.write_text('{"POPC": {"C2": 3}}')
            self.assertEqual(mf.read_isomer_review(log, chirality)["chiral_well_formed"], 3)

    def test_malformed_chirality_definitions_are_found(self):
        from membraneforger.mstool_worker import malformed_chirals
        bonds = [("C2", "C1"), ("C2", "H2"), ("C2", "O2"), ("C2", "C3"), ("C3", "C4")]
        good, bad = ["O2", "C2", "C1", "H2", "C3"], ["O2", "C2", "C1", "H2", "C4"]  # C4 is not bonded to C2
        self.assertEqual(malformed_chirals([good], bonds), [])
        self.assertEqual(malformed_chirals([good, bad], bonds), ["C2"])

    def test_gm3_chirality_definitions_are_well_formed_and_correct(self):
        """Every GM3 definition names four neighbours of its centre and, on the sugars, reads correct on a reference."""
        import json
        import xml.etree.ElementTree as ET

        from membraneforger.mstool_worker import malformed_chirals
        gm3_block = (DATA / "map.dat").read_text().split("RESI GM3")[1].split("RESI ")[0]
        block = gm3_block.split("[ chiral ]")[1].split("[")[0]
        chirals = [line.split() for line in block.splitlines() if line.strip()]
        gm3 = next(r for r in ET.parse(DATA / "GM3.xml").getroot().iter("Residue") if r.get("name") == "GM3")
        bonds = [(b.get("atomName1"), b.get("atomName2")) for b in gm3.iter("Bond")]
        self.assertEqual(len(chirals), 19)
        self.assertEqual(malformed_chirals(chirals, bonds), [])
        reference = json.loads((Path(__file__).parent / "data/gm3_sugar_reference.json").read_text())["sugars"]
        checked = 0
        for target, centre, c, d, e in chirals:
            sugar = next((m for m in reference.values() if all(n in m for n in (target, centre, c, d, e))), None)
            if sugar is None:  # ceramide centres and the C22 methylene (prochiral H names) are not on the sugar models
                continue
            t, x, p, q, r = (np.array(sugar[n]) for n in (target, centre, c, d, e))
            self.assertGreater((t - x) @ np.cross(q - p, r - q), 0, f"GM3 {centre} reads flipped on the reference")
            checked += 1
        self.assertEqual(checked, 16)  # 5 glucose, 5 galactose, 6 sialic acid centres
        # GM3 carries the same ceramide as PSM: its C4=C5 double bond and amide need the same trans definitions
        trans = lambda blk: sorted(tuple(l.split()) for l in blk.partition("[ trans ]")[2].split("[")[0].splitlines() if l.strip())
        psm_block = (DATA / "map.dat").read_text().split("RESI PSM")[1].split("RESI ")[0]
        self.assertEqual(trans(gm3_block), trans(psm_block))
        for a, b, c, d in trans(gm3_block):
            self.assertTrue({frozenset((a, b)), frozenset((b, c)), frozenset((c, d))} <= {frozenset(x) for x in bonds})

    def test_gm3_linkage_beads_hold_the_linked_atoms(self):
        """The Martini 3 DPG3 bonds between residues must join map.dat beads that contain the glycosidic bond."""
        import re
        import xml.etree.ElementTree as ET
        block = (DATA / "map.dat").read_text().split("RESI GM3")[1].split("RESI ")[0].split("[ chiral ]")[0]
        beads, current = {}, None
        for line in block.splitlines():
            found = re.match(r"\[\s*(\S+)\s*\]", line.strip())
            if found:
                current = beads.setdefault(found.group(1), [])
            elif current is not None:
                current += line.split()
        gm3 = next(r for r in ET.parse(DATA / "GM3.xml").getroot().iter("Residue") if r.get("name") == "GM3")
        bonds = {frozenset((b.get("atomName1"), b.get("atomName2"))) for b in gm3.iter("Bond")}
        rename = {(res, cg): aa for res, names in mf.MARTINI3_GLYCOLIPIDS["GM3"] for cg, aa in names.items()}
        # gm3_final.itp (DPG3, martini3001 v1.0): bonds GLC A-CER AM1, GLC B-GAL A, GAL B-NMC A
        for (r1, b1), (r2, b2), link in ((("GLC", "A"), ("CER", "AM1"), ("O1", "C1S")),
                                         (("GLC", "B"), ("GAL", "A"), ("O4", "C7")),
                                         (("GAL", "B"), ("NMC", "A"), ("O10", "C15"))):
            first, second = beads[rename[(r1, b1)]], beads[rename[(r2, b2)]]
            self.assertIn(frozenset(link), bonds)
            self.assertTrue(link[0] in first and link[1] in second, f"{r1} {b1}-{r2} {b2} does not hold {'-'.join(link)}")

    def test_worker_declares_tested_versions(self):
        self.assertIn("0.3.9", mf.backmapping.TESTED_MSTOOL)


class RingPiercing(unittest.TestCase):
    """A bond through the inside of a ring must be found; a bond beside it, or the ring's own bonds, must not."""

    def setUp(self):
        angle = np.arange(6) * np.pi / 3
        self.ring_xyz = np.c_[0.14 * np.cos(angle), 0.14 * np.sin(angle), np.zeros(6)] + 2.0
        self.rings, self.sizes, self.cell = np.array([[0, 1, 2, 3, 4, 5]]), np.array([6]), np.array([5.0, 5.0, 5.0])

    def piercings(self, start, end, shift=0.0):
        xyz = np.vstack([self.ring_xyz, [start, end]]) + shift
        return mf.ring_piercings(np.mod(xyz, self.cell), self.cell, self.rings, self.sizes, np.array([[6, 7], [0, 1]]))

    def test_bond_through_ring_centre(self):
        self.assertEqual(self.piercings([2.0, 2.0, 1.93], [2.0, 2.0, 2.08]), [(0, 0)])

    def test_bond_beside_ring(self):
        self.assertEqual(self.piercings([2.3, 2.0, 1.93], [2.3, 2.0, 2.08]), [])

    def test_bond_above_ring_does_not_cross(self):
        self.assertEqual(self.piercings([2.0, 2.0, 2.05], [2.0, 2.0, 2.20]), [])

    def test_detected_across_the_periodic_boundary(self):
        self.assertEqual(self.piercings([2.0, 2.0, 1.93], [2.0, 2.0, 2.08], shift=np.array([3.0, 0.0, 3.0])), [(0, 0)])

    def verdict(self, pierced, contact=0.2, review=None, last=False):
        staged = {"rings": {"pierced": pierced}, "solute_atoms": 8,
                  "contact": {"distance_nm": contact, "atoms": ["POPS:H13Y(atom 20)", "ILE:CD(atom 3)"]}}
        return mf.backmap_verdict(staged, review or {}, mf.Settings(), last)

    def test_backmap_verdict(self):
        poke = {"ring": ["PHE:CG(atom 5)"], "bond": ["DOPE:C25(atom 9)", "DOPE:H5R(atom 10)"], "ring_atoms": [5, 6], "bond_atoms": [9, 10]}
        heavy = {"ring": ["PHE:CG(atom 5)"], "bond": ["DOPE:C25(atom 9)", "DOPE:C26(atom 12)"], "ring_atoms": [5, 6], "bond_atoms": [9, 12]}
        internal = {"ring": ["TRP:CE2(atom 2)"], "bond": ["TRP:NE1(atom 6)", "TRP:CZ2(atom 7)"], "ring_atoms": [2, 3], "bond_atoms": [6, 7]}
        self.assertTrue(self.verdict([])["accept"])
        self.assertIn("threads the ring", self.verdict([heavy])["problem"])
        self.assertIn("threads the ring", self.verdict([heavy], last=True)["problem"])       # never tolerated
        self.assertIn("pokes through", self.verdict([poke])["problem"])
        self.assertTrue(self.verdict([poke], last=True)["accept"])                           # left to EM and the audit
        self.assertTrue(self.verdict([internal])["accept"])                                  # re-backmapping cannot change the input
        self.assertEqual(self.verdict([internal])["counts"]["solute_internal_ring_threadings"], 1)
        self.assertIn("starts 0.08 A from", self.verdict([], contact=0.008)["problem"])
        self.assertIn("wrong configuration", self.verdict([], review={"chiral_well_formed": 1})["problem"])
        self.assertIn("after 5 backmapping attempt(s)", failure(mf.pipeline.refuse_backmap, "x", 5))

    def test_only_heavy_atom_threading_triggers_a_new_backmap(self):
        report = {"pierced": [{"ring": ["PHE:CG(atom 5)"], "bond": ["DOPE:C25(atom 9)", "DOPE:H5R(atom 10)"]},
                              {"ring": ["PHE:CG(atom 5)"], "bond": ["DOPE:C25(atom 9)", "DOPE:C26(atom 12)"]}]}
        heavy = mf.heavy_atom_piercings(report)
        self.assertEqual([p["bond"][1] for p in heavy], ["DOPE:C26(atom 12)"])
        self.assertIn("minimization cannot undo", failure(mf.pipeline.refuse_threaded_rings, heavy, 3))
        self.assertIsNone(failure(mf.pipeline.refuse_threaded_rings, [], 3))

    def test_five_membered_ring_and_stretched_bond(self):
        angle = np.arange(5) * 2 * np.pi / 5
        ring = np.c_[0.12 * np.cos(angle), 0.12 * np.sin(angle), np.zeros(5)] + 2.0
        xyz = np.vstack([ring, [[2.0, 2.0, 1.65], [2.0, 2.0, 2.35]]])  # a 0.7 nm bond straight through the centre
        found = mf.ring_piercings(xyz, self.cell, np.array([[0, 1, 2, 3, 4, -1]]), np.array([5]), np.array([[5, 6]]))
        self.assertEqual(found, [(0, 0)])

    def test_rings_from_bond_graph(self):
        fused = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (5, 6), (6, 7), (7, 8), (8, 4), (8, 9)]  # indole-like 6+5
        self.assertEqual(sorted(len(r) for r in mf.molecule_rings(fused)), [5, 6])
        self.assertEqual(len(mf.molecule_rings(mf.itp_bonds([FORCEFIELD / "toppar/CHL1.itp"])["CHL1"])), 4)


class AuditGeometry(unittest.TestCase):
    """Regression: a chain written in another periodic image was reported as 4.6 nm of drift (wave-2 audit failures)."""

    def test_periodic_image_is_not_drift(self):
        rng = np.random.default_rng(1)
        before = rng.uniform(1.0, 4.0, size=(50, 3))
        cell = np.array([5.0, 5.0, 8.0])
        after = before + rng.normal(scale=0.01, size=before.shape)
        after[:20] += cell * np.array([1, 0, -1])  # one "chain" wrapped by the writer
        self.assertGreater(mf.superposed_rmsd(before, after), 1.0)
        self.assertLess(mf.superposed_rmsd(before, mf.same_image(before, after, cell)), 0.03)

    def test_superposition_removes_rigid_motion(self):
        rng = np.random.default_rng(2)
        P = rng.normal(size=(30, 3))
        R = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        self.assertLess(mf.superposed_rmsd(P, P @ R.T + 3.0), 1e-9)


class AuditTopology(unittest.TestCase):
    def test_includes_follow_ifdef_and_else(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "top.top").write_text('#include "a.itp"\n#ifdef POSRES\n#include "p.itp"\n#else\n#include "q.itp"\n#endif\n')
            for name in ("a.itp", "p.itp", "q.itp"):
                (root / name).write_text("; empty\n")
            self.assertEqual([p.name for p in mf.topology_includes(root / "top.top", set())], ["top.top", "a.itp", "q.itp"])
            self.assertEqual([p.name for p in mf.topology_includes(root / "top.top", {"POSRES"})], ["top.top", "a.itp", "p.itp"])
            (root / "q.itp").unlink()
            self.assertIn("does not exist", failure(mf.topology_includes, root / "top.top", set()))

    def test_molecules_and_atoms_are_parsed(self):
        atoms = mf.itp_atoms([FORCEFIELD / "toppar/POPC.itp"])["POPC"]
        self.assertEqual(len(atoms), 134)
        self.assertAlmostEqual(sum(q for _, q in atoms), 0.0, places=4)
        self.assertEqual(len(mf.itp_bonds([FORCEFIELD / "toppar/POPC.itp"])["POPC"]), 133)


class AuditDihedralRestraints(unittest.TestCase):
    """CHARMM-GUI dihedral restraints: configurations and double bonds fail the audit, rings out of their chair are reported."""

    def test_glpa_restraints_are_parsed(self):
        rows = mf.itp_dihedral_restraints([FORCEFIELD / "toppar/GLPA.itp"])["GLPA"]
        targets = [t for _, t in rows]
        self.assertEqual(len(rows), 23)
        self.assertEqual(sum(abs(abs(t) - 120) < 1 for t in targets), 4)  # ceramide C2S, C3S; Neu5Ac C7, C8
        self.assertEqual(sum(abs(abs(t) - 60) < 1 for t in targets), 18)  # three pyranose chairs
        self.assertEqual(sum(abs(t) == 180 for t in targets), 1)  # sphingosine C4=C5 trans
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "toppar").mkdir()
            (Path(tmp) / "toppar/POPC.itp").write_text((FORCEFIELD / "toppar/POPC.itp").read_text())
            self.assertTrue(mf.has_dihedral_restraints(Path(tmp)))
            (Path(tmp) / "toppar/POPC.itp").write_text("[ moleculetype ]\nX 1\n")
            self.assertFalse(mf.has_dihedral_restraints(Path(tmp)))

    def test_dihedrals_use_the_minimum_image(self):
        quad = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 1.0]]) + 2.0
        cell = np.array([4.0, 4.0, 4.0])
        wrapped = quad.copy()
        wrapped[3] -= cell * np.array([0, 1, 0])  # written in another periodic image
        for xyz in (quad, wrapped):
            self.assertAlmostEqual(float(mf.dihedrals_pbc(xyz, cell, np.array([[0, 1, 2, 3]]))[0]), -90.0, places=6)

    def test_wrong_configuration_fails_and_ring_is_reported(self):
        itp = ("[ moleculetype ]\nX 1\n[ atoms ]\n" + "".join(f"{i} C 1 X C{i} {i} 0.0 12.0\n" for i in range(1, 5))
               + "[ dihedral_restraints ]\n1 2 3 4 1 -120.0 2.5 1000\n1 2 3 4 1 60.0 2.5 1000\n")
        good = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 1.0]]  # -90: the side of the -120 target
        bad = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, -1.0]]  # +90: the mirror image
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            (run / "x.itp").write_text(itp)
            (run / "topol.top").write_text('#include "x.itp"\n[ system ]\nt\n[ molecules ]\nX 2\n')
            lines = [f"{1 + k // 4:5d}X    {'C' + str(k % 4 + 1):>5s}{k + 1:5d}{x + 1:8.3f}{y + 1:8.3f}{z + 2 + 3 * (k // 4):8.3f}"
                     for k, (x, y, z) in enumerate(good + bad)]
            (run / "em.gro").write_text("t\n8\n" + "\n".join(lines) + "\n  10.0 10.0 10.0\n")
            report, fails = mf.check_dihedral_restraints(run, "em.gro")
        self.assertEqual(report["X"]["wrong_configurations"], 1)  # -120 row: only the mirror image is wrong
        self.assertEqual(len(fails), 1)
        self.assertIn("X: 1 configuration restraint(s) on the wrong side", fails[0])
        self.assertEqual((report["X"]["rings"], report["X"]["rings_out_of_chair"]), (2, 1))  # +60 row: reported only


class StageRunner(unittest.TestCase):
    def setUp(self):
        self.out = Path(tempfile.mkdtemp())
        (self.out / mf.LOG_NAME).write_text("")
        self.session = mf.Session(out=self.out, name="t", gmx="gmx", forcefield=FORCEFIELD, data=None, python=sys.executable,
                                  ntomp=1, nsteps=1)

    def test_failures_become_stage_failures_with_a_hint(self):
        def broken():
            raise SystemExit("nope")
        with self.assertRaises(mf.StageFailure) as caught:
            mf.run_stage(self.session, "topology", broken)
        self.assertEqual((caught.exception.stage, caught.exception.what), ("topology", "nope"))
        self.assertIn(mf.LOG_NAME, caught.exception.inspect)
        self.assertIn("topology", self.session.timings)

    def test_unexpected_exceptions_are_named(self):
        with self.assertRaises(mf.StageFailure) as caught:
            mf.run_stage(self.session, "box", lambda: 1 / 0)
        self.assertIn("ZeroDivisionError", caught.exception.what)

    def test_zero_backmap_attempts_is_refused(self):
        self.session.settings = mf.Settings(backmap_attempts=0)
        with self.assertRaises(mf.StageFailure):
            mf.build_from_two_inputs(self.session, AA_PDB, CG_GRO, {})

    def test_wrong_stereochemistry_is_refused(self):
        self.assertIn("wrong configuration", failure(mf.pipeline.refuse_wrong_stereochemistry, 2, 3))
        self.assertIsNone(failure(mf.pipeline.refuse_wrong_stereochemistry, 0, 3))


class StructureIO(unittest.TestCase):
    def test_pdb_round_trip(self):
        atoms, _ = mf.read_pdb(AA_PDB)
        with tempfile.TemporaryDirectory() as tmp:
            mf.write_pdb(atoms[:500], Path(tmp) / "x.pdb", "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1")
            again, cryst = mf.read_pdb(Path(tmp) / "x.pdb")
        self.assertTrue(cryst.startswith("CRYST1"))
        self.assertEqual([(a["atom"], a["resname"], a["resid"], a["chain"]) for a in again],
                         [(a["atom"], a["resname"], a["resid"], a["chain"]) for a in atoms[:500]])
        self.assertTrue(np.allclose(mf.xyz_nm(again), mf.xyz_nm(atoms[:500]), atol=1e-3))

    def test_cg_pdb_is_read_in_nm(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cg.pdb"
            path.write_text("CRYST1   50.000   50.000   60.000  90.00  90.00  90.00 P 1           1\n"
                            "ATOM      1  BB  GLY A   1      10.000  10.000  10.000  1.00  0.00\n")
            atoms, box = mf.read_cg(path)
        self.assertEqual(box, [5.0, 5.0, 6.0])
        self.assertAlmostEqual(atoms[0]["x"], 1.0)

    def test_itp_parsing(self):
        itp = mf.read_itp(FORCEFIELD / "toppar/POPC.itp")
        self.assertEqual((itp["mol"], len(itp["atoms"])), ("POPC", 134))
        self.assertAlmostEqual(sum(a["charge"] for a in itp["atoms"]), 0.0, places=4)
        self.assertEqual(len(mf.read_itp(FORCEFIELD / "toppar/GLPA.itp")["atoms"]), len(mf.GM3_XML_TO_GLPA))

    def test_gm3_plan_follows_glpa_itp(self):
        order = [a["atom"] for a in mf.read_itp(FORCEFIELD / "toppar/GLPA.itp")["atoms"]]
        self.assertEqual([dst for _, dst in mf.GM3_XML_TO_GLPA], order)
        self.assertEqual(sorted(src for src, _ in mf.GM3_XML_TO_GLPA), sorted(MAPPING["GM3"]["atoms"]))

    def test_mapping_atoms_match_gromacs_topologies(self):
        for aa, itp_name in (("CHL1", "CHL1"), ("POPC", "POPC"), ("SAPI25", "SAPI25"), ("PSM", "PSM"), ("DOPS", "DOPS")):
            itp = [a["atom"] for a in mf.read_itp(FORCEFIELD / f"toppar/{itp_name}.itp")["atoms"]]
            self.assertEqual(sorted(MAPPING[aa]["atoms"]), sorted(itp), aa)


class Alignment(unittest.TestCase):
    def test_semi_global_alignment(self):
        pairs = mf.align_sequences("AAAKLMNPQRSTVWY", "KLMNPQRXTVWYGG")
        self.assertEqual(pairs[0], (3, 0))
        self.assertEqual(len(pairs), 12)

    def test_kabsch_recovers_a_known_transform(self):
        rng = np.random.default_rng(0)
        P = rng.normal(size=(40, 3))
        angle = 0.7
        R = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        Q = P @ R.T + np.array([3.0, -2.0, 5.0])
        found_R, found_t = mf.kabsch(P, Q)
        self.assertTrue(np.allclose(found_R, R, atol=1e-9))
        self.assertTrue(np.allclose(found_t, [3.0, -2.0, 5.0], atol=1e-9))
        self.assertAlmostEqual(np.linalg.det(found_R), 1.0)

    def test_segments_are_unwrapped_across_the_boundary(self):
        cell = np.array([50.0, 50.0, 50.0])
        chain = [residue(i + 1, "GLY", ["BB"]) for i in range(6)]
        for i, beads in enumerate(chain):
            beads[0].update(x=(4.8 + 0.35 * i) % 5.0, y=1.0, z=1.0)
        segments = mf.cg_protein_segments(chain, cell, 1.0)
        self.assertEqual(len(segments), 1)
        steps = np.linalg.norm(np.diff(segments[0]["xyz"], axis=0), axis=1)
        self.assertTrue(np.allclose(steps, 3.5, atol=1e-6))

    def test_chain_breaks_split_segments(self):
        chain = [residue(i + 1, "GLY", ["BB"]) for i in range(4)] + [residue(1, "ALA", ["BB", "SC1"])]
        for i, beads in enumerate(chain):
            beads[0].update(x=1.0 + 0.35 * i, y=1.0, z=1.0)
        self.assertEqual(len(mf.cg_protein_segments(chain, np.array([90.0] * 3), 1.0)), 2)


class Disulfides(unittest.TestCase):
    def setUp(self):
        atoms, _ = mf.read_pdb(AA_PDB)
        self.chains = {}
        for a in atoms:
            if a["resname"] in mf.AMINO:
                self.chains.setdefault(a["chain"], []).append(a)

    def test_detected_from_geometry(self):
        self.assertEqual(sorted(mf.find_disulfides(self.chains)["R"]), [(43, 67), (58, 100), (81, 121), (224, 294)])

    def test_inter_chain_disulfide_is_refused(self):
        chains = copy.deepcopy(self.chains)
        sg = next(a for a in chains["R"] if a["resname"] == "CYS" and a["atom"] == "SG")
        other = next(a for a in chains["B"] if a["resname"] == "CYS" and a["atom"] == "SG")
        other.update(x=sg["x"] + 2.0, y=sg["y"], z=sg["z"])
        self.assertIn("inter-chain disulfide", failure(mf.find_disulfides, chains))

    def test_three_way_contact_is_refused(self):
        chains = copy.deepcopy(self.chains)
        sgs = [a for a in chains["R"] if a["resname"] == "CYS" and a["atom"] == "SG"]
        free = next(a for a in sgs if a["resid"] not in (43, 67, 58, 100, 81, 121, 224, 294))
        bonded = next(a for a in sgs if a["resid"] == 43)
        free.update(x=bonded["x"] + 2.0, y=bonded["y"], z=bonded["z"])
        self.assertIn("ambiguous", failure(mf.find_disulfides, chains))


class Subprocesses(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        (self.tmp / mf.LOG_NAME).write_text("")

    def test_non_zero_exit(self):
        message = failure(mf.run_command, self.tmp, [sys.executable, "-c", "import sys; open('p.out','w').write('x'); sys.exit(3)"],
                          produces=("p.out",))
        self.assertIn("failed (exit 3)", message)

    def test_stale_output_cannot_satisfy_a_step(self):
        (self.tmp / "stale.out").write_text("old")
        message = failure(mf.run_command, self.tmp, [sys.executable, "-c", "pass"], produces=("stale.out",))
        self.assertIn("exited 0 but did not write", message)
        self.assertFalse((self.tmp / "stale.out").exists())

    def test_empty_output_is_rejected(self):
        self.assertIn("did not write", failure(mf.run_command, self.tmp, [sys.executable, "-c", "open('e.out','w')"],
                                               produces=("e.out",)))

    def test_missing_executable(self):
        self.assertIn("cannot run", failure(mf.run_command, self.tmp, ["/nonexistent/gmx", "x"]))
        self.assertIn("GROMACS not found", failure(mf.find_gromacs, "/nonexistent/gmx"))

    def test_python_traceback_is_summarised(self):
        message = failure(mf.run_command, self.tmp, [sys.executable, "-c", "import module_that_is_not_installed"])
        self.assertIn("ModuleNotFoundError", message)

    def test_success_returns_output_and_logs_debug(self):
        output = mf.run_command(self.tmp, [sys.executable, "-c", "print('hello'); open('g.out','w').write('1')"], produces=("g.out",))
        self.assertIn("hello", output)
        self.assertIn("DEBUG: $", (self.tmp / mf.LOG_NAME).read_text())

    def test_log_levels(self):
        mf.log(self.tmp, "x", "CHECK")
        self.assertIn("CHECK: x", (self.tmp / mf.LOG_NAME).read_text())
        with self.assertRaises(ValueError):
            mf.log(self.tmp, "x", "LOUD")


class OutputProtection(unittest.TestCase):
    def test_input_inside_output_directory_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / "membrane.pdb").write_text("x")
            self.assertIn("would be overwritten", failure(mf.protect_inputs, [out / "membrane.pdb"], out))
            mf.protect_inputs([out / "membrane.pdb"], out, keep=out / "membrane.pdb")  # the --membrane route keeps its source
            mf.protect_inputs([out / "my_inputs.pdb"], out)

    def test_input_in_a_generated_subdirectory_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            self.assertIn("would be overwritten", failure(mf.protect_inputs, [out / "work" / "aa.pdb"], out))
            self.assertIn("would be overwritten", failure(mf.protect_inputs, [out / "toppar" / "x" / "cg.gro"], out))

    def test_membrane_exemption_covers_only_membrane_pdb(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            mf.protect_inputs([out / "membrane.pdb"], out, keep=out / "membrane.pdb")
            self.assertIn("would be overwritten", failure(mf.protect_inputs, [out / "prot-memb.pdb"], out, keep=out / "prot-memb.pdb"))

    def test_output_inside_an_installation_directory_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            forcefield = Path(tmp) / "ff"
            for out in (forcefield, forcefield / "runs" / "a", Path(tmp)):
                self.assertIn("overlaps the installation directory", failure(mf.protect_inputs, [], out, resources=(forcefield, None)))
            mf.protect_inputs([], Path(tmp) / "elsewhere", resources=(forcefield, None))

    def test_stale_outputs_are_cleared(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / "em.gro").write_text("old")
            (out / "toppar").mkdir()
            (out / "toppar" / "old.itp").write_text("old")
            (out / "notes.txt").write_text("mine")
            removed = mf.clear_stale_outputs(out)
            self.assertEqual(sorted(removed), ["em.gro", "toppar"])
            self.assertEqual([p.name for p in out.iterdir()], ["notes.txt"])

    def test_changed_input_is_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "a.pdb"
            path.write_text("one")
            hashes = {path: mf.sha256(path)}
            mf.check_inputs_unchanged(hashes)
            path.write_text("two")
            self.assertIn("changed during the build", failure(mf.check_inputs_unchanged, hashes))


if __name__ == "__main__":
    unittest.main()
