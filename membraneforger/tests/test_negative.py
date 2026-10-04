"""Negative tests: every invalid state must be refused, and a refusal must never leave a verified em.gro."""
import copy
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from .common import AA_PDB, CG_GRO, CHOL3, FORCEFIELD, MAPPING, POPC, REPO, failure, mf, residue, small_system

AA_LINES = [l for l in AA_PDB.read_text().splitlines() if l[:6] in ("ATOM  ", "HETATM")]
GRO_LINES = CG_GRO.read_text().splitlines()


class Martini2AndUnknownContent(unittest.TestCase):
    def test_martini2_alanine(self):
        self.assertIn("Martini 2 protein bead inventory", failure(mf.classify_cg, residue(1, "ALA", ["BB"]) + small_system()[2:], MAPPING))

    def test_martini2_tryptophan(self):
        bad = residue(1, "TRP", ["BB", "SC1", "SC2", "SC3", "SC4"]) + small_system()[8:]
        self.assertIn("Martini 2 protein", failure(mf.classify_cg, bad, MAPPING))

    def test_martini2_cholesterol(self):
        bad = small_system()[:21] + residue(5, "CHOL", CHOL3[:6] + CHOL3[7:])
        message = failure(mf.classify_cg, bad, MAPPING)
        self.assertIn("not the Martini 3 bead set", message)
        self.assertIn("Martini 2 inventory", message)

    def test_martini2_water_and_ions(self):
        self.assertIn("Martini 2 polarizable water", failure(mf.classify_cg, small_system() + residue(8, "PW", ["W", "WP", "WM"]), MAPPING))
        self.assertIn("Martini 2 ion naming", failure(mf.classify_cg, small_system() + residue(8, "ION", ["NA+"]), MAPPING))

    def test_unknown_residue_is_counted_not_dropped(self):
        bad = small_system() + residue(8, "XYZ", ["Q1", "Q2"]) + residue(9, "XYZ", ["Q1", "Q2"])
        self.assertIn("2 x XYZ [Q1 Q2]", failure(mf.classify_cg, bad, MAPPING))

    def test_unexpected_bead_inventory(self):
        self.assertIn("unknown residue or bead set", failure(mf.classify_cg, small_system() + residue(8, "POPC", POPC[:-1] + ["ZZ9"]), MAPPING))
        self.assertIn("not the Martini 3 beads BB SC1 SC2",
                      failure(mf.classify_cg, residue(1, "LYS", ["BB", "SC1", "SC2", "CA"]) + small_system(), MAPPING))

    def test_incomplete_glycolipid(self):
        bad = small_system() + residue(8, "GLC", ["A", "B", "C", "V"]) + residue(9, "GAL", ["A", "B", "C", "V"])
        self.assertIn("not a complete Martini 3 GM3", failure(mf.classify_cg, bad, MAPPING))

    def test_unsupported_mapping(self):
        without = {k: v for k, v in MAPPING.items() if k != "CHL1"}
        self.assertIn("installed mapping for CHL1", failure(mf.classify_cg, small_system(), without))

    def test_membrane_residue_without_gromacs_topology(self):
        self.assertIn("no GROMACS topology", failure(mf.pipeline.require_membrane_topologies, [{"aa": "DPPC"}], FORCEFIELD))

    def test_no_protein_or_no_membrane(self):
        self.assertIn("no Martini 3 protein", failure(mf.classify_cg, small_system()[9:], MAPPING))
        self.assertIn("no membrane lipids", failure(mf.classify_cg, small_system()[:9] + small_system()[-2:], MAPPING))


class InvalidCoarseGrainedInput(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def read(self, name, text):
        path = self.tmp / name
        path.write_text(text)
        return failure(mf.read_cg, path)

    def test_truncated_gro(self):
        self.assertIn("header declares", self.read("t.gro", "\n".join(GRO_LINES[:500]) + "\n"))

    def test_wrong_extension(self):
        self.assertIn("expected a .pdb or .gro", self.read("x.xyz", "hi"))

    def test_pdb_without_box(self):
        self.assertIn("no CRYST1", self.read("n.pdb", "ATOM      1  BB  GLY A   1       1.000   1.000   1.000  1.00  0.00\n"))

    def test_triclinic_box_is_refused(self):
        self.assertIn("not orthorhombic", self.read("t.gro", "t\n1\n    1GLY     BB    1   1.000   1.000   1.000\n 5 5 5 0 0 1 0 0 0\n"))
        self.assertIn("not orthorhombic", self.read("t.pdb", "CRYST1   50.000   50.000   60.000  90.00  90.00 120.00 P 1           1\n"
                                                             "ATOM      1  BB  GLY A   1      10.000  10.000  10.000  1.00  0.00\n"))

    def test_non_finite_coordinate(self):
        self.assertIn("non-finite", self.read("n.gro", "t\n1\n    1GLY     BB    1     nan   1.000   1.000\n 5 5 5\n"))

    def test_multiple_bilayers_are_refused(self):
        atoms, box = mf.read_cg(CG_GRO)
        membrane = copy.deepcopy(mf.classify_cg(atoms, MAPPING)["membrane"])
        for mol in membrane[::2]:  # move every other lipid half a box up: two stacked bilayers
            for b in mol["beads"]:
                b["z"] = (b["z"] + box[2] / 2) % box[2]
        self.assertIn("expected one planar bilayer", failure(mf.make_membrane_whole, membrane, box))


class InvalidAllAtomInput(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def read(self, name, lines):
        path = self.tmp / name
        path.write_text("\n".join(lines) + "\n")
        return failure(mf.read_all_atom, path, FORCEFIELD)

    def test_empty(self):
        self.assertIn("no ATOM/HETATM", self.read("e.pdb", ["REMARK"]))

    def test_malformed_record(self):
        self.assertIn("malformed", self.read("b.pdb", AA_LINES[:5] + [AA_LINES[5][:30] + "   abc  " + AA_LINES[5][38:]]))

    def test_no_protein(self):
        self.assertIn("no protein residues", self.read("l.pdb", [l for l in AA_LINES if l[17:20] == "GDP"]))

    def test_duplicate_atoms(self):
        self.assertIn("duplicate atoms", self.read("d.pdb", AA_LINES[:50] + AA_LINES[:1]))

    def test_solvent(self):
        water = "HETATM 9999  OH2 TIP3W   1       0.000   0.000   0.000  1.00  0.00"
        self.assertIn("solvent or bulk ions", self.read("w.pdb", AA_LINES[:50] + [water]))

    def test_ligand_without_topology(self):
        message = self.read("u.pdb", AA_LINES[:50] + ["HETATM 9999  C1  ZZZ X   1       0.000   0.000   0.000  1.00  0.00"])
        self.assertIn("no topology", message)
        self.assertIn("ZZZ", message)

    def test_insertion_code_is_refused(self):
        self.assertIn("insertion codes", self.read("i.pdb", AA_LINES[:50] + [AA_LINES[50][:26] + "A" + AA_LINES[50][27:]]))

    def test_non_finite_coordinate(self):
        self.assertIn("non-finite", self.read("n.pdb", AA_LINES[:50] + [AA_LINES[50][:30] + "     nan" + AA_LINES[50][38:]]))


class MappingRefusals(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.aa, _ = mf.read_all_atom(AA_PDB, FORCEFIELD)
        atoms, cls.box = mf.read_cg(CG_GRO)
        cls.raw = atoms
        cls.cg = mf.classify_cg(atoms, MAPPING)
        cls.slab = mf.make_membrane_whole(cls.cg["membrane"], cls.box)

    def test_cg_protein_without_all_atom_counterpart(self):
        other = [dict(a, resname="GLY") if a["resname"] in mf.MARTINI3_PROTEIN else a for a in self.raw]
        other = [a for a in other if not (a["resname"] == "GLY" and a["atom"] != "BB")]
        protein = mf.classify_cg(other, MAPPING)["protein"]
        self.assertIn("no all-atom counterpart", failure(mf.map_all_atom_to_cg, self.aa, protein, self.box, self.slab))

    def test_unacceptable_fit(self):
        bent = [dict(a, z=a["z"] * 0.3) if a["chain"] == "R" else a for a in self.aa]
        self.assertIn("do not superpose", failure(mf.map_all_atom_to_cg, bent, self.cg["protein"], self.box, self.slab))

    def test_fit_limits_are_configurable(self):
        strict = mf.Settings(fit_max_core_rmsd_a=0.5)
        self.assertIn("limits 0.5 A", failure(mf.map_all_atom_to_cg, self.aa, self.cg["protein"], self.box, self.slab, strict))

    def test_a_displaced_chain_is_refused_even_when_the_complex_core_fits(self):
        moved = [dict(a, x=a["x"] + 60.0) if a["chain"] == "P" else a for a in self.aa]
        message = failure(mf.map_all_atom_to_cg, moved, self.cg["protein"], self.box, self.slab)
        self.assertIn("chain P does not sit where the coarse-grained frame has it", message)

    def test_a_moderately_shifted_chain_is_kept_with_a_warning(self):
        shifted = [dict(a, x=a["x"] + 7.0) if a["chain"] == "P" else a for a in self.aa]
        fit = mf.map_all_atom_to_cg(shifted, self.cg["protein"], self.box, self.slab)
        self.assertTrue(any("chain P deviates" in note and "keeps its all-atom pose" in note for note in fit["notes"]))

    def test_ambiguous_identical_chains(self):
        twins = [a for a in self.aa if a["chain"] == "P"] + [dict(a, chain="Y") for a in self.aa if a["chain"] == "P"]
        peptide = self.cg["protein"][-27:]
        both = [b for r in peptide for b in r] + [dict(b) for r in peptide for b in r]
        self.assertIn("ambiguous", failure(mf.map_all_atom_to_cg, twins, mf.cg_residues(both), self.box, self.slab))


class BackmappingRefusals(unittest.TestCase):
    def test_lost_molecule(self):
        membrane = [{"aa": "POPC"}, {"aa": "POPC"}]
        one = [[{"atom": n, "resname": "POPC", "x": 0.0, "y": 0.0, "z": 0.0} for n in MAPPING["POPC"]["atoms"]]]
        self.assertIn("changed the membrane inventory: POPC 1/2", failure(mf.check_backmapped_inventory, one, membrane, MAPPING))

    def test_incomplete_molecule(self):
        short = [[{"atom": n, "resname": "POPC", "x": 0.0, "y": 0.0, "z": 0.0} for n in MAPPING["POPC"]["atoms"][:-3]]]
        self.assertIn("incomplete", failure(mf.check_backmapped_inventory, short, [{"aa": "POPC"}], MAPPING))

    def test_non_finite_coordinates(self):
        bad = [[{"atom": n, "resname": "POPC", "x": float("nan"), "y": 0.0, "z": 0.0} for n in MAPPING["POPC"]["atoms"]]]
        self.assertIn("non-finite", failure(mf.check_backmapped_inventory, bad, [{"aa": "POPC"}], MAPPING))


class CommandLineRefusals(unittest.TestCase):
    """End-to-end refusals through the real entry point; none may leave em.gro behind."""

    def run_cli(self, *args):
        result = subprocess.run([sys.executable, "-m", "membraneforger", *map(str, args)], cwd=REPO, text=True, capture_output=True)
        return result.returncode, result.stdout + result.stderr

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.gmx = os.environ.get("MEMBRANEFORGER_GMX") or shutil.which("gmx") or "gmx"

    def test_missing_input_and_missing_arguments(self):
        self.assertIn("missing input file", self.run_cli("--all-atom", "nope.pdb", "--coarse-grain", CG_GRO)[1])
        self.assertIn("missing input file", self.run_cli("--aa", AA_PDB, "--cg", "nope.gro")[1])
        self.assertIn("--aa is required", self.run_cli("--cg", CG_GRO)[1])
        self.assertIn("--aa is required", self.run_cli("--cg", "1")[1])

    def test_bundled_membrane_is_the_default_and_needs_no_cg(self):
        out = self.tmp / "bundled"
        code, output = self.run_cli("--aa", AA_PDB, "--out", out, "--gmx", self.gmx, "--mstool-python", "/usr/bin/python3")
        self.assertEqual(code, 1)
        self.assertIn("KOR1_cg_cellmem.gro", output)  # the bundled frame was picked up as input
        self.assertIn("ERROR: mstool:", output)  # and the build stopped at the mstool stage, not at the parser
        code, output = self.run_cli("--aa", AA_PDB, "--cg", "2", "--out", out, "--gmx", self.gmx, "--mstool-python", "/usr/bin/python3")
        self.assertIn("GPR1_cg_cellmem.gro", output)

    def test_lipid_and_box_options_are_checked(self):
        self.assertIn("unknown lipid XYZ", self.run_cli("--aa", AA_PDB, "--dellipid", "XYZ")[1])
        self.assertIn("unknown lipid", self.run_cli("--aa", AA_PDB, "--addlipid", "DPPC")[1])
        self.assertIn("three positive edge lengths", self.run_cli("--aa", AA_PDB, "--box", "80", "80", "0")[1])
        self.assertIn("need --aa and --cg", self.run_cli("--membrane", AA_PDB, "--box", "80", "80", "100")[1])

    def test_wrong_executable(self):
        code, output = self.run_cli("--all-atom", AA_PDB, "--coarse-grain", CG_GRO, "--out", self.tmp / "o", "--gmx", "/no/gmx")
        self.assertNotEqual(code, 0)
        self.assertIn("GROMACS not found", output)

    def test_missing_mapping_data(self):
        code, output = self.run_cli("--all-atom", AA_PDB, "--coarse-grain", CG_GRO, "--out", self.tmp / "o", "--gmx", self.gmx,
                                    "--data", self.tmp)
        self.assertNotEqual(code, 0)
        self.assertIn("has no map.dat", output)

    def test_interpreter_without_mstool(self):
        out = self.tmp / "nomstool"
        code, output = self.run_cli("--all-atom", AA_PDB, "--coarse-grain", CG_GRO, "--out", out, "--gmx", self.gmx,
                                    "--mstool-python", "/usr/bin/python3")
        self.assertEqual(code, 1)
        self.assertIn("ERROR: mstool:", output)
        self.assertIn("Inspect:", output)
        self.assertFalse((out / "em.gro").exists())
        self.assertTrue((out / "run_manifest.json").exists())

    def test_stale_outputs_do_not_survive_a_failed_run(self):
        out = self.tmp / "stale"
        out.mkdir()
        (out / "em.gro").write_text("old valid-looking")
        (out / "topol.top").write_text("old")
        bad = self.tmp / "bad.gro"
        bad.write_text("t\n1\n    1XYZ     Q1    1   1.000   1.000   1.000\n 5 5 5\n")
        code, output = self.run_cli("--all-atom", AA_PDB, "--coarse-grain", bad, "--out", out, "--gmx", self.gmx)
        self.assertEqual(code, 1)
        self.assertIn("ERROR: classify:", output)
        self.assertFalse((out / "em.gro").exists())
        self.assertFalse((out / "topol.top").exists())

    def test_input_overwrite_is_refused_and_input_untouched(self):
        out = self.tmp / "clash"
        out.mkdir()
        shutil.copy(AA_PDB, out / "membrane.pdb")
        before = mf.sha256(out / "membrane.pdb")
        code, output = self.run_cli("--all-atom", out / "membrane.pdb", "--coarse-grain", CG_GRO, "--out", out, "--gmx", self.gmx)
        self.assertEqual(code, 1)
        self.assertIn("would be overwritten", output)
        self.assertEqual(mf.sha256(out / "membrane.pdb"), before)


if __name__ == "__main__":
    unittest.main()
