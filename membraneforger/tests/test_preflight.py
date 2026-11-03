"""--check (preflight.py), the dependency guard of __main__.py and the equilibration mdp templates."""
import contextlib
import io
import re
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from .. import __main__ as entry
from .. import preflight
from .common import AA_PDB, CG_GRO, DATA, FORCEFIELD, REPO, mf

MDP_DIR = REPO / "membraneforger" / "mdp"


def statuses(checks):
    """{item: status} of a check list (the last line of an item wins)."""
    return {c.item: c.status for c in checks}


def check(all_atom=AA_PDB, coarse_grain=CG_GRO, out=None, gmx_problem="GROMACS not found", orientation=None, **kwargs):
    out = out or Path(tempfile.mkdtemp()) / "new_build"
    orientation = orientation or mf.OrientationRequest(mode="none")
    return preflight.check_build(all_atom=all_atom, coarse_grain=coarse_grain, membrane=None, out=out, forcefield=FORCEFIELD,
                                 data=DATA, gmx=None, gmx_problem=gmx_problem, python=sys.executable, orientation=orientation,
                                 ntomp=4, **kwargs), out


class PreflightOnTheBundledExample(unittest.TestCase):
    def test_reports_every_item_and_writes_nothing(self):
        checks, out = check()
        found = statuses(checks)
        self.assertEqual(found["GROMACS"], preflight.FAIL)                 # no gmx in the test environment
        self.assertEqual(found["Python requirements"], preflight.OK)
        self.assertEqual(found["all-atom input"], preflight.OK)
        self.assertEqual(found["coarse-grained frame"], preflight.OK)
        self.assertEqual(found["backmapping data"], preflight.OK)
        self.assertEqual(found["output directory"], preflight.OK)
        self.assertEqual(found["histidines"], preflight.WARN)              # GCGR has HIS residues that become HSD
        self.assertEqual(found["chain R termini"], preflight.INFO)
        self.assertIn("ACE", next(c.text for c in checks if c.item == "chain R termini"))
        self.assertIn("R:43-67", next(c.text for c in checks if c.item == "disulfides"))
        self.assertFalse(out.exists())                                     # nothing was created
        self.assertFalse(any(c.item == "PPM 3.0" for c in checks))         # not with --orientation none

    def test_orientation_needs_a_chain_when_there_are_several(self):
        checks, _ = check(orientation=mf.OrientationRequest(mode="auto"))
        orientation = [c for c in checks if c.item == "orientation"]
        self.assertTrue(any(c.status == preflight.FAIL and "--orient-chain" in c.text for c in orientation))
        self.assertTrue(any(c.status == preflight.WARN and "--nterm-side" in c.text for c in orientation))
        checks, _ = check(orientation=mf.OrientationRequest(mode="auto", chains=("R",), nterm_side="out"))
        self.assertFalse([c for c in checks if c.item == "orientation" and c.status == preflight.FAIL])
        checks, _ = check(orientation=mf.OrientationRequest(mode="auto", chains=("Z",), nterm_side="out"))
        self.assertTrue(any(c.item == "orientation" and c.status == preflight.FAIL and "Z" in c.text for c in checks))

    def test_ppm_is_a_warning_in_auto_mode_and_a_failure_in_ppm_mode(self):
        with mock.patch.dict("os.environ", {"MEMBRANEFORGER_PPM": ""}), mock.patch("shutil.which", return_value=None):
            auto, _ = check(orientation=mf.OrientationRequest(mode="auto", chains=("R",), nterm_side="out"))
            ppm, _ = check(orientation=mf.OrientationRequest(mode="ppm", chains=("R",), nterm_side="out"))
        self.assertEqual(statuses(auto)["PPM 3.0"], preflight.WARN)
        self.assertEqual(statuses(ppm)["PPM 3.0"], preflight.FAIL)

    def test_a_chain_without_a_terminal_oxygen_fails_before_the_build(self):
        lines = [l for l in AA_PDB.read_text().splitlines() if not (l.startswith(("ATOM", "HETATM")) and l[21] == "A" and l[12:16].strip() == "OXT")]
        tmp = Path(tempfile.mkdtemp())
        (tmp / "no_oxt.pdb").write_text("\n".join(lines) + "\n")
        checks, _ = check(all_atom=tmp / "no_oxt.pdb")
        self.assertEqual(statuses(checks)["chain A termini"], preflight.FAIL)
        self.assertIn("OXT", next(c.text for c in checks if c.item == "chain A termini"))

    def test_an_input_the_build_refuses_is_a_failure_with_the_same_message(self):
        tmp = Path(tempfile.mkdtemp())
        (tmp / "water.pdb").write_text("\n".join(AA_PDB.read_text().splitlines()[:60]
                                                 + ["HETATM 9999  OH2 TIP3W   1       0.000   0.000   0.000  1.00  0.00"]) + "\n")
        checks, _ = check(all_atom=tmp / "water.pdb")
        self.assertEqual(statuses(checks)["all-atom input"], preflight.FAIL)
        self.assertIn("solvent or bulk ions", next(c.text for c in checks if c.item == "all-atom input"))

    def test_an_earlier_build_in_the_output_directory_is_a_failure_without_overwrite(self):
        out = Path(tempfile.mkdtemp()) / "taken"
        out.mkdir()
        (out / "em.gro").write_text("x\n")
        checks, _ = check(out=out)
        self.assertEqual(statuses(checks)["output directory"], preflight.FAIL)
        checks, _ = check(out=out, overwrite=True)
        self.assertEqual(statuses(checks)["output directory"], preflight.OK)

    def test_print_checks_exit_code_and_verdict_line(self):
        checks, _ = check()
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = preflight.print_checks(checks, "demo")
        self.assertEqual(code, 1)
        self.assertIn("demo\tCHECK FAIL", buffer.getvalue())
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = preflight.print_checks([c for c in checks if c.status != preflight.FAIL], "demo")
        self.assertEqual(code, 0)
        self.assertIn("demo\tCHECK PASS", buffer.getvalue())

    def test_the_command_line_runs_the_check_without_gromacs(self):
        out = Path(tempfile.mkdtemp()) / "cli_check"
        buffer = io.StringIO()
        with mock.patch("shutil.which", return_value=None), contextlib.redirect_stdout(buffer):
            code = mf.main(["--check", "--aa", str(AA_PDB), "--cg", str(CG_GRO), "--orientation", "none", "--out", str(out)])
        self.assertEqual(code, 1)
        self.assertIn("FAIL: GROMACS", buffer.getvalue())
        self.assertIn("OK: all-atom input", buffer.getvalue())
        self.assertFalse(out.exists())


class DependencyGuard(unittest.TestCase):
    def test_nothing_is_missing_here(self):
        self.assertEqual(entry.missing_dependencies(), [])

    def test_a_missing_package_is_named_by_its_pip_name(self):
        self.assertEqual(entry.missing_dependencies({"numpy": "numpy", "no_such_module_xyz": "pandas<3"}), ["pandas<3"])

    def test_run_prints_the_install_line_instead_of_a_traceback(self):
        err = io.StringIO()
        with mock.patch.object(entry, "missing_dependencies", return_value=["openmm", "pandas<3"]), contextlib.redirect_stderr(err):
            code = entry.run(["--list-membranes"])
        self.assertEqual(code, 2)
        self.assertIn('pip install openmm "pandas<3"', err.getvalue())


class EquilibrationTemplates(unittest.TestCase):
    STEPS = [f"step6.{i}_equilibration.mdp" for i in range(1, 7)]

    def read(self, name):
        text = (MDP_DIR / name).read_text()
        return {k.strip(): v.split(";")[0].strip() for k, v in
                (line.split("=", 1) for line in text.splitlines() if "=" in line and not line.lstrip().startswith(";"))}

    def test_every_equilibration_step_sets_the_restraint_macros_the_topology_expects(self):
        for name in self.STEPS:
            define = self.read(name)["define"]
            for macro in ("-DPOSRES ", "-DPOSRES_FC_BB=", "-DPOSRES_FC_SC=", "-DPOSRES_FC_LIPID=", "-DDIHRES ", "-DDIHRES_FC="):
                self.assertIn(macro, define + " ", name)
            self.assertNotIn("BBDIHRES", define, f"{name} must never define the EM-only backbone restraints")

    def test_the_force_constants_decrease_and_end_at_the_backbone_only(self):
        bb = [float(re.search(r"POSRES_FC_BB=([\d.]+)", self.read(n)["define"]).group(1)) for n in self.STEPS]
        self.assertEqual(bb, sorted(bb, reverse=True))
        last = self.read(self.STEPS[-1])["define"]
        self.assertIn("-DPOSRES_FC_SC=0.0", last)
        self.assertIn("-DPOSRES_FC_LIPID=0.0", last)
        self.assertIn("-DDIHRES_FC=0.0", last)

    def test_production_has_no_define_and_the_groups_are_the_index_groups(self):
        production = self.read("step7_production.mdp")
        self.assertNotIn("define", production)
        self.assertEqual(production["pcoupl"], "Parrinello-Rahman")
        for name in self.STEPS + ["step7_production.mdp"]:
            values = self.read(name)
            self.assertEqual(values["tc_grps"].split(), list(mf.INDEX_GROUPS[1:]), name)
            self.assertEqual(values["comm_grps"].split(), list(mf.INDEX_GROUPS[1:]), name)
            self.assertEqual((values["rvdw"], values["rcoulomb"], values["vdw-modifier"]), ("1.2", "1.2", "Force-switch"), name)

    def test_the_schedule_is_the_charmm_gui_length(self):
        total_ps = sum(float(self.read(n)["dt"]) * int(self.read(n)["nsteps"]) for n in self.STEPS)
        self.assertAlmostEqual(total_ps, 1875.0)
        self.assertEqual(self.read(self.STEPS[0])["gen-vel"], "yes")
        self.assertEqual(self.read(self.STEPS[0])["pcoupl"], "no")
        self.assertEqual(self.read(self.STEPS[2])["pcoupl"], "C-rescale")


if __name__ == "__main__":
    unittest.main()
