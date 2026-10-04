from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
MODULE_PATH = ROOT / "membraneforger" / "stages" / "legacy_impl.py"
SPEC = importlib.util.spec_from_file_location("membraneforger_stage_impl", MODULE_PATH)
workflow = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
sys.modules["membraneforger_stage_impl"] = workflow
SPEC.loader.exec_module(workflow)


def _write_mdp_files(root: Path) -> None:
    mdp_dir = root / "resources" / "cg_mdp"
    mdp_dir.mkdir(parents=True)
    for name in ("step6.0_minimization.mdp", *workflow.CG_EQUILIBRATION_MDPS):
        (mdp_dir / name).write_text("integrator = steep\n", encoding="utf-8")


def _write_stage2_inputs(root: Path) -> Path:
    out = root / "outputs" / "run" / "stage2"
    toppar = out / "toppar"
    toppar.mkdir(parents=True)
    (out / "cg_topology.top").write_text(
        "\n".join(
            [
                '#include "toppar/Protein_0.itp"',
                '#include "toppar/POPC.itp"',
                '#include "toppar/W.itp"',
                "",
                "[ molecules ]",
                "Protein_0 1",
                "POPC 6",
                "W 2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    (toppar / "Protein_0.itp").write_text("[ moleculetype ]\nProtein_0 1\n[ atoms ]\n1 P1 1 PROT BB 1 0\n", encoding="utf-8")
    (toppar / "POPC.itp").write_text("[ moleculetype ]\nPOPC 1\n[ atoms ]\n1 P4 1 POPC PO4 1 0\n", encoding="utf-8")
    (toppar / "W.itp").write_text("[ moleculetype ]\nW 1\n[ atoms ]\n1 P4 1 W W 1 0\n", encoding="utf-8")
    atoms = [
        {"resid": 1, "resname": "PROT", "name": "BB", "x": 1.0, "y": 1.0, "z": 3.5},
        *[
            {"resid": 2 + i, "resname": "POPC", "name": "PO4", "x": 2.0 + i, "y": 1.0, "z": z}
            for i, z in enumerate((2.0, 2.1, 2.2, 5.0, 5.1, 5.2))
        ],
        {"resid": 8, "resname": "W", "name": "W", "x": 1.0, "y": 2.0, "z": 1.0},
        {"resid": 9, "resname": "W", "name": "W", "x": 1.0, "y": 3.0, "z": 6.0},
    ]
    workflow.write_gro(out / "cg_scaffold_input.gro", "INITIAL", atoms, "10.00000 10.00000 10.00000")
    return out


def test_stage2_uses_final_cg_equilibration_for_backmap_handoff(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    _write_mdp_files(tmp_path)
    out = _write_stage2_inputs(tmp_path)
    calls: list[tuple[str, list[str]]] = []

    def fake_run_logged(cmd, cwd, log_dir, label, timeout=None, input_text=None, env=None):
        calls.append((label, cmd))
        if cmd[1] == "grompp":
            (cwd / cmd[cmd.index("-o") + 1]).write_text("tpr\n", encoding="utf-8")
        if cmd[1] == "mdrun":
            _, atoms, box_line, _ = workflow.parse_gro(out / "cg_scaffold_input.gro")
            output = cmd[cmd.index("-c") + 1]
            title = "FINAL_EQ" if "step6_6" in output else ("MINIMIZED" if output == "cg_minimized.gro" else "EQUIL")
            workflow.write_gro(cwd / output, title, atoms, box_line)
            deffnm = cmd[cmd.index("-deffnm") + 1]
            (cwd / f"{deffnm}.cpt").write_text("cpt\n", encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(workflow, "gmx_command", lambda: ["gmx"])
    monkeypatch.setattr(workflow, "gmx_mdrun_command", lambda: ["gmx", "mdrun"])
    monkeypatch.setattr(workflow, "run_logged", fake_run_logged)

    report = workflow.run_cg_minimize_equilibrate(
        tmp_path, out, tmp_path / "logs", {}, out / "cg_scaffold_input.gro", out / "cg_topology.top",
        "stage2_scaffold_equil", "stage2_scaffold_grompp", "stage2_scaffold_mdrun",
    )

    assert report["sequence"][-1]["mdp"] == "step6.6_equilibration.mdp"
    assert "FINAL_EQ" in (out / "cg_equilibrated_scaffold.gro").read_text(encoding="utf-8").splitlines()[0]
    assert (out / "cg_backmap_input.gro").read_text(encoding="utf-8") == (
        out / "cg_equilibrated_scaffold.gro"
    ).read_text(encoding="utf-8")
    assert (out / "cg_backmap_input.gro").read_text(encoding="utf-8") != (out / "cg_minimized.gro").read_text(encoding="utf-8")
    assert any(label == "stage2_scaffold_equil_step6_6_equilibration_mdrun" for label, _ in calls)


def test_stage2_failed_cg_equilibration_blocks_handoff(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    _write_mdp_files(tmp_path)
    out = _write_stage2_inputs(tmp_path)

    def fake_run_logged(cmd, cwd, log_dir, label, timeout=None, input_text=None, env=None):
        if label == "stage2_scaffold_equil_step6_3_equilibration_mdrun":
            return subprocess.CompletedProcess(cmd, 1, "", "failed")
        if cmd[1] == "mdrun":
            _, atoms, box_line, _ = workflow.parse_gro(out / "cg_scaffold_input.gro")
            output = cmd[cmd.index("-c") + 1]
            workflow.write_gro(cwd / output, "PARTIAL", atoms, box_line)
            deffnm = cmd[cmd.index("-deffnm") + 1]
            (cwd / f"{deffnm}.cpt").write_text("cpt\n", encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(workflow, "gmx_command", lambda: ["gmx"])
    monkeypatch.setattr(workflow, "gmx_mdrun_command", lambda: ["gmx", "mdrun"])
    monkeypatch.setattr(workflow, "run_logged", fake_run_logged)

    with pytest.raises(workflow.ValidationError, match="step6.3_equilibration.mdp"):
        workflow.run_cg_minimize_equilibrate(
            tmp_path, out, tmp_path / "logs", {}, out / "cg_scaffold_input.gro", out / "cg_topology.top",
            "stage2_scaffold_equil", "stage2_scaffold_grompp", "stage2_scaffold_mdrun",
        )
    assert not (out / "cg_backmap_input.gro").exists()
    assert not (out / "cg_equilibrated_scaffold.gro").exists()


def test_finalize_insane_topology_uses_stage1_moleculetype_names(tmp_path: Path) -> None:
    martini = tmp_path / "resources" / "forcefields" / "martini"
    martini.mkdir(parents=True)
    (martini / "martini_v3.0.0.itp").write_text("; base Martini FF\n", encoding="utf-8")
    (martini / "martini_v3.0.0_phospholipids_v1.itp").write_text(
        "[ moleculetype ]\nPOPC 1\n[ atoms ]\n1 P4 1 POPC PO4 1 0\n",
        encoding="utf-8",
    )
    (martini / "martini_v3.0.0_solvents_v1.itp").write_text(
        "[ moleculetype ]\nW 1\n[ atoms ]\n1 P4 1 W W 1 0\n",
        encoding="utf-8",
    )
    stage1 = tmp_path / "outputs" / "run" / "stage1"
    stage1.mkdir(parents=True)
    (stage1 / "ReplacementProtein_A.itp").write_text(
        "[ moleculetype ]\nActualProteinA 1\n[ atoms ]\n1 P1 1 PROT BB 1 0\n",
        encoding="utf-8",
    )
    work_stage1 = tmp_path / "work" / "run" / "stage1"
    work_stage1.mkdir(parents=True)
    (work_stage1 / "replacement_protein.top").write_text("[ molecules ]\nActualProteinA 1\n", encoding="utf-8")
    out = tmp_path / "outputs" / "run" / "stage2"
    toppar = out / "toppar"
    toppar.mkdir(parents=True)
    for src in martini.glob("*.itp"):
        (toppar / src.name).write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
    (out / "cg_topology.top").write_text(
        '#include "martini.itp"\n\n[ molecules ]\nProtein 1\nPOPC 1\nW 1\n',
        encoding="utf-8",
    )
    atoms = [
        {"resid": 1, "resname": "PROT", "name": "BB", "x": 1.0, "y": 1.0, "z": 1.0},
        {"resid": 2, "resname": "POPC", "name": "PO4", "x": 2.0, "y": 1.0, "z": 1.0},
        {"resid": 3, "resname": "W", "name": "W", "x": 3.0, "y": 1.0, "z": 1.0},
    ]
    workflow.write_gro(out / "cg_system.gro", "CG", atoms, "5.00000 5.00000 5.00000")

    report = workflow.finalize_insane_topology(tmp_path, {"coarse_grained": {"martini_forcefield": "martini3001"}}, "run", out / "cg_topology.top", stage1, out)
    text = (out / "cg_topology.top").read_text(encoding="utf-8")
    assert '#include "toppar/martini_v3.0.0.itp"' in text
    assert '#include "toppar/ReplacementProtein_A.itp"' in text
    assert "ActualProteinA" in text
    assert "Protein 1" not in text
    assert report["status"] == "PASS"
    assert report["protein_moleculetype_names"] == ["ActualProteinA"]


def test_validate_topology_contract_rejects_unresolved_include(tmp_path: Path) -> None:
    top = tmp_path / "cg_topology.top"
    top.write_text('#include "toppar/missing.itp"\n\n[ molecules ]\nPOPC 1\n', encoding="utf-8")
    with pytest.raises(workflow.ValidationError, match="unresolved include"):
        workflow.validate_topology_contract(top)


def test_validate_topology_contract_rejects_unknown_molecule(tmp_path: Path) -> None:
    toppar = tmp_path / "toppar"
    toppar.mkdir()
    (toppar / "POPC.itp").write_text("[ moleculetype ]\nPOPC 1\n[ atoms ]\n1 P4 1 POPC PO4 1 0\n", encoding="utf-8")
    top = tmp_path / "cg_topology.top"
    top.write_text('#include "toppar/POPC.itp"\n\n[ molecules ]\nDOES_NOT_EXIST 1\n', encoding="utf-8")
    with pytest.raises(workflow.ValidationError, match="lack matching"):
        workflow.validate_topology_contract(top)


def test_local_thread_count_respects_slurm_allocation(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("MEMBRANEFORGER_LOCAL_THREADS", raising=False)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "6")
    assert workflow.local_thread_count() == 6
