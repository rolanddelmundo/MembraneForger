from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "membraneforger_workflow.py"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location("membraneforger_workflow_for_config", MODULE_PATH)
workflow = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
sys.modules["membraneforger_workflow_for_config"] = workflow
SPEC.loader.exec_module(workflow)


def test_public_config_keys_are_authoritative(tmp_path: Path) -> None:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "\n".join(
            [
                "input:",
                "  pdb: inputs/public_input.pdb",
                "coarse_grained:",
                "  martini_forcefield: martini3001",
                "backmapping:",
                "  mstool_root: resources/vendor/mstool",
                "all_atom:",
                "  salt_concentration_molar: 0.2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    loaded = workflow.load_config(cfg)
    assert loaded["stage1"]["input_pdb"] == "inputs/public_input.pdb"
    assert loaded["stage3"]["input_aa_protlig"] == "inputs/public_input.pdb"
    assert loaded["stage4"]["salt_concentration_molar"] == 0.2
    assert loaded["_meta"]["normalized_config_version"] == 1


def test_conflicting_public_and_legacy_values_are_rejected(tmp_path: Path) -> None:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "coarse_grained:\n  martini_forcefield: martini3001\nstage1:\n  martinize_forcefield: martini22\n",
        encoding="utf-8",
    )
    with pytest.raises(workflow.ContractError, match="conflicting duplicate"):
        workflow.load_config(cfg)


def test_unknown_scientific_config_key_is_rejected(tmp_path: Path) -> None:
    cfg = tmp_path / "config.yaml"
    cfg.write_text("simulation:\n  silently_ignored_option: true\n", encoding="utf-8")
    with pytest.raises(workflow.ContractError, match="unknown configuration key simulation.silently_ignored_option"):
        workflow.load_config(cfg)


def test_unimplemented_cg_production_option_is_rejected(tmp_path: Path) -> None:
    cfg = tmp_path / "config.yaml"
    cfg.write_text("simulation:\n  run_cg_production: true\n", encoding="utf-8")
    with pytest.raises(workflow.ContractError, match="run_cg_production is not implemented"):
        workflow.load_config(cfg)
