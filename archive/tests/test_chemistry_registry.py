from __future__ import annotations

import pytest

from membraneforger.chemistry import ChemistryError, component_for, validate_membrane_composition


def _config(component: str, value: float = 1.0, *, plugin: bool = False, all_atom: bool = True) -> dict:
    return {
        "membrane": {"composition": {"upper": {component: value}, "lower": {"POPC": 1}}},
        "backmapping": {"enabled": True},
        "all_atom": {"enabled": all_atom},
        "stage3": {"plugins": {"glycolipid_template_finalize": {"enabled": plugin}}},
    }


def test_registry_resolves_aliases() -> None:
    assert component_for("CHL1").name == "CHOL"
    assert component_for("GM3").name == "DPG3"


def test_validate_membrane_composition_accepts_supported_lipid() -> None:
    report = validate_membrane_composition(_config("POPC"))
    assert report["status"] == "PASS"
    assert report["normalized_insane_composition"]["upper"] == {"POPC": 1.0}


def test_validate_membrane_composition_rejects_nonpositive_values() -> None:
    with pytest.raises(ChemistryError, match="must be positive"):
        validate_membrane_composition(_config("POPC", 0.0))


def test_validate_membrane_composition_rejects_unknown_component() -> None:
    with pytest.raises(ChemistryError, match="unsupported membrane component"):
        validate_membrane_composition(_config("NOTALIPID"))


def test_validate_membrane_composition_requires_glycolipid_adapter_for_all_atom() -> None:
    with pytest.raises(ChemistryError, match="lacks all-atom topology support"):
        validate_membrane_composition(_config("GM3", plugin=False, all_atom=True))
    report = validate_membrane_composition(_config("GM3", plugin=True, all_atom=True))
    assert report["normalized_insane_composition"]["upper"] == {"DPG3": 1.0}
