from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "membraneforger" / "stages" / "legacy_impl.py"
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location("membraneforger_stage_impl", MODULE_PATH)
workflow = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
sys.modules["membraneforger_stage_impl"] = workflow
SPEC.loader.exec_module(workflow)


ALA_EXPECTED = ["N", "HN", "CA", "HA", "C", "O", "CB", "HB1", "HB2", "HB3"]


def ala_duplicate_residue(second_sidechain_offset: float = 0.0) -> pd.DataFrame:
    names = ["N", "HN", "CA", "HA", "C", "O", "CB", "HB1", "HB2", "HB3", "CB", "HB1", "HB2", "HB3"]
    coords = [
        (0.0, 0.0, 0.0),
        (0.0, 0.1, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 0.1, 0.0),
        (2.0, 0.0, 0.0),
        (2.1, 0.0, 0.0),
        (1.1, 2.7, 0.0),
        (1.2, 3.7, 0.0),
        (1.3, 3.6, 0.0),
        (1.4, 3.5, 0.0),
        (3.0 + second_sidechain_offset, 0.0, 0.0),
        (3.8 + second_sidechain_offset, 0.0, 0.0),
        (3.7 + second_sidechain_offset, 0.2, 0.0),
        (3.7 + second_sidechain_offset, -0.2, 0.0),
    ]
    return pd.DataFrame(
        {
            "id": list(range(len(names))),
            "chain": ["A"] * len(names),
            "resid": [1] * len(names),
            "resname": ["ALA"] * len(names),
            "name": names,
            "x": [c[0] for c in coords],
            "y": [c[1] for c in coords],
            "z": [c[2] for c in coords],
        }
    )


def test_ala_duplicate_orphaned_fragment_keeps_sc1_connected_set() -> None:
    residue = ala_duplicate_residue()
    drop, reason = workflow.select_duplicate_repair_indices(residue, ALA_EXPECTED, np.array([3.05, 0.0, 0.0]))
    kept_names = residue[~residue.index.isin(drop)]["name"].tolist()
    assert kept_names.count("CB") == 1
    assert 10 not in drop
    assert "CG SC1 provenance" in reason


def test_ambiguous_duplicate_coordinates_refuse_to_guess() -> None:
    residue = ala_duplicate_residue(second_sidechain_offset=-1.85)
    with pytest.raises(workflow.ValidationError):
        workflow.select_duplicate_repair_indices(residue, ALA_EXPECTED, np.array([1.125, 1.35, 0.0]))


def test_duplicate_residue_ids_in_different_chains_are_not_duplicates() -> None:
    residue = ala_duplicate_residue().iloc[:10].copy()
    two_chains = pd.concat([residue, residue.assign(chain="B")], ignore_index=True)
    rows = workflow.duplicate_atom_rows(two_chains, [])
    assert rows == [
        "chain\tresid\tinsertion_code\tresname\tatom_name\tduplicate_count\tatom_indices\tcoordinates\tbonded_neighbors\tselected_action\treason"
    ]


def test_protein_fragment_break_is_not_joined() -> None:
    atoms = pd.DataFrame(
        {
            "chain": ["A", "A", "A", "A", "A", "A"],
            "resid": [1, 1, 1, 2, 2, 2],
            "resname": ["GLY"] * 6,
            "name": ["N", "CA", "C", "N", "CA", "C"],
            "x": [0.0, 1.0, 2.0, 10.0, 11.0, 12.0],
            "y": [0.0] * 6,
            "z": [0.0] * 6,
        }
    )
    fragments = workflow.protein_fragments(atoms)
    assert len(fragments) == 2
    assert fragments[0]["break_after"].startswith("1->2 C-N")


def test_clean_residue_has_no_duplicate_repair() -> None:
    residue = ala_duplicate_residue().iloc[:10].copy()
    drop, reason = workflow.select_duplicate_repair_indices(residue, ALA_EXPECTED)
    assert drop == set()
    assert reason == "clean"


def test_membrane_template_finalizer_is_not_protein_terminal_repair() -> None:
    source = MODULE_PATH.read_text()
    assert "finalize_stage3_membrane_templates" in source
    assert hasattr(workflow, "terminalize_protein_fragments")


def test_mstool_ungroup_kwargs_match_pinned_api_without_ss(tmp_path: Path) -> None:
    def pinned_ungroup(
        structure,
        out,
        mapping,
        mapping_add,
        backbone,
        water_resname,
        water_number,
        water_chain_dms,
        sort,
        use_AA_structure,
        AA_shrink_factor,
    ):
        return None

    kwargs = workflow.mstool_ungroup_kwargs(tmp_path / "out.dms", ["protein.dat"], ["extra.dat"])
    assert "ss" not in kwargs
    report = workflow.preflight_mstool_ungroup_api(pinned_ungroup, kwargs)
    assert report["status"] == "PASS"
    with pytest.raises(workflow.ValidationError, match="unsupported argument"):
        workflow.preflight_mstool_ungroup_api(pinned_ungroup, {**kwargs, "ss": 3.5})


class DummyUniverse:
    def __init__(self, atoms: pd.DataFrame) -> None:
        self.atoms = atoms
        self.written: list[Path] = []

    def write(self, path: str) -> None:
        output = Path(path)
        output.write_text("derived DPG3-to-GM3 input\n", encoding="utf-8")
        self.written.append(output)


def dpg3_atoms() -> pd.DataFrame:
    rows = [{"resname": "ALA", "name": "BB", "resid": 1, "resn": 1}]
    rows.extend(
        {"resname": resname, "name": name, "resid": 2, "resn": 2}
        for resname, name in workflow.DPG3_SOURCE_PATTERN
    )
    return pd.DataFrame(rows)


def test_dpg3_normalizer_reconstructs_one_canonical_gm3(tmp_path: Path) -> None:
    universe = DummyUniverse(dpg3_atoms())
    report = workflow.normalize_dpg3_to_gm3(universe, tmp_path)
    gm3 = universe.atoms[universe.atoms["resname"] == "GM3"]
    assert report["status"] == "PASS"
    assert report["gm3_molecules"] == 1
    assert report["retained_beads"] == 20
    assert report["omitted_virtual_sites"] == 3
    assert gm3["name"].tolist() == list(workflow.GM3_CANONICAL_BEADS)
    assert gm3["resid"].nunique() == 1
    assert not universe.atoms["resname"].isin(workflow.DPG3_COMPONENT_RESNAMES).any()
    assert universe.written == [tmp_path / "cg_backmap_input_chained_gm3.dms"]


def test_dpg3_normalizer_rejects_partial_component_window(tmp_path: Path) -> None:
    universe = DummyUniverse(dpg3_atoms().iloc[:-1].copy())
    with pytest.raises(workflow.ValidationError, match="malformed DPG3 component window"):
        workflow.normalize_dpg3_to_gm3(universe, tmp_path)


def test_stage1_input_preprocessing_rejects_ligand_in_selected_chain(tmp_path: Path) -> None:
    source = tmp_path / "input.pdb"
    source.write_text(
        "\n".join(
            [
                "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N",
                "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C",
                "ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00  0.00           C",
                "ATOM      4  O   ALA A   1       3.000   0.000   0.000  1.00  0.00           O",
                "HETATM    5  C1  LIG A 900       4.000   0.000   0.000  1.00  0.00           C",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(workflow.ValidationError, match="nonprotein components"):
        workflow.write_protein_only_pdb(source, tmp_path / "protein_only.pdb", ["A"])


def test_stage1_input_preprocessing_can_exclude_nonprotein_chain(tmp_path: Path) -> None:
    source = tmp_path / "input.pdb"
    source.write_text(
        "\n".join(
            [
                "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N",
                "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C",
                "ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00  0.00           C",
                "ATOM      4  O   ALA A   1       3.000   0.000   0.000  1.00  0.00           O",
                "HETATM    5  C1  LIG B 900       4.000   0.000   0.000  1.00  0.00           C",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    report = workflow.write_protein_only_pdb(source, tmp_path / "protein_only.pdb", ["A"])
    assert report["kept_atom_records"] == 4
    assert report["skipped"]["unselected_chain"] == 1


def test_martinize2_command_does_not_force_coil_by_default() -> None:
    cmd, effective = workflow.martinize2_command("martinize2", "in.pdb", "out.pdb", "out.top", {})
    assert "-ss" not in cmd
    assert effective["secondary_structure_policy"] == "auto"
    fallback_cmd, fallback_effective = workflow.martinize2_command(
        "martinize2",
        "in.pdb",
        "out.pdb",
        "out.top",
        {"secondary_structure": {"policy": "fallback", "fallback": "C"}},
    )
    assert fallback_cmd[-2:] == ["-ss", "C"]
    assert fallback_effective["warnings"]
