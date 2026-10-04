from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import aa_stage4
from membraneforger.config import normalize_config
from membraneforger.stages import legacy_impl


def write_protein_pdb(path: Path, chains: tuple[str, ...] = ("A",)) -> None:
    atoms = []
    serial = 1
    for chain in chains:
        for resid, resname, offset in ((1, "ALA", 0.0), (2, "GLY", 4.0)):
            for name, dx, element in (("N", 0.0, "N"), ("CA", 1.4, "C"), ("C", 2.8, "C"), ("O", 3.6, "O")):
                atoms.append({
                    "record": "ATOM",
                    "name": name,
                    "resname": resname,
                    "chain": chain,
                    "resid": resid,
                    "icode": "",
                    "x": offset + dx,
                    "y": float(ord(chain) - ord("A")) * 5.0,
                    "z": 0.0,
                    "element": element,
                    "serial": serial,
                })
                serial += 1
    aa_stage4._write_pdb(path, [], atoms, np.array([50.0, 50.0, 50.0]))


def test_charmm_terminal_patch_answers_use_ace_and_ct3(tmp_path: Path) -> None:
    pdb = tmp_path / "protein.pdb"
    write_protein_pdb(pdb)
    answers = legacy_impl.stage4_terminal_patch_answers({}, pdb)
    assert answers == ["ACE", "CT3"]


def test_stage4_pdb2gmx_input_puts_disulfides_before_terminal_patches(tmp_path: Path) -> None:
    pdb = tmp_path / "protein.pdb"
    write_protein_pdb(pdb)
    answers = legacy_impl.stage4_pdb2gmx_input(
        {"pdb2gmx_disulfide_answers": ["y", "n"], "n_terminal_patch": "ACE", "c_terminal_patch": "CT3"},
        pdb,
    )
    assert answers.splitlines() == ["y", "n", "ACE", "CT3"]


def test_inferred_disulfide_answers_do_not_overrun_terminal_patch_prompts(tmp_path: Path) -> None:
    pdb = tmp_path / "protein.pdb"
    atoms = []
    serial = 1
    for resid, x in ((1, 0.0), (2, 2.05), (3, 12.0)):
        for name, dx, element in (("N", -1.0, "N"), ("CA", -0.5, "C"), ("C", 0.5, "C"), ("O", 1.0, "O"), ("SG", 0.0, "S")):
            atoms.append({
                "record": "ATOM",
                "name": name,
                "resname": "CYS",
                "chain": "A",
                "resid": resid,
                "icode": "",
                "x": x + dx,
                "y": 0.0,
                "z": 0.0,
                "element": element,
                "serial": serial,
            })
            serial += 1
    aa_stage4._write_pdb(pdb, [], atoms, np.array([50.0, 50.0, 50.0]))
    assert legacy_impl.inferred_disulfide_prompt_count(pdb) == 1
    answers = legacy_impl.stage4_pdb2gmx_input({"n_terminal_patch": "ACE", "c_terminal_patch": "CT3"}, pdb)
    assert answers.splitlines() == ["y", "ACE", "CT3"]


def test_unavailable_charmm_terminal_patch_fails_before_pdb2gmx(tmp_path: Path) -> None:
    charmm = tmp_path / "charmm36.ff"
    charmm.mkdir()
    (charmm / "merged.n.tdb").write_text("[ None ]\n[ NH3+ ]\n", encoding="utf-8")
    (charmm / "merged.c.tdb").write_text("[ None ]\n[ COO- ]\n", encoding="utf-8")
    with pytest.raises(legacy_impl.ValidationError, match="n_terminal_patch=ACE"):
        legacy_impl.validate_charmm_terminal_patch_resources(charmm, {"n_terminal_patch": "ACE", "c_terminal_patch": "CT3"})


def test_available_charmm_terminal_patches_pass_resource_validation(tmp_path: Path) -> None:
    charmm = tmp_path / "charmm36.ff"
    charmm.mkdir()
    (charmm / "merged.n.tdb").write_text("[ None ]\n[ ACE ]\n", encoding="utf-8")
    (charmm / "merged.c.tdb").write_text("[ None ]\n[ CT3 ]\n", encoding="utf-8")
    legacy_impl.validate_charmm_terminal_patch_resources(charmm, {"n_terminal_patch": "ACE", "c_terminal_patch": "CT3"})


def test_legacy_c_terminal_cap_nme_normalizes_to_ct3() -> None:
    cfg = normalize_config({"stage4": {"c_terminal_cap": "NME"}})
    assert cfg["stage4"]["c_terminal_patch"] == "CT3"
    assert "c_terminal_cap" not in cfg["stage4"]


def test_conflicting_cap_and_patch_settings_fail() -> None:
    with pytest.raises(Exception, match="conflicting terminal configuration"):
        normalize_config({"stage4": {"c_terminal_patch": "none", "c_terminal_cap": "NME"}})


def test_stage4_rejects_explicit_synthetic_nme_residue(tmp_path: Path) -> None:
    path = tmp_path / "nme.pdb"
    aa_stage4._write_pdb(
        path,
        [],
        [
            {"record": "ATOM", "name": "CA", "resname": "ALA", "chain": "A", "resid": 1, "icode": "", "x": 1.0, "y": 1.0, "z": 1.0, "element": "C"},
            {"record": "ATOM", "name": "CH3", "resname": "NME", "chain": "A", "resid": 2, "icode": "", "x": 2.0, "y": 1.0, "z": 1.0, "element": "C"},
        ],
        np.array([50.0, 50.0, 50.0]),
    )
    with pytest.raises(aa_stage4.Stage4Error, match="explicit terminal cap residues"):
        aa_stage4.reject_explicit_charmm_terminal_cap_residues(path, "Stage 3 canonical final_all_atom.pdb")


def test_synthetic_explicit_ace_residue_is_not_required_for_patch_selection(tmp_path: Path) -> None:
    pdb = tmp_path / "protein.pdb"
    write_protein_pdb(pdb)
    labels = aa_stage4.explicit_charmm_terminal_cap_residues(pdb)
    answers = legacy_impl.stage4_terminal_patch_answers({"n_terminal_patch": "ACE", "c_terminal_patch": "CT3"}, pdb)
    assert labels == []
    assert answers == ["ACE", "CT3"]


def test_stage4_rejects_pdb2gmx_explicit_nme_residue(tmp_path: Path) -> None:
    pdb = tmp_path / "protein.pdb"
    write_protein_pdb(pdb)
    gro = tmp_path / "protein.gro"
    gro.write_text(
        "\n".join([
            "protein",
            "3",
            "    1NME    CH3    1   0.000   0.000   0.000",
            "    1NME      N    2   0.100   0.000   0.000",
            "    1NME      C    3   0.200   0.000   0.000",
            "   5.00000   5.00000   5.00000",
        ]) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(aa_stage4.Stage4Error, match="pdb2gmx generated explicit terminal cap residues"):
        aa_stage4.replace_protein_with_pdb2gmx_coordinates(pdb, gro, tmp_path / "out.pdb")


def test_terminal_residue_identity_survives_stage4_replacement(tmp_path: Path) -> None:
    pdb = tmp_path / "protein.pdb"
    write_protein_pdb(pdb)
    gro = tmp_path / "protein.gro"
    gro.write_text(
        "\n".join([
            "protein",
            "9",
            "    1ALA      N    1   0.000   0.000   0.000",
            "    1ALA     CA    2   0.140   0.000   0.000",
            "    1ALA      C    3   0.280   0.000   0.000",
            "    1ALA    CAY    4  -0.100   0.000   0.000",
            "    2GLY      N    5   0.400   0.000   0.000",
            "    2GLY     CA    6   0.540   0.000   0.000",
            "    2GLY      C    7   0.680   0.000   0.000",
            "    2GLY     OT1   8   0.780   0.000   0.000",
            "    2GLY    CAT    9   0.860   0.000   0.000",
            "   5.00000   5.00000   5.00000",
        ]) + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.pdb"
    aa_stage4.replace_protein_with_pdb2gmx_coordinates(pdb, gro, out)
    _, atoms, _ = aa_stage4._pdb_atoms(out)
    assert {atom["resname"] for atom in atoms if atom["resid"] == 1} == {"ALA"}
    assert {atom["resname"] for atom in atoms if atom["resid"] == 2} == {"GLY"}
    assert "NME" not in {atom["resname"] for atom in atoms}
    assert "ACE" not in {atom["resname"] for atom in atoms}


def test_multichain_terminal_patches_apply_per_chain(tmp_path: Path) -> None:
    pdb = tmp_path / "multichain.pdb"
    write_protein_pdb(pdb, ("A", "B"))
    answers = legacy_impl.stage4_terminal_patch_answers({}, pdb)
    assert answers == ["ACE", "CT3", "ACE", "CT3"]


def test_requested_absent_terminal_patch_chain_fails(tmp_path: Path) -> None:
    pdb = tmp_path / "protein.pdb"
    write_protein_pdb(pdb, ("A",))
    with pytest.raises(legacy_impl.ValidationError, match="absent protein chains"):
        legacy_impl.stage4_terminal_patch_answers({"terminal_patch_chains": ["Z"]}, pdb)


def test_multichain_duplicate_resid_ambiguity_remains_active(tmp_path: Path) -> None:
    pdb = tmp_path / "multichain.pdb"
    write_protein_pdb(pdb, ("A", "B"))
    gro = tmp_path / "protein.gro"
    gro.write_text(
        "\n".join([
            "protein",
            "6",
            "    1ALA      N    1   0.000   0.000   0.000",
            "    1ALA     CA    2   0.100   0.000   0.000",
            "    1ALA      C    3   0.200   0.000   0.000",
            "    1ALA      N    4   0.000   0.300   0.000",
            "    1ALA     CA    5   0.100   0.300   0.000",
            "    1ALA      C    6   0.200   0.300   0.000",
            "   5.00000   5.00000   5.00000",
        ]) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(aa_stage4.Stage4Error, match="duplicate residue numbers across chains"):
        aa_stage4.replace_protein_with_pdb2gmx_coordinates(pdb, gro, tmp_path / "out.pdb")
