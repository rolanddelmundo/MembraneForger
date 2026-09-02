from __future__ import annotations

from copy import deepcopy
from typing import Any


class ConfigError(RuntimeError):
    pass


PUBLIC_TO_LEGACY = {
    ("input", "pdb"): (("stage1", "input_pdb"), ("stage3", "input_aa_protlig")),
    ("coarse_grained", "martini_forcefield"): (("stage1", "martinize_forcefield"),),
    ("coarse_grained", "protein_mode"): (("stage1", "membrane_protein_mode"),),
    ("membrane", "build_mode"): (("stage2", "build_mode"),),
    ("backmapping", "mstool_root"): (("stage3", "resources", "vendor_mstool"),),
    ("all_atom", "charmm36_forcefield"): (("stage4", "resources", "charmm36_forcefield"),),
    ("all_atom", "membrane_toppar"): (("stage4", "resources", "membrane_toppar"),),
    ("all_atom", "ligand_params"): (("stage4", "ligand_params"),),
    ("all_atom", "salt_concentration_molar"): (("stage4", "salt_concentration_molar"),),
    ("all_atom", "n_terminal_patch"): (("stage4", "n_terminal_patch"),),
    ("all_atom", "c_terminal_patch"): (("stage4", "c_terminal_patch"),),
}


ALLOWED_KEYS = {
    (): {"run", "input", "membrane", "coarse_grained", "simulation", "backmapping", "all_atom", "compute", "review", "advanced", "stage1", "stage2", "stage3", "stage4", "_meta"},
    ("input",): {"pdb", "orientation"},
    ("membrane",): {"build_mode", "composition", "box", "protein_image_distance_nm", "z_distance_nm", "salt_concentration_molar"},
    ("membrane", "composition"): {"upper", "lower", "units"},
    ("membrane", "box"): {"shape", "x_nm", "y_nm", "z_nm"},
    ("coarse_grained",): {"martini_forcefield", "protein_mode"},
    ("simulation",): {"run_cg_minimization", "run_cg_equilibration", "run_cg_production"},
    ("backmapping",): {"enabled", "mstool_root"},
    ("all_atom",): {"enabled", "charmm36_forcefield", "membrane_toppar", "ligand_params", "salt_concentration_molar", "n_terminal_patch", "c_terminal_patch", "n_terminal_cap", "c_terminal_cap"},
    ("compute",): {"threads", "gmx", "ntmpi", "ntomp", "gpu"},
    ("review",): {"checkpoints", "require_input_review_on_warnings"},
    ("advanced",): {"legacy_scaffold_replace"},
    ("run",): {"id", "overwrite", "update_latest", "require_inputs_for_check"},
    ("stage1",): {"input_pdb", "orientation_enabled", "orientation_chain", "protein_chains", "membrane_protein_mode", "martinize_forcefield", "insane_input", "resources", "role", "secondary_structure", "position_restraints", "elastic_network"},
    ("stage1", "resources"): {"cg_mdp", "forcefields"},
    ("stage1", "secondary_structure"): {"policy", "assignment", "fallback"},
    ("stage1", "position_restraints"): {"enabled", "selection", "force_constant"},
    ("stage1", "elastic_network"): {"enabled", "options"},
    ("stage2",): {"mode", "production", "build_mode", "resources", "handoff", "scaffold_coordinates", "scaffold_topology", "scaffold_toppar", "salt_concentration_molar", "insane_timeout_seconds", "cg_grompp_timeout_seconds", "cg_minimization_timeout_seconds", "cg_equilibration_timeout_seconds", "max_cg_midplane_gap_nm"},
    ("stage2", "resources"): {"cg_mdp"},
    ("stage2", "handoff"): {"mode"},
    ("stage3",): {"input_aa_protlig", "input_opm_reference", "resources", "placement", "plugins"},
    ("stage3", "resources"): {"mappings", "templates", "membrane_templates", "membrane_toppar", "forcefields", "vendor_mstool"},
    ("stage3", "placement"): {"mode", "use_opm", "fixed_selection", "moving_selection", "fixed_chain", "moving_chain", "z_translation_nm", "tilt_degrees", "rotation_degrees"},
    ("stage3", "plugins"): {"glycolipid_template_finalize"},
    ("stage3", "plugins", "glycolipid_template_finalize"): {"enabled", "source_residue", "target_molecule_type", "target_itp"},
    ("stage4",): {"resources", "ligand_params", "dssp_required", "charmm36_required", "membrane_toppar_required", "cgenff_required", "ligand_params_required", "parameter_generation_enabled", "pyrosetta_plugins_enabled", "slice_buffer_nm", "minimum_slice_buffer_nm", "z_slab_padding_nm", "salt_concentration_molar", "pdb2gmx_timeout_seconds", "pdb2gmx_disulfide_answers", "pdb2gmx_disulfide_answer_count", "n_terminal_patch", "c_terminal_patch", "n_terminal_cap", "c_terminal_cap", "terminal_patch_chains", "solvent_membrane_z_clearance_nm", "hard_clash_cutoff_nm"},
    ("stage4", "resources"): {"aa_mdp", "charmm36_forcefield", "membrane_toppar", "ligand_params_root_env", "charmm36_root_env", "cgenff_root_env", "cgenff_root"},
}


def _has(data: dict[str, Any], path: tuple[str, ...]) -> bool:
    current: Any = data
    for key in path:
        if not isinstance(current, dict) or key not in current:
            return False
        current = current[key]
    return True


def _get(data: dict[str, Any], path: tuple[str, ...]) -> Any:
    current: Any = data
    for key in path:
        current = current[key]
    return current


def _set(data: dict[str, Any], path: tuple[str, ...], value: Any) -> None:
    current = data
    for key in path[:-1]:
        current = current.setdefault(key, {})
    current[path[-1]] = value


def _canonical_terminal_patch(value: Any, *, terminus: str, legacy_cap: bool) -> str | None:
    if value is None:
        return None
    text = str(value).strip().upper()
    if text in {"", "NONE", "NULL", "FALSE"}:
        return "none"
    if legacy_cap and terminus == "C" and text == "NME":
        return "CT3"
    return text


def _normalize_terminal_caps(normalized: dict[str, Any], explicit: dict[str, Any], errors: list[str]) -> None:
    stage4 = normalized.setdefault("stage4", {})
    if not isinstance(stage4, dict):
        return
    all_atom = normalized.get("all_atom")
    public_all_atom = all_atom if isinstance(all_atom, dict) else {}
    for terminus, patch_key, cap_key, default in (
        ("N", "n_terminal_patch", "n_terminal_cap", "ACE"),
        ("C", "c_terminal_patch", "c_terminal_cap", "CT3"),
    ):
        if patch_key not in stage4 and patch_key in public_all_atom:
            stage4[patch_key] = public_all_atom[patch_key]
        if cap_key not in stage4 and cap_key in public_all_atom:
            stage4[cap_key] = public_all_atom[cap_key]
        patch = _canonical_terminal_patch(stage4.get(patch_key, default), terminus=terminus, legacy_cap=False)
        cap = _canonical_terminal_patch(stage4.get(cap_key), terminus=terminus, legacy_cap=True)
        if patch is not None and cap is not None and patch != cap:
            errors.append(f"conflicting terminal configuration: stage4.{patch_key} != stage4.{cap_key}")
        selected = cap if patch is None else patch
        if selected is None:
            selected = default
        allowed = {"ACE", "none"} if terminus == "N" else {"CT3", "none"}
        if selected not in allowed:
            errors.append(
                f"stage4.{patch_key} must be {'ACE or none' if terminus == 'N' else 'CT3 or none'}; "
                f"use stage4.{cap_key}: NME only as legacy input for CHARMM CT3"
            )
        stage4[patch_key] = selected
        if cap_key in stage4:
            del stage4[cap_key]
    if isinstance(all_atom, dict):
        all_atom["n_terminal_patch"] = stage4.get("n_terminal_patch", "ACE")
        all_atom["c_terminal_patch"] = stage4.get("c_terminal_patch", "CT3")
        all_atom.pop("n_terminal_cap", None)
        all_atom.pop("c_terminal_cap", None)


def _validate_unknown_keys(data: dict[str, Any], errors: list[str], path: tuple[str, ...] = ()) -> None:
    allowed = ALLOWED_KEYS.get(path)
    if allowed is not None:
        for key in data:
            if key not in allowed:
                errors.append(f"unknown configuration key {'.'.join((*path, key))}")
    for key, value in data.items():
        child_path = (*path, key)
        if isinstance(value, dict) and (child_path in ALLOWED_KEYS or len(child_path) <= 3):
            _validate_unknown_keys(value, errors, child_path)


def normalize_config(config: dict[str, Any], explicit_config: dict[str, Any] | None = None) -> dict[str, Any]:
    normalized = deepcopy(config)
    explicit = explicit_config or {}
    errors: list[str] = []
    warnings: list[str] = []
    _validate_unknown_keys(explicit, errors)
    for public_path, legacy_paths in PUBLIC_TO_LEGACY.items():
        if not _has(normalized, public_path):
            continue
        public_value = _get(normalized, public_path)
        for legacy_path in legacy_paths:
            if _has(explicit, legacy_path) and _get(explicit, legacy_path) != public_value:
                errors.append(
                    "conflicting duplicate configuration values: "
                    f"{'.'.join(public_path)} != {'.'.join(legacy_path)}"
                )
            _set(normalized, legacy_path, public_value)
    for section in ("stage1", "stage2", "stage3", "stage4"):
        if section in explicit:
            warnings.append(f"{section}.* is legacy compatibility configuration; prefer public sections when available")
    simulation = normalized.get("simulation", {})
    if isinstance(simulation, dict) and simulation.get("run_cg_production"):
        errors.append("simulation.run_cg_production is not implemented; set it false")
    _normalize_terminal_caps(normalized, explicit, errors)
    stage4 = normalized.get("stage4", {})
    if isinstance(stage4, dict):
        if stage4.get("parameter_generation_enabled"):
            errors.append("stage4.parameter_generation_enabled is not implemented in this workflow")
        if stage4.get("pyrosetta_plugins_enabled"):
            errors.append("stage4.pyrosetta_plugins_enabled is not implemented in this workflow")
    if errors:
        raise ConfigError("; ".join(errors))
    normalized.setdefault("_meta", {})["normalized_config_version"] = 1
    normalized["_meta"]["normalization_warnings"] = warnings
    return normalized
