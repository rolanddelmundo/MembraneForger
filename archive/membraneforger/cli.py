#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

try:
    import yaml
except Exception:  # pragma: no cover - diagnosed by doctor
    yaml = None

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from membraneforger_paths import PathResolutionError, display_path, gmx_command, resolve_cli_path, resolve_mstool  # noqa: E402
from external_dependencies import DependencyError, preflight_stage4_charmm_resources, resolve_stage4_resource_paths  # noqa: E402
import membraneforger_workflow as legacy_workflow  # noqa: E402

STAGES = ("stage1", "stage2", "stage3", "stage4")
RUN_MARKER = ".membraneforger-run.json"
STAGE_LABELS = {
    "stage1": "stage1_cg_setup",
    "stage2": "stage2_cg_simulation",
    "stage3": "stage3_backmapping",
    "stage4": "stage4_all_atom",
}
TERMINAL_SUCCESS = {"COMPLETE", "REVIEW_REQUIRED"}
ATOM_RECORDS = ("ATOM  ", "HETATM")
LOCAL_EXTERNAL_RESOURCE_PREFIXES = (
    "resources/external/",
    "resources/forcefields/charmm36/toppar/",
)


class CliError(RuntimeError):
    pass


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_yaml(path: Path) -> dict[str, Any]:
    if yaml is None:
        raise CliError("PyYAML is required. Fix: ./setup.sh --auto")
    data = yaml.safe_load(path.read_text(encoding="utf-8")) if path.is_file() else {}
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise CliError(f"configuration root must be a mapping: {path}")
    return data


def write_yaml(path: Path, data: dict[str, Any], header: str | None = None) -> None:
    if yaml is None:
        raise CliError("PyYAML is required. Fix: ./setup.sh --auto")
    path.parent.mkdir(parents=True, exist_ok=True)
    body = yaml.safe_dump(data, sort_keys=False)
    prefix = f"{header.rstrip()}\n" if header else ""
    path.write_text(prefix + body, encoding="utf-8")


def run_id_from_path(path: Path) -> str:
    text = path.name.strip() or "run"
    run_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    return run_id.strip("._-") or "run"


def read_status(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "status.json"
    if path.is_file():
        return json.loads(path.read_text(encoding="utf-8"))
    stages = {stage: {"state": "NOT_STARTED"} for stage in STAGES}
    return {"schema_version": 1, "run_dir": str(run_dir), "stages": stages, "created_at": now(), "updated_at": now()}


def write_status(run_dir: Path, status: dict[str, Any]) -> None:
    status["updated_at"] = now()
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "status.json").write_text(json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def update_stage(run_dir: Path, stage: str, state: str, **extra: Any) -> None:
    status = read_status(run_dir)
    item = status.setdefault("stages", {}).setdefault(stage, {})
    item.update(extra)
    item["state"] = state
    item["updated_at"] = now()
    write_status(run_dir, status)


def repo_commit() -> str:
    result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    return result.stdout.strip() if result.returncode == 0 else "not-a-git-checkout"


def pdb_report(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise CliError(f"input PDB is not readable: {path}")
    atoms = het = altloc = waters = 0
    chains: set[str] = set()
    residues: set[tuple[str, str, str]] = set()
    missing_backbone: list[str] = []
    by_residue: dict[tuple[str, str, str], set[str]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(ATOM_RECORDS):
                continue
            if line.startswith("ATOM  "):
                atoms += 1
            else:
                het += 1
            chain = line[21:22].strip() or "_"
            resname = line[17:20].strip()
            resid = line[22:26].strip()
            atom = line[12:16].strip()
            if line[16:17].strip():
                altloc += 1
            if resname in {"HOH", "WAT", "TIP3", "SOL"}:
                waters += 1
            chains.add(chain)
            key = (chain, resid, resname)
            residues.add(key)
            by_residue.setdefault(key, set()).add(atom)
    for (chain, resid, resname), names in sorted(by_residue.items())[:10000]:
        if resname not in {"HOH", "WAT", "TIP3", "SOL"} and {"N", "CA", "C", "O"} - names:
            missing_backbone.append(f"{chain}:{resname}{resid}")
    return {
        "path": str(path),
        "sha256": file_sha256(path),
        "atom_records": atoms,
        "hetatm_records": het,
        "chains": sorted(chains),
        "residue_count": len(residues),
        "alternate_location_records": altloc,
        "waters": waters,
        "missing_backbone_residues": missing_backbone[:50],
    }


def generated_config(pdb: Path, run_dir: Path, run_id: str) -> dict[str, Any]:
    staged = run_dir / "input" / "original.pdb"
    return {
        "input": {"pdb": str(staged), "orientation": "user_review_required"},
        "membrane": {
            "build_mode": "de_novo_insane",
            "composition": {"upper": {}, "lower": {}},
            "box": {"shape": "rectangular", "x_nm": None, "y_nm": None, "z_nm": None},
            "salt_concentration_molar": 0.15,
        },
        "coarse_grained": {"martini_forcefield": "martini3001", "protein_mode": "whole"},
        "simulation": {"run_cg_minimization": True, "run_cg_equilibration": True, "run_cg_production": False},
        "backmapping": {"enabled": True, "mstool_root": "resources/vendor/mstool"},
        "all_atom": {
            "enabled": True,
            "charmm36_forcefield": "resources/external/charmm36.ff",
            "membrane_toppar": "resources/external/charmm36_membrane",
            "ligand_params": [],
            "salt_concentration_molar": 0.15,
            "n_terminal_patch": "ACE",
            "c_terminal_patch": "CT3",
        },
        "compute": {"threads": "$MEMBRANEFORGER_LOCAL_THREADS", "gmx": "$GMX_BIN"},
        "review": {"checkpoints": ["stage2"], "require_input_review_on_warnings": True},
        "advanced": {"legacy_scaffold_replace": False},
        "run": {"id": run_id, "overwrite": False, "update_latest": False, "require_inputs_for_check": False},
        "stage1": {
            "input_pdb": str(staged),
            "orientation_enabled": False,
            "orientation_chain": None,
            "membrane_protein_mode": "whole",
            "martinize_forcefield": "martini3001",
        },
        "stage2": {"mode": "production", "production": False, "build_mode": "de_novo_insane"},
        "stage3": {
            "input_aa_protlig": str(staged),
            "input_opm_reference": None,
            "resources": {"membrane_templates": "resources/templates", "membrane_toppar": "resources/external/charmm36_membrane"},
            "placement": {
                "mode": "aa_reference_to_backmapped",
                "use_opm": False,
                "fixed_selection": "protein",
                "moving_selection": "protein",
                "fixed_chain": None,
                "moving_chain": None,
            },
            "plugins": {"glycolipid_template_finalize": {"enabled": False}},
        },
        "stage4": {
            "resources": {
                "aa_mdp": "resources/aa_mdp",
                "charmm36_forcefield": "resources/external/charmm36.ff",
                "membrane_toppar": "resources/external/charmm36_membrane",
            },
            "ligand_params": [],
            "dssp_required": False,
            "charmm36_required": True,
            "cgenff_required": False,
            "ligand_params_required": False,
            "parameter_generation_enabled": False,
            "pyrosetta_plugins_enabled": False,
            "n_terminal_patch": "ACE",
            "c_terminal_patch": "CT3",
            "slice_buffer_nm": 1.0,
            "z_slab_padding_nm": 1.5,
            "salt_concentration_molar": 0.15,
        },
    }


def cmd_init(args: argparse.Namespace) -> int:
    pdb = resolve_cli_path(args.pdb)
    run_dir = resolve_cli_path(args.output_dir)
    if run_dir.exists() and any(run_dir.iterdir()) and not args.overwrite:
        raise CliError(f"output directory already exists and is not empty: {run_dir}; use --overwrite")
    run_id = args.run_id or run_id_from_path(run_dir)
    staged = run_dir / "input" / "original.pdb"
    staged.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(pdb, staged)
    for name in ("config", "logs", "provenance", "final/cg", "final/all_atom"):
        (run_dir / name).mkdir(parents=True, exist_ok=True)
    report = pdb_report(staged)
    config = generated_config(pdb, run_dir, run_id)
    write_yaml(
        run_dir / "config" / "config.yaml",
        config,
        "# MembraneForger canonical run config. Fill membrane.box/composition and licensed AA paths before a full run.",
    )
    provenance = {
        "schema_version": 1,
        "membraneforger_commit": repo_commit(),
        "command": " ".join(sys.argv),
        "created_at": now(),
        "input": report,
    }
    (run_dir / "provenance" / "run.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    status = read_status(run_dir)
    status.update({"run_id": run_id, "config": str(run_dir / "config" / "config.yaml")})
    status["stages"]["stage1"]["state"] = "READY"
    write_status(run_dir, status)
    (run_dir / RUN_MARKER).write_text(
        json.dumps({"schema_version": 1, "run_id": run_id, "created_by": "MembraneForger"}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"PASS init run_dir={run_dir}")
    print(f"INFO config={run_dir / 'config' / 'config.yaml'}")
    print("INFO next: ./run_pipeline.sh doctor --config " + str(run_dir / "config" / "config.yaml"))
    return 0


def config_path_from_args(args: argparse.Namespace) -> Path:
    if getattr(args, "config", None):
        return resolve_cli_path(args.config)
    if not getattr(args, "run_dir", None):
        raise CliError("either --config or --run-dir is required")
    run_dir = resolve_cli_path(args.run_dir)
    return run_dir / "config" / "config.yaml"


def run_dir_from_config(config_path: Path, explicit: str | None = None) -> Path:
    if explicit:
        return resolve_cli_path(explicit)
    if config_path.name == "config.yaml" and config_path.parent.name == "config":
        return config_path.parent.parent.resolve()
    return (ROOT / "runs" / run_id_from_path(config_path.parent)).resolve()


def stage_range(start: str | None, through: str | None, only: str | None = None) -> list[str]:
    if only:
        if only not in STAGES:
            raise CliError(f"unknown stage: {only}")
        return [only]
    start_i = STAGES.index(start) if start else 0
    end_i = STAGES.index(through) if through else len(STAGES) - 1
    if start_i > end_i:
        raise CliError("--from stage must not be after --through stage")
    return list(STAGES[start_i : end_i + 1])


def legacy_config(config_path: Path) -> tuple[dict[str, Any], str]:
    config = legacy_workflow.load_config(config_path)
    run_id = legacy_workflow.run_id_from(config, None)
    return config, run_id


def ensure_ready_for_stage(run_dir: Path, stage: str, accept_review: bool = False) -> bool:
    status = read_status(run_dir)
    state = status.get("stages", {}).get(stage, {}).get("state", "NOT_STARTED")
    if state == "REVIEW_REQUIRED" and not accept_review:
        print(f"BLOCKED {stage} is waiting for researcher review")
        print(f"Fix: inspect {run_dir / STAGE_LABELS[stage]} then rerun with --accept-review")
        return False
    return True


def configured_review_checkpoints(config_path: Path) -> set[str]:
    try:
        config = load_yaml(config_path)
    except Exception:
        return set()
    review = config.get("review", {})
    checkpoints = review.get("checkpoints", []) if isinstance(review, dict) else []
    return {str(item) for item in checkpoints or []}


def link_or_copy(src: Path, dst: Path) -> None:
    if not src.exists():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        if dst.is_dir() and not dst.is_symlink():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    try:
        dst.symlink_to(src.resolve(), target_is_directory=src.is_dir())
    except OSError:
        if src.is_dir():
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)


def sync_outputs(run_dir: Path, run_id: str, stage: str) -> None:
    src = ROOT / "outputs" / run_id / stage
    link_or_copy(src, run_dir / STAGE_LABELS[stage])
    final = run_dir / "final"
    if stage == "stage2":
        link_or_copy(src / "cg_backmap_input.gro", final / "cg" / "backmap_input.gro")
        link_or_copy(src / "cg_topology.top", final / "cg" / "topology.top")
        link_or_copy(src / "cg_equilibrated_scaffold.gro", final / "cg" / "structure.gro")
    if stage == "stage3":
        link_or_copy(src / "final_all_atom.pdb", final / "all_atom" / "backmapped.pdb")
        link_or_copy(src / "aa_backmapped.dms", final / "all_atom" / "backmapped.dms")
    if stage == "stage4":
        link_or_copy(src / "minimized_all_atom.gro", final / "all_atom" / "minimized_all_atom.gro")
        link_or_copy(src / "minimized_all_atom.top", final / "all_atom" / "topology.top")


def run_stage(config_path: Path, run_dir: Path, stage: str, accept_review: bool = False, mode: str | None = None) -> int:
    if not ensure_ready_for_stage(run_dir, stage, accept_review):
        return 3
    config, run_id = legacy_config(config_path)
    update_stage(run_dir, stage, "RUNNING", started_at=now())
    argv = ["--config", str(config_path), "--run-id", run_id]
    if mode:
        argv += ["--mode", mode]
    rc = legacy_workflow.stage_main(stage, argv)
    if rc == 0:
        sync_outputs(run_dir, run_id, stage)
        checkpoints = configured_review_checkpoints(config_path)
        state = "REVIEW_REQUIRED" if stage in checkpoints and not accept_review else "COMPLETE"
        update_stage(run_dir, stage, state, completed_at=now(), root_output=str(ROOT / "outputs" / run_id / stage))
        next_index = STAGES.index(stage) + 1
        if next_index < len(STAGES):
            update_stage(run_dir, STAGES[next_index], "READY")
        if state == "REVIEW_REQUIRED":
            print(f"REVIEW_REQUIRED {stage}: inspect {run_dir / STAGE_LABELS[stage]}")
            print(f"Continue: ./run_pipeline.sh resume --run-dir {run_dir} --accept-review")
    else:
        update_stage(run_dir, stage, "FAILED", failed_at=now(), exit_code=rc)
    return rc


def cmd_run(args: argparse.Namespace) -> int:
    config_path = config_path_from_args(args)
    run_dir = run_dir_from_config(config_path, getattr(args, "run_dir", None))
    selected = stage_range(args.from_stage, args.through, None)
    for stage in selected:
        rc = run_stage(config_path, run_dir, stage, accept_review=args.accept_review)
        if rc != 0:
            return rc
        state = read_status(run_dir)["stages"][stage]["state"]
        if state == "REVIEW_REQUIRED":
            return 3
    return 0


def cmd_stage(args: argparse.Namespace) -> int:
    config_path = config_path_from_args(args)
    run_dir = run_dir_from_config(config_path, args.run_dir)
    return run_stage(config_path, run_dir, args.stage, accept_review=args.accept_review, mode=args.mode)


def cmd_resume(args: argparse.Namespace) -> int:
    run_dir = resolve_cli_path(args.run_dir)
    status = read_status(run_dir)
    config_path = Path(status.get("config") or run_dir / "config" / "config.yaml")
    for stage in STAGES:
        state = status.get("stages", {}).get(stage, {}).get("state", "NOT_STARTED")
        if state == "REVIEW_REQUIRED" and not args.accept_review:
            print(f"BLOCKED {stage} review required")
            print(f"Continue: ./run_pipeline.sh resume --run-dir {run_dir} --accept-review")
            return 3
        if state not in TERMINAL_SUCCESS:
            return cmd_run(argparse.Namespace(config=str(config_path), run_dir=str(run_dir), from_stage=stage, through=args.through, accept_review=args.accept_review))
    print("PASS all stages complete")
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    run_dir = resolve_cli_path(args.run_dir)
    status = read_status(run_dir)
    if args.json:
        print(json.dumps(status, indent=2, sort_keys=True))
    else:
        for stage in STAGES:
            print(f"{stage:<6} {status.get('stages', {}).get(stage, {}).get('state', 'NOT_STARTED')}")
    return 0


def python_module_status(module: str, required: bool) -> dict[str, Any]:
    try:
        importlib.import_module(module)
        try:
            version = importlib.metadata.version(module)
        except Exception:
            version = "import-ok"
        return {"name": module, "status": "PASS", "detail": version, "required": required}
    except Exception as exc:
        return {"name": module, "status": "FAIL" if required else "SKIP", "detail": str(exc), "required": required}


def executable_status(name: str, command: list[str], required: bool, fix: str) -> dict[str, Any]:
    exe = shutil.which(command[0]) if len(command) == 1 or "/" not in command[0] else command[0]
    if not exe:
        return {"name": name, "status": "FAIL" if required else "SKIP", "detail": "not found", "fix": fix, "required": required}
    try:
        result = subprocess.run([exe, *(command[1:] or ["--version"])], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=10)
        first = (result.stdout or "").splitlines()[0] if (result.stdout or "").splitlines() else f"exit {result.returncode}"
        ok = result.returncode == 0
        return {"name": name, "status": "PASS" if ok else ("FAIL" if required else "WARN"), "detail": f"{exe} {first}", "fix": fix, "required": required}
    except Exception as exc:
        return {"name": name, "status": "FAIL" if required else "WARN", "detail": str(exc), "fix": fix, "required": required}


def resource_status(path: Path, name: str, required: bool, fix: str) -> dict[str, Any]:
    ok = path.exists()
    return {
        "name": name,
        "status": "PASS" if ok else ("FAIL" if required else "SKIP"),
        "detail": str(path),
        "fix": "" if ok else fix,
        "required": required,
    }


def stage4_resource_status(path: Path, required_path: Path, name: str, required: bool, source: str, fix: str) -> dict[str, Any]:
    ok = required_path.exists()
    status = "PASS" if ok else ("BLOCKED" if required else "SKIP")
    detail = f"{path} via {source}"
    if not ok and required_path != path:
        detail = f"missing {required_path}; configured root {path} via {source}"
    return {
        "name": name,
        "status": status,
        "detail": detail,
        "fix": "" if ok else fix,
        "required": required,
    }


def manifest_status() -> dict[str, Any]:
    manifest = ROOT / "resources" / "RESOURCE_MANIFEST.tsv"
    if not manifest.is_file():
        return {
            "name": "Resource manifest",
            "status": "FAIL",
            "detail": "missing resources/RESOURCE_MANIFEST.tsv",
            "fix": "restore resources/RESOURCE_MANIFEST.tsv",
            "required": True,
        }
    failures: list[str] = []
    listed: set[str] = set()
    for line in manifest.read_text(encoding="utf-8").splitlines()[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 7:
            failures.append("malformed manifest row")
            continue
        relpath, expected = parts[0], parts[5]
        listed.add(relpath)
        path = ROOT / relpath
        if not path.is_file():
            failures.append(f"missing {relpath}")
        elif file_sha256(path) != expected:
            failures.append(f"checksum {relpath}")
    for path in (ROOT / "resources").rglob("*"):
        if path.name == ".DS_Store":
            continue
        if path.is_file() and path.name != "RESOURCE_MANIFEST.tsv":
            if path.relative_to(ROOT).as_posix().startswith("resources/vendor/mstool/"):
                continue
            if path.relative_to(ROOT).as_posix().startswith(LOCAL_EXTERNAL_RESOURCE_PREFIXES) and path.name != "README.md":
                continue
            rel = path.relative_to(ROOT).as_posix()
            if rel not in listed:
                failures.append(f"unlisted {rel}")
    if failures:
        return {
            "name": "Resource manifest",
            "status": "FAIL",
            "detail": "; ".join(failures[:5]),
            "fix": "restore the resource or regenerate the manifest after auditing provenance",
            "required": True,
        }
    return {"name": "Resource manifest", "status": "PASS", "detail": "checksums verified", "required": True}


def doctor_rows(config_path: Path | None) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    rows.append({"name": "Python", "status": "PASS", "detail": sys.version.split()[0], "required": True})
    for module in ("yaml", "numpy", "Cython", "openmm"):
        rows.append(python_module_status(module, required=module in {"yaml", "numpy", "Cython"}))
    rows.append(python_module_status("vermouth", required=True))
    rows.append(executable_status("martinize2", [os.environ.get("MARTINIZE2_BIN", "martinize2"), "--version"], True, "conda env update -f environments/environment.yml"))
    gmx = gmx_command()
    rows.append(executable_status("GROMACS", gmx, True, "install GROMACS or set GMX_BIN"))
    insane = ROOT / "scripts" / "insane_M3_lipids_new.py"
    if insane.is_file():
        rows.append({
            "name": "INSANE bundled exact script",
            "status": "PASS",
            "detail": f"{display_path(insane, ROOT)} sha256={file_sha256(insane)}",
            "required": True,
        })
    else:
        rows.append({
            "name": "INSANE bundled exact script",
            "status": "FAIL",
            "detail": "missing scripts/insane_M3_lipids_new.py",
            "fix": "restore scripts/insane_M3_lipids_new.py",
            "required": True,
        })
    try:
        resolution = resolve_mstool(root=ROOT, import_module=False)
        rows.append({"name": "mstool", "status": "PASS", "detail": resolution.provenance(ROOT), "required": True})
    except Exception as exc:
        rows.append({"name": "mstool", "status": "FAIL", "detail": str(exc), "fix": "./setup.sh --bootstrap mstool", "required": True})
    rows.append(manifest_status())
    rows.append(resource_status(ROOT / "resources" / "forcefields" / "martini" / "martini_v3.0.0.itp", "Martini 3 FF", True, "./setup.sh --bootstrap martini"))
    rows.append(resource_status(ROOT / "resources" / "aa_mdp" / "minim.mdp", "AA minimization MDP", True, "restore resources/aa_mdp/minim.mdp"))
    if config_path and config_path.is_file():
        try:
            config = legacy_workflow.load_config(config_path)
            aa_enabled = bool(config.get("all_atom", {}).get("enabled", True))
            stage4 = config.get("stage4", {})
            charmm_required = aa_enabled and bool(stage4.get("charmm36_required", False))
            membrane_required = aa_enabled and bool(stage4.get("membrane_toppar_required", stage4.get("charmm36_required", False)))
            paths = resolve_stage4_resource_paths(config, root=ROOT)
            rows.append(stage4_resource_status(
                paths.charmm36_forcefield,
                paths.charmm36_forcefield / "forcefield.itp",
                "CHARMM36 GROMACS force field",
                charmm_required,
                paths.charmm36_source,
                "Stage 4 needs a CHARMM36 molecular-parameter package formatted for GROMACS. "
                "MembraneForger does not redistribute the local package used during development.\n"
                "       Configure:\n"
                "         MEMBRANEFORGER_CHARMM36_ROOT=/path/to/charmm36.ff",
            ))
            rows.append(stage4_resource_status(
                paths.membrane_toppar,
                paths.membrane_toppar,
                "CHARMM membrane parameter files",
                membrane_required,
                paths.membrane_toppar_source,
                "Stage 4 also requires the compatible CHARMM topology/parameter resource tree.\n"
                "       Configure:\n"
                "         MEMBRANEFORGER_MEMBRANE_TOPPAR=/path/to/toppar",
            ))
            if charmm_required or membrane_required:
                preflight = preflight_stage4_charmm_resources(config, root=ROOT, required=True)
                if preflight.errors:
                    rows.append({
                        "name": "Stage 4 CHARMM resource preflight",
                        "status": "BLOCKED",
                        "detail": " | ".join(preflight.errors[:3]),
                        "fix": (
                            "Provide paired, compatible CHARMM36 GROMACS and toppar resources. "
                            "This blocks Stage 4 for the requested configuration; generic public checks can still pass."
                        ),
                        "required": True,
                    })
                else:
                    metadata = ", ".join(f"{key}={value}" for key, value in sorted(preflight.metadata.items())) or "metadata not declared by resources"
                    rows.append({
                        "name": "Stage 4 CHARMM resource preflight",
                        "status": "PASS",
                        "detail": metadata,
                        "required": True,
                    })
            plugin = config.get("stage3", {}).get("plugins", {}).get("glycolipid_template_finalize", {})
            if plugin.get("enabled"):
                target = Path(os.path.expandvars(os.path.expanduser(str(plugin.get("target_itp", "")))))
                rows.append(resource_status(target, "glycolipid finalizer target ITP", True, "set stage3.plugins.glycolipid_template_finalize.target_itp"))
        except legacy_workflow.ContractError as exc:
            rows.append({"name": "configuration", "status": "FAIL", "detail": str(exc), "required": True})
        except DependencyError as exc:
            rows.append({"name": "Stage 4 external resources", "status": "BLOCKED", "detail": str(exc), "required": True})
    else:
        rows.append({
            "name": "Stage 4 external resources",
            "status": "INFO",
            "detail": "generic doctor checks public dependencies only; use --config to check CHARMM36/toppar readiness",
            "required": False,
        })
    return rows


def cmd_doctor(args: argparse.Namespace) -> int:
    config_path = resolve_cli_path(args.config) if args.config else None
    rows = doctor_rows(config_path)
    if args.json:
        print(json.dumps(rows, indent=2, sort_keys=True))
    else:
        for row in rows:
            status = row["status"]
            print(f"[{status:<7}] {row['name']:<30} {row.get('detail', '')}")
            if status in {"FAIL", "BLOCKED"} and row.get("fix"):
                print(f"       Fix: {row['fix']}")
    return 1 if any(row["status"] in {"FAIL", "BLOCKED"} and row.get("required") for row in rows) else 0


def gromacs_version_report() -> dict[str, Any]:
    try:
        command = gmx_command()
    except PathResolutionError as exc:
        return {"command": [], "version": f"unavailable: {exc}", "returncode": 1}
    try:
        result = subprocess.run([*command, "--version"], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=15)
        first = (result.stdout or "").splitlines()[0] if (result.stdout or "").splitlines() else f"exit {result.returncode}"
        return {"command": command, "version": first, "returncode": result.returncode}
    except Exception as exc:
        return {"command": command, "version": f"unavailable: {exc}", "returncode": 1}


def cmd_scientific_integration(args: argparse.Namespace) -> int:
    config_path = resolve_cli_path(args.config)
    config, run_id = legacy_config(config_path)
    run_dir = run_dir_from_config(config_path, args.run_dir)
    provenance_dir = run_dir / "provenance"
    provenance_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "schema_version": 1,
        "created_at": now(),
        "command": " ".join(sys.argv),
        "membraneforger_commit": repo_commit(),
        "config": str(config_path),
        "run_id": run_id,
        "gromacs": gromacs_version_report(),
        "stage4_resource_preflight": preflight_stage4_charmm_resources(
            config,
            root=ROOT,
            required=True,
            include_file_hashes=True,
        ).as_dict(),
        "stage4_construction": "NOT RUN",
        "full_scientific_integration": "NOT RUN",
    }
    (provenance_dir / "stage4_resource_preflight.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    preflight = report["stage4_resource_preflight"]
    if preflight["errors"]:
        print("BLOCKED Stage 4 resource preflight")
        for item in preflight["errors"]:
            print(f"  {item}")
        print(f"INFO wrote provenance: {provenance_dir / 'stage4_resource_preflight.json'}")
        return 1
    print("PASS Stage 4 resource preflight")
    print(f"INFO wrote provenance: {provenance_dir / 'stage4_resource_preflight.json'}")
    if args.preflight_only:
        return 0
    print("INFO running requested stage-gated scientific workflow")
    rc = cmd_run(argparse.Namespace(
        config=str(config_path),
        run_dir=str(run_dir),
        from_stage=args.from_stage,
        through=args.through,
        accept_review=args.accept_review,
    ))
    status = read_status(run_dir)
    if status.get("stages", {}).get("stage4", {}).get("state") in TERMINAL_SUCCESS:
        print("PASS Stage 4 construction")
    if rc == 0:
        print("PASS full scientific integration workflow completed")
    else:
        print("BLOCKED full scientific integration workflow did not complete")
    return rc


def cmd_validate(args: argparse.Namespace) -> int:
    config_path = resolve_cli_path(args.config)
    config, run_id = legacy_config(config_path)
    return legacy_workflow.check_contract(config, run_id, STAGES)


def cmd_dry_run(args: argparse.Namespace) -> int:
    config_path = resolve_cli_path(args.config)
    config, run_id = legacy_config(config_path)
    rc = legacy_workflow.check_contract(config, run_id, STAGES)
    if rc:
        return rc
    run_dir = run_dir_from_config(config_path, getattr(args, "run_dir", None))
    print(f"DRY RUN run_dir={run_dir}")
    return legacy_workflow.dry_run(config, run_id, STAGES)


def cmd_clean(args: argparse.Namespace) -> int:
    run_dir = resolve_cli_path(args.run_dir)
    root_real = ROOT.resolve()
    runs_root = (ROOT / "runs").resolve()
    run_real = run_dir.resolve(strict=False)
    if run_dir.is_symlink():
        raise CliError(f"refusing to clean symlink path: {run_dir}")
    if run_real in {root_real, runs_root} or root_real in run_real.parents and run_real.parent != runs_root:
        raise CliError(f"refusing to clean non-run directory: {run_dir}")
    if run_real.parent != runs_root:
        raise CliError(f"refusing to clean outside runs root {runs_root}: {run_dir}")
    marker = run_real / RUN_MARKER
    if args.confirm:
        if not marker.is_file():
            raise CliError(f"refusing to clean run directory without {RUN_MARKER}: {run_dir}")
        try:
            marker_data = json.loads(marker.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise CliError(f"invalid run marker {marker}: {exc}") from exc
        if marker_data.get("created_by") != "MembraneForger" or not marker_data.get("run_id"):
            raise CliError(f"refusing to clean unrecognized run marker: {marker}")
    targets = [run_dir]
    for target in targets:
        print(f"REMOVE {target}")
    if args.confirm:
        for target in targets:
            if target.exists():
                shutil.rmtree(target)
        print("PASS clean")
    else:
        print("DRY RUN only; rerun with --confirm")
    return 0


def legacy_passthrough(argv: list[str]) -> int:
    return legacy_workflow.main(argv)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="run_pipeline.sh", description="MembraneForger pipeline controller")
    sub = parser.add_subparsers(dest="command")
    p = sub.add_parser("doctor", help="check dependencies and resources")
    p.add_argument("--config")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_doctor)
    p = sub.add_parser("scientific-integration", help="local maintainer preflight and full scientific workflow")
    p.add_argument("--config", default=str(ROOT / "config" / "workflow.yaml"))
    p.add_argument("--run-dir")
    p.add_argument("--from", dest="from_stage", choices=STAGES)
    p.add_argument("--through", choices=STAGES, default="stage4")
    p.add_argument("--accept-review", action="store_true")
    p.add_argument("--preflight-only", action="store_true")
    p.set_defaults(func=cmd_scientific_integration)
    p = sub.add_parser("init", help="create a run directory and commented config")
    p.add_argument("--pdb", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--run-id")
    p.add_argument("--overwrite", action="store_true")
    p.set_defaults(func=cmd_init)
    p = sub.add_parser("validate", help="validate config and stage contracts")
    p.add_argument("--config", required=True)
    p.set_defaults(func=cmd_validate)
    p = sub.add_parser("dry-run", help="show planned stage commands")
    p.add_argument("--config", required=True)
    p.add_argument("--run-dir")
    p.set_defaults(func=cmd_dry_run)
    p = sub.add_parser("run", help="run stages")
    p.add_argument("--config", required=True)
    p.add_argument("--run-dir")
    p.add_argument("--from", dest="from_stage", choices=STAGES)
    p.add_argument("--through", choices=STAGES)
    p.add_argument("--accept-review", action="store_true")
    p.set_defaults(func=cmd_run)
    p = sub.add_parser("resume", help="continue a run from status.json")
    p.add_argument("--run-dir", required=True)
    p.add_argument("--through", choices=STAGES)
    p.add_argument("--accept-review", action="store_true")
    p.set_defaults(func=cmd_resume)
    p = sub.add_parser("stage", help="run one stage")
    p.add_argument("stage", choices=STAGES)
    p.add_argument("--config")
    p.add_argument("--run-dir", required=True)
    p.add_argument("--mode")
    p.add_argument("--accept-review", action="store_true")
    p.set_defaults(func=cmd_stage)
    p = sub.add_parser("status", help="show run status")
    p.add_argument("--run-dir", required=True)
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_status)
    p = sub.add_parser("clean", help="safely remove a run directory")
    p.add_argument("--run-dir", required=True)
    p.add_argument("--confirm", action="store_true")
    p.set_defaults(func=cmd_clean)
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] not in {"doctor", "scientific-integration", "init", "validate", "dry-run", "run", "resume", "stage", "status", "clean", "-h", "--help"}:
        return legacy_passthrough(argv)
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    try:
        return int(args.func(args))
    except (CliError, PathResolutionError, legacy_workflow.ContractError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
