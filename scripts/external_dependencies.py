from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from membraneforger_paths import REPO_ROOT, resolve_config_path


class DependencyError(RuntimeError):
    pass


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _resolve_executable(value: str | None, default_name: str | None = None) -> Path | None:
    if value:
        p = Path(os.path.expanduser(os.path.expandvars(value)))
        if p.is_absolute() or "/" in value:
            return p.resolve()
        found = shutil.which(value)
        return Path(found).resolve() if found else None
    if default_name:
        found = shutil.which(default_name)
        return Path(found).resolve() if found else None
    return None


def _version(cmd: list[str], timeout: int = 10) -> str:
    try:
        result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout)
    except Exception as exc:
        return f"version unavailable: {exc}"
    return (result.stdout or f"exit {result.returncode}").splitlines()[0]


@dataclass(frozen=True)
class DependencyRecord:
    name: str
    status: str
    path: str
    version: str = ""
    sha256: str = ""
    detail: str = ""

    def as_dict(self) -> dict[str, str]:
        return {
            "name": self.name,
            "status": self.status,
            "path": self.path,
            "version": self.version,
            "sha256": self.sha256,
            "detail": self.detail,
        }


@dataclass(frozen=True)
class Stage4ResourcePaths:
    charmm36_forcefield: Path
    membrane_toppar: Path
    aa_mdp: Path
    charmm36_source: str
    membrane_toppar_source: str
    aa_mdp_source: str


@dataclass(frozen=True)
class Stage4ResourcePreflight:
    paths: Stage4ResourcePaths
    status: str
    errors: tuple[str, ...]
    metadata: dict[str, str]
    tree_hashes: dict[str, Any]

    def as_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "errors": list(self.errors),
            "paths": {
                "charmm36_forcefield": str(self.paths.charmm36_forcefield),
                "membrane_toppar": str(self.paths.membrane_toppar),
                "aa_mdp": str(self.paths.aa_mdp),
                "charmm36_source": self.paths.charmm36_source,
                "membrane_toppar_source": self.paths.membrane_toppar_source,
                "aa_mdp_source": self.paths.aa_mdp_source,
            },
            "metadata": self.metadata,
            "tree_hashes": self.tree_hashes,
        }


CHARMM36_REQUIRED_FILES = (
    "forcefield.itp",
    "ffbonded.itp",
    "ffnonbonded.itp",
    "merged.rtp",
    "merged.n.tdb",
    "merged.c.tdb",
    "tip3p.itp",
    "ions.itp",
)
TOPPAR_REQUIRED_FILES = (
    "forcefield.itp",
    "top_all36_prot.rtf",
    "par_all36m_prot.prm",
)


def _config_path(config: dict) -> Path | None:
    raw = config.get("_meta", {}).get("config_path") if isinstance(config.get("_meta"), dict) else None
    return Path(raw) if raw else None


def _stage4_resources(config: dict) -> dict:
    stage4 = config.get("stage4", {}) if isinstance(config.get("stage4", {}), dict) else {}
    resources = stage4.get("resources", {}) if isinstance(stage4.get("resources", {}), dict) else {}
    return resources


def resolve_config_or_env_path(
    config: dict,
    key: str,
    env: str,
    default: str,
    *,
    root: Path = REPO_ROOT,
) -> tuple[Path, str]:
    raw_env = os.environ.get(env)
    if raw_env:
        expanded = os.path.expanduser(os.path.expandvars(raw_env))
        if "$" in expanded:
            raise DependencyError(f"{env} contains an unresolved environment variable: {raw_env!r}")
        return Path(expanded).resolve(), env
    value = _stage4_resources(config).get(key, default)
    return resolve_config_path(value, config_path=_config_path(config), root=root), "config"


def resolve_stage4_resource_paths(config: dict, *, root: Path = REPO_ROOT) -> Stage4ResourcePaths:
    charmm36, charmm36_source = resolve_config_or_env_path(
        config,
        "charmm36_forcefield",
        "MEMBRANEFORGER_CHARMM36_ROOT",
        "resources/external/charmm36.ff",
        root=root,
    )
    membrane, membrane_source = resolve_config_or_env_path(
        config,
        "membrane_toppar",
        "MEMBRANEFORGER_MEMBRANE_TOPPAR",
        "resources/external/charmm36_membrane",
        root=root,
    )
    aa_mdp = resolve_config_path(
        _stage4_resources(config).get("aa_mdp", "resources/aa_mdp"),
        config_path=_config_path(config),
        root=root,
    )
    return Stage4ResourcePaths(charmm36, membrane, aa_mdp, charmm36_source, membrane_source, "config")


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace") if path.is_file() else ""


def _first_pattern(text: str, patterns: tuple[str, ...]) -> str:
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return " ".join(group for group in match.groups() if group).strip()
    return ""


def _month_year_tokens(text: str) -> set[str]:
    months = "Jan(?:uary)?|Feb(?:ruary)?|Mar(?:ch)?|Apr(?:il)?|May|Jun(?:e)?|Jul(?:y)?|Aug(?:ust)?|Sep(?:t(?:ember)?)?|Oct(?:ober)?|Nov(?:ember)?|Dec(?:ember)?"
    tokens: set[str] = set()
    for match in re.finditer(rf"\b({months})\.?\s+([12][0-9]{{3}})\b", text, flags=re.IGNORECASE):
        month = match.group(1).rstrip(".").lower()[:3]
        tokens.add(f"{month}-{match.group(2)}")
    return tokens


def _terminal_database_entries(charmm36: Path, terminus: str) -> set[str]:
    filename = "merged.n.tdb" if terminus == "N" else "merged.c.tdb"
    path = charmm36 / filename
    entries: set[str] = set()
    for raw in _read_text(path).splitlines():
        stripped = raw.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            entries.add(stripped[1:-1].strip().upper())
    return entries


def _has_charmm_patch(path: Path, patch: str) -> bool:
    return re.search(rf"^PRES\s+{re.escape(patch)}\b", _read_text(path), flags=re.IGNORECASE | re.MULTILINE) is not None


def _tree_digest(path: Path, *, include_files: bool = False) -> dict[str, Any]:
    if not path.is_dir():
        return {"exists": False, "file_count": 0, "size_bytes": 0, "sha256": ""}
    files = sorted(p for p in path.rglob("*") if p.is_file())
    digest = hashlib.sha256()
    file_rows: list[dict[str, str | int]] = []
    size = 0
    for file_path in files:
        rel = file_path.relative_to(path).as_posix()
        file_hash = _sha256(file_path)
        file_size = file_path.stat().st_size
        size += file_size
        digest.update(rel.encode("utf-8"))
        digest.update(b"\0")
        digest.update(str(file_size).encode("ascii"))
        digest.update(b"\0")
        digest.update(file_hash.encode("ascii"))
        digest.update(b"\0")
        if include_files:
            file_rows.append({"path": rel, "size_bytes": file_size, "sha256": file_hash})
    result: dict[str, Any] = {
        "exists": True,
        "file_count": len(files),
        "size_bytes": size,
        "sha256": digest.hexdigest(),
    }
    if include_files:
        result["files"] = file_rows
    return result


def stage4_resource_metadata(paths: Stage4ResourcePaths) -> dict[str, str]:
    charmm_doc = _read_text(paths.charmm36_forcefield / "forcefield.doc")
    charmm_itp = _read_text(paths.charmm36_forcefield / "forcefield.itp")
    toppar_forcefield = _read_text(paths.membrane_toppar / "forcefield.itp")
    top_prot = _read_text(paths.membrane_toppar / "top_all36_prot.rtf")
    par_prot = _read_text(paths.membrane_toppar / "par_all36m_prot.prm")
    combined = "\n".join([charmm_doc, charmm_itp, toppar_forcefield, top_prot, par_prot])
    metadata = {
        "charmm36_release": _first_pattern(
            charmm_doc,
            (
                r"CHARMM36\s+all-atom\s+force\s+field\s*\(([^)]+)\)",
                r"CHARMM36.*?\b((?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Sept|Oct|Nov|Dec)[a-z]*\.?\s+[12][0-9]{3})\b",
            ),
        ),
        "cgenff_version": _first_pattern(charmm_doc + "\n" + charmm_itp, (r"CGenFF\s+([0-9]+(?:\.[0-9]+)*)",)),
        "protein_forcefield_revision": "CHARMM36m" if re.search(r"CHARMM36m", combined, flags=re.IGNORECASE) else "",
        "gromacs_conversion": "charmm2gmx.py" if "charmm2gmx.py" in charmm_itp else "",
        "toppar_conversion": "CHARMM-GUI FF-Converter" if "CHARMM-GUI FF-Converter" in toppar_forcefield else "",
        "toppar_protein_topology_date": _first_pattern(top_prot, (r">+\s*([^<>\n]*[12][0-9]{3}[^<>\n]*)\s*<+",)),
        "toppar_protein_parameter_date": _first_pattern(par_prot, (r">+\s*([^<>\n]*[12][0-9]{3}[^<>\n]*)\s*<+",)),
        "ljpme_variant": "LJPME" if re.search(r"\bljpme\b", combined, flags=re.IGNORECASE) else "",
    }
    return {key: value for key, value in metadata.items() if value}


def _requested_terminal_patches(config: dict[str, Any]) -> dict[str, str]:
    stage4 = config.get("stage4", {}) if isinstance(config.get("stage4", {}), dict) else {}
    return {
        "n_terminal_patch": str(stage4.get("n_terminal_patch", "ACE")).strip().upper(),
        "c_terminal_patch": str(stage4.get("c_terminal_patch", "CT3")).strip().upper(),
    }


def preflight_stage4_charmm_resources(
    config: dict[str, Any],
    *,
    root: Path = REPO_ROOT,
    required: bool = True,
    include_file_hashes: bool = False,
) -> Stage4ResourcePreflight:
    paths = resolve_stage4_resource_paths(config, root=root)
    errors: list[str] = []
    if not required:
        return Stage4ResourcePreflight(
            paths=paths,
            status="SKIP",
            errors=(),
            metadata=stage4_resource_metadata(paths),
            tree_hashes={
                "charmm36_forcefield": _tree_digest(paths.charmm36_forcefield, include_files=include_file_hashes),
                "membrane_toppar": _tree_digest(paths.membrane_toppar, include_files=include_file_hashes),
            },
        )

    if not paths.charmm36_forcefield.is_dir():
        errors.append(
            "Stage 4 is blocked: CHARMM36 GROMACS force-field directory is missing. "
            f"Resolved {paths.charmm36_forcefield} via {paths.charmm36_source}. "
            "Configure MEMBRANEFORGER_CHARMM36_ROOT=/path/to/charmm36.ff or stage4.resources.charmm36_forcefield."
        )
    else:
        for rel in CHARMM36_REQUIRED_FILES:
            if not (paths.charmm36_forcefield / rel).is_file():
                errors.append(f"STRUCTURAL: Stage 4 CHARMM36 GROMACS force field is incomplete: missing {paths.charmm36_forcefield / rel}")

    if not paths.membrane_toppar.is_dir():
        errors.append(
            "Stage 4 is blocked: CHARMM membrane topology/parameter directory is missing. "
            f"Resolved {paths.membrane_toppar} via {paths.membrane_toppar_source}. "
            "Configure MEMBRANEFORGER_MEMBRANE_TOPPAR=/path/to/toppar or stage4.resources.membrane_toppar."
        )
    else:
        for rel in TOPPAR_REQUIRED_FILES:
            if not (paths.membrane_toppar / rel).is_file():
                errors.append(f"STRUCTURAL: Stage 4 CHARMM membrane parameter tree is incomplete: missing {paths.membrane_toppar / rel}")
        protein_rtf = paths.membrane_toppar / "top_all36_prot.rtf"
        if protein_rtf.is_file():
            if not _has_charmm_patch(protein_rtf, "ACE"):
                errors.append(f"STRUCTURAL: Stage 4 CHARMM protein topology lacks PRES ACE: {protein_rtf}")
            if not _has_charmm_patch(protein_rtf, "CT3"):
                errors.append(f"STRUCTURAL: Stage 4 CHARMM protein topology lacks PRES CT3: {protein_rtf}")

    patches = _requested_terminal_patches(config)
    if paths.charmm36_forcefield.is_dir():
        if patches["n_terminal_patch"] != "NONE" and patches["n_terminal_patch"] not in _terminal_database_entries(paths.charmm36_forcefield, "N"):
            errors.append(
                f"STRUCTURAL: Stage 4 uses gmx pdb2gmx -ter to apply terminal patches, but n_terminal_patch={patches['n_terminal_patch']} is not available in "
                f"{paths.charmm36_forcefield / 'merged.n.tdb'}; do not use explicit ACE/NME residues as a workaround."
            )
        if patches["c_terminal_patch"] != "NONE" and patches["c_terminal_patch"] not in _terminal_database_entries(paths.charmm36_forcefield, "C"):
            errors.append(
                f"STRUCTURAL: Stage 4 uses gmx pdb2gmx -ter to apply terminal patches, but c_terminal_patch={patches['c_terminal_patch']} is not available in "
                f"{paths.charmm36_forcefield / 'merged.c.tdb'}; CHARMM methylamide C-terminal patch is CT3, not NME."
            )

    metadata = stage4_resource_metadata(paths)
    if paths.charmm36_forcefield.is_dir() and paths.membrane_toppar.is_dir():
        charmm_text = "\n".join([
            _read_text(paths.charmm36_forcefield / "forcefield.doc"),
            _read_text(paths.charmm36_forcefield / "forcefield.itp"),
        ])
        toppar_text = "\n".join([
            _read_text(paths.membrane_toppar / "top_all36_prot.rtf"),
            _read_text(paths.membrane_toppar / "par_all36m_prot.prm"),
            _read_text(paths.membrane_toppar / "toppar_history"),
            _read_text(paths.membrane_toppar / "toppar_all.history"),
        ])
        common_release_tokens = sorted(_month_year_tokens(charmm_text) & _month_year_tokens(toppar_text))
        if common_release_tokens:
            metadata["compatibility_release_token"] = ",".join(common_release_tokens)
        else:
            errors.append(
                "PROVENANCE: Stage 4 CHARMM resource compatibility is not established: the GROMACS force-field tree and "
                "toppar tree do not expose a shared release/date marker. Use a paired CHARMM36 GROMACS force field "
                "and matching toppar resources from an established distribution."
            )

    return Stage4ResourcePreflight(
        paths=paths,
        status="PASS" if not errors else "FAIL",
        errors=tuple(errors),
        metadata=metadata,
        tree_hashes={
            "charmm36_forcefield": _tree_digest(paths.charmm36_forcefield, include_files=include_file_hashes),
            "membrane_toppar": _tree_digest(paths.membrane_toppar, include_files=include_file_hashes),
        },
    )


def resolve_dssp(*, required: bool = False) -> DependencyRecord:
    path = _resolve_executable(os.environ.get("DSSP_BIN"), "mkdssp")
    if not path or not path.is_file():
        if required:
            raise DependencyError("DSSP/mkdssp is required for this mode; install mkdssp or set DSSP_BIN")
        return DependencyRecord("DSSP", "OPTIONAL_MISSING", os.environ.get("DSSP_BIN", "mkdssp"))
    return DependencyRecord("DSSP", "PASS", str(path), _version([str(path), "--version"]), _sha256(path))


def resolve_pyrosetta(*, required: bool = False) -> DependencyRecord:
    configured = os.environ.get("PYROSETTA_PYTHON")
    python = _resolve_executable(configured, None)
    if not python or not python.is_file():
        if required:
            raise DependencyError("PyRosetta is required for this mode; set PYROSETTA_PYTHON to a licensed Python environment")
        return DependencyRecord("PyRosetta", "OPTIONAL_MISSING", configured or "PYROSETTA_PYTHON")
    code = "import json, pyrosetta; print(json.dumps({'version': getattr(pyrosetta, '__version__', 'import-ok')}))"
    result = subprocess.run([str(python), "-c", code], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=20)
    if result.returncode != 0:
        if required:
            raise DependencyError("PyRosetta is required for this mode; configured PYROSETTA_PYTHON cannot import pyrosetta")
        return DependencyRecord("PyRosetta", "OPTIONAL_UNAVAILABLE", str(python), detail=(result.stderr or result.stdout).strip())
    try:
        version = json.loads(result.stdout.splitlines()[-1]).get("version", "import-ok")
    except Exception:
        version = "import-ok"
    return DependencyRecord("PyRosetta", "PASS", str(python), version, _sha256(python))


def resolve_rosetta_bin(*, required: bool = False) -> DependencyRecord:
    raw = os.environ.get("ROSETTA_BIN")
    path = Path(os.path.expanduser(os.path.expandvars(raw))).resolve() if raw else None
    if not path or not path.exists():
        if required:
            raise DependencyError("Rosetta utilities are required for this mode; set ROSETTA_BIN to a licensed Rosetta bin directory")
        return DependencyRecord("Rosetta utilities", "OPTIONAL_MISSING", raw or "ROSETTA_BIN")
    if path.is_file():
        return DependencyRecord("Rosetta utilities", "PASS", str(path), _version([str(path), "--help"]), _sha256(path))
    candidates = sorted(p for p in path.iterdir() if p.is_file() and os.access(p, os.X_OK))
    detail = f"{len(candidates)} executable files detected"
    return DependencyRecord("Rosetta utilities", "PASS", str(path), detail=detail)


def resolve_rosetta_database(*, required: bool = False) -> DependencyRecord:
    raw = os.environ.get("ROSETTA_DATABASE")
    path = Path(os.path.expanduser(os.path.expandvars(raw))).resolve() if raw else None
    if not path or not path.is_dir():
        if required:
            raise DependencyError("Rosetta database is required for this mode; set ROSETTA_DATABASE to a licensed database directory")
        return DependencyRecord("Rosetta database", "OPTIONAL_MISSING", raw or "ROSETTA_DATABASE")
    return DependencyRecord("Rosetta database", "PASS", str(path), detail="directory present")


def resolve_molfile_to_params(*, required: bool = False) -> DependencyRecord:
    path = _resolve_executable(os.environ.get("MOLFILE_TO_PARAMS"), "molfile_to_params.py")
    if not path or not path.is_file():
        if required:
            raise DependencyError("molfile_to_params.py is required for this mode; set MOLFILE_TO_PARAMS to the licensed Rosetta utility")
        return DependencyRecord("molfile_to_params", "OPTIONAL_MISSING", os.environ.get("MOLFILE_TO_PARAMS", "molfile_to_params.py"))
    return DependencyRecord("molfile_to_params", "PASS", str(path), _version([str(path), "--help"]), _sha256(path))


def resolve_directory_env(name: str, env: str, *, required: bool = False) -> DependencyRecord:
    raw = os.environ.get(env, "")
    path = Path(os.path.expanduser(os.path.expandvars(raw))).resolve() if raw else None
    if not path or not path.is_dir():
        if required:
            raise DependencyError(f"{name} is required for this mode; set {env} to a readable licensed directory")
        return DependencyRecord(name, "OPTIONAL_MISSING", raw or env)
    files = sorted(p for p in path.rglob("*") if p.is_file())
    digest = ""
    if files:
        h = hashlib.sha256()
        for p in files[:200]:
            h.update(str(p.relative_to(path)).encode())
            h.update(_sha256(p).encode())
        digest = h.hexdigest()
    return DependencyRecord(name, "PASS", str(path), sha256=digest, detail=f"{len(files)} files")


def resolve_directory_config_or_env(name: str, value: str | None, env: str, *, required: bool = False) -> DependencyRecord:
    candidates: list[tuple[str, Path]] = []
    raw = os.environ.get(env, "")
    if raw:
        expanded = os.path.expanduser(os.path.expandvars(raw))
        if "$" not in expanded:
            candidates.append((env, Path(expanded).resolve()))
    if value:
        expanded = os.path.expanduser(os.path.expandvars(value))
        if "$" not in expanded:
            candidates.append(("config", Path(expanded).resolve()))
    for source, path in candidates:
        if path.is_dir():
            files = sorted(p for p in path.rglob("*") if p.is_file())
            digest = ""
            if files:
                h = hashlib.sha256()
                for p in files[:200]:
                    h.update(str(p.relative_to(path)).encode())
                    h.update(_sha256(p).encode())
                digest = h.hexdigest()
            return DependencyRecord(name, "PASS", str(path), sha256=digest, detail=f"{len(files)} files via {source}")
    if required:
        raise DependencyError(f"{name} is required for this mode; set {env} or the corresponding stage4.resources path")
    return DependencyRecord(name, "OPTIONAL_MISSING", value or raw or env)


def resolve_stage4_dependencies(config: dict, *, required: bool = False) -> list[DependencyRecord]:
    stage4 = config.get("stage4", {})
    resources = stage4.get("resources", {}) if isinstance(stage4.get("resources", {}), dict) else {}
    records: list[DependencyRecord] = []
    records.append(resolve_dssp(required=bool(stage4.get("dssp_required", False)) and required))
    records.append(resolve_directory_config_or_env("CHARMM36", resources.get("charmm36_forcefield"), "MEMBRANEFORGER_CHARMM36_ROOT", required=bool(stage4.get("charmm36_required", False)) and required))
    records.append(resolve_directory_config_or_env("CGenFF", resources.get("cgenff_root"), "MEMBRANEFORGER_CGENFF_ROOT", required=bool(stage4.get("cgenff_required", False)) and required))
    records.append(resolve_directory_env("Ligand params", "MEMBRANEFORGER_LIGAND_PARAMS_ROOT", required=bool(stage4.get("ligand_params_required", False)) and required))
    records.append(resolve_pyrosetta(required=bool(stage4.get("pyrosetta_plugins_enabled", False)) and required))
    rosetta_required = bool(stage4.get("parameter_generation_enabled", False)) and required
    records.append(resolve_rosetta_bin(required=rosetta_required))
    records.append(resolve_rosetta_database(required=rosetta_required))
    records.append(resolve_molfile_to_params(required=rosetta_required))
    return records
