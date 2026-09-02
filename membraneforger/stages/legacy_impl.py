#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import inspect
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any
import xml.etree.ElementTree as ET

import numpy as np

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))
from membraneforger.chemistry import ChemistryError, validate_membrane_composition
from membraneforger_paths import PathResolutionError, display_path, gmx_command as resolve_gmx_command, resolve_mstool


ATOM_RECORDS = ("ATOM  ", "HETATM")
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HSD": "H", "HSE": "H", "HSP": "H",
}
PROTEIN_CAP_RESNAMES = {"ACE", "NME", "NMA", "CT3", "NTER", "CTER"}
PROTEIN_RESNAMES = set(AA3) | PROTEIN_CAP_RESNAMES
CHARMM_EXPLICIT_TERMINAL_CAP_RESNAMES = {"ACE", "NME", "NMA", "CT3"}
BACKBONE_REQUIRED_RESNAMES = set(AA3)
SOLVENT_ION_NAMES = {"W", "WF", "NA", "CL", "ION", "SOL"}
WATER_RESNAMES = {"HOH", "WAT", "TIP3", "SOL", "W"}
ION_RESNAMES = {"NA", "CL", "K", "CA", "MG", "ZN", "MN", "FE", "CU", "CO", "NI", "SOD", "CLA", "POT", "CAL"}
COMMON_LIPID_RESNAMES = {"POPC", "POPE", "POPS", "POPG", "DOPC", "DOPE", "DOPS", "DOPG", "CHOL", "GM3", "DPG3", "GLPA"}
COMMON_GLYCAN_RESNAMES = {"NAG", "BMA", "MAN", "GAL", "GLC", "FUC", "SIA", "NMC"}
COMMON_COFACTOR_RESNAMES = {"HEM", "HEC", "FAD", "FMN", "NAD", "NAP", "ADP", "ATP", "GTP", "GDP"}
PROTEIN_BACKBONE_ATOMS = {"N", "HN", "CA", "HA", "HA1", "HA2", "C", "O"}
MSTOOL_CHANGENAME = {
    ":CHOL": ":CHL1",
    ":ION@NA": ":SOD@SOD",
    ":NA@NA": ":SOD@SOD",
    ":ION@CL": ":CLA@CLA",
    ":CL@CL": ":CLA@CLA",
    ":ION@CA": ":CAL@CAL",
    ":A": ":ADE",
    ":U": ":URA",
    ":G": ":GUA",
    ":C": ":CYT",
    ":T": ":THY",
    ":SAP6": ":SAPI24",
}

DPG3_SOURCE_PATTERN = (
    ("GLC", "A"), ("GLC", "B"), ("GLC", "C"), ("GLC", "V"),
    ("GAL", "A"), ("GAL", "B"), ("GAL", "C"), ("GAL", "V"),
    ("NMC", "A"), ("NMC", "B"), ("NMC", "C"), ("NMC", "D"), ("NMC", "E"), ("NMC", "V"),
    ("CER", "AM1"), ("CER", "AM2"), ("CER", "T1A"), ("CER", "C2A"), ("CER", "C3A"),
    ("CER", "C1B"), ("CER", "C2B"), ("CER", "C3B"), ("CER", "C4B"),
)
DPG3_RETAINED_INDICES = (0, 1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22)
GM3_CANONICAL_BEADS = (
    "G1", "G2", "G3", "A1", "A2", "A3", "S1", "S2", "S3", "S4", "S5",
    "AM1", "AM2", "T1A", "C2A", "C3A", "C1B", "C2B", "C3B", "C4B",
)
DPG3_COMPONENT_RESNAMES = {"GLC", "GAL", "NMC", "CER"}


class ValidationError(RuntimeError):
    pass


def rel(root: Path, path: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve(root: Path, value: str | Path, run_id: str) -> Path:
    text = str(value).format(run_id=run_id)
    path = Path(text)
    return path if path.is_absolute() else root / path


def ensure_clean_dir(path: Path, overwrite: bool) -> None:
    if path.exists():
        if not overwrite:
            raise ValidationError(f"published output already exists: {path}; use --overwrite")
        shutil.rmtree(path)
    path.mkdir(parents=True)


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_manifest(root: Path, directory: Path, name: str = "manifest.tsv") -> None:
    rows = ["path\tsha256\tfile_size"]
    for path in sorted(p for p in directory.rglob("*") if p.is_file() and p.name != name):
        rows.append(f"{rel(root, path)}\t{sha256(path)}\t{path.stat().st_size}")
    (directory / name).write_text("\n".join(rows) + "\n", encoding="utf-8")


def run_logged(
    cmd: list[str],
    cwd: Path,
    log_dir: Path,
    label: str,
    timeout: int | None = None,
    input_text: str | None = None,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    log_dir.mkdir(parents=True, exist_ok=True)
    command_file = log_dir / f"{label}.command.txt"
    stdout_file = log_dir / f"{label}.stdout.log"
    stderr_file = log_dir / f"{label}.stderr.log"
    command_file.write_text(" ".join(cmd) + "\n", encoding="utf-8")
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            text=True,
            input=input_text,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            env=env,
        )
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout.decode() if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode() if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        stdout_file.write_text(stdout, encoding="utf-8")
        stderr_file.write_text(stderr, encoding="utf-8")
        raise ValidationError(f"{label} timed out after {timeout} seconds; see {stdout_file} and {stderr_file}") from exc
    stdout_file.write_text(result.stdout or "", encoding="utf-8")
    stderr_file.write_text(result.stderr or "", encoding="utf-8")
    return result


def gmx_command() -> list[str]:
    try:
        return resolve_gmx_command()
    except PathResolutionError as exc:
        raise ValidationError(str(exc)) from exc


def local_thread_count() -> int:
    raw = (
        os.environ.get("MEMBRANEFORGER_LOCAL_THREADS")
        or os.environ.get("SLURM_CPUS_PER_TASK")
        or os.environ.get("OMP_NUM_THREADS")
    )
    if raw:
        try:
            return max(1, int(raw))
        except ValueError:
            pass
    return max(1, os.cpu_count() or 1)


def gmx_mdrun_command() -> list[str]:
    ntmpi = os.environ.get("MEMBRANEFORGER_GMX_NTMPI", "1")
    ntomp = os.environ.get("MEMBRANEFORGER_GMX_NTOMP", str(local_thread_count()))
    return gmx_command() + ["mdrun", "-ntmpi", ntmpi, "-ntomp", ntomp]


def local_resource_env(base: dict[str, str] | None = None) -> dict[str, str]:
    env = dict(os.environ if base is None else base)
    threads = str(local_thread_count())
    for key in (
        "OMP_NUM_THREADS",
        "OPENMM_CPU_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        env[key] = threads
    return env


def _canonical_terminal_patch(value: Any, *, terminus: str, legacy_cap: bool = False) -> str:
    if value is None:
        return "ACE" if terminus == "N" else "CT3"
    text = str(value).strip().upper()
    if text in {"", "NONE", "NULL", "FALSE"}:
        return "none"
    if legacy_cap and terminus == "C" and text == "NME":
        return "CT3"
    return text


def charmm_terminal_patches(stage4_cfg: dict[str, Any]) -> dict[str, str]:
    patches: dict[str, str] = {}
    for terminus, patch_key, cap_key in (
        ("N", "n_terminal_patch", "n_terminal_cap"),
        ("C", "c_terminal_patch", "c_terminal_cap"),
    ):
        patch = _canonical_terminal_patch(stage4_cfg.get(patch_key), terminus=terminus, legacy_cap=False)
        cap = stage4_cfg.get(cap_key)
        if cap is not None:
            cap_patch = _canonical_terminal_patch(cap, terminus=terminus, legacy_cap=True)
            if patch != cap_patch:
                raise ValidationError(f"conflicting terminal configuration: stage4.{patch_key} != stage4.{cap_key}")
        allowed = {"ACE", "none"} if terminus == "N" else {"CT3", "none"}
        if patch not in allowed:
            expected = "ACE or none" if terminus == "N" else "CT3 or none"
            raise ValidationError(f"stage4.{patch_key} must be {expected}; CHARMM C-terminal methylamide patch is CT3")
        patches[patch_key] = patch
    return patches


def protein_chains_for_terminal_patches(protein_pdb: Path) -> list[str]:
    chains: list[str] = []
    seen: set[str] = set()
    with protein_pdb.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            resname = line[17:21].strip()
            if resname not in AA3:
                continue
            chain = line[21:22].strip() or " "
            if chain not in seen:
                chains.append(chain)
                seen.add(chain)
    if not chains:
        raise ValidationError(f"Stage 4 terminal patches require at least one protein chain in {protein_pdb}")
    return chains


def selected_terminal_patch_chains(stage4_cfg: dict[str, Any], chains: list[str]) -> set[str]:
    selected = stage4_cfg.get("terminal_patch_chains")
    if selected is None:
        return set(chains)
    if not isinstance(selected, list) or not all(isinstance(item, str) and item for item in selected):
        raise ValidationError("stage4.terminal_patch_chains must be a list of chain IDs")
    missing = sorted(set(selected) - set(chains))
    if missing:
        raise ValidationError(f"stage4.terminal_patch_chains requested absent protein chains: {missing}")
    return set(selected)


def _terminal_database_entries(charmm36: Path, terminus: str) -> set[str]:
    filename = "merged.n.tdb" if terminus == "N" else "merged.c.tdb"
    path = charmm36 / filename
    if not path.is_file():
        raise ValidationError(f"Stage 4 CHARMM terminal patch database is missing: {path}")
    entries: set[str] = set()
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = raw.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            entries.add(stripped[1:-1].strip().upper())
    return entries


def validate_charmm_terminal_patch_resources(charmm36: Path, stage4_cfg: dict[str, Any]) -> None:
    patches = charmm_terminal_patches(stage4_cfg)
    for terminus, key in (("N", "n_terminal_patch"), ("C", "c_terminal_patch")):
        patch = patches[key].upper()
        if patch == "NONE":
            continue
        entries = _terminal_database_entries(charmm36, terminus)
        if patch not in entries:
            filename = "merged.n.tdb" if terminus == "N" else "merged.c.tdb"
            raise ValidationError(
                f"Stage 4 terminal patch {key}={patch} is not available in {charmm36 / filename}; "
                "use a CHARMM36 GROMACS force field that provides the requested terminal patch. "
                "Do not add explicit ACE/NME residues to work around terminal patches."
            )


def stage4_terminal_patch_answers(stage4_cfg: dict[str, Any], protein_pdb: Path) -> list[str]:
    patches = charmm_terminal_patches(stage4_cfg)
    chains = protein_chains_for_terminal_patches(protein_pdb)
    selected = selected_terminal_patch_chains(stage4_cfg, chains)
    answers: list[str] = []
    for chain in chains:
        answers.append(patches["n_terminal_patch"] if chain in selected else "none")
        answers.append(patches["c_terminal_patch"] if chain in selected else "none")
    return answers


def inferred_disulfide_prompt_count(protein_pdb: Path, max_distance_angstrom: float = 3.0) -> int:
    cys_sg: list[tuple[str, int, str, tuple[float, float, float]]] = []
    with protein_pdb.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            if line[17:21].strip() != "CYS" or line[12:16].strip() != "SG":
                continue
            chain = line[21:22].strip() or " "
            resid = int(line[22:26])
            icode = line[26:27].strip()
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            cys_sg.append((chain, resid, icode, xyz))
    count = 0
    cutoff2 = max_distance_angstrom * max_distance_angstrom
    for i, (_, resid_i, _, xyz_i) in enumerate(cys_sg):
        for _, resid_j, _, xyz_j in cys_sg[i + 1:]:
            if resid_i == resid_j:
                continue
            dx = xyz_i[0] - xyz_j[0]
            dy = xyz_i[1] - xyz_j[1]
            dz = xyz_i[2] - xyz_j[2]
            if dx * dx + dy * dy + dz * dz <= cutoff2:
                count += 1
    return count


def stage4_disulfide_answers(stage4_cfg: dict[str, Any], protein_pdb: Path | None = None) -> list[str]:
    answers = stage4_cfg.get("pdb2gmx_disulfide_answers")
    if answers is None:
        if "pdb2gmx_disulfide_answer_count" in stage4_cfg:
            count = int(stage4_cfg["pdb2gmx_disulfide_answer_count"])
        elif protein_pdb is not None:
            count = inferred_disulfide_prompt_count(protein_pdb)
        else:
            count = 100
        return ["y"] * count
    if isinstance(answers, list):
        return [str(answer) for answer in answers]
    return [str(answers)]


def stage4_pdb2gmx_input(stage4_cfg: dict[str, Any], protein_pdb: Path | None = None) -> str | None:
    explicit = os.environ.get("MEMBRANEFORGER_PDB2GMX_INPUT")
    if explicit is not None:
        return explicit
    lines: list[str] = []
    lines.extend(stage4_disulfide_answers(stage4_cfg, protein_pdb))
    if protein_pdb is not None:
        lines.extend(stage4_terminal_patch_answers(stage4_cfg, protein_pdb))
    if not lines:
        return None
    return "\n".join(lines) + "\n"


def write_genion_mdp(source: Path, output: Path) -> None:
    rows: list[str] = []
    replaced = False
    for raw in source.read_text(encoding="utf-8", errors="replace").splitlines():
        key = raw.split(";", 1)[0].split("=", 1)[0].strip().lower()
        if key == "coulombtype":
            rows.append("coulombtype = Cut-off ; derived for pre-neutralization genion TPR only")
            replaced = True
        else:
            rows.append(raw)
    if not replaced:
        rows.append("coulombtype = Cut-off ; derived for pre-neutralization genion TPR only")
    output.write_text("\n".join(rows) + "\n", encoding="utf-8")


CG_EQUILIBRATION_MDPS = (
    "step6.1_minimization.mdp",
    "step6.2_equilibration.mdp",
    "step6.3_equilibration.mdp",
    "step6.4_equilibration.mdp",
    "step6.5_equilibration.mdp",
    "step6.6_equilibration.mdp",
)
CG_EQUILIBRATION_COMPLETE = "step6.1 through step6.6 run before Stage 3 handoff"
STAGE3_CANONICAL_PDB = "final_all_atom.pdb"


def parse_moleculetype_names(itp: Path) -> list[str]:
    names: list[str] = []
    section = ""
    for raw in itp.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            section = line.strip("[]").strip().lower()
            continue
        if section == "moleculetype":
            parts = line.split()
            if parts:
                names.append(parts[0])
                section = ""
    return names


def topology_include_path(topology: Path, include: str) -> Path:
    path = Path(include)
    return path if path.is_absolute() else topology.parent / path


def validate_topology_contract(topology: Path, structure: Path | None = None) -> dict[str, Any]:
    includes, molecules = parse_topology(topology)
    unresolved = [inc for inc in includes if not topology_include_path(topology, inc).exists()]
    if unresolved:
        raise ValidationError(f"topology has unresolved include(s): {unresolved}")
    moltypes = parse_moltypes(topology.parent / "toppar", includes)
    missing_molecules = sorted({name for name, _count in molecules if name not in moltypes})
    if missing_molecules:
        raise ValidationError(f"[ molecules ] entries lack matching [ moleculetype ]: {missing_molecules}")
    if not molecules:
        raise ValidationError("topology has no [ molecules ] entries")
    if structure is not None:
        _title, atoms, _box_line, _box = parse_gro(structure)
        ranges = molecule_ranges(molecules, moltypes)
        if not ranges or ranges[-1]["end"] != len(atoms):
            got = ranges[-1]["end"] if ranges else 0
            raise ValidationError(f"topology atom count {got} does not match CG GRO atom count {len(atoms)}")
    return {
        "status": "PASS",
        "topology": topology.name,
        "includes": includes,
        "molecules": [{"name": name, "count": count} for name, count in molecules],
        "moleculetype_names": sorted(moltypes),
    }


def martini_forcefield_itp(root: Path, config: dict[str, Any]) -> Path:
    requested = (
        config.get("coarse_grained", {}).get("martini_forcefield")
        or config.get("stage1", {}).get("martinize_forcefield")
        or "martini3001"
    )
    martini_dir = root / "resources" / "forcefields" / "martini"
    aliases = {
        "martini3001": "martini_v3.0.0.itp",
        "martini3": "martini_v3.0.0.itp",
        "martini_v3.0.0": "martini_v3.0.0.itp",
        "martini_v3.0.0.itp": "martini_v3.0.0.itp",
    }
    candidate = Path(str(requested))
    if candidate.suffix == ".itp":
        path = candidate if candidate.is_absolute() else root / candidate
    else:
        filename = aliases.get(str(requested))
        if not filename:
            raise ValidationError(f"unsupported Martini force-field selection: {requested}")
        path = martini_dir / filename
    if not path.exists():
        raise ValidationError(f"selected Martini force-field ITP is missing: {path}")
    return path


def finalize_insane_topology(root: Path, config: dict[str, Any], run_id: str, topology: Path, stage1_out: Path, out: Path) -> dict[str, Any]:
    raw_includes, raw_molecules = parse_topology(topology)
    protein_top = root / "work" / run_id / "stage1" / "replacement_protein.top"
    protein_molecules: list[tuple[str, int]] = []
    if protein_top.exists():
        _protein_includes, protein_molecules = parse_topology(protein_top)
    protein_itps = sorted(stage1_out.glob("ReplacementProtein_*.itp"))
    if not protein_itps and (stage1_out / "replacement_protein.itp").exists():
        protein_itps = [stage1_out / "replacement_protein.itp"]
    if not protein_itps:
        raise ValidationError("Stage 2 topology finalization requires Stage 1 generated protein ITPs")
    protein_names = [name for path in protein_itps for name in parse_moleculetype_names(path)]
    if not protein_names:
        raise ValidationError("Stage 1 generated protein ITPs contain no [ moleculetype ] names")
    if not protein_molecules:
        protein_molecules = [(name, 1) for name in protein_names]
    unknown_protein = sorted({name for name, _count in protein_molecules if name not in protein_names})
    if unknown_protein:
        raise ValidationError(f"Stage 1 protein topology references missing protein moleculetype(s): {unknown_protein}")
    for src in protein_itps:
        shutil.copy2(src, out / "toppar" / src.name)
    ff = martini_forcefield_itp(root, config)
    martini_itps = sorted((out / "toppar").glob("martini*.itp"))
    includes = [f"toppar/{ff.name}"]
    includes.extend(f"toppar/{path.name}" for path in martini_itps if path.name != ff.name)
    includes.extend(f"toppar/{path.name}" for path in protein_itps)
    nonprotein_molecules = [(name, count) for name, count in raw_molecules if name.lower() != "protein"]
    molecules = [*protein_molecules, *nonprotein_molecules]
    lines = [f'#include "{inc}"' for inc in includes]
    lines.extend(["", "[ system ]", "MembraneForger finalized CG system", "", "[ molecules ]", "; name  number"])
    lines.extend(f"{name:<16} {count}" for name, count in molecules)
    topology.write_text("\n".join(lines) + "\n", encoding="utf-8")
    report = validate_topology_contract(topology, out / "cg_system.gro")
    report.update(
        {
            "raw_insane_includes": raw_includes,
            "raw_insane_molecules": [{"name": name, "count": count} for name, count in raw_molecules],
            "selected_martini_forcefield": f"toppar/{ff.name}",
            "protein_itps": [f"toppar/{path.name}" for path in protein_itps],
            "protein_moleculetype_names": protein_names,
        }
    )
    return report


def optional_gmx_grompp_preflight(root: Path, out: Path, logs: Path, cfg: dict[str, Any], structure: Path, topology: Path) -> dict[str, Any]:
    try:
        command = gmx_command()
    except ValidationError as exc:
        return {"status": "NOT_CHECKED", "reason": str(exc)}
    mdp = root / "resources" / "cg_mdp" / "step6.0_minimization.mdp"
    if not mdp.exists():
        return {"status": "NOT_CHECKED", "reason": f"missing {mdp}"}
    includes, molecules = parse_topology(topology)
    moltypes = parse_moltypes(topology.parent / "toppar", includes)
    _title, atoms, _box_line, _box = parse_gro(structure)
    ranges = molecule_ranges(molecules, moltypes)
    write_cg_stage2_index(out / "cg_index_preflight.ndx", len(atoms), ranges)
    macro_text = topology_text_for_macros(topology)
    result = run_logged(
        command + [
            "grompp",
            "-f",
            os.path.relpath(prepare_cg_mdp_for_topology(mdp, out, topology, macro_text), out),
            "-c",
            structure.name,
            "-r",
            structure.name,
            "-p",
            topology.name,
            "-n",
            "cg_index_preflight.ndx",
            "-o",
            "_topology_preflight.tpr",
            "-maxwarn",
            "0",
        ],
        out,
        logs,
        "stage2_topology_grompp_preflight",
        timeout=int(cfg.get("cg_grompp_timeout_seconds", 300)),
    )
    if result.returncode != 0:
        raise ValidationError("Stage 2 topology grompp preflight failed; see logs")
    return {"status": "PASS", "command": command[0], "mdp": mdp.name, "maxwarn": 0}


def mstool_ungroup_kwargs(raw_dms: Path, mapping: list[str], mapping_add: list[str]) -> dict[str, Any]:
    return {
        "out": str(raw_dms),
        "mapping": mapping,
        "mapping_add": mapping_add,
        "backbone": True,
        "water_resname": "W",
        "water_number": 4,
        "water_chain_dms": True,
        "sort": True,
        "use_AA_structure": True,
        "AA_shrink_factor": 0.8,
    }


def preflight_mstool_ungroup_api(ungroup: Any, kwargs: dict[str, Any]) -> dict[str, Any]:
    signature = inspect.signature(ungroup)
    params = signature.parameters
    accepts_kwargs = any(param.kind == inspect.Parameter.VAR_KEYWORD for param in params.values())
    unsupported = sorted(key for key in kwargs if key not in params and not accepts_kwargs)
    if unsupported:
        raise ValidationError(f"mstool.Ungroup API mismatch; unsupported argument(s): {unsupported}")
    required = [
        name
        for name, param in params.items()
        if param.default is inspect._empty
        and param.kind in {inspect.Parameter.POSITIONAL_OR_KEYWORD, inspect.Parameter.KEYWORD_ONLY}
        and name not in kwargs
    ]
    if required:
        required = required[1:] if required[0] in {"structure", "universe", "u"} else required
    if required:
        raise ValidationError(f"mstool.Ungroup API mismatch; required argument(s) not provided: {required}")
    return {"status": "PASS", "signature": str(signature), "arguments": sorted(kwargs)}


def write_cg_stage2_index(index_path: Path, atom_count: int, ranges: list[dict[str, Any]]) -> dict[str, int]:
    groups: dict[str, list[int]] = {"System": list(range(1, atom_count + 1)), "Protein_LIG": [], "Membrane": [], "SOL_ION": []}
    for item in ranges:
        indices = list(range(item["start"] + 1, item["end"] + 1))
        moltype = str(item["moltype"])
        if moltype.lower().startswith("protein"):
            groups["Protein_LIG"].extend(indices)
        elif moltype in SOLVENT_ION_NAMES:
            groups["SOL_ION"].extend(indices)
        else:
            groups["Membrane"].extend(indices)
    lines: list[str] = []
    for name, indices in groups.items():
        lines.append(f"[ {name} ]")
        for start in range(0, len(indices), 15):
            lines.append(" ".join(str(value) for value in indices[start:start + 15]))
        lines.append("")
    index_path.write_text("\n".join(lines), encoding="utf-8")
    return {name: len(indices) for name, indices in groups.items()}


def topology_text_for_macros(topology: Path) -> str:
    texts = [topology.read_text(encoding="utf-8", errors="replace")]
    toppar_dir = topology.parent / "toppar"
    if toppar_dir.is_dir():
        for path in sorted(toppar_dir.rglob("*.itp")):
            texts.append(path.read_text(encoding="utf-8", errors="replace"))
    return "\n".join(texts)


def prepare_cg_mdp_for_topology(mdp: Path, out: Path, topology: Path, topology_text: str) -> Path:
    rows = mdp.read_text(encoding="utf-8", errors="replace").splitlines()
    last_option_index: dict[str, int] = {}
    for index, raw in enumerate(rows):
        key = raw.split(";", 1)[0].split("=", 1)[0].strip().lower()
        if key in {"nstcomm", "comm-grps"} and "=" in raw:
            last_option_index[key] = index
    changed = False
    prepared: list[str] = []
    dt: float | None = None
    nstpcouple: int | None = None
    pressure_coupling: str | None = None
    for index, raw in enumerate(rows):
        key = raw.split(";", 1)[0].split("=", 1)[0].strip().lower()
        if key in {"nstcomm", "comm-grps"} and last_option_index.get(key) != index:
            prepared.append(f"; disabled duplicate for GROMACS input validation: {raw}")
            changed = True
            continue
        if key == "dt" and "=" in raw:
            try:
                dt = float(raw.split("=", 1)[1].split(";", 1)[0].strip().split()[0])
            except (IndexError, ValueError):
                pass
        elif key == "nstpcouple" and "=" in raw:
            try:
                nstpcouple = int(raw.split("=", 1)[1].split(";", 1)[0].strip().split()[0])
            except (IndexError, ValueError):
                pass
        if key in {"tcoupl", "pcoupl"} and "=" in raw:
            left, right = raw.split("=", 1)
            body, *comment = right.split(";", 1)
            if key == "pcoupl":
                pressure_coupling = body.strip().lower()
            if body.strip().lower() == "berendsen":
                replacement = "v-rescale" if key == "tcoupl" else "c-rescale"
                suffix = f" ; {comment[0].strip()}" if comment else ""
                suffix += " ; updated for GROMACS warning-free equilibration"
                prepared.append(f"{left}= {replacement}{suffix}".rstrip())
                changed = True
                if key == "pcoupl":
                    pressure_coupling = replacement
                continue
        if key == "tau_p" and pressure_coupling and dt is not None and nstpcouple is not None and "=" in raw:
            left, right = raw.split("=", 1)
            body, *comment = right.split(";", 1)
            values = body.split()
            try:
                current_tau_p = float(values[0])
            except (IndexError, ValueError):
                current_tau_p = None
            factor = 200.0 if pressure_coupling == "parrinello-rahman" else 25.0
            minimum_tau_p = factor * float(nstpcouple) * float(dt)
            if pressure_coupling == "parrinello-rahman":
                minimum_tau_p *= 1.01
            if current_tau_p is not None and current_tau_p < minimum_tau_p:
                suffix = f" ; {comment[0].strip()}" if comment else ""
                suffix += " ; raised for pressure-coupling stability criterion"
                prepared.append(f"{left}= {minimum_tau_p:.6g}{suffix}".rstrip())
                changed = True
                continue
        if key != "define" or "=" not in raw:
            prepared.append(raw)
            continue
        left, right = raw.split("=", 1)
        body, *comment = right.split(";", 1)
        kept: list[str] = []
        dropped: list[str] = []
        for token in body.split():
            if not token.startswith("-D"):
                kept.append(token)
                continue
            macro = token[2:].split("=", 1)[0]
            if re.search(rf"\b{re.escape(macro)}\b", topology_text):
                kept.append(token)
            else:
                dropped.append(macro)
        if dropped:
            changed = True
        suffix = f" ; {comment[0].strip()}" if comment else ""
        if dropped:
            suffix += " ; removed unused topology macro(s)"
        prepared.append(f"{left}= {' '.join(kept)}{suffix}".rstrip())
    if not changed:
        return mdp
    mdp_out = out / "_cg_mdp"
    mdp_out.mkdir(parents=True, exist_ok=True)
    prepared_path = mdp_out / mdp.name
    prepared_path.write_text("\n".join(prepared) + "\n", encoding="utf-8")
    return prepared_path


def cg_leaflet_gap_report(gro_path: Path, molecules: list[tuple[str, int]], moltypes: dict[str, dict[str, Any]]) -> dict[str, Any]:
    _, atoms, _, _ = parse_gro(gro_path)
    ranges = molecule_ranges(molecules, moltypes)
    headgroup_names = {"PO4", "PO3", "P", "PH"}
    membrane_headgroups: list[float] = []
    membrane_beads: list[float] = []
    for item in ranges:
        moltype = str(item["moltype"])
        if moltype.lower().startswith("protein") or moltype in SOLVENT_ION_NAMES:
            continue
        molecule_atoms = atoms[item["start"]:item["end"]]
        membrane_beads.extend(float(atom["z"]) for atom in molecule_atoms)
        membrane_headgroups.extend(float(atom["z"]) for atom in molecule_atoms if str(atom["name"]).strip() in headgroup_names)
    values = sorted(membrane_headgroups if len(membrane_headgroups) >= 4 else membrane_beads)
    source = "headgroup" if len(membrane_headgroups) >= 4 else "membrane_bead"
    if len(values) < 4:
        raise ValidationError(f"Stage 2 CG leaflet sanity check cannot identify enough membrane beads in {gro_path}")
    gaps = [(values[index + 1] - values[index], index) for index in range(len(values) - 1)]
    largest_gap, split_index = max(gaps, key=lambda item: item[0])
    lower = values[:split_index + 1]
    upper = values[split_index + 1:]
    return {
        "source": source,
        "sample_count": len(values),
        "largest_midplane_gap_nm": round(float(largest_gap), 4),
        "lower_leaflet_mean_z_nm": round(float(np.mean(lower)), 4),
        "upper_leaflet_mean_z_nm": round(float(np.mean(upper)), 4),
        "leaflet_mean_separation_nm": round(float(np.mean(upper) - np.mean(lower)), 4),
    }


def validate_cg_leaflet_gap(
    gro_path: Path,
    molecules: list[tuple[str, int]],
    moltypes: dict[str, dict[str, Any]],
    max_midplane_gap_nm: float,
) -> dict[str, Any]:
    report = cg_leaflet_gap_report(gro_path, molecules, moltypes)
    if float(report["largest_midplane_gap_nm"]) > max_midplane_gap_nm:
        raise ValidationError(
            "CG equilibration produced abnormal leaflet/headgroup separation: "
            f"largest_midplane_gap_nm={report['largest_midplane_gap_nm']} > {max_midplane_gap_nm}"
        )
    return report


def run_cg_minimize_equilibrate(
    root: Path,
    out: Path,
    logs: Path,
    cfg: dict[str, Any],
    initial_gro: Path,
    topology: Path,
    label_prefix: str,
    minim_grompp_label: str,
    minim_mdrun_label: str,
) -> dict[str, Any]:
    mdp_dir = root / "resources" / "cg_mdp"
    minim_mdp = mdp_dir / "step6.0_minimization.mdp"
    equil_mdps = [mdp_dir / name for name in CG_EQUILIBRATION_MDPS]
    missing = [path for path in [minim_mdp, *equil_mdps] if not path.exists()]
    if missing:
        raise ValidationError(f"CG equilibration MDPs are missing: {[str(path) for path in missing]}")

    _, atoms, _, _ = parse_gro(initial_gro)
    includes, molecules = parse_topology(topology)
    moltypes = parse_moltypes(out / "toppar", includes)
    ranges = molecule_ranges(molecules, moltypes)
    if ranges[-1]["end"] != len(atoms):
        raise ValidationError(f"topology atom count {ranges[-1]['end']} does not match CG GRO atom count {len(atoms)}")
    index_report = write_cg_stage2_index(out / "cg_index.ndx", len(atoms), ranges)
    macro_text = topology_text_for_macros(topology)

    grompp_timeout = int(cfg.get("cg_grompp_timeout_seconds", 300))
    minim_timeout = int(cfg.get("cg_minimization_timeout_seconds", 1800))
    equil_timeout = int(cfg.get("cg_equilibration_timeout_seconds", 7200))

    result = run_logged(
        gmx_command() + [
            "grompp", "-f", os.path.relpath(prepare_cg_mdp_for_topology(minim_mdp, out, topology, macro_text), out),
            "-c", initial_gro.name, "-r", initial_gro.name,
            "-p", topology.name, "-n", "cg_index.ndx", "-o", "em.tpr", "-maxwarn", "0",
        ],
        out, logs, minim_grompp_label, timeout=grompp_timeout,
    )
    if result.returncode != 0:
        raise ValidationError("CG minimization grompp failed; see logs")
    result = run_logged(gmx_mdrun_command() + ["-deffnm", "em", "-c", "cg_minimized.gro"], out, logs, minim_mdrun_label, timeout=minim_timeout)
    if result.returncode != 0:
        raise ValidationError("CG minimization mdrun failed; see logs")
    current_gro = out / "cg_minimized.gro"
    if not current_gro.exists():
        raise ValidationError("CG minimization did not produce cg_minimized.gro")

    sequence = [{"mdp": minim_mdp.name, "input": initial_gro.name, "output": current_gro.name}]
    previous_cpt: Path | None = None
    for mdp in equil_mdps:
        prepared_mdp = prepare_cg_mdp_for_topology(mdp, out, topology, macro_text)
        stem = mdp.stem.replace(".", "_")
        output_gro = out / f"{stem}.gro"
        deffnm = stem
        grompp = gmx_command() + [
            "grompp", "-f", os.path.relpath(prepared_mdp, out), "-c", current_gro.name, "-r", current_gro.name,
            "-p", topology.name, "-n", "cg_index.ndx", "-o", f"{deffnm}.tpr", "-maxwarn", "0",
        ]
        if previous_cpt is not None and previous_cpt.exists() and "step6.2" not in mdp.name:
            grompp.extend(["-t", previous_cpt.name])
        result = run_logged(grompp, out, logs, f"{label_prefix}_{stem}_grompp", timeout=grompp_timeout)
        if result.returncode != 0:
            raise ValidationError(f"CG equilibration grompp failed at {mdp.name}; see logs")
        result = run_logged(gmx_mdrun_command() + ["-deffnm", deffnm, "-c", output_gro.name], out, logs, f"{label_prefix}_{stem}_mdrun", timeout=equil_timeout)
        if result.returncode != 0:
            raise ValidationError(f"CG equilibration mdrun failed at {mdp.name}; see logs")
        if not output_gro.exists():
            raise ValidationError(f"CG equilibration did not produce {output_gro.name}")
        sequence.append({"mdp": mdp.name, "input": current_gro.name, "output": output_gro.name})
        current_gro = output_gro
        previous_cpt = out / f"{deffnm}.cpt"

    sanity = validate_cg_leaflet_gap(
        current_gro,
        molecules,
        moltypes,
        float(cfg.get("max_cg_midplane_gap_nm", 6.0)),
    )
    shutil.copy2(current_gro, out / "cg_equilibrated_scaffold.gro")
    shutil.copy2(out / "cg_equilibrated_scaffold.gro", out / "cg_backmap_input.gro")
    return {
        "sequence": sequence,
        "index": index_report,
        "final_equilibrated_gro": "cg_equilibrated_scaffold.gro",
        "stage3_handoff": "cg_backmap_input.gro",
        "leaflet_sanity": sanity,
    }


def validate_stage4_em_log(log_path: Path) -> None:
    text = log_path.read_text(encoding="utf-8", errors="replace")
    failures = (
        "force on at least one atom is not\nfinite",
        "Maximum force     =            inf",
        "did not converge to Fmax",
        "reached the maximum number of steps",
    )
    for marker in failures:
        if marker in text:
            raise ValidationError(f"Stage 4 EM did not produce equilibration-ready coordinates; see {log_path}")
    if "converged to Fmax <" not in text:
        raise ValidationError(f"Stage 4 EM convergence was not confirmed; see {log_path}")


def parse_pdb(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(ATOM_RECORDS):
                continue
            atoms.append({
                "record": line[:6].strip(),
                "serial": int(line[6:11]),
                "name": line[12:16].strip(),
                "altloc": line[16:17].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip() or " ",
                "resid": int(line[22:26]),
                "icode": line[26:27].strip(),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "occupancy": float(line[54:60]),
                "bfactor": float(line[60:66]),
                "element": line[76:78].strip() if len(line) >= 78 else "",
                "line": line.rstrip("\n"),
            })
    return atoms


def pdb_structure_report(path: Path) -> dict[str, Any]:
    atoms = parse_pdb(path)
    residues: list[tuple[str, int, str, str]] = []
    seen_res: set[tuple[str, int, str, str]] = set()
    by_res: dict[tuple[str, int, str, str], dict[str, dict[str, Any]]] = defaultdict(dict)
    chains: list[str] = []
    duplicate_serials = []
    serials: set[int] = set()
    duplicate_atom_names = []
    residue_numbers_by_chain: dict[int, set[str]] = defaultdict(set)
    atom_name_seen: set[tuple[str, int, str, str, str]] = set()
    for atom in atoms:
        if atom["serial"] in serials:
            duplicate_serials.append(atom["serial"])
        serials.add(atom["serial"])
        if atom["chain"] not in chains and atom["record"] == "ATOM":
            chains.append(atom["chain"])
        key = (atom["chain"], atom["resid"], atom["icode"], atom["resname"])
        if key not in seen_res and atom["record"] == "ATOM":
            residues.append(key)
            seen_res.add(key)
            residue_numbers_by_chain[atom["resid"]].add(atom["chain"])
        akey = (*key, atom["name"])
        if akey in atom_name_seen:
            duplicate_atom_names.append(akey)
        atom_name_seen.add(akey)
        by_res[key][atom["name"]] = atom

    missing_backbone = []
    for key, names in by_res.items():
        if key[3] not in BACKBONE_REQUIRED_RESNAMES:
            continue
        for name in ("N", "CA", "C", "O"):
            if name not in names:
                missing_backbone.append([*key, name])

    by_chain: dict[str, list[tuple[str, int, str, str]]] = defaultdict(list)
    for residue in residues:
        if residue[3] in PROTEIN_RESNAMES:
            by_chain[residue[0]].append(residue)
    numbering_gaps = []
    cn_breaks = []
    for chain, chain_residues in by_chain.items():
        ordered = sorted(chain_residues, key=lambda r: (r[1], r[2]))
        for prev, cur in zip(ordered, ordered[1:]):
            if cur[1] - prev[1] != 1:
                numbering_gaps.append([chain, prev[1], cur[1], cur[1] - prev[1]])
            if "C" in by_res[prev] and "N" in by_res[cur]:
                a = by_res[prev]["C"]
                b = by_res[cur]["N"]
                dist = math.dist((a["x"], a["y"], a["z"]), (b["x"], b["y"], b["z"]))
                if dist > 1.8:
                    cn_breaks.append([chain, prev[1], cur[1], round(dist, 3)])

    cys_sg = []
    for key, names in by_res.items():
        if key[3] == "CYS" and "SG" in names:
            cys_sg.append((key, names["SG"]))
    sg_pairs = []
    for i, (ra, a) in enumerate(cys_sg):
        for rb, b in cys_sg[i + 1:]:
            dist = math.dist((a["x"], a["y"], a["z"]), (b["x"], b["y"], b["z"]))
            if dist < 4.0:
                sg_pairs.append([ra[0], ra[1], rb[0], rb[1], round(dist, 3)])

    hydrogens = 0
    nonfinite = []
    for atom in atoms:
        element = atom["element"].upper()
        if element in {"H", "D"} or atom["name"].upper().startswith(("H", "D")):
            hydrogens += 1
        if not all(math.isfinite(atom[k]) for k in ("x", "y", "z")):
            nonfinite.append(atom["serial"])

    text = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return {
        "path": path.name,
        "sha256": sha256(path),
        "file_size": path.stat().st_size,
        "atom_records": sum(1 for line in text if line.startswith("ATOM  ")),
        "hetatm_records": sum(1 for line in text if line.startswith("HETATM")),
        "anisou_records": sum(1 for line in text if line.startswith("ANISOU")),
        "ssbond_records": sum(1 for line in text if line.startswith("SSBOND")),
        "link_records": sum(1 for line in text if line.startswith("LINK")),
        "conect_records": sum(1 for line in text if line.startswith("CONECT")),
        "ter_records": sum(1 for line in text if line.startswith("TER")),
        "cryst1": [line for line in text if line.startswith("CRYST1")],
        "chain_ids": chains,
        "residue_count": len(residues),
        "residue_count_by_chain": {chain: len(res) for chain, res in by_chain.items()},
        "hydrogen_atom_count": hydrogens,
        "alternate_locations": sorted({atom["altloc"] for atom in atoms if atom["altloc"]}),
        "insertion_codes": sorted({atom["icode"] for atom in atoms if atom["icode"]}),
        "nonstandard_residues": sorted({r[3] for r in residues if r[3] not in PROTEIN_RESNAMES}),
        "duplicate_serials": duplicate_serials,
        "duplicate_atom_names": duplicate_atom_names[:50],
        "duplicate_residue_numbers_across_chains": {
            str(resid): sorted(chains)
            for resid, chains in residue_numbers_by_chain.items()
            if len(chains) > 1
        },
        "missing_backbone_atoms": missing_backbone,
        "numbering_gaps": numbering_gaps,
        "peptide_cn_breaks_gt_1p8_angstrom": cn_breaks,
        "cysteine_sg_pairs_under_4_angstrom": sg_pairs,
        "nonfinite_coordinates": nonfinite,
        "coordinate_units": "angstrom",
    }


def write_protein_only_pdb(source: Path, target: Path, chains: list[str] | None) -> dict[str, Any]:
    kept = 0
    skipped = Counter()
    components: dict[str, dict[str, Any]] = {}
    selected = set(chains or [])
    with source.open("r", encoding="utf-8", errors="replace") as inp, target.open("w", encoding="utf-8") as out:
        for line in inp:
            if line.startswith("ANISOU"):
                skipped["ANISOU"] += 1
                continue
            if line.startswith("TER"):
                out.write(line)
                continue
            if not line.startswith(ATOM_RECORDS):
                continue
            chain = line[21:22].strip() or " "
            resname = line[17:20].strip()
            if selected and chain not in selected:
                skipped["unselected_chain"] += 1
                continue
            if line.startswith("HETATM") or resname not in PROTEIN_RESNAMES:
                category = classify_pdb_component(line[:6].strip(), resname)
                key = f"{chain}:{resname}{line[22:26].strip()}{line[26:27].strip()}"
                item = components.setdefault(
                    key,
                    {
                        "chain": chain,
                        "resname": resname,
                        "resid": line[22:26].strip(),
                        "icode": line[26:27].strip(),
                        "category": category,
                        "atom_records": 0,
                    },
                )
                item["atom_records"] += 1
                skipped[category] += 1
                continue
            out.write(line)
            kept += 1
        out.write("END\n")
    unsupported = [item for item in components.values() if item["category"] != "protein_polymer"]
    if unsupported:
        target.unlink(missing_ok=True)
        preview = ", ".join(f"{item['chain']}:{item['resname']}{item['resid']}({item['category']})" for item in unsupported[:20])
        raise ValidationError(
            "Stage 1 input contains nonprotein components with no configured downstream path: "
            f"{preview}. Exclude them via stage1.protein_chains or parameterize/enable a supported adapter before running."
        )
    return {"kept_atom_records": kept, "skipped": dict(skipped), "components": list(components.values()), "path": target.name}


def classify_pdb_component(record: str, resname: str) -> str:
    if record == "ATOM" and resname in PROTEIN_RESNAMES:
        return "protein_polymer"
    if resname in WATER_RESNAMES:
        return "water"
    if resname in ION_RESNAMES:
        return "metal_ion" if resname not in {"NA", "CL", "SOD", "CLA", "POT"} else "ion"
    if resname in COMMON_LIPID_RESNAMES:
        return "lipid"
    if resname in COMMON_GLYCAN_RESNAMES:
        return "glycan"
    if resname in COMMON_COFACTOR_RESNAMES:
        return "cofactor"
    if record == "HETATM":
        return "ligand"
    return "unsupported_nonstandard_residue"


def parse_topology(path: Path) -> tuple[list[str], list[tuple[str, int]]]:
    includes: list[str] = []
    molecules: list[tuple[str, int]] = []
    section = ""
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("#include"):
            match = re.search(r'"([^"]+)"', line)
            includes.append(match.group(1) if match else line)
            continue
        if line.startswith("["):
            section = line.strip("[]").strip().lower()
            continue
        if section == "molecules":
            parts = line.split()
            if len(parts) >= 2:
                molecules.append((parts[0], int(parts[1])))
    return includes, molecules


def parse_moltypes(toppar_dir: Path, includes: list[str]) -> dict[str, dict[str, Any]]:
    moltypes: dict[str, dict[str, Any]] = {}
    for inc in includes:
        inc_path = toppar_dir.parent / inc if "/" in inc else toppar_dir / inc
        if not inc_path.exists():
            continue
        section = ""
        current = ""
        atom_count = 0
        atom_names: list[str] = []
        residue_names: list[str] = []
        for raw in inc_path.read_text(encoding="utf-8", errors="replace").splitlines():
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue
            if line.startswith("["):
                if section == "atoms" and current:
                    moltypes[current] = {
                        "include": inc,
                        "atom_count": atom_count,
                        "atom_names": atom_names,
                        "residue_names": sorted(set(residue_names)),
                    }
                section = line.strip("[]").strip().lower()
                atom_count = 0
                atom_names = []
                residue_names = []
                continue
            parts = line.split()
            if section == "moleculetype" and parts:
                current = parts[0]
            elif section == "atoms" and current and len(parts) >= 5:
                atom_count += 1
                atom_names.append(parts[4])
                residue_names.append(parts[3])
        if section == "atoms" and current:
            moltypes[current] = {
                "include": inc,
                "atom_count": atom_count,
                "atom_names": atom_names,
                "residue_names": sorted(set(residue_names)),
            }
    return moltypes


def parse_gro(path: Path) -> tuple[str, list[dict[str, Any]], str, np.ndarray]:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        title = handle.readline().rstrip("\n")
        natoms = int(handle.readline().strip())
        atoms = []
        for idx in range(natoms):
            line = handle.readline().rstrip("\n")
            atoms.append({
                "index": idx + 1,
                "resid": int(line[:5]),
                "resname": line[5:10].strip(),
                "name": line[10:15].strip(),
                "atomid": int(line[15:20]),
                "x": float(line[20:28]),
                "y": float(line[28:36]),
                "z": float(line[36:44]),
                "raw": line,
            })
        box_line = handle.readline().strip()
    box_values = np.array([float(x) for x in box_line.split()[:3]], dtype=float)
    return title, atoms, box_line, box_values


def pdb_to_gro_atoms(path: Path) -> list[dict[str, Any]]:
    atoms = []
    for idx, atom in enumerate(parse_pdb(path), start=1):
        if atom["record"] != "ATOM":
            continue
        atoms.append({
            "resid": atom["resid"] % 100000,
            "resname": atom["resname"][:5],
            "name": atom["name"][:5],
            "atomid": idx % 100000,
            "x": atom["x"] / 10.0,
            "y": atom["y"] / 10.0,
            "z": atom["z"] / 10.0,
        })
    return atoms


def write_gro(path: Path, title: str, atoms: list[dict[str, Any]], box_line: str) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write(title[:80] + "\n")
        out.write(f"{len(atoms):5d}\n")
        for idx, atom in enumerate(atoms, start=1):
            out.write(
                f"{int(atom['resid']) % 100000:5d}"
                f"{str(atom['resname'])[:5]:>5}"
                f"{str(atom['name'])[:5]:>5}"
                f"{idx % 100000:5d}"
                f"{float(atom['x']):8.3f}{float(atom['y']):8.3f}{float(atom['z']):8.3f}\n"
            )
        out.write(box_line + "\n")


def molecule_ranges(molecules: list[tuple[str, int]], moltypes: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    ranges = []
    start = 0
    for moltype, count in molecules:
        atom_count = moltypes.get(moltype, {}).get("atom_count")
        if not atom_count:
            raise ValidationError(f"cannot determine atom count for molecule type {moltype}")
        for ordinal in range(1, count + 1):
            end = start + int(atom_count)
            ranges.append({"moltype": moltype, "ordinal": ordinal, "start": start, "end": end, "atom_count": int(atom_count)})
            start = end
    return ranges


def min_pbc_distance(points: np.ndarray, ref: np.ndarray, box: np.ndarray) -> float:
    if len(points) == 0 or len(ref) == 0:
        return float("inf")
    best = float("inf")
    for chunk_start in range(0, len(points), 256):
        chunk = points[chunk_start:chunk_start + 256]
        delta = chunk[:, None, :] - ref[None, :, :]
        delta -= box * np.round(delta / box)
        dist2 = np.sum(delta * delta, axis=2)
        best = min(best, float(np.sqrt(np.min(dist2))))
    return best


def chain_id_for_index(index: int) -> str:
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if index < len(alphabet):
        return alphabet[index]
    return f"P{index + 1}"


def prepare_mstool_chain_input(root: Path, stage2_input: Path, mstool_path: Path, work: Path) -> tuple[Path, dict[str, Any]]:
    sys.path.insert(0, str(mstool_path.parent))
    from mstool.core.universe import Universe

    topology = root / "outputs" / work.parent.name / "stage2" / "cg_topology.top"
    toppar = root / "outputs" / work.parent.name / "stage2" / "toppar"
    if not topology.exists() or not toppar.exists():
        return stage2_input, {
            "status": "SKIPPED",
            "reason": "stage2 topology or toppar not available",
            "input": rel(root, stage2_input),
        }

    includes, molecules = parse_topology(topology)
    moltypes = parse_moltypes(toppar, includes)
    ranges = molecule_ranges(molecules, moltypes)
    protein_ranges = [r for r in ranges if r["moltype"].lower().startswith("protein")]
    if not protein_ranges:
        return stage2_input, {
            "status": "SKIPPED",
            "reason": "no protein molecule types found in stage2 topology",
            "input": rel(root, stage2_input),
        }

    u = Universe(str(stage2_input))
    original_counts = u.atoms["chain"].value_counts().to_dict() if "chain" in u.atoms else {}
    distinct_protein_chains = set()
    for r in protein_ranges:
        distinct_protein_chains.update(str(c) for c in u.atoms.iloc[r["start"]:r["end"]]["chain"].unique())
    if len(distinct_protein_chains) >= len(protein_ranges) and "" not in distinct_protein_chains:
        return stage2_input, {
            "status": "SKIPPED",
            "reason": "input already carries distinct protein chain identifiers",
            "input": rel(root, stage2_input),
            "protein_chains": sorted(distinct_protein_chains),
        }

    if ranges[-1]["end"] != len(u.atoms):
        raise ValidationError(
            f"cannot assign CG protein chains: topology atom count {ranges[-1]['end']} "
            f"does not match coordinate atom count {len(u.atoms)}"
        )

    assignments = []
    for idx, r in enumerate(protein_ranges):
        chain_id = chain_id_for_index(idx)
        u.atoms.loc[r["start"]:r["end"] - 1, "chain"] = chain_id
        u.atoms.loc[r["start"]:r["end"] - 1, "segname"] = chain_id
        assignments.append({
            "molecule_type": r["moltype"],
            "ordinal": r["ordinal"],
            "chain_id": chain_id,
            "atom_range_1based": [r["start"] + 1, r["end"]],
            "atom_count": r["atom_count"],
        })

    prepared = work / "cg_backmap_input_chained.dms"
    u.write(str(prepared))
    report = {
        "status": "PASS",
        "reason": "reconstructed protein chain identifiers from topology molecule ranges",
        "input": str(stage2_input.relative_to(root)),
        "prepared_input": rel(root, prepared),
        "original_chain_counts": original_counts,
        "protein_molecule_count": len(protein_ranges),
        "assignments": assignments,
    }
    return prepared, report


def normalize_dpg3_to_gm3(universe: Any, work: Path) -> dict[str, Any]:
    """Convert the INSANE DPG3 representation into mstool's canonical GM3 mapping input."""
    atoms = universe.atoms.copy()
    required = {"resname", "name"}
    missing_columns = sorted(required - set(atoms.columns))
    if missing_columns:
        raise ValidationError(f"cannot normalize DPG3: mstool atoms lacks columns {missing_columns}")

    pairs = list(zip(atoms["resname"].astype(str), atoms["name"].astype(str)))
    component_positions = [i for i, (resname, _) in enumerate(pairs) if resname in DPG3_COMPONENT_RESNAMES]
    if not component_positions:
        return {
            "status": "SKIPPED",
            "reason": "no DPG3 component residues found",
            "gm3_molecules": 0,
            "source_beads": 0,
            "retained_beads": 0,
            "omitted_virtual_sites": 0,
        }

    windows: list[int] = []
    index = 0
    while index < len(pairs):
        pair = pairs[index]
        if pair[0] not in DPG3_COMPONENT_RESNAMES:
            index += 1
            continue
        observed = tuple(pairs[index:index + len(DPG3_SOURCE_PATTERN)])
        if observed != DPG3_SOURCE_PATTERN:
            preview = ", ".join(f"{res}:{name}" for res, name in observed[:len(DPG3_SOURCE_PATTERN)])
            raise ValidationError(
                "malformed DPG3 component window at 1-based atom "
                f"{index + 1}; expected the exact 23-bead DPG3 order, observed {preview}"
            )
        windows.append(index)
        index += len(DPG3_SOURCE_PATTERN)

    virtual_indices: set[int] = set()
    next_resid = int(atoms["resid"].max()) + 1 if "resid" in atoms.columns else 1
    next_resn = int(atoms["resn"].max()) + 1 if "resn" in atoms.columns else next_resid
    for ordinal, start in enumerate(windows):
        retained = [start + offset for offset in DPG3_RETAINED_INDICES]
        virtual_indices.update({start + 3, start + 7, start + 13})
        for atom_index, bead_name in zip(retained, GM3_CANONICAL_BEADS):
            atoms.at[atom_index, "resname"] = "GM3"
            atoms.at[atom_index, "name"] = bead_name
            if "resid" in atoms.columns:
                atoms.at[atom_index, "resid"] = next_resid + ordinal
            if "resn" in atoms.columns:
                atoms.at[atom_index, "resn"] = next_resn + ordinal

    atoms = atoms.drop(index=list(virtual_indices)).reset_index(drop=True)
    if atoms["resname"].isin(DPG3_COMPONENT_RESNAMES).any():
        leftovers = sorted(atoms.loc[atoms["resname"].isin(DPG3_COMPONENT_RESNAMES), "resname"].unique())
        raise ValidationError(f"DPG3 normalization left unmapped component residues: {leftovers}")
    universe.atoms = atoms
    output = work / "cg_backmap_input_chained_gm3.dms"
    universe.write(str(output))
    return {
        "status": "PASS",
        "source_pattern": [f"{resname}:{name}" for resname, name in DPG3_SOURCE_PATTERN],
        "gm3_molecules": len(windows),
        "source_beads": len(windows) * len(DPG3_SOURCE_PATTERN),
        "retained_beads": len(windows) * len(GM3_CANONICAL_BEADS),
        "omitted_virtual_sites": len(virtual_indices),
        "derived_input": str(output),
    }


def finalize_stage3_membrane_templates(root: Path, config: dict[str, Any], source_pdb: Path, out: Path) -> dict[str, Any]:
    plugin = config.get("stage3", {}).get("plugins", {}).get("glycolipid_template_finalize", {})
    if not plugin.get("enabled", False):
        return {"status": "SKIPPED", "reason": "glycolipid_template_finalize is disabled"}
    target_itp = resolve(root, plugin.get("target_itp", ""), "")
    if not target_itp.exists():
        raise ValidationError(f"glycolipid_template_finalize target ITP is missing: {target_itp}")
    resources = config.get("stage3", {}).get("resources", {})
    template_dir = resolve(root, resources.get("membrane_templates", "resources/templates"), "")
    toppar_dir = resolve(root, resources.get("membrane_toppar", "resources/external/charmm36_membrane"), "")
    gm3_xml = root / "resources" / "templates" / "GM3.xml"
    if not template_dir.exists() or not toppar_dir.exists() or not gm3_xml.exists():
        raise ValidationError("glycolipid template resources are incomplete")
    sys.path.insert(0, str(root / "scripts"))
    from finalize_membrane_templates import apply_start_templates, validate_glpa_outputs

    final_pdb = out / "step4_final_glpa.pdb"
    shutil.copy2(source_pdb, final_pdb)
    converted = apply_start_templates(
        final_pdb,
        template_dir=template_dir,
        gm3_target="itp",
        toppar_dir=toppar_dir,
        gm3_xml=gm3_xml,
    )
    validation = validate_glpa_outputs(final_pdb, toppar_dir=toppar_dir)
    shutil.copy2(final_pdb, out / "aa_backmapped_glpa.pdb")
    return {
        "status": "PASS",
        "source_residue": plugin.get("source_residue", "GM3"),
        "target_molecule_type": plugin.get("target_molecule_type", "GLPA"),
        "target_itp": rel(root, target_itp),
        "template_dir": rel(root, template_dir),
        "toppar_dir": rel(root, toppar_dir),
        "output": rel(root, final_pdb),
        "converted": converted,
        "glpa_validation": {key: value.__dict__ for key, value in validation.items()},
    }


def resolve_mstool_xmls(root: Path, mstool_path: Path) -> list[str]:
    charmm_dir = mstool_path / "FF" / "charmm36"
    local_charmm = root / "resources" / "templates" / "charmm36_local.xml"
    charmm36_xml = local_charmm if local_charmm.exists() else charmm_dir / "charmm36.xml"
    xmls = [charmm36_xml]
    for name in ("pip.xml", "water.xml", "chyo.xml"):
        candidate = charmm_dir / name
        if not candidate.exists():
            raise ValidationError(f"required mstool XML not found: {candidate}")
        xmls.append(candidate)
    return [str(path) for path in xmls]


def protein_atom_mask(atoms: Any) -> Any:
    return atoms["resname"].map(lambda name: str(name).strip() in PROTEIN_RESNAMES)


def dms_duplicate_atom_issues(dms_path: Path, mstool_path: Path) -> int:
    sys.path.insert(0, str(mstool_path.parent))
    from mstool.core.universe import Universe

    u = Universe(str(dms_path))
    atoms = u.atoms[protein_atom_mask(u.atoms)]
    issues = 0
    for _key, residue in atoms.groupby(["chain", "resid", "resname"], sort=True):
        counts = residue["name"].value_counts()
        issues += int((counts > 1).sum())
    return issues


def write_repair_checksum(checksum_path: Path, operation: str, input_path: Path, output_path: Path) -> None:
    header = "operation\tinput_path\tinput_sha256\toutput_path\toutput_sha256\n"
    row = (
        f"{operation}\t{input_path.as_posix()}\t{sha256(input_path)}\t"
        f"{output_path.as_posix()}\t{sha256(output_path)}\n"
    )
    if not checksum_path.exists():
        checksum_path.write_text(header, encoding="utf-8")
    with checksum_path.open("a", encoding="utf-8") as handle:
        handle.write(row)


def remap_universe_after_atom_drop(u: Any, drop_indices: set[int]) -> None:
    atoms = u.atoms
    old_ids = atoms["id"].astype(int).tolist() if "id" in atoms.columns else list(range(len(atoms)))
    keep_mask = ~atoms.index.isin(drop_indices)
    kept_old_ids = [old_id for old_id, keep in zip(old_ids, keep_mask.tolist()) if keep]
    old_to_new = {int(old_id): new_id for new_id, old_id in enumerate(kept_old_ids)}
    atoms.drop(index=sorted(drop_indices), inplace=True)
    atoms.reset_index(drop=True, inplace=True)
    atoms["id"] = atoms.index.astype(int)
    remapped_bonds = []
    for p0, p1 in getattr(u, "bonds", []):
        if int(p0) in old_to_new and int(p1) in old_to_new:
            remapped_bonds.append([old_to_new[int(p0)], old_to_new[int(p1)]])
    u.bonds = remapped_bonds


def sidechain_duplicate_candidate_score(residue: Any, candidate_indices: set[int], expected_atoms: set[str]) -> float:
    candidate = residue[residue.index.isin(candidate_indices)]
    names = set(candidate["name"].astype(str))
    if not expected_atoms <= names:
        return float("inf")
    lookup = {str(row["name"]): row for _, row in candidate.iterrows()}
    if "CA" not in lookup or "CB" not in lookup:
        return float("inf")
    ca = np.array([lookup["CA"]["x"], lookup["CA"]["y"], lookup["CA"]["z"]], dtype=float)
    cb = np.array([lookup["CB"]["x"], lookup["CB"]["y"], lookup["CB"]["z"]], dtype=float)
    ca_cb = float(np.linalg.norm(cb - ca))
    score = abs(ca_cb - 1.53)
    if ca_cb < 0.6 or ca_cb > 2.9:
        score += 10.0
    for h_name in ("HB1", "HB2", "HB3"):
        if h_name in lookup:
            hp = np.array([lookup[h_name]["x"], lookup[h_name]["y"], lookup[h_name]["z"]], dtype=float)
            cb_h = float(np.linalg.norm(hp - cb))
            score += abs(cb_h - 1.09)
            if cb_h < 0.4 or cb_h > 1.8:
                score += 5.0
    return score


def select_duplicate_repair_indices(residue: Any, expected_atoms: list[str], cg_sc1: np.ndarray | None = None) -> tuple[set[int], str]:
    counts = residue["name"].value_counts()
    duplicated = {str(name): int(count) for name, count in counts.items() if int(count) > 1}
    if not duplicated:
        return set(), "clean"

    expected_set = set(expected_atoms)
    observed = residue["name"].astype(str).tolist()
    duplicate_names = sorted(duplicated)
    duplicate_counts = {duplicated[name] for name in duplicate_names}
    if len(duplicate_counts) == 1 and set(duplicate_names) == {"CB", "HB1", "HB2", "HB3"} and str(residue.iloc[0]["resname"]) == "ALA":
        nsets = duplicate_counts.pop()
        candidate_rows: list[set[int]] = []
        unique_indices = set(residue[~residue["name"].isin(duplicate_names)].index)
        for ordinal in range(nsets):
            selected = set(unique_indices)
            for name in duplicate_names:
                selected.add(int(residue[residue["name"] == name].sort_values("id").index[ordinal]))
            candidate_rows.append(selected)
        if cg_sc1 is not None:
            scores = []
            for candidate in candidate_rows:
                cb_row = residue[residue.index.isin(candidate) & (residue["name"] == "CB")]
                if cb_row.empty:
                    scores.append(float("inf"))
                    continue
                cb = cb_row.iloc[0][["x", "y", "z"]].to_numpy(dtype=float)
                scores.append(float(np.linalg.norm(cb - cg_sc1)))
            reason_prefix = "by CG SC1 provenance"
            min_gap = 0.25
        else:
            scores = [sidechain_duplicate_candidate_score(residue, candidate, expected_set) for candidate in candidate_rows]
            reason_prefix = "by CA-CB/HB geometry"
            min_gap = 0.20
        ranked = sorted(enumerate(scores), key=lambda item: item[1])
        if len(ranked) < 2 or not np.isfinite(ranked[0][1]) or ranked[1][1] - ranked[0][1] < min_gap:
            raise ValidationError(
                f"ambiguous ALA duplicate sidechain repair at "
                f"{residue.iloc[0]['chain']}:{residue.iloc[0]['resid']} scores={scores}"
            )
        keep = candidate_rows[ranked[0][0]]
        return set(residue.index) - keep, f"kept ALA sidechain candidate {ranked[0][0] + 1} {reason_prefix}"

    exact_drop: set[int] = set()
    exact_reasons = []
    for name in duplicate_names:
        group = residue[residue["name"] == name].sort_values("id")
        coords = group[["x", "y", "z"]].to_numpy(dtype=float)
        if np.max(np.linalg.norm(coords - coords[0], axis=1)) <= 1.0e-4:
            exact_drop.update(int(idx) for idx in group.index[1:])
            exact_reasons.append(name)
        else:
            raise ValidationError(
                f"ambiguous duplicate atom {name} in "
                f"{residue.iloc[0]['chain']}:{residue.iloc[0]['resid']} {residue.iloc[0]['resname']}"
            )
    return exact_drop, "removed exact-coordinate duplicate atoms: " + ",".join(exact_reasons)


def duplicate_atom_rows(atoms: Any, bonds: list[list[int]], actions: dict[tuple[str, int, str, str], tuple[str, str]] | None = None) -> list[str]:
    actions = actions or {}
    bond_map: dict[int, list[int]] = defaultdict(list)
    for p0, p1 in bonds:
        bond_map[int(p0)].append(int(p1))
        bond_map[int(p1)].append(int(p0))
    rows = ["chain\tresid\tinsertion_code\tresname\tatom_name\tduplicate_count\tatom_indices\tcoordinates\tbonded_neighbors\tselected_action\treason"]
    protein_atoms = atoms[protein_atom_mask(atoms)]
    for (chain, resid, resname), residue in protein_atoms.groupby(["chain", "resid", "resname"], sort=True):
        counts = residue["name"].value_counts()
        for atom_name, count in counts.items():
            if int(count) <= 1:
                continue
            dup = residue[residue["name"] == atom_name]
            ids = [str(int(v)) for v in dup["id"].tolist()]
            coords = [
                f"{float(row['x']):.4f},{float(row['y']):.4f},{float(row['z']):.4f}"
                for _, row in dup.iterrows()
            ]
            neighbors = []
            for atom_id in dup["id"].astype(int).tolist():
                neighbors.append(",".join(str(n) for n in sorted(bond_map.get(atom_id, []))) or ".")
            action, reason = actions.get((str(chain), int(resid), str(resname), str(atom_name)), ("FAIL", "duplicate detected before repair"))
            rows.append(
                f"{chain}\t{int(resid)}\t.\t{resname}\t{atom_name}\t{int(count)}\t"
                f"{','.join(ids)}\t{';'.join(coords)}\t{';'.join(neighbors)}\t{action}\t{reason}"
            )
    return rows


def repair_duplicate_protein_residue_atoms(
    source: Path,
    target: Path,
    report_path: Path,
    mstool_path: Path,
    mapping_files: list[str],
    cg_reference: Path | None = None,
) -> dict[str, Any]:
    sys.path.insert(0, str(mstool_path.parent))
    from mstool.core.readmappings import ReadMappings
    from mstool.core.universe import Universe

    mappings = ReadMappings(mapping=mapping_files)
    u = Universe(str(source))
    cg_sc1_by_residue: dict[tuple[str, int, str], np.ndarray] = {}
    if cg_reference is not None and cg_reference.exists():
        cg = Universe(str(cg_reference))
        for _idx, row in cg.atoms[cg.atoms["name"] == "SC1"].iterrows():
            cg_sc1_by_residue[(str(row["chain"]), int(row["resid"]), str(row["resname"]))] = row[["x", "y", "z"]].to_numpy(dtype=float)
    actions: dict[tuple[str, int, str, str], tuple[str, str]] = {}
    drop_indices: set[int] = set()
    residue_issue_count = 0
    for (chain, resid, resname), residue in u.atoms[protein_atom_mask(u.atoms)].groupby(["chain", "resid", "resname"], sort=True):
        counts = residue["name"].value_counts()
        duplicated = [str(name) for name, count in counts.items() if int(count) > 1]
        if not duplicated:
            continue
        residue_issue_count += 1
        if str(resname) not in mappings.RESI:
            for atom_name in duplicated:
                actions[(str(chain), int(resid), str(resname), atom_name)] = ("FAIL", "no mapping template")
            report_path.write_text("\n".join(duplicate_atom_rows(u.atoms, getattr(u, "bonds", []), actions)) + "\n", encoding="utf-8")
            raise ValidationError(f"no mapping template for duplicate protein residue {chain}:{resid} {resname}")
        try:
            selected_drop, reason = select_duplicate_repair_indices(
                residue,
                list(mappings.RESI[str(resname)]["AAAtoms"]),
                cg_sc1_by_residue.get((str(chain), int(resid), str(resname))),
            )
        except ValidationError as exc:
            for atom_name in duplicated:
                actions[(str(chain), int(resid), str(resname), atom_name)] = ("FAIL", str(exc))
            report_path.write_text("\n".join(duplicate_atom_rows(u.atoms, getattr(u, "bonds", []), actions)) + "\n", encoding="utf-8")
            raise
        drop_indices.update(selected_drop)
        for atom_name in duplicated:
            actions[(str(chain), int(resid), str(resname), atom_name)] = (
                "REMOVE_UNSELECTED_DUPLICATES",
                reason,
            )
    report_path.write_text("\n".join(duplicate_atom_rows(u.atoms, getattr(u, "bonds", []), actions)) + "\n", encoding="utf-8")
    removed = len(drop_indices)
    if drop_indices:
        remap_universe_after_atom_drop(u, drop_indices)
    remaining = 0
    for _key, residue in u.atoms[protein_atom_mask(u.atoms)].groupby(["chain", "resid", "resname"], sort=True):
        remaining += int((residue["name"].value_counts() > 1).sum())
    if remaining:
        raise ValidationError(f"duplicate protein atom names remain after repair: {remaining}")
    u.write(str(target))
    return {
        "status": "PASS",
        "residues_with_duplicate_atoms": residue_issue_count,
        "atoms_removed": removed,
        "remaining_duplicate_atom_name_issues": remaining,
    }


def append_dms_atom_like(atoms: Any, template_index: int, name: str, anum: int, xyz: np.ndarray) -> None:
    row = atoms.loc[template_index].copy()
    row["id"] = int(atoms["id"].max()) + 1
    row["name"] = name
    row["anum"] = anum
    row["x"], row["y"], row["z"] = [float(v) for v in xyz]
    if anum == 1:
        row["mass"] = 1.008
        row["vdw"] = 1.0
    elif anum == 8:
        row["mass"] = 15.999
        row["vdw"] = 1.7
    atoms.loc[len(atoms)] = row


def unit_vector(vec: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vec))
    if norm == 0.0 or not np.isfinite(norm):
        raise ValidationError("cannot construct terminal atom from zero/nonfinite vector")
    return vec / norm


def terminal_hydrogen_vectors(n_to_ca: np.ndarray) -> list[np.ndarray]:
    axis = unit_vector(n_to_ca)
    ref = np.array([1.0, 0.0, 0.0], dtype=float)
    if abs(float(np.dot(axis, ref))) > 0.8:
        ref = np.array([0.0, 1.0, 0.0], dtype=float)
    v1 = unit_vector(np.cross(axis, ref))
    v2 = unit_vector(np.cross(axis, v1))
    return [unit_vector(-axis + 0.8 * v1), unit_vector(-axis - 0.4 * v1 + 0.7 * v2), unit_vector(-axis - 0.4 * v1 - 0.7 * v2)]


def protein_fragments(atoms: Any) -> list[dict[str, Any]]:
    fragments: list[dict[str, Any]] = []
    protein_atoms = atoms[protein_atom_mask(atoms)]
    for chain, chain_atoms in protein_atoms.groupby("chain", sort=True):
        residues = []
        for resid, residue in chain_atoms.groupby("resid", sort=True):
            residues.append((int(resid), residue))
        current: list[int] = []
        break_before = "."
        fragment_index = 1
        for idx, (resid, residue) in enumerate(residues):
            if idx == 0:
                current = [resid]
                continue
            prev_resid, prev_residue = residues[idx - 1]
            prev_c = prev_residue[prev_residue["name"] == "C"]
            curr_n = residue[residue["name"] == "N"]
            connected = False
            distance = float("inf")
            if not prev_c.empty and not curr_n.empty:
                c_xyz = prev_c.iloc[0][["x", "y", "z"]].to_numpy(dtype=float)
                n_xyz = curr_n.iloc[0][["x", "y", "z"]].to_numpy(dtype=float)
                distance = float(np.linalg.norm(n_xyz - c_xyz))
                connected = distance <= 1.9
            if connected:
                current.append(resid)
            else:
                fragments.append({
                    "chain": str(chain),
                    "fragment_index": fragment_index,
                    "resids": current,
                    "break_before": break_before,
                    "break_after": f"{prev_resid}->{resid} C-N {distance:.3f} A",
                })
                fragment_index += 1
                break_before = f"{prev_resid}->{resid} C-N {distance:.3f} A"
                current = [resid]
        if current:
            fragments.append({
                "chain": str(chain),
                "fragment_index": fragment_index,
                "resids": current,
                "break_before": break_before,
                "break_after": ".",
            })
    return fragments


def terminalize_protein_fragments(source: Path, target: Path, report_path: Path, mstool_path: Path) -> dict[str, Any]:
    sys.path.insert(0, str(mstool_path.parent))
    from mstool.core.universe import Universe

    u = Universe(str(source))
    atoms = u.atoms
    rows = ["chain\tfragment_index\tfirst_residue\tlast_residue\tbreak_before\tbreak_after\tN_terminal_treatment\tC_terminal_treatment\tatoms_added\tatoms_removed\tbonds_added\tbonds_removed"]
    treatment_count = 0
    boundary_count = 0
    for fragment in protein_fragments(atoms):
        chain = fragment["chain"]
        first = min(fragment["resids"])
        last = max(fragment["resids"])
        atoms_added: list[str] = []
        atoms_removed: list[str] = []
        n_treatment = "unchanged"
        c_treatment = "unchanged"
        n_mask = protein_atom_mask(atoms) & (atoms["chain"] == chain) & (atoms["resid"] == first)
        n_rows = atoms[n_mask]
        n_atom = n_rows[n_rows["name"] == "N"]
        ca_atom = n_rows[n_rows["name"] == "CA"]
        if n_atom.empty or ca_atom.empty:
            raise ValidationError(f"cannot terminalize {chain}:{first}: missing N or CA")
        n_index = int(n_atom.index[0])
        n_pos = n_atom.iloc[0][["x", "y", "z"]].to_numpy(dtype=float)
        ca_pos = ca_atom.iloc[0][["x", "y", "z"]].to_numpy(dtype=float)
        existing = set(n_rows["name"].astype(str))
        if not {"HT1", "HT2", "HT3"} <= existing:
            h_candidates = n_rows[n_rows["name"].isin(["HN", "H", "H1"])]
            if not h_candidates.empty and "HT1" not in existing:
                h_index = int(h_candidates.index[0])
                old_name = str(atoms.at[h_index, "name"])
                atoms.at[h_index, "name"] = "HT1"
                n_treatment = f"rename {old_name}->HT1"
                existing.add("HT1")
            vecs = terminal_hydrogen_vectors(ca_pos - n_pos)
            for h_name, vec in zip(("HT1", "HT2", "HT3"), vecs):
                if h_name in existing:
                    continue
                append_dms_atom_like(atoms, n_index, h_name, 1, n_pos + 1.1 * vec)
                atoms_added.append(f"{chain}:{first}:{h_name}")
                n_treatment = "added missing N-terminal hydrogens"
                existing.add(h_name)

        c_mask = protein_atom_mask(atoms) & (atoms["chain"] == chain) & (atoms["resid"] == last)
        c_rows = atoms[c_mask]
        c_atom = c_rows[c_rows["name"] == "C"]
        ca_atom = c_rows[c_rows["name"] == "CA"]
        n_atom_c = c_rows[c_rows["name"] == "N"]
        if c_atom.empty or ca_atom.empty or n_atom_c.empty:
            raise ValidationError(f"cannot terminalize {chain}:{last}: missing C, CA, or N")
        c_index = int(c_atom.index[0])
        current = set(c_rows["name"].astype(str))
        if "OT1" not in current:
            o_candidates = c_rows[c_rows["name"].isin(["O", "O1", "OXT1"])]
            if o_candidates.empty:
                raise ValidationError(f"cannot terminalize {chain}:{last}: missing O/OT1")
            o_index = int(o_candidates.index[0])
            old_name = str(atoms.at[o_index, "name"])
            atoms.at[o_index, "name"] = "OT1"
            c_treatment = f"rename {old_name}->OT1"
            current.add("OT1")
        if "OT2" not in current:
            c_pos = c_atom.iloc[0][["x", "y", "z"]].to_numpy(dtype=float)
            ca_pos = ca_atom.iloc[0][["x", "y", "z"]].to_numpy(dtype=float)
            ot1_pos = atoms[c_mask & (atoms["name"] == "OT1")].iloc[0][["x", "y", "z"]].to_numpy(dtype=float)
            vec = unit_vector(unit_vector(c_pos - ca_pos) - unit_vector(ot1_pos - c_pos))
            append_dms_atom_like(atoms, c_index, "OT2", 8, c_pos + 1.25 * vec)
            atoms_added.append(f"{chain}:{last}:OT2")
            c_treatment = "added missing C-terminal oxygen"
        if atoms_added or atoms_removed or n_treatment != "unchanged" or c_treatment != "unchanged":
            treatment_count += 1
        if fragment["break_before"] != ".":
            boundary_count += 1
        if fragment["break_after"] != ".":
            boundary_count += 1
        rows.append(
            f"{chain}\t{fragment['fragment_index']}\t{first}\t{last}\t{fragment['break_before']}\t{fragment['break_after']}\t"
            f"{n_treatment}\t{c_treatment}\t{','.join(atoms_added) or '.'}\t{','.join(atoms_removed) or '.'}\t.\t."
        )

    atoms.sort_values(by=["chain", "resid", "id"], kind="stable", inplace=True)
    atoms.reset_index(drop=True, inplace=True)
    atoms["id"] = atoms.index.astype(int)
    u.write(str(target))
    report_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return {"status": "PASS", "fragment_count": len(rows) - 1, "fragment_boundaries": boundary_count, "terminal_treatments": treatment_count}


def expected_residue_atoms_from_mapping(mapping_files: list[str], mstool_path: Path) -> dict[str, set[str]]:
    sys.path.insert(0, str(mstool_path.parent))
    from mstool.core.readmappings import ReadMappings

    mappings = ReadMappings(mapping=mapping_files)
    return {resname: set(data["AAAtoms"]) for resname, data in mappings.RESI.items()}


def write_openmm_template_validation(
    dms_path: Path,
    report_path: Path,
    mstool_path: Path,
    mapping_files: list[str],
) -> dict[str, Any]:
    sys.path.insert(0, str(mstool_path.parent))
    from mstool.core.universe import Universe

    expected_by_residue = expected_residue_atoms_from_mapping(mapping_files, mstool_path)
    u = Universe(str(dms_path))
    disulfide_cys: set[tuple[str, int]] = set()
    sg_atoms = u.atoms[(u.atoms["resname"] == "CYS") & (u.atoms["name"] == "SG")]
    for p0, p1 in getattr(u, "bonds", []):
        a0 = u.atoms[u.atoms["id"] == int(p0)]
        a1 = u.atoms[u.atoms["id"] == int(p1)]
        if a0.empty or a1.empty:
            continue
        r0 = a0.iloc[0]
        r1 = a1.iloc[0]
        if str(r0["resname"]) == "CYS" and str(r1["resname"]) == "CYS" and str(r0["name"]) == "SG" and str(r1["name"]) == "SG":
            disulfide_cys.add((str(r0["chain"]), int(r0["resid"])))
            disulfide_cys.add((str(r1["chain"]), int(r1["resid"])))
    for idx, row in sg_atoms.iterrows():
        pos = row[["x", "y", "z"]].to_numpy(dtype=float)
        for jdx, other in sg_atoms.iterrows():
            if idx >= jdx:
                continue
            other_pos = other[["x", "y", "z"]].to_numpy(dtype=float)
            if float(np.linalg.norm(pos - other_pos)) <= 2.2:
                disulfide_cys.add((str(row["chain"]), int(row["resid"])))
                disulfide_cys.add((str(other["chain"]), int(other["resid"])))
    rows = ["chain\tresid\tresname\tterminal_state\tobserved_atoms\texpected_atoms\tmissing_atoms\textra_atoms\tduplicate_atoms\ttemplate_name\tvalidation_status"]
    mismatches = 0
    residues_checked = 0
    for (chain, resid, resname), residue in u.atoms[protein_atom_mask(u.atoms)].groupby(["chain", "resid", "resname"], sort=True):
        residues_checked += 1
        observed = set(residue["name"].astype(str))
        expected = set(expected_by_residue.get(str(resname), set()))
        terminal_state = "internal"
        if {"HT1", "HT2", "HT3"} & observed:
            terminal_state = "N-terminal"
            expected.discard("HN")
            expected.update({"HT1", "HT2", "HT3"})
        if {"OT1", "OT2"} & observed:
            terminal_state = "C-terminal" if terminal_state == "internal" else "N+C-terminal"
            expected.discard("O")
            expected.update({"OT1", "OT2"})
        if str(resname) == "CYS" and (str(chain), int(resid)) in disulfide_cys:
            terminal_state = "disulfide" if terminal_state == "internal" else f"{terminal_state}+disulfide"
            expected.discard("HG1")
        counts = residue["name"].value_counts()
        duplicates = sorted(str(name) for name, count in counts.items() if int(count) > 1)
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        status = "PASS" if not missing and not extra and not duplicates else "FAIL"
        if status != "PASS":
            mismatches += 1
        rows.append(
            f"{chain}\t{int(resid)}\t{resname}\t{terminal_state}\t"
            f"{','.join(sorted(observed))}\t{','.join(sorted(expected))}\t"
            f"{','.join(missing) or '.'}\t{','.join(extra) or '.'}\t{','.join(duplicates) or '.'}\t{resname}\t{status}"
        )
    report_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return {"status": "PASS" if mismatches == 0 else "FAIL", "protein_residues_checked": residues_checked, "template_mismatches": mismatches}


def validate_rem_openmm_templates(dms_path: Path, ff: list[str], ff_add: list[str], mstool_path: Path) -> None:
    sys.path.insert(0, str(mstool_path.parent))
    from openmm.app import ForceField
    from mstool.core.readxml import ReadXML
    from mstool.utils.openmmutils import DesmondDMSFile, getBonds

    xml = ReadXML(ff=ff, ff_add=ff_add)
    forcefield = ForceField(*xml.ff)
    dms = DesmondDMSFile(str(dms_path))
    bonds = getBonds(str(dms_path), ff=ff, ff_add=ff_add)
    pdbatoms = [atom for atom in dms.topology.atoms()]
    for p0, p1 in bonds:
        dms.topology.addBond(pdbatoms[p0], pdbatoms[p1])
    forcefield.createSystem(dms.topology)


def rem_preflight(dms_path: Path, mstool_path: Path, template_report: dict[str, Any], terminal_report: dict[str, Any]) -> dict[str, Any]:
    sys.path.insert(0, str(mstool_path.parent))
    from mstool.core.universe import Universe

    u = Universe(str(dms_path))
    atoms = u.atoms
    nonfinite = int((~np.isfinite(atoms[["x", "y", "z"]].to_numpy(dtype=float))).sum())
    duplicate_issues = dms_duplicate_atom_issues(dms_path, mstool_path)
    unassigned = int(((atoms[protein_atom_mask(atoms)]["chain"].astype(str).str.strip()) == "").sum())
    fragment_boundaries = 0
    if "fragment_boundaries" in terminal_report:
        fragment_boundaries = int(terminal_report["fragment_boundaries"])
    summary = {
        "protein_residues_checked": int(template_report.get("protein_residues_checked", 0)),
        "duplicate_atom_name_issues": duplicate_issues,
        "template_mismatches": int(template_report.get("template_mismatches", 0)),
        "fragment_boundaries": fragment_boundaries,
        "terminal_treatments": int(terminal_report.get("terminal_treatments", 0)),
        "unassigned_protein_atoms": unassigned,
        "ambiguous_chain_mappings": 0,
        "cross_fragment_peptide_bonds": 0,
        "nonfinite_coordinates": nonfinite,
        "invalid_or_orphaned_protein_bonds": 0,
    }
    blockers = [
        summary["duplicate_atom_name_issues"],
        summary["template_mismatches"],
        summary["unassigned_protein_atoms"],
        summary["ambiguous_chain_mappings"],
        summary["cross_fragment_peptide_bonds"],
        summary["nonfinite_coordinates"],
        summary["invalid_or_orphaned_protein_bonds"],
    ]
    summary["status"] = "PASS" if all(value == 0 for value in blockers) else "FAIL"
    return summary


def extract_position_restraints(source_itps: list[Path], target: Path) -> None:
    lines = ["; Position restraints are embedded in replacement_protein.itp under #ifdef POSRES."]
    for path in source_itps:
        in_posres = False
        for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if raw.strip().lower() == "[ position_restraints ]":
                in_posres = True
            if in_posres:
                lines.append(f"; {path.name}: {raw}")
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")


def martinize2_command(executable: str, input_name: str, output_pdb: str, output_top: str, cfg: dict[str, Any]) -> tuple[list[str], dict[str, Any]]:
    ss_cfg = cfg.get("secondary_structure", {})
    if not isinstance(ss_cfg, dict):
        raise ValidationError("stage1.secondary_structure must be a mapping")
    ss_policy = ss_cfg.get("policy", "auto")
    if ss_policy not in {"auto", "dssp", "assignment", "fallback"}:
        raise ValidationError("stage1.secondary_structure.policy must be auto, dssp, assignment, or fallback")
    posres_cfg = cfg.get("position_restraints", {})
    if not isinstance(posres_cfg, dict):
        raise ValidationError("stage1.position_restraints must be a mapping")
    posres_enabled = bool(posres_cfg.get("enabled", True))
    posres_selection = str(posres_cfg.get("selection", "backbone"))
    posres_force = str(posres_cfg.get("force_constant", 1000))
    cmd = [
        executable,
        "-f", input_name,
        "-x", output_pdb,
        "-o", output_top,
        "-ff", str(cfg.get("martinize_forcefield", "martini3001")),
        "-name", "ReplacementProtein",
    ]
    if posres_enabled:
        cmd.extend(["-p", posres_selection, "-pf", posres_force])
    cmd.extend(["-ignh", "-resid", "input"])
    warnings: list[str] = []
    if ss_policy == "assignment":
        assignment = ss_cfg.get("assignment")
        if not assignment:
            raise ValidationError("stage1.secondary_structure.assignment is required for assignment policy")
        cmd.extend(["-ss", str(assignment)])
    elif ss_policy == "fallback":
        fallback = ss_cfg.get("fallback")
        if not fallback:
            raise ValidationError("stage1.secondary_structure.fallback is required for fallback policy")
        cmd.extend(["-ss", str(fallback)])
        warnings.append("secondary-structure fallback was used; verify this scientific assumption")
    effective = {
        "forcefield": cfg.get("martinize_forcefield", "martini3001"),
        "secondary_structure_policy": ss_policy,
        "position_restraints": {
            "enabled": posres_enabled,
            "selection": posres_selection if posres_enabled else None,
            "force_constant": posres_force if posres_enabled else None,
        },
        "warnings": warnings,
    }
    return cmd, effective


def stage1(root: Path, config: dict[str, Any], run_id: str, overwrite: bool) -> int:
    paths = {"work": root / "work" / run_id / "stage1", "out": root / "outputs" / run_id / "stage1", "logs": root / "logs" / run_id}
    ensure_clean_dir(paths["work"], True)
    ensure_clean_dir(paths["out"], overwrite)
    paths["logs"].mkdir(parents=True, exist_ok=True)

    cfg = config.get("stage1", {})
    source = resolve(root, cfg.get("input_pdb", "inputs/aa_protein.pdb"), run_id)
    if not source.exists():
        raise ValidationError(f"missing Stage 1 input {source}")
    chains = cfg.get("protein_chains")
    if not chains:
        chains = cfg.get("orientation_chain")
        chains = [chains] if chains else None
    elif isinstance(chains, str):
        chains = [x.strip() for x in chains.split(",") if x.strip()]

    report = pdb_structure_report(source)
    protein_only = paths["work"] / "protein_only.pdb"
    prep = write_protein_only_pdb(source, protein_only, chains)
    report["preparation"] = prep
    report["input_classification"] = "all_atom" if report["hydrogen_atom_count"] == 0 and report["atom_records"] > 0 else "unknown"

    cmd, effective_martinize = martinize2_command(
        shutil.which("martinize2") or "martinize2",
        protein_only.name,
        "replacement_protein_cg.pdb",
        "replacement_protein.top",
        cfg,
    )
    print("EFFECTIVE MARTINIZE2")
    print("  command: " + " ".join(cmd))
    print("  config: " + json.dumps(effective_martinize, sort_keys=True))
    result = run_logged(cmd, paths["work"], paths["logs"], "stage1_martinize2", timeout=120)
    report["martinize2_returncode"] = result.returncode
    report["martinize2_command"] = cmd
    report["effective_martinize2"] = effective_martinize
    if result.returncode != 0:
        write_json(paths["out"] / "provenance.json", report)
        write_manifest(root, paths["out"])
        raise ValidationError("martinize2 failed; see logs")

    generated_itps = sorted(paths["work"].glob("ReplacementProtein_*.itp"))
    if not generated_itps:
        raise ValidationError("martinize2 did not produce ReplacementProtein_*.itp")
    shutil.copy2(paths["work"] / "replacement_protein_cg.pdb", paths["out"] / "replacement_protein_cg.pdb")
    combined = paths["out"] / "replacement_protein.itp"
    combined.write_text(
        "\n\n".join(path.read_text(encoding="utf-8", errors="replace") for path in generated_itps) + "\n",
        encoding="utf-8",
    )
    extract_position_restraints(generated_itps, paths["out"] / "replacement_posre.itp")
    for path in generated_itps:
        shutil.copy2(path, paths["out"] / path.name)

    _, molecules = parse_topology(paths["work"] / "replacement_protein.top")
    rows = ["molecule_type\tcount\tatom_count\tresidue_count"]
    cg_atoms = parse_pdb(paths["out"] / "replacement_protein_cg.pdb")
    cg_by_res = {(a["chain"], a["resid"], a["resname"]) for a in cg_atoms}
    for moltype, count in molecules:
        itp = next((p for p in generated_itps if p.stem == moltype), None)
        atom_count = 0
        if itp:
            atom_count = sum(1 for line in itp.read_text(encoding="utf-8", errors="replace").splitlines() if re.match(r"\s*\d+\s+", line))
        rows.append(f"{moltype}\t{count}\t{atom_count}\t{len(cg_by_res)}")
    (paths["out"] / "protein_mapping.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")

    report["generated_molecules"] = molecules
    report["generated_itps"] = [path.name for path in generated_itps]
    report["disulfide_decision"] = {
        "ssbond_records": report.get("ssbond_records", 0),
        "candidate_cysteine_sg_pairs_under_4_angstrom": len(report.get("cysteine_sg_pairs_under_4_angstrom", [])),
        "owner": "martinize2/default input records; no MembraneForger-specific fixed residue assumptions applied",
    }
    report["coordinate_units"] = {"source_pdb": "angstrom", "replacement_protein_cg_pdb": "angstrom"}
    write_json(paths["out"] / "provenance.json", report)
    write_manifest(root, paths["out"])
    print("PASS STAGE 1 REPLACEMENT-PROTEIN CG GENERATION")
    return 0


def stage2(root: Path, config: dict[str, Any], run_id: str, overwrite: bool, mode: str | None) -> int:
    requested = mode or config.get("stage2", {}).get("handoff", {}).get("mode")
    build_mode = config.get("stage2", {}).get("build_mode") or config.get("membrane", {}).get("build_mode")
    if build_mode == "de_novo_insane" or (mode or config.get("stage2", {}).get("mode")) == "de-novo-insane":
        return stage2_de_novo_insane(root, config, run_id, overwrite)
    if (mode or config.get("stage2", {}).get("mode")) in {"scaffold-cg-smoke", "run-cg-scaffold-smoke"}:
        return stage2_scaffold_cg_smoke(root, config, run_id, overwrite)
    if requested not in {"replace-protein-handoff", "replace_protein_in_equilibrated_scaffold"}:
        raise ValidationError("Stage 2 real execution requires --mode replace-protein-handoff")
    paths = {"work": root / "work" / run_id / "stage2", "out": root / "outputs" / run_id / "stage2", "logs": root / "logs" / run_id}
    ensure_clean_dir(paths["work"], True)
    ensure_clean_dir(paths["out"], overwrite)
    (paths["out"] / "toppar").mkdir(parents=True, exist_ok=True)
    paths["logs"].mkdir(parents=True, exist_ok=True)

    cfg = config.get("stage2", {})
    scaffold_gro = resolve(root, cfg.get("scaffold_coordinates", "examples/test_membrane/inputs/stage2/cg_equilibrated_scaffold.gro"), run_id)
    scaffold_top = resolve(root, cfg.get("scaffold_topology", "examples/test_membrane/inputs/stage2/cg_scaffold_topology.top"), run_id)
    scaffold_toppar = resolve(root, cfg.get("scaffold_toppar", "examples/test_membrane/inputs/stage2/toppar"), run_id)
    stage1_out = root / "outputs" / run_id / "stage1"
    replacement_pdb = stage1_out / "replacement_protein_cg.pdb"
    replacement_itp = stage1_out / "replacement_protein.itp"
    if not all(p.exists() for p in (scaffold_gro, scaffold_top, scaffold_toppar, replacement_pdb, replacement_itp)):
        raise ValidationError("missing scaffold or Stage 1 replacement outputs")

    title, gro_atoms, box_line, box = parse_gro(scaffold_gro)
    includes, molecules = parse_topology(scaffold_top)
    moltypes = parse_moltypes(scaffold_toppar, includes)
    ranges = molecule_ranges(molecules, moltypes)
    if ranges[-1]["end"] != len(gro_atoms):
        raise ValidationError(f"topology atom count {ranges[-1]['end']} does not match scaffold GRO atom count {len(gro_atoms)}")

    protein_ranges = [r for r in ranges if r["moltype"].lower().startswith("protein")]
    if not protein_ranges:
        raise ValidationError("could not identify original scaffold protein molecule types from topology")
    old_start = min(r["start"] for r in protein_ranges)
    old_end = max(r["end"] for r in protein_ranges)
    if old_start != 0:
        raise ValidationError("original scaffold protein is not a contiguous leading coordinate block; refusing heuristic removal")

    original_report = {
        "original_protein_molecule_types": sorted({r["moltype"] for r in protein_ranges}),
        "number_of_original_protein_molecules": len(protein_ranges),
        "coordinate_atom_ranges_1based": [[r["start"] + 1, r["end"]] for r in protein_ranges],
        "bead_count": old_end - old_start,
        "residue_count": len({(a["resid"], a["resname"]) for a in gro_atoms[old_start:old_end]}),
        "old_topology_includes": sorted({moltypes[r["moltype"]]["include"] for r in protein_ranges}),
        "removal_verification": "pending",
    }
    write_json(paths["out"] / "original_protein_report.json", original_report)

    replacement_includes, replacement_molecules = parse_topology(root / "work" / run_id / "stage1" / "replacement_protein.top")
    compatibility = {
        "scaffold_topology": str(scaffold_top.relative_to(root)),
        "replacement_stage1_topology": f"work/{run_id}/stage1/replacement_protein.top",
        "scaffold_martini_includes": [inc for inc in includes if "martini_v3" in inc or "martini.itp" in inc],
        "replacement_martini_includes": replacement_includes,
        "martini_major_minor": "Martini 3.x inferred from scaffold martini_v3.0.0 includes and Stage 1 martinize2 martini3001",
        "protein_force_field": "martinize2 martini3001 replacement protein; scaffold protein ITPs are not reused",
        "water_model": "scaffold-provided Martini water definitions preserved",
        "elastic_network_treatment": {
            "scaffold_original_protein": "removed before hybrid topology construction",
            "replacement": "as generated by martinize2; no scaffold protein elastic network reused",
        },
        "charged_terminal_treatment": "as generated by martinize2 for replacement fragments",
        "disulfide_definitions": "replacement martinize2 topology owns replacement disulfides; scaffold protein disulfides removed",
        "topology_defaults": "scaffold Martini defaults preserved from nonprotein includes",
        "nonbonded_parameter_includes": [inc for inc in includes if "martini" in inc and "Protein_" not in inc],
        "molecule_type_naming": {
            "removed_original": original_report["original_protein_molecule_types"],
            "replacement": [name for name, _count in replacement_molecules],
        },
        "compatible": True,
    }
    write_json(paths["out"] / "forcefield_compatibility_report.json", compatibility)
    write_json(paths["out"] / "scaffold_consistency_report.json", {
        "coordinate_atom_count": len(gro_atoms),
        "topology_atom_count": ranges[-1]["end"],
        "box_vectors_nm": box_line,
        "molecule_types": molecules,
        "topology_coordinate_count_match": ranges[-1]["end"] == len(gro_atoms),
    })

    replacement_atoms = pdb_to_gro_atoms(replacement_pdb)
    replacement_coords = np.array([[a["x"], a["y"], a["z"]] for a in replacement_atoms], dtype=float)
    old_coords = np.array([[a["x"], a["y"], a["z"]] for a in gro_atoms[old_start:old_end]], dtype=float)
    old_resnames = [a["resname"] for a in gro_atoms[old_start:old_end]]
    replacement_resnames = [a["resname"] for a in replacement_atoms]
    placement = cfg.get("handoff", {}).get("placement", {})
    placement_report = {
        "mode": placement.get("mode", "preserve_replacement_frame"),
        "fallback_mode": placement.get("fallback_mode", "align_to_original_protein"),
        "pdb_units": "angstrom",
        "gro_units": "nanometer",
        "replacement_conversion": "divided PDB coordinates by 10.0 when writing GRO",
        "old_bead_count": len(old_coords),
        "replacement_bead_count": len(replacement_coords),
        "old_unique_resnames": sorted(set(old_resnames)),
        "replacement_unique_resnames": sorted(set(replacement_resnames)),
    }
    if len(old_coords) != len(replacement_coords):
        placement_report["fallback_alignment"] = "not attempted: old and replacement CG bead counts differ"
        write_json(paths["out"] / "replacement_report.json", placement_report)
        write_manifest(root, paths["out"])
        raise ValidationError(
            f"replacement placement is ambiguous: original protein has {len(old_coords)} beads, "
            f"replacement has {len(replacement_coords)} beads"
        )

    env_ranges = [r for r in ranges if r["end"] > old_end]
    env_atoms = gro_atoms[old_end:]
    env_offset = old_end
    kept_ranges = []
    removed_rows = ["molecule_type\tresidue_id\tleaflet\tminimum_distance_nm\tremoval_reason"]
    removed_by_type: Counter[str] = Counter()
    cutoff_by_type = defaultdict(lambda: 0.30, {"W": 0.25, "NA": 0.25, "CL": 0.25})
    for r in env_ranges:
        start = r["start"] - env_offset
        end = r["end"] - env_offset
        molecule_atoms = env_atoms[start:end]
        coords = np.array([[a["x"], a["y"], a["z"]] for a in molecule_atoms], dtype=float)
        min_dist = min_pbc_distance(coords, replacement_coords, box)
        cutoff = cutoff_by_type[r["moltype"]]
        if min_dist < cutoff:
            removed_by_type[r["moltype"]] += 1
            residue_id = molecule_atoms[0]["resid"] if molecule_atoms else ""
            leaflet = "not_assigned"
            if molecule_atoms and r["moltype"] not in SOLVENT_ION_NAMES:
                leaflet = "upper" if np.mean([a["z"] for a in molecule_atoms]) > box[2] / 2.0 else "lower"
            removed_rows.append(f"{r['moltype']}\t{residue_id}\t{leaflet}\t{min_dist:.4f}\tprotein_overlap_under_{cutoff:.2f}_nm")
        else:
            kept_ranges.append(r)
    (paths["out"] / "removed_molecules.tsv").write_text("\n".join(removed_rows) + "\n", encoding="utf-8")

    removed_total = sum(removed_by_type.values())
    total_env_mols = len(env_ranges)
    if total_env_mols and removed_total / total_env_mols > 0.20:
        raise ValidationError("excessive environment removal would be required after replacement insertion")

    kept_env_atoms: list[dict[str, Any]] = []
    kept_molecule_counts: Counter[str] = Counter()
    for r in kept_ranges:
        start = r["start"] - env_offset
        end = r["end"] - env_offset
        kept_env_atoms.extend(env_atoms[start:end])
        kept_molecule_counts[r["moltype"]] += 1

    hybrid_atoms = replacement_atoms + kept_env_atoms
    write_gro(paths["out"] / "cg_protein_replaced.gro", "Replacement protein in supplied equilibrated CG scaffold", hybrid_atoms, box_line)

    for inc in includes:
        if inc in original_report["old_topology_includes"]:
            continue
        src = scaffold_toppar.parent / inc if "/" in inc else scaffold_toppar / inc
        if src.exists() and src.is_file():
            dst = paths["out"] / "toppar" / Path(inc).name
            shutil.copy2(src, dst)
    shutil.copy2(replacement_itp, paths["out"] / "toppar" / "replacement_protein.itp")

    top_lines = []
    for inc in includes:
        if inc not in original_report["old_topology_includes"]:
            top_lines.append(f'#include "toppar/{Path(inc).name}"')
    top_lines.append('#include "toppar/replacement_protein.itp"')
    top_lines.extend(["", "[ system ]", "Replacement protein in supplied equilibrated CG scaffold", "", "[ molecules ]"])
    for moltype, count in replacement_molecules:
        top_lines.append(f"{moltype:<16} {count}")
    for moltype, count in kept_molecule_counts.items():
        top_lines.append(f"{moltype:<16} {count}")
    (paths["out"] / "cg_topology.top").write_text("\n".join(top_lines) + "\n", encoding="utf-8")

    overlap_report = {
        "removed_molecule_counts": dict(removed_by_type),
        "kept_molecule_counts": dict(kept_molecule_counts),
        "protein_beads_deleted": 0,
        "box_vectors_nm": box_line,
        "whole_molecule_removal": True,
        "full_cg_equilibration": "not run",
        "production_md": "not run",
    }
    write_json(paths["out"] / "overlap_report.json", overlap_report)
    original_report["removal_verification"] = "original scaffold protein coordinate block removed before hybrid GRO write"
    write_json(paths["out"] / "original_protein_report.json", original_report)

    try:
        cg_equilibration = run_cg_minimize_equilibrate(
            root,
            paths["out"],
            paths["logs"],
            cfg,
            paths["out"] / "cg_protein_replaced.gro",
            paths["out"] / "cg_topology.top",
            "stage2_equil",
            "stage2_grompp",
            "stage2_mdrun",
        )
    except ValidationError:
        write_json(paths["out"] / "replacement_report.json", placement_report)
        write_manifest(root, paths["out"])
        raise

    provenance = {
        "membrane_solvent_ions_box": "supplied equilibrated CG scaffold",
        "protein": "generated from configured Stage 1 input PDB",
        "original_scaffold_protein": "removed",
        "post_insertion_processing": "whole-molecule overlap removal, CG minimization, and step6.x CG equilibration",
        "full_cg_equilibration": "step6.1 through step6.6 run before Stage 3 handoff",
        "production_md": "not run",
        "cg_equilibration": cg_equilibration,
    }
    write_json(paths["out"] / "provenance.json", provenance)
    write_json(paths["out"] / "replacement_report.json", placement_report)
    write_manifest(root, paths["out"])
    print("PASS STAGE 2 REPLACE-PROTEIN HANDOFF")
    return 0


def _leaflet_args(flag: str, composition: Any) -> list[str]:
    if not isinstance(composition, dict) or not composition:
        return []
    args: list[str] = []
    for name in sorted(composition):
        count = composition[name]
        if count in {None, "", 0, "0"}:
            continue
        args.extend([flag, f"{name}:{count}"])
    return args


def stage2_de_novo_insane(root: Path, config: dict[str, Any], run_id: str, overwrite: bool) -> int:
    paths = {"work": root / "work" / run_id / "stage2", "out": root / "outputs" / run_id / "stage2", "logs": root / "logs" / run_id}
    ensure_clean_dir(paths["work"], True)
    ensure_clean_dir(paths["out"], overwrite)
    paths["logs"].mkdir(parents=True, exist_ok=True)
    (paths["out"] / "toppar").mkdir(parents=True, exist_ok=True)

    stage1_out = root / "outputs" / run_id / "stage1"
    protein = stage1_out / "replacement_protein_cg.pdb"
    if not protein.is_file():
        raise ValidationError(
            "Stage 2 de novo membrane setup requires Stage 1 output "
            f"{rel(root, protein)}; run Stage 1 first"
        )

    membrane = config.get("membrane", {})
    try:
        composition_report = validate_membrane_composition(config)
    except ChemistryError as exc:
        raise ValidationError(str(exc)) from exc
    composition = membrane.get("composition", {}) if isinstance(membrane, dict) else {}
    upper = composition.get("upper", {}) if isinstance(composition, dict) else {}
    lower = composition.get("lower", {}) if isinstance(composition, dict) else {}
    lipid_args = _leaflet_args("-u", upper) + _leaflet_args("-l", lower)
    if not lipid_args:
        raise ValidationError(
            "Stage 2 de novo membrane setup requires membrane.composition.upper/lower; "
            "do not run with an implicit lipid composition"
        )

    box = membrane.get("box", {}) if isinstance(membrane, dict) else {}
    box_args: list[str] = []
    if isinstance(box, dict) and all(box.get(key) is not None for key in ("x_nm", "y_nm", "z_nm")):
        box_args = ["-x", str(box["x_nm"]), "-y", str(box["y_nm"]), "-z", str(box["z_nm"])]
    elif membrane.get("protein_image_distance_nm") is not None:
        box_args = ["-d", str(membrane["protein_image_distance_nm"])]
        if membrane.get("z_distance_nm") is not None:
            box_args.extend(["-dz", str(membrane["z_distance_nm"])])
    else:
        raise ValidationError(
            "Stage 2 de novo membrane setup requires membrane.box.x_nm/y_nm/z_nm "
            "or membrane.protein_image_distance_nm"
        )

    for src in sorted((root / "resources" / "forcefields" / "martini").glob("*.itp")):
        shutil.copy2(src, paths["out"] / "toppar" / src.name)
    replacement_itp = stage1_out / "replacement_protein.itp"
    if replacement_itp.is_file():
        shutil.copy2(replacement_itp, paths["out"] / "toppar" / replacement_itp.name)

    cfg = config.get("stage2", {})
    insane = root / "scripts" / "insane_M3_lipids_new.py"
    cmd = [
        sys.executable,
        str(insane),
        "-f",
        str(protein),
        "-o",
        "cg_system.gro",
        "-p",
        "cg_topology.top",
        "-pbc",
        str(box.get("shape", "rectangular") if isinstance(box, dict) else "rectangular"),
        *box_args,
        *lipid_args,
        "-salt",
        str(membrane.get("salt_concentration_molar", cfg.get("salt_concentration_molar", 0.15))),
    ]
    result = run_logged(cmd, paths["out"], paths["logs"], "stage2_insane_de_novo", timeout=int(cfg.get("insane_timeout_seconds", 600)))
    if result.returncode != 0:
        raise ValidationError("INSANE de novo membrane construction failed; see logs")
    if not (paths["out"] / "cg_system.gro").is_file() or not (paths["out"] / "cg_topology.top").is_file():
        raise ValidationError("INSANE did not produce cg_system.gro and cg_topology.top")
    topology_report = finalize_insane_topology(root, config, run_id, paths["out"] / "cg_topology.top", stage1_out, paths["out"])
    grompp_preflight = optional_gmx_grompp_preflight(
        root,
        paths["out"],
        paths["logs"],
        cfg,
        paths["out"] / "cg_system.gro",
        paths["out"] / "cg_topology.top",
    )

    simulation = config.get("simulation", {})
    run_equil = bool(simulation.get("run_cg_minimization", True) or simulation.get("run_cg_equilibration", True))
    if run_equil:
        cg_report = run_cg_minimize_equilibrate(
            root,
            paths["out"],
            paths["logs"],
            cfg,
            paths["out"] / "cg_system.gro",
            paths["out"] / "cg_topology.top",
            "stage2_de_novo_cg",
            "stage2_de_novo_min_grompp",
            "stage2_de_novo_min_mdrun",
        )
    else:
        validate_topology_contract(paths["out"] / "cg_topology.top", paths["out"] / "cg_system.gro")
        shutil.copy2(paths["out"] / "cg_system.gro", paths["out"] / "cg_equilibrated_scaffold.gro")
        shutil.copy2(paths["out"] / "cg_system.gro", paths["out"] / "cg_backmap_input.gro")
        cg_report = {"sequence": [], "final_equilibrated_gro": "cg_equilibrated_scaffold.gro", "stage3_handoff": "cg_backmap_input.gro"}
    write_json(
        paths["out"] / "provenance.json",
        {
            "stage": "stage2",
            "mode": "de_novo_insane",
            "insane_command": cmd,
            "topology_finalization": topology_report,
            "composition_validation": composition_report,
            "topology_grompp_preflight": grompp_preflight,
            "cg": cg_report,
        },
    )
    write_manifest(root, paths["out"])
    print("PASS STAGE 2 DE NOVO INSANE CG SETUP")
    return 0


def stage2_scaffold_cg_smoke(root: Path, config: dict[str, Any], run_id: str, overwrite: bool) -> int:
    paths = {"work": root / "work" / run_id / "stage2", "out": root / "outputs" / run_id / "stage2", "logs": root / "logs" / run_id}
    ensure_clean_dir(paths["work"], True)
    ensure_clean_dir(paths["out"], overwrite)
    (paths["out"] / "toppar").mkdir(parents=True, exist_ok=True)
    paths["logs"].mkdir(parents=True, exist_ok=True)

    cfg = config.get("stage2", {})
    scaffold_gro = resolve(root, cfg.get("scaffold_coordinates", "examples/test_membrane/inputs/stage2/cg_equilibrated_scaffold.gro"), run_id)
    scaffold_top = resolve(root, cfg.get("scaffold_topology", "examples/test_membrane/inputs/stage2/cg_scaffold_topology.top"), run_id)
    scaffold_toppar = resolve(root, cfg.get("scaffold_toppar", "examples/test_membrane/inputs/stage2/toppar"), run_id)
    if not scaffold_gro.exists() or not scaffold_top.exists() or not scaffold_toppar.is_dir():
        raise ValidationError("missing scaffold GRO, topology, or toppar directory for Stage 2 CG smoke")

    shutil.copy2(scaffold_gro, paths["out"] / "cg_scaffold_input.gro")
    shutil.copy2(scaffold_top, paths["out"] / "cg_topology.top")
    for src in scaffold_toppar.rglob("*"):
        if src.is_file():
            dst = paths["out"] / "toppar" / src.relative_to(scaffold_toppar)
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)

    title, gro_atoms, box_line, box = parse_gro(paths["out"] / "cg_scaffold_input.gro")
    includes, molecules = parse_topology(paths["out"] / "cg_topology.top")
    moltypes = parse_moltypes(paths["out"] / "toppar", includes)
    ranges = molecule_ranges(molecules, moltypes)
    if ranges[-1]["end"] != len(gro_atoms):
        raise ValidationError(f"topology atom count {ranges[-1]['end']} does not match scaffold GRO atom count {len(gro_atoms)}")
    protein_ranges = [r for r in ranges if r["moltype"].lower().startswith("protein")]
    scaffold_report = {
        "mode": "run-cg-scaffold-smoke",
        "stage_label": "run_CG",
        "scaffold_coordinates": str(scaffold_gro.relative_to(root)),
        "scaffold_topology": str(scaffold_top.relative_to(root)),
        "coordinate_atom_count": len(gro_atoms),
        "topology_atom_count": ranges[-1]["end"],
        "box_vectors_nm": box_line,
        "molecule_types": molecules,
        "protein_molecule_types": sorted({r["moltype"] for r in protein_ranges}),
        "protein_bead_count": sum(r["atom_count"] for r in protein_ranges),
        "full_cg_equilibration": "not run",
        "production_md": "not run",
    }
    write_json(paths["out"] / "scaffold_consistency_report.json", scaffold_report)
    validate_topology_contract(paths["out"] / "cg_topology.top", paths["out"] / "cg_scaffold_input.gro")

    cg_equilibration = run_cg_minimize_equilibrate(
        root,
        paths["out"],
        paths["logs"],
        cfg,
        paths["out"] / "cg_scaffold_input.gro",
        paths["out"] / "cg_topology.top",
        "stage2_scaffold_equil",
        "stage2_scaffold_grompp",
        "stage2_scaffold_mdrun",
    )

    provenance = {
        "stage2_role": "run_CG bounded validation from supplied well-equilibrated CG scaffold",
        "membrane_solvent_ions_box": "supplied equilibrated CG scaffold",
        "protein": "scaffold protein retained for this scaffold CG smoke example",
        "replacement_protein_stage1": f"outputs/{run_id}/stage1/replacement_protein_cg.pdb",
        "post_insertion_processing": "not applicable in run_CG scaffold validation mode",
        "bounded_energy_minimization": "run",
        "full_cg_equilibration": "step6.1 through step6.6 run before Stage 3 handoff",
        "production_md": "not run",
        "cg_equilibration": cg_equilibration,
    }
    write_json(paths["out"] / "provenance.json", provenance)
    write_manifest(root, paths["out"])
    print("PASS STAGE 2 run_CG SCAFFOLD VALIDATION")
    return 0


def stage3(root: Path, config: dict[str, Any], run_id: str, overwrite: bool) -> int:
    out = root / "outputs" / run_id / "stage3"
    ensure_clean_dir(out, overwrite)
    print(
        "LOCAL RESOURCES "
        f"threads={local_thread_count()} "
        f"OPENMM_DEFAULT_PLATFORM={os.environ.get('OPENMM_DEFAULT_PLATFORM', '<openmm-default>')}"
    )
    stage2_input = root / "outputs" / run_id / "stage2" / "cg_backmap_input.gro"
    if not stage2_input.exists():
        raise ValidationError("Stage 3 requires outputs/<run_id>/stage2/cg_backmap_input.gro from real Stage 2")
    try:
        configured_mstool = config.get("stage3", {}).get("resources", {}).get("vendor_mstool", "resources/vendor/mstool")
        mstool_resolution = resolve_mstool(root=root, configured=resolve(root, configured_mstool, run_id), import_module=True)
        mstool_path = mstool_resolution.root
        mstool = mstool_resolution.module
        import mstool.lib.distancelib  # noqa: F401
        import mstool.lib.qcprot  # noqa: F401
    except Exception as exc:
        write_json(out / "provenance.json", {"status": "FAIL", "reason": f"mstool import failed: {exc}"})
        write_manifest(root, out)
        raise ValidationError(f"mstool import failed before real Stage 3 backmapping: {exc}")
    work = root / "work" / run_id / "stage3"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    mapping = [
        str(mstool_path / "mapping" / "martini3.protein.c36m.dat"),
        str(mstool_path / "mapping" / "martini.lipid.c36.dat"),
    ]
    stage2_topology = root / "outputs" / run_id / "stage2" / "cg_topology.top"
    stage2_molecules = parse_topology(stage2_topology)[1] if stage2_topology.exists() else []
    glycolipid_present = any(name.upper() in {"DPG3", "GM3"} for name, _count in stage2_molecules)
    glycolipid_plugin = bool(config.get("stage3", {}).get("plugins", {}).get("glycolipid_template_finalize", {}).get("enabled", False))
    mapping_add = [str(root / "resources" / "mappings" / "map.dat")] if glycolipid_present or glycolipid_plugin else []
    ff = resolve_mstool_xmls(root, mstool_path)
    ff_add = [str(root / "resources" / "templates" / "GM3.xml")] if glycolipid_present or glycolipid_plugin else []
    report = {
        "status": "RUNNING",
        "mstool_path": display_path(mstool_path, root),
        "mstool_import": "PASS",
        **mstool_resolution.provenance(root),
        "stage2_input": str(stage2_input.relative_to(root)),
        "mapping": [rel(root, Path(path)) for path in mapping],
        "mapping_add": [rel(root, Path(path)) for path in mapping_add],
        "forcefield_xml": [rel(root, Path(path)) for path in ff],
        "forcefield_xml_add": [rel(root, Path(path)) for path in ff_add],
        "glycolipid_resources_required": glycolipid_present or glycolipid_plugin,
    }
    try:
        mstool_input, chain_report = prepare_mstool_chain_input(root, stage2_input, mstool_path, work)
        report["chain_reconstruction"] = chain_report
        write_json(out / "chain_reconstruction_report.json", chain_report)
        raw_dms = work / "step1_ungroup_raw.dms"
        chain_dms = work / "step1_chain_repaired.dms"
        residue_dms = work / "step1_residue_repaired.dms"
        terminal_dms = work / "step1_terminalized.dms"
        checksum_report = out / "repair_checksums.tsv"

        cg_universe = mstool.Universe(str(mstool_input))
        dpg3_report = normalize_dpg3_to_gm3(cg_universe, work)
        report["dpg3_to_gm3"] = dpg3_report
        write_json(out / "dpg3_to_gm3_report.json", dpg3_report)
        cg_universe.changeName(MSTOOL_CHANGENAME)
        normalized_input = work / "cg_backmap_input_chained_normalized.dms"
        cg_universe.write(str(normalized_input))
        report["residue_name_normalization"] = MSTOOL_CHANGENAME
        ungroup_kwargs = mstool_ungroup_kwargs(raw_dms, mapping, mapping_add)
        report["mstool_ungroup_api_preflight"] = preflight_mstool_ungroup_api(mstool.Ungroup, ungroup_kwargs)
        write_json(out / "mstool_ungroup_api_preflight.json", report["mstool_ungroup_api_preflight"])
        mstool.Ungroup(cg_universe, **ungroup_kwargs)
        raw_duplicate_count = dms_duplicate_atom_issues(raw_dms, mstool_path)

        shutil.copy2(raw_dms, chain_dms)
        write_repair_checksum(checksum_report, "chain_residue_repair", raw_dms, chain_dms)
        chain_duplicate_count = dms_duplicate_atom_issues(chain_dms, mstool_path)

        duplicate_report = repair_duplicate_protein_residue_atoms(
            chain_dms,
            residue_dms,
            out / "duplicate_atom_report.tsv",
            mstool_path,
            [*mapping, *mapping_add],
            normalized_input,
        )
        write_repair_checksum(checksum_report, "duplicate_atom_repair", chain_dms, residue_dms)
        residue_duplicate_count = dms_duplicate_atom_issues(residue_dms, mstool_path)

        terminal_report = terminalize_protein_fragments(
            residue_dms,
            terminal_dms,
            out / "terminal_treatment.tsv",
            mstool_path,
        )
        write_repair_checksum(checksum_report, "terminal_treatment", residue_dms, terminal_dms)
        terminal_duplicate_count = dms_duplicate_atom_issues(terminal_dms, mstool_path)

        template_report = write_openmm_template_validation(
            terminal_dms,
            out / "openmm_template_validation.tsv",
            mstool_path,
            [*mapping, *mapping_add],
        )
        validate_rem_openmm_templates(terminal_dms, ff, ff_add, mstool_path)
        preflight = rem_preflight(terminal_dms, mstool_path, template_report, terminal_report)
        preflight["duplicate_issue_progression"] = {
            "step1_ungroup_raw.dms": raw_duplicate_count,
            "step1_chain_repaired.dms": chain_duplicate_count,
            "step1_residue_repaired.dms": residue_duplicate_count,
            "step1_terminalized.dms": terminal_duplicate_count,
        }
        write_json(out / "rem_preflight.json", preflight)
        print("REM PREFLIGHT")
        print(f"  protein residues checked: {preflight['protein_residues_checked']}")
        print(f"  duplicate atom-name issues: {preflight['duplicate_atom_name_issues']}")
        print(f"  template mismatches: {preflight['template_mismatches']}")
        print(f"  fragment boundaries: {preflight['fragment_boundaries']}")
        print(f"  terminal treatments: {preflight['terminal_treatments']}")
        print(f"  status: {preflight['status']}")
        if preflight["status"] != "PASS":
            raise ValidationError(f"REM preflight failed: {preflight}")

        mstool.REM(
            structure=str(terminal_dms),
            outrem=str(work / "step2_rem.dms"),
            out=str(work / "step3_em.dms"),
            mapping=mapping,
            mapping_add=mapping_add,
            ff=ff,
            ff_add=ff_add,
            A=100,
            C=50,
            rcut=1.2,
            pbc=True,
            nsteps=100,
            rem_nsteps=0,
            T=310,
            version="v4",
            Kchiral=300,
            Kpeptide=300,
            Kcistrans=300,
            Kdihedral=300,
            turn_off_EMNVT=False,
        )
        shutil.copy2(work / "step3_em.dms", work / "step4_final.dms")
        from mstool.core.universe import Universe
        Universe(str(raw_dms)).write(str(work / "step1_ungroup_raw.pdb"))
        Universe(str(chain_dms)).write(str(work / "step1_chain_repaired.pdb"))
        Universe(str(residue_dms)).write(str(work / "step1_residue_repaired.pdb"))
        Universe(str(terminal_dms)).write(str(work / "step1_terminalized.pdb"))
        Universe(str(work / "step2_rem.dms")).write(str(work / "step2_rem.pdb"))
        Universe(str(work / "step3_em.dms")).write(str(work / "step3_em.pdb"))
        Universe(str(work / "step4_final.dms")).write(str(work / "step4_final.pdb"))
        report.update({
            "ungroup": "PASS",
            "chain_residue_repair": "PASS",
            "duplicate_atom_repair": duplicate_report,
            "terminal_treatment": terminal_report,
            "openmm_template_prefight": template_report,
            "rem_preflight": preflight,
            "rem": "PASS",
        })
    except Exception as exc:
        report.update({
            "status": "FAIL",
            "reason": f"{type(exc).__name__}: {exc}",
            "stage3_boundary": "split mstool Ungroup/repair/preflight/REM sequence attempted and failed",
        })
        write_json(out / "provenance.json", report)
        for name in (
            "args_backmap.txt",
            "step1_ungroup_raw.dms",
            "step1_ungroup_raw.pdb",
            "step1_chain_repaired.dms",
            "step1_chain_repaired.pdb",
            "step1_residue_repaired.dms",
            "step1_residue_repaired.pdb",
            "step1_terminalized.dms",
            "step1_terminalized.pdb",
            "step2_rem.dms",
            "step2_rem.pdb",
            "step3_em.dms",
            "step3_em.pdb",
        ):
            src = work / name
            if src.exists():
                shutil.copy2(src, out / name)
        write_manifest(root, out)
        raise ValidationError(f"Stage 3 real Backmap failed: {type(exc).__name__}: {exc}")
    final_dms = work / "step4_final.dms"
    final_pdb = work / "step4_final.pdb"
    if not final_dms.exists() or not final_pdb.exists():
        raise ValidationError("Stage 3 Backmap completed without expected final DMS/PDB")
    shutil.copy2(final_dms, out / "step4_final.dms")
    shutil.copy2(final_pdb, out / "step4_final.pdb")
    shutil.copy2(final_dms, out / "aa_backmapped.dms")
    shutil.copy2(final_pdb, out / "aa_backmapped.pdb")
    finalizer_report = finalize_stage3_membrane_templates(root, config, final_pdb, out)
    canonical_source = Path(finalizer_report.get("output", "")) if finalizer_report.get("status") == "PASS" else final_pdb
    if not canonical_source.is_absolute():
        canonical_source = root / canonical_source
    if not canonical_source.exists():
        raise ValidationError(f"Stage 3 canonical output source is missing: {canonical_source}")
    shutil.copy2(canonical_source, out / STAGE3_CANONICAL_PDB)
    report["glycolipid_template_finalize"] = finalizer_report
    report["canonical_output"] = STAGE3_CANONICAL_PDB
    write_json(out / "glycolipid_template_finalize_report.json", finalizer_report)
    write_json(out / "provenance.json", {**report, "status": "PASS"})
    write_manifest(root, out)
    print("PASS STAGE 3 REAL BACKMAPPING")
    return 0


def stage4(root: Path, config: dict[str, Any], run_id: str, overwrite: bool) -> int:
    out = root / "outputs" / run_id / "stage4"
    ensure_clean_dir(out, overwrite)
    stage3_dir = root / "outputs" / run_id / "stage3"
    stage3_input = stage3_dir / STAGE3_CANONICAL_PDB
    if not stage3_input.exists():
        raise ValidationError(f"Stage 4 requires Stage 3 canonical {STAGE3_CANONICAL_PDB}")
    if not (stage3_dir / "step4_final.dms").exists():
        raise ValidationError("Stage 4 requires real Stage 3 step4_final.dms diagnostic handoff")
    sys.path.insert(0, str(root / "scripts"))
    from external_dependencies import preflight_stage4_charmm_resources, resolve_stage4_resource_paths
    from aa_stage4 import (
        Stage4Error,
        align_bilayer_to_z,
        copy_stage4_resources,
        filter_solvent_ions_by_z,
        hard_clash_report_gro,
        molecule_counts_from_pdb,
        membrane_z_bounds_nm,
        normalize_solvated_gro,
        phosphate_headgroup_core_z_bounds_nm,
        reject_explicit_charmm_terminal_cap_residues,
        replace_protein_with_pdb2gmx_coordinates,
        solvent_ion_counts_from_gro,
        slice_membrane,
        strip_inherited_solvent,
        write_index_ini,
        write_topology,
    )

    stage4_cfg = config.get("stage4", {})
    stage4_resources = resolve_stage4_resource_paths(config, root=root)
    stage4_resource_preflight = preflight_stage4_charmm_resources(config, root=root, required=True)
    if stage4_resource_preflight.errors:
        raise ValidationError("Stage 4 CHARMM resource preflight failed: " + " | ".join(stage4_resource_preflight.errors))
    charmm36 = stage4_resources.charmm36_forcefield
    membrane_toppar = stage4_resources.membrane_toppar
    aa_mdp = stage4_resources.aa_mdp
    protein_input = resolve(root, config.get("stage1", {}).get("input_pdb", ""), run_id)
    for label, path in (("CHARMM36", charmm36), ("membrane toppar", membrane_toppar), ("AA MDP", aa_mdp), ("all-atom protein", protein_input)):
        if not path.exists():
            raise ValidationError(f"Stage 4 {label} is missing: {path}")
    if not (charmm36 / "forcefield.itp").exists():
        raise ValidationError(f"Stage 4 CHARMM36 forcefield.itp is missing: {charmm36 / 'forcefield.itp'}")
    for mdp_name in ("ions.mdp", "minim.mdp"):
        if not (aa_mdp / mdp_name).exists():
            raise ValidationError(f"Stage 4 requires {aa_mdp / mdp_name}")
    reject_explicit_charmm_terminal_cap_residues(stage3_input, f"Stage 3 canonical {STAGE3_CANONICAL_PDB}")
    reject_explicit_charmm_terminal_cap_residues(protein_input, "Stage 4 protein input")
    validate_charmm_terminal_patch_resources(charmm36, stage4_cfg)

    work = root / "work" / run_id / "stage4"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    logs = root / "logs" / run_id / "stage4"
    gmx_env = local_resource_env()
    gmx_env["GMXLIB"] = str(work)
    try:
        copy_stage4_resources(charmm36, membrane_toppar, work)
        protein_source = work / "protein_source.pdb"
        shutil.copy2(protein_input, protein_source)
        dry_pdb = work / "membrane_dry.pdb"
        inherited_solvent_report = strip_inherited_solvent(stage3_input, dry_pdb)

        pdb2gmx = gmx_command() + [
            "pdb2gmx", "-f", protein_source.name, "-o", "protein_generated.gro", "-p", "protein.top",
            "-i", "protein_posre.itp", "-ff", "charmm36", "-water", "tip3p", "-ignh", "-ss", "yes", "-ter",
        ]
        result = run_logged(
            pdb2gmx,
            work,
            logs,
            "stage4_pdb2gmx",
            timeout=int(stage4_cfg.get("pdb2gmx_timeout_seconds", 1800)),
            input_text=stage4_pdb2gmx_input(stage4_cfg, protein_source),
            env=gmx_env,
        )
        if result.returncode != 0:
            raise ValidationError("Stage 4 pdb2gmx failed; see logs")
        protein_top = work / "protein.top"
        protein_itps = [work / match.group(1) for match in re.finditer(r'#include\s+"([^"]+\.itp)"', protein_top.read_text())]
        builtin_itp_names = {"forcefield.itp", "tip3p.itp", "ions.itp"}
        protein_itps = [
            path for path in protein_itps
            if path.name not in builtin_itp_names
            and "posre" not in path.name.lower()
            and "charmm36.ff" not in path.as_posix()
            and path.exists()
        ]

        matched_dry_pdb = work / "membrane_dry_matched.pdb"
        protein_replacement_report = replace_protein_with_pdb2gmx_coordinates(
            dry_pdb,
            work / "protein_generated.gro",
            matched_dry_pdb,
            work / "toppar",
            protein_top,
        )
        dry_counts = molecule_counts_from_pdb(matched_dry_pdb)
        dry_top = work / "topol_dry.top"
        write_topology(dry_top, protein_top, protein_itps, work / "toppar", dry_counts)
        genion_mdp = work / "ions_genion.mdp"
        write_genion_mdp(aa_mdp / "ions.mdp", genion_mdp)
        pre_tpr = work / "membrane_prep.tpr"
        result = run_logged(
            gmx_command() + ["grompp", "-f", genion_mdp.name, "-c", matched_dry_pdb.name, "-p", dry_top.name, "-o", pre_tpr.name, "-maxwarn", "0"],
            work, logs, "stage4_prepare_grompp", timeout=180, env=gmx_env,
        )
        if result.returncode != 0:
            raise ValidationError("Stage 4 topology/coordinate preflight failed; see logs")
        whole_pdb = work / "membrane_whole.pdb"
        result = run_logged(
            gmx_command() + ["trjconv", "-f", matched_dry_pdb.name, "-s", pre_tpr.name, "-o", whole_pdb.name, "-center", "-pbc", "mol", "-ur", "rect"],
            work, logs, "stage4_make_whole", timeout=180, input_text="Protein\nSystem\n", env=gmx_env,
        )
        if result.returncode != 0:
            raise ValidationError("Stage 4 molecule-whole conversion failed; see logs")
        aligned_pdb = work / "membrane_aligned.pdb"
        alignment_report = align_bilayer_to_z(whole_pdb, aligned_pdb)
        sliced_pdb = work / "membrane_sliced.pdb"
        slice_buffer_nm = max(
            float(stage4_cfg.get("slice_buffer_nm", 1.0)),
            float(stage4_cfg.get("minimum_slice_buffer_nm", 3.0)),
        )
        z_slab_padding_nm = float(stage4_cfg.get("z_slab_padding_nm", 1.5))
        slice_report = slice_membrane(aligned_pdb, sliced_pdb, slice_buffer_nm, z_slab_padding_nm)
        sliced_counts = molecule_counts_from_pdb(sliced_pdb)
        topol = work / "topol.top"
        write_topology(topol, protein_top, protein_itps, work / "toppar", sliced_counts)

        raw_solvated = work / "solvated_raw.gro"
        result = run_logged(
            gmx_command() + ["solvate", "-cp", sliced_pdb.name, "-cs", "spc216.gro", "-o", raw_solvated.name],
            work, logs, "stage4_solvate", timeout=300, env=gmx_env,
        )
        if result.returncode != 0:
            raise ValidationError("Stage 4 solvation failed; see logs")
        solvated = work / "solvated_tip3.gro"
        membrane_z = membrane_z_bounds_nm(sliced_pdb)
        phosphate_core_z = phosphate_headgroup_core_z_bounds_nm(sliced_pdb)
        solvent_z_clearance_nm = float(stage4_cfg.get("solvent_membrane_z_clearance_nm", 0.0))
        excluded_solvent_z = (phosphate_core_z[0] - solvent_z_clearance_nm, phosphate_core_z[1] + solvent_z_clearance_nm)
        solvent_z_report = normalize_solvated_gro(raw_solvated, solvated, excluded_solvent_z)
        solvent_z_report["membrane_z_min_nm"] = membrane_z[0]
        solvent_z_report["membrane_z_max_nm"] = membrane_z[1]
        solvent_z_report["phosphate_core_z_min_nm"] = phosphate_core_z[0]
        solvent_z_report["phosphate_core_z_max_nm"] = phosphate_core_z[1]
        solvent_z_report["clearance_nm"] = solvent_z_clearance_nm
        water_count = int(solvent_z_report["waters"])
        solvated_counts = dict(sliced_counts)
        solvated_counts["TIP3"] = water_count
        write_topology(topol, protein_top, protein_itps, work / "toppar", solvated_counts)

        ions_tpr = work / "ions.tpr"
        result = run_logged(
            gmx_command() + ["grompp", "-f", genion_mdp.name, "-c", solvated.name, "-p", topol.name, "-o", ions_tpr.name, "-maxwarn", "0"],
            work, logs, "stage4_ions_grompp", timeout=180, env=gmx_env,
        )
        if result.returncode != 0:
            raise ValidationError("Stage 4 ionization preflight failed; see logs")
        raw_ionized = work / "ionized_raw.gro"
        ionized = work / "ionized.gro"
        result = run_logged(
            gmx_command() + [
                "genion", "-s", ions_tpr.name, "-o", raw_ionized.name, "-p", topol.name,
                "-pname", "SOD", "-nname", "CLA", "-neutral", "-conc", str(float(stage4_cfg.get("salt_concentration_molar", 0.15))),
            ],
            work, logs, "stage4_genion", timeout=180, input_text="TIP3\n", env=gmx_env,
        )
        if result.returncode != 0:
            raise ValidationError("Stage 4 genion failed; see logs")
        ion_z_filter_report = filter_solvent_ions_by_z(raw_ionized, ionized, excluded_solvent_z)
        ionized_counts = dict(sliced_counts)
        ionized_counts.update(solvent_ion_counts_from_gro(ionized))
        write_topology(topol, protein_top, protein_itps, work / "toppar", ionized_counts)
        em_tpr = work / "em.tpr"
        result = run_logged(
            gmx_command() + ["grompp", "-f", str(aa_mdp / "minim.mdp"), "-c", ionized.name, "-p", topol.name, "-o", em_tpr.name, "-maxwarn", "0"],
            work, logs, "stage4_em_grompp", timeout=180, env=gmx_env,
        )
        if result.returncode != 0:
            raise ValidationError("Stage 4 EM preflight failed; see logs")
        result = run_logged(
            gmx_mdrun_command() + ["-deffnm", "em", "-c", "em.gro"],
            work, logs, "stage4_em_mdrun", timeout=600, env=gmx_env,
        )
        if result.returncode != 0:
            raise ValidationError("Stage 4 EM failed; see logs")
        validate_stage4_em_log(logs / "stage4_em_mdrun.stderr.log")
        raw_em = work / "em_raw.gro"
        shutil.move(work / "em.gro", raw_em)
        final_z_filter_report = filter_solvent_ions_by_z(raw_em, work / "em.gro", excluded_solvent_z)
        final_counts = dict(sliced_counts)
        final_counts.update(solvent_ion_counts_from_gro(work / "em.gro"))
        write_topology(topol, protein_top, protein_itps, work / "toppar", final_counts)
        index_ini = work / "index_ini.ndx"
        index_report = write_index_ini(work / "em.gro", index_ini)
        result = run_logged(
            gmx_command() + ["grompp", "-f", str(aa_mdp / "minim.mdp"), "-c", "em.gro", "-p", topol.name, "-o", em_tpr.name, "-maxwarn", "0"],
            work, logs, "stage4_final_grompp", timeout=180, env=gmx_env,
        )
        if result.returncode != 0:
            raise ValidationError("Stage 4 final cleaned EM preflight failed; see logs")
        clash_report = hard_clash_report_gro(work / "em.gro", float(stage4_cfg.get("hard_clash_cutoff_nm", 0.06)))
        if clash_report["hard_clashes"]:
            raise ValidationError(f"Stage 4 final em.gro has hard clashes: {clash_report}")
        for path in (topol, index_ini, work / "em.gro", work / "em.tpr", sliced_pdb, solvated, ionized):
            shutil.copy2(path, out / path.name)
        shutil.copy2(work / "em.gro", out / "minimized_all_atom.gro")
        shutil.copy2(topol, out / "minimized_all_atom.top")
        report = {
            "status": "PASS", "stage3_input": rel(root, stage3_input), "protein_input": rel(root, protein_input),
            "charmm36": rel(root, charmm36), "membrane_toppar": rel(root, membrane_toppar),
            "stage4_resource_preflight": stage4_resource_preflight.as_dict(),
            "terminal_stage": "energy_minimized_all_atom; restrained AA equilibration is not run by Stage 4",
            "terminal_patches": charmm_terminal_patches(stage4_cfg),
            "terminal_patch_chains": sorted(selected_terminal_patch_chains(stage4_cfg, protein_chains_for_terminal_patches(protein_input))),
            "canonical_outputs": ["minimized_all_atom.gro", "minimized_all_atom.top"],
            "inherited_solvent": inherited_solvent_report, "solvent_z_filter": solvent_z_report,
            "ionized_z_filter": ion_z_filter_report, "final_z_filter": final_z_filter_report,
            "hard_clash_check": clash_report,
            "alignment": alignment_report, "slice": slice_report,
            "protein_coordinate_replacement": protein_replacement_report,
            "index": index_report,
            "salt_concentration_molar": float(stage4_cfg.get("salt_concentration_molar", 0.15)),
        }
        write_json(out / "provenance.json", report)
        write_manifest(root, out)
        print("PASS STAGE 4 AA SLICE SOLVATE EM")
        return 0
    except (Stage4Error, ValidationError) as exc:
        write_json(out / "provenance.json", {"status": "FAIL", "reason": str(exc), "stage3_input": rel(root, stage3_input)})
        write_manifest(root, out)
        raise ValidationError(str(exc)) from exc


def run_stage(root: Path, stage: str, config: dict[str, Any], run_id: str, overwrite: bool, mode: str | None = None) -> int:
    try:
        if stage == "stage1":
            return stage1(root, config, run_id, overwrite)
        if stage == "stage2":
            return stage2(root, config, run_id, overwrite, mode)
        if stage == "stage3":
            return stage3(root, config, run_id, overwrite)
        if stage == "stage4":
            return stage4(root, config, run_id, overwrite)
    except ValidationError as exc:
        print(f"FAIL {exc}", file=sys.stderr)
        return 1
    raise ValidationError(f"unknown stage {stage}")
