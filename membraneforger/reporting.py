#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// run_manifest.json: what ran, on what, with what result.
#//=============================================================
"""Write the machine-readable record of one build."""
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

__all__ = ['SOURCE_AT_IMPORT', 'source_state', 'software_versions', 'artifact_hashes', 'orientation_summary', 'write_run_manifest']

PACKAGE = Path(__file__).resolve().parent
ENVIRONMENT = ("MEMBRANEFORGER_TOPPAR", "MEMBRANEFORGER_DATA", "OMP_NUM_THREADS", "GMXLIB", "CONDA_PREFIX", "TMPDIR",
               "SLURM_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_CPUS_PER_TASK", "HOSTNAME")


def source_state() -> dict:
    """Identify the code that ran: the git commit when there is one, and always a hash of every package file."""
    files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(PACKAGE.glob("*.py"))}
    git = subprocess.run(["git", "-C", str(PACKAGE), "rev-parse", "HEAD"], text=True, capture_output=True)
    return {"git_commit": git.stdout.strip() if git.returncode == 0 else None, "package_dir": str(PACKAGE),
            "package_sha256": files}


SOURCE_AT_IMPORT = source_state()  # the code that is actually running, even if the files change during the build


def software_versions(session) -> dict:
    """Collect the versions of the interpreters and libraries the build used."""
    import networkx
    import numpy
    import scipy
    return {"python": sys.version.split()[0], "python_executable": sys.executable, "platform": platform.platform(),
            "numpy": numpy.__version__, "scipy": scipy.__version__, "networkx": networkx.__version__,
            "gromacs": session.record.get("gromacs_version"), "gromacs_command": session.gmx,
            "mstool": session.record.get("mstool", {}).get("mstool_version"),
            "mstool_path": session.record.get("mstool", {}).get("mstool_path"), "mstool_python": session.python}


def artifact_hashes(out: Path) -> dict:
    """Hash every regular file directly in the output directory (and toppar/) so outputs can be compared later."""
    files = [p for p in sorted(out.iterdir()) if p.is_file() and p.name != "run_manifest.json" and not p.name.startswith("#")]
    files += [p for p in sorted((out / "toppar").rglob("*")) if p.is_file()] if (out / "toppar").is_dir() else []
    return {str(p.relative_to(out)): {"sha256": hashlib.sha256(p.read_bytes()).hexdigest(), "bytes": p.stat().st_size}
            for p in files}


def orientation_summary(report: dict | None) -> dict | None:
    """The orientation record for the manifest: the full report plus flat keys for the fields a reader looks up first."""
    if not report:
        return report
    reference, ppm = report.get("reference") or {}, report.get("ppm") or {}
    fit = report.get("fit") or {}
    return {**report,
            "reference_source": reference.get("source") or (ppm.get("executable") and f"PPM run {ppm['executable']}") or None,
            "reference_sha256": reference.get("sha256") or ppm.get("output_sha256"),
            "ppm_version": ppm.get("version"),
            "ppm_parameters": {k: ppm.get(k) for k in ("membrane_code", "curvature", "nterm_side", "heteroatoms_submitted",
                                                     "input", "executable_sha256", "res_lib_sha256")} if ppm else None,
            "fit_rmsd": fit.get("core_rmsd_A", fit.get("rmsd_A")),
            "translation_vector": report.get("translation_A"),
            "membrane_center": report.get("membrane_center_A")}


def write_run_manifest(session, status: str, error: dict | None, started: str, hashes: dict, command: list) -> Path:
    """Write run_manifest.json for a finished (passed or failed) build and return its path."""
    record = session.record
    manifest = {
        "status": status, "error": error, "started": started, "finished": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "name": session.name, "command": command, "working_directory": os.getcwd(),
        "output_directory": str(session.out),
        "work_directory": str(session.work) if session.work else None,  # kept only when the build failed
        "source": SOURCE_AT_IMPORT, "software": software_versions(session),
        "environment": {k: os.environ[k] for k in ENVIRONMENT if k in os.environ},
        "inputs": {str(p): {"sha256": digest, "bytes": p.stat().st_size} for p, digest in hashes.items()},
        "resources": {"forcefield": str(session.forcefield), "data": str(session.data) if session.data else None,
                      "threads": session.ntomp, "max_em_steps": session.nsteps},
        "settings": vars(session.settings) if hasattr(session.settings, "__dict__") else
                    {f: getattr(session.settings, f) for f in session.settings.__dataclass_fields__},
        "stage_seconds": session.timings,
        "orientation": orientation_summary(record.get("orientation")), "slice": record.get("slice"), "slice_check": record.get("slice_check"),
        "membrane_validation": record.get("membrane_validation"), "box": record.get("box"),
        "coarse_grain": record.get("coarse_grain"), "aa_cg_mapping": record.get("aa_cg_mapping"),
        "embedding": record.get("embedding"), "lipid_edits": record.get("lipid_edits"), "box_trim": record.get("box_trim"),
        "backmap_isomer_review": record.get("backmap_isomer_review"),
        "backmap_attempts": record.get("backmap_attempts"), "mstool": record.get("mstool"),
        "ring_piercing_before_em": record.get("ring_piercing_before_em"),
        "closest_membrane_solute_contact_before_em": record.get("closest_membrane_solute_contact_before_em"),
        "disulfides": record.get("disulfides"), "removed_lipids": record.get("removed_lipids"), "box_nm": record.get("box_nm"),
        "topology": record.get("topology"),
        "em": record.get("em"), "audit": record.get("audit"),
        "artifacts": artifact_hashes(session.out), "screenshots": [],
    }
    path = session.out / "run_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    return path
