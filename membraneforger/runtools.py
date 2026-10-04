#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Logging, hashing and the one subprocess wrapper.
#//=============================================================
"""Levelled logging, file hashing and the single wrapper every external command goes through."""
import hashlib
import os
import re
import shutil
import subprocess
from pathlib import Path

__all__ = ['LOG_NAME', 'LEVELS', 'log', 'sha256', 'run_command', 'find_gromacs', 'gromacs_version']

LOG_NAME = "membranebuilder.log"
LEVELS = ("INFO", "CHECK", "PASS", "WARN", "ERROR", "DEBUG")


def log(out: Path, message: str, level: str = "INFO") -> None:
    """Print one levelled message and append it to the build log."""
    if level not in LEVELS:
        raise ValueError(f"unknown log level {level}")
    line = f"{level}: {message}"
    if level != "DEBUG":
        print(f"[{out.name}] {line}", flush=True)
    with (out / LOG_NAME).open("a") as fh:
        fh.write(line + "\n")


def sha256(path: Path) -> str:
    """Hash a file so we can tell if it changed mid-build."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_command(out: Path, cmd: list, stdin: str | None = None, env: dict | None = None, produces: tuple = (),
                cwd: Path | None = None) -> str:
    """Run an external command without a shell, log its transcript, and fail unless it exits 0 and writes its outputs."""
    # Files listed in `produces` are deleted first and must exist and be non-empty afterwards, so a stale or
    # partial file from an earlier run can never stand in for the output of a command that did not finish.
    cmd = [str(c) for c in cmd]
    targets = [(cwd or out) / p for p in produces]
    for target in targets:
        target.unlink(missing_ok=True)
    env = dict(env if env is not None else os.environ, GMX_MAXBACKUP="-1")  # never leave #file.N# backups behind
    try:
        result = subprocess.run(cmd, cwd=cwd or out, text=True, capture_output=True, input=stdin, env=env)
    except OSError as exc:
        raise SystemExit(f"cannot run {cmd[0]}: {exc}") from None
    output = result.stdout + result.stderr
    log(out, f"$ {' '.join(cmd)}\n{output}", "DEBUG")
    step = next((c for c in cmd[1:] if not c.startswith("-") and "/" not in c), Path(cmd[0]).name)
    if result.returncode:
        errors = re.findall(r"(?ms)^(?:ERROR|WARNING|Fatal error).*?(?=^\s*$)", output)
        if "Traceback (most recent call last)" in output:  # a python step: its last line is the error
            errors = [output.strip().splitlines()[-1]]
        raise SystemExit(f"{step} failed (exit {result.returncode}): "
                         f"{' | '.join(' '.join(e.split()) for e in errors[:4]) or output[-800:]}")
    missing = [t.name for t in targets if not t.is_file() or t.stat().st_size == 0]
    if missing:
        raise SystemExit(f"{step} exited 0 but did not write {missing}")
    return output


def find_gromacs(requested: str | None) -> str:
    """Return the GROMACS command to use, refusing to guess when the requested one is absent."""
    command = requested or next((g for g in ("gmx", "gmx_mpi", "gmx_d") if shutil.which(g)), None)
    if command is None or not shutil.which(command.split()[0]):
        raise SystemExit(f"GROMACS not found ({command or 'gmx, gmx_mpi, gmx_d'}); pass --gmx")
    return command


def gromacs_version(gmx: str) -> str:
    """Report the version string of the GROMACS command in use."""
    result = subprocess.run(gmx.split() + ["--version"], text=True, capture_output=True)
    found = re.search(r"GROMACS version:\s*(\S+)", result.stdout + result.stderr)
    if result.returncode or not found:
        raise SystemExit(f"{gmx} --version did not report a GROMACS version")
    return found.group(1)
