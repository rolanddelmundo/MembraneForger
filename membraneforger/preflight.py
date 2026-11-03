#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// --check: the fast pre-flight run before a build.
#//=============================================================
"""Check the inputs, programs and output directory of a build in seconds, without building anything.

`python -m membraneforger --check ...` takes the same options as a build and runs what the first minutes of one would
run anyway, in the order a new user hits them: GROMACS, the Python requirements, the mstool interpreter, the all-atom
input (with the conventions pdb2gmx will apply to it), the coarse-grained frame, PPM or OPM for the orientation, the
backmapping data and the output directory. Nothing is written. Every line is OK, WARN, FAIL or INFO, and the command
exits 1 when anything failed, so a wrong input or a missing program is found before the job is queued."""
import importlib
import subprocess
import sys
from pathlib import Path

from .config import AMINO, RENAME_RESIDUE
from .mstool_worker import TESTED_MSTOOL
from .orientation import OrientationFailure, OrientationRequest, find_ppm_executable
from .runtools import gromacs_version
from .structio import read_pdb, residues_in_order
from .topology import chain_termini, find_disulfides
from .validation import protect_inputs, read_all_atom, read_cg, refuse_overwrite

__all__ = ['Check', 'check_build', 'print_checks', 'run_check']

OK, WARN, FAIL, INFO = "OK", "WARN", "FAIL", "INFO"
PYTHON_REQUIREMENTS = ("numpy", "scipy", "networkx", "pandas", "openmm")


class Check:
    """One line of the report: a status, the item checked and what was found."""

    def __init__(self, status: str, item: str, text: str):
        self.status, self.item, self.text = status, item, text

    def __str__(self) -> str:
        return f"{self.status}: {self.item}: {self.text}"


def _mstool_version(python: str) -> tuple[str | None, str]:
    """Import mstool in the backmapping interpreter and return (version, error text)."""
    try:
        result = subprocess.run([python, "-c", "import mstool; print(getattr(mstool, '__version__', 'unknown'))"],
                                text=True, capture_output=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return None, str(exc)
    if result.returncode:
        last = (result.stderr.strip().splitlines() or [result.stdout.strip() or "no output"])[-1]
        return None, last
    return result.stdout.strip(), ""


def check_programs(gmx: str | None, gmx_problem: str | None, python: str) -> list[Check]:
    """GROMACS, the Python requirements of this interpreter, and mstool in the backmapping interpreter."""
    checks = []
    if gmx_problem:
        checks.append(Check(FAIL, "GROMACS", f"{gmx_problem}; the build runs grompp, solvate, genion and mdrun through it"))
    else:
        try:
            checks.append(Check(OK, "GROMACS", f"{gmx} is GROMACS {gromacs_version(gmx)}"))
        except SystemExit as exc:
            checks.append(Check(FAIL, "GROMACS", str(exc)))
    missing = [m for m in PYTHON_REQUIREMENTS if importlib.util.find_spec(m) is None]
    if missing:
        checks.append(Check(FAIL, "Python requirements", f"{sys.executable} cannot import {', '.join(missing)} "
                                                         "(pip install . from the repository installs them)"))
    else:
        checks.append(Check(OK, "Python requirements", f"numpy, scipy, networkx, pandas and openmm import in {sys.executable}"))
    version, error = _mstool_version(python)
    if version is None:
        checks.append(Check(FAIL, "mstool", f"{python} cannot import mstool ({error}); install mstool 0.3.9 or 0.3.10 there, or "
                                            "point --mstool-python at an interpreter that has it"))
    elif version not in TESTED_MSTOOL:
        checks.append(Check(WARN, "mstool", f"version {version} in {python}; tested with {' and '.join(TESTED_MSTOOL)}"))
    else:
        checks.append(Check(OK, "mstool", f"version {version} in {python}"))
    return checks


def check_all_atom(path: Path, forcefield: Path, orientation: OrientationRequest, exclude_disulfides: tuple = ()) -> list[Check]:
    """Read the all-atom input as the build will, and say what pdb2gmx and the orientation will do with it."""
    try:
        atoms, notes = read_all_atom(path, forcefield)
    except SystemExit as exc:
        return [Check(FAIL, "all-atom input", f"{exc} (see docs/preparing_inputs.md)")]
    checks = [Check(OK, "all-atom input", notes[0])]
    chains: dict[str, list[dict]] = {}
    for a in atoms:
        if a["resname"] in AMINO:
            chains.setdefault(a["chain"] or a["segid"] or "A", []).append(a)
    sizes = ", ".join(f"{label} ({len(residues_in_order(chain))} residues)" for label, chain in chains.items())
    checks.append(Check(INFO, "protein chains", sizes))
    for label, chain in chains.items():
        residues = residues_in_order(chain)
        try:
            nterm, cterm = chain_termini(residues)
            checks.append(Check(INFO, f"chain {label} termini", f"N {nterm}, C {cterm}"))
        except SystemExit as exc:
            checks.append(Check(FAIL, f"chain {label} termini", f"{exc}; a chain must end in OXT/OT1/OT2 or a CT2/CT3 cap "
                                                                "(see docs/preparing_inputs.md)"))
    renamed = {}
    for a in atoms:
        if a["resname"] in RENAME_RESIDUE and a["atom"] == "CA":
            renamed[a["resname"]] = renamed.get(a["resname"], 0) + 1
    if renamed.get("HIS"):
        checks.append(Check(WARN, "histidines", f"{renamed['HIS']} HIS residues will be built as HSD (neutral, H on ND1); name "
                                                "them HID/HSD, HIE/HSE or HIP/HSP in the input to choose the protonation"))
    try:
        found = find_disulfides(chains, exclude_disulfides)
        pairs = [(label, a, b) for label, bonds in found.items() for a, b in bonds]
        checks.append(Check(INFO, "disulfides", (", ".join(f"{label}:{a}-{b}" for label, a, b in pairs) or "none") +
                            (f"; excluded {', '.join(f'{c}:{a}-{c}:{b}' for c, a, b in exclude_disulfides)}" if exclude_disulfides else "")))
    except SystemExit as exc:
        checks.append(Check(FAIL, "disulfides", str(exc)))
    if orientation.mode == "none":
        checks.append(Check(WARN, "orientation", "none: the coordinates are used as given, so the membrane normal must already be along z "
                                                 "with the extracellular side at +z"))
    else:
        if len(chains) > 1 and not orientation.chains:
            checks.append(Check(FAIL, "orientation", f"{len(chains)} protein chains and no --orient-chain: name the chain that spans "
                                                     f"the membrane (one of {', '.join(chains)})"))
        elif orientation.chains and (set(orientation.chains) - set(chains)):
            checks.append(Check(FAIL, "orientation", f"--orient-chain {','.join(orientation.chains)} is not a protein chain of the "
                                                     f"input ({', '.join(chains)})"))
        if orientation.nterm_side == "auto" and orientation.mode in ("auto", "ppm"):
            checks.append(Check(WARN, "orientation", "no --nterm-side: PPM needs the side (in or out) of residue 1 of the anchor chain as it "
                                                     "is in the file, unless an exact OPM entry gives it"))
    return checks


def check_coarse_grain(path: Path | None, membrane: Path | None) -> list[Check]:
    """Read the coarse-grained frame (or the --membrane PDB) as the build will."""
    if membrane is not None:
        try:
            atoms, cryst1 = read_pdb(membrane)
        except (ValueError, IndexError) as exc:
            return [Check(FAIL, "--membrane", f"{membrane.name}: {exc}")]
        if not cryst1:
            return [Check(FAIL, "--membrane", f"{membrane.name} has no CRYST1 box")]
        return [Check(OK, "--membrane", f"{membrane.name}: {len(atoms)} atoms, CRYST1 present")]
    try:
        atoms, box = read_cg(path)
    except SystemExit as exc:
        return [Check(FAIL, "coarse-grained frame", str(exc))]
    return [Check(OK, "coarse-grained frame", f"{path.name}: {len(atoms)} beads, box {box[0]:.2f} x {box[1]:.2f} x {box[2]:.2f} nm")]


def check_orientation_programs(orientation: OrientationRequest) -> list[Check]:
    """PPM 3.0 and the OPM reference, as far as they can be checked without running them."""
    if orientation.mode == "none":
        return []
    checks = []
    if orientation.opm_file:
        checks.append(Check(OK if orientation.opm_file.is_file() else FAIL, "OPM file", str(orientation.opm_file)))
    elif orientation.mode in ("auto", "opm"):
        checks.append(Check(INFO, "OPM", "the exact OPM entry is downloaded when a PDB ID is known (--pdb-id or the HEADER record), "
                                         "cached in ~/.cache/membraneforger/opm; offline, pass --opm-file"))
    if orientation.mode in ("auto", "ppm"):
        try:
            exe = find_ppm_executable(orientation.ppm_exe)
            reslib = exe.parent / "res.lib"
            checks.append(Check(OK if reslib.is_file() else FAIL, "PPM 3.0",
                                f"{exe}" + ("" if reslib.is_file() else f"; res.lib is missing next to it ({reslib})")))
        except OrientationFailure as exc:
            status = FAIL if orientation.mode == "ppm" else WARN
            checks.append(Check(status, "PPM 3.0", f"{exc}" + ("" if status == FAIL else
                                                                "; with --orientation auto it is needed unless an exact OPM entry matches")))
    return checks


def check_output(out: Path, inputs: list[Path], membrane: Path | None, forcefield: Path, data: Path,
                 overwrite: bool, overwrite_restart: bool) -> list[Check]:
    """The backmapping data and the output directory, exactly as the build refuses or accepts them."""
    has_map = (data / "map.dat").is_file()
    checks = [Check(OK if has_map else FAIL, "backmapping data", f"{data} ({'has' if has_map else 'no'} map.dat)")]
    try:
        protect_inputs(inputs, out, keep=membrane, resources=(forcefield, data))
        refuse_overwrite(out, keep=membrane, overwrite=overwrite, overwrite_restart=overwrite_restart)
        state = "exists, no earlier build" if out.exists() else "will be created"
        checks.append(Check(OK, "output directory", f"{out} ({state}{'; --overwrite given' if overwrite else ''})"))
    except SystemExit as exc:
        checks.append(Check(FAIL, "output directory", str(exc)))
    return checks


def check_build(all_atom: Path | None, coarse_grain: Path | None, membrane: Path | None, out: Path, forcefield: Path, data: Path,
                gmx: str | None, gmx_problem: str | None, python: str, orientation: OrientationRequest, ntomp: int,
                overwrite: bool = False, overwrite_restart: bool = False, exclude_disulfides: tuple = ()) -> list[Check]:
    """Every pre-flight check of one build, in the order the build would meet them."""
    checks = check_programs(gmx, gmx_problem, python)
    if all_atom is not None:
        checks += check_all_atom(all_atom, forcefield, orientation, exclude_disulfides)
    checks += check_coarse_grain(coarse_grain, membrane)
    if membrane is None:
        checks += check_orientation_programs(orientation)
    inputs = [p for p in (all_atom, coarse_grain, membrane) if p]
    checks += check_output(out, inputs, membrane, forcefield, data, overwrite, overwrite_restart)
    checks.append(Check(INFO, "resources", f"{ntomp} threads; a GPCR-sized complex builds in about 15-30 minutes on 8 cores with "
                                           "about 1.5 GB of memory (--ntomp sets the threads)"))
    return checks


def print_checks(checks: list[Check], name: str) -> int:
    """Print the report and the verdict; returns the exit code (1 when anything failed)."""
    for check in checks:
        print(f"[{name}] {check}")
    failed = [c for c in checks if c.status == FAIL]
    warned = [c for c in checks if c.status == WARN]
    verdict = "FAIL" if failed else "PASS"
    print(f"{name}\tCHECK {verdict}\t{len(failed)} failed, {len(warned)} warning(s); nothing was built")
    return int(bool(failed))


def run_check(name: str, **kwargs) -> int:
    """check_build followed by print_checks."""
    return print_checks(check_build(**kwargs), name)
