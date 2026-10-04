#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Input validation and output-directory protection.
#//=============================================================
"""Read and validate both inputs, and keep a build from overwriting inputs or reusing stale outputs."""
import math
import shutil
from collections import Counter
from pathlib import Path

import numpy as np

from .config import AMINO, GENERATED, SOLVENT
from .runtools import LOG_NAME, sha256
from .structio import read_pdb, xyz_nm

__all__ = ['read_cg', 'read_all_atom', 'protect_inputs', 'clear_stale_outputs', 'check_inputs_unchanged']

def read_cg(path: Path) -> tuple[list[dict], list[float]]:
    """Read a coarse-grained .gro or .pdb into beads (nm) and an orthorhombic box (nm)."""
    try:
        if path.suffix.lower() == ".gro":
            lines = path.read_text(errors="replace").splitlines()
            n = int(lines[1])
            if n <= 0 or len(lines) < n + 3:
                raise ValueError(f"header declares {n} atoms but the file has {max(len(lines) - 3, 0)} atom lines")
            atoms = [{"resid": int(l[0:5]), "resname": l[5:10].strip(), "atom": l[10:15].strip(), "chain": "",
                      "x": float(l[20:28]), "y": float(l[28:36]), "z": float(l[36:44])} for l in lines[2:2 + n]]
            cell = [float(v) for v in lines[2 + n].split()]
            if len(cell) not in (3, 9) or any(abs(v) > 1e-6 for v in cell[3:]):
                raise ValueError("box is missing or not orthorhombic")
            box = cell[:3]
        elif path.suffix.lower() == ".pdb":
            if sum(l.startswith("MODEL ") for l in path.read_text(errors="replace").splitlines()) > 1:
                raise ValueError("multiple MODEL records")
            raw, cryst1 = read_pdb(path)
            if not cryst1:
                raise ValueError("no CRYST1 box")
            cryst = cryst1.ljust(54)
            if not all(abs(float(cryst[i:i + 7]) - 90.0) < 0.02 for i in (33, 40, 47)):
                raise ValueError("box is not orthorhombic")
            box = [float(cryst[6:15]) / 10.0, float(cryst[15:24]) / 10.0, float(cryst[24:33]) / 10.0]
            atoms = [dict(a, x=a["x"] / 10.0, y=a["y"] / 10.0, z=a["z"] / 10.0) for a in raw if a["altloc"] in " A"]
        else:
            raise ValueError("expected a .pdb or .gro file")
        if not atoms:
            raise ValueError("no atoms")
        if not np.isfinite(xyz_nm(atoms)).all() or not all(math.isfinite(v) and v > 0 for v in box):
            raise ValueError("non-finite coordinates or a non-positive box")
    except (ValueError, IndexError) as exc:
        raise SystemExit(f"coarse-grained input {path.name}: {exc}")
    return atoms, box


def read_all_atom(path: Path, forcefield: Path) -> tuple[list[dict], list[str]]:
    """Read and validate the all-atom protein/ligand PDB."""
    # Returns its atoms (altloc A only) and notes for the log.
    text = path.read_text(errors="replace").splitlines()
    try:
        atoms, _ = read_pdb(path)
    except (ValueError, IndexError) as exc:
        raise SystemExit(f"all-atom input {path.name}: malformed ATOM/HETATM record ({exc})")
    if not atoms:
        raise SystemExit(f"all-atom input {path.name}: no ATOM/HETATM records")
    if sum(l.startswith("MODEL ") for l in text) > 1:
        raise SystemExit(f"all-atom input {path.name}: multiple MODEL records")
    if any(l[:6] in ("ATOM  ", "HETATM") and l[26:27].strip() for l in text):
        raise SystemExit(f"all-atom input {path.name}: residue insertion codes are not supported")
    if not np.isfinite(xyz_nm(atoms)).all():
        raise SystemExit(f"all-atom input {path.name}: non-finite coordinates")
    alternates = sorted({f"{a['chain']}:{a['resname']}{a['resid']}" for a in atoms if a["altloc"] not in " A"})
    atoms = [a for a in atoms if a["altloc"] in " A"]
    if not any(a["resname"] in AMINO for a in atoms):
        raise SystemExit(f"all-atom input {path.name}: no protein residues")
    seen = Counter((a["chain"], a["segid"], a["resid"], a["resname"], a["atom"]) for a in atoms)
    repeats = [k for k, n in seen.items() if n > 1]
    if repeats:
        raise SystemExit(f"all-atom input {path.name}: {len(repeats)} duplicate atoms, e.g. "
                         f"{repeats[0][0]}:{repeats[0][3]}{repeats[0][2]}:{repeats[0][4]}")
    others = sorted({a["resname"] for a in atoms if a["resname"] not in AMINO})
    bulk = [n for n in others if n in SOLVENT | {"HOH", "SOL", "WAT", "TIP", "NA", "CL", "K"}]
    if bulk:
        raise SystemExit(f"all-atom input {path.name}: contains solvent or bulk ions {bulk}; "
                         "supply only the protein/ligand complex")
    known = {p.stem for p in (forcefield / "toppar").glob("*.itp")} | ({"V6G"} if (forcefield / "V6G_gromacs").is_dir() else set())
    missing = [n for n in others if n not in known]
    if missing:
        raise SystemExit(f"all-atom input {path.name}: no topology in {forcefield} for residue(s) {missing}")
    return atoms, [f"{path.name}: {len(atoms)} atoms, {sum(a['resname'] in AMINO for a in atoms)} protein, "
                   f"ligands {', '.join(others) or 'none'}; dropped altloc B of {', '.join(alternates) or 'none'}"]


def protect_inputs(inputs: list[Path], out: Path, keep: Path | None = None, resources: tuple = ()) -> None:
    """Refuse an output directory in which the build would overwrite or delete an input or an installation resource."""
    # `keep` is the --membrane source: it may live in the output directory only under the name membrane.pdb,
    # the one generated file that route never writes.
    for path in inputs:
        if out not in path.parents:
            continue
        first = path.relative_to(out).parts[0]
        if (first in GENERATED or first == LOG_NAME) and not (path == keep and path.name == "membrane.pdb" and path.parent == out):
            raise SystemExit(f"input {path} would be overwritten by the build; choose another --out")
    for resource in resources:
        if resource is not None and (out == resource or resource in out.parents or out in resource.parents):
            raise SystemExit(f"--out {out} overlaps the installation directory {resource}; choose another --out")


def clear_stale_outputs(out: Path, keep: Path | None = None) -> list[str]:
    """Delete every file a previous build could have left in the output directory; returns what was removed."""
    removed = []
    for name in GENERATED:
        target = out / name
        if target == keep or not (target.exists() or target.is_symlink()):
            continue
        if target.is_dir() and not target.is_symlink():
            shutil.rmtree(target)
        else:
            target.unlink()
        removed.append(name)
    return removed


def check_inputs_unchanged(hashes: dict) -> None:
    """Fail if any input file no longer has the hash recorded at the start of the build."""
    for path, digest in hashes.items():
        if sha256(path) != digest:
            raise SystemExit(f"{path.name} changed during the build")
