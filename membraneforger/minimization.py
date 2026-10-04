#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Energy minimization and its validation.
#//=============================================================
"""Run steepest-descent minimization into em.unverified.gro and validate the result before it may be published."""
import math
import re
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .config import EM_MDP, FMAX_TARGET
from .runtools import run_command
from .structio import read_gro, wrap, xyz_nm
from .validation import check_inputs_unchanged

__all__ = ['run_em', 'closest_contact', 'validate_em']

def run_em(out: Path, topology: dict, gmx: str, ntomp: int, nsteps: int) -> None:
    """Run steepest-descent minimization (PME when neutral)."""
    coulomb = "PME" if abs(topology["charge"]) < 1e-3 else "Cut-off"
    (out / "em.mdp").write_text(EM_MDP.format(emtol=FMAX_TARGET, nsteps=nsteps, coulombtype=coulomb))
    for stale in ("em.gro", "em.unverified.gro", "em.log", "em.edr", "em.trr"):
        (out / stale).unlink(missing_ok=True)
    run_command(out, gmx.split() + ["grompp", "-f", "em.mdp", "-c", "solv_ions.gro", "-p", "topol.top",
                                "-n", "index_ini.ndx", "-o", "em.tpr", "-po", "mdout.mdp", "-maxwarn", "0"],
            produces=("em.tpr",))
    (out / "mdout.mdp").unlink(missing_ok=True)
    run_command(out, gmx.split() + ["mdrun", "-deffnm", "em", "-c", "em.unverified.gro", "-ntmpi", "1", "-ntomp", str(ntomp)],
            produces=("em.unverified.gro", "em.log"))


def closest_contact(out: Path, number: int) -> str:
    """Name the atom closest to the one EM blew up on, to point at the overlap."""
    atoms, box = read_gro(out / "solv_ions.gro")
    label = lambda a: f"{a['resname']}{a['resid']}:{a['atom']}"
    xyz, target = wrap(xyz_nm(atoms), box), atoms[number - 1]
    others = [i for i, a in enumerate(atoms) if (a["resid"], a["resname"]) != (target["resid"], target["resname"])]
    d, k = cKDTree(xyz[others], boxsize=box[:3]).query(xyz[number - 1])
    return f"atom {number} {label(target)} starts {10 * d:.2f} A from {label(atoms[others[k]])} before EM"


def validate_em(system: dict, topology: dict, index: dict) -> dict:
    """Check EM converged and that the minimized structure, topology, index and box are all consistent."""
    out = system["out"]
    check_inputs_unchanged(system["inputs"])
    text = (out / "em.log").read_text(errors="replace")
    fmax, energy = re.findall(r"Maximum force\s*=\s*(\S+)", text), re.findall(r"Potential Energy\s*=\s*(\S+)", text)
    if not fmax or not energy or not math.isfinite(float(energy[-1])) or not math.isfinite(float(fmax[-1])):
        culprit = re.findall(r"on atom (\d+)", text)
        raise SystemExit(f"non-finite EM force: {closest_contact(out, int(culprit[-1])) if culprit else 'no force report in em.log'}")
    if not float(fmax[-1]) < FMAX_TARGET:
        raise SystemExit(f"EM final Fmax {float(fmax[-1]):.1f} >= {FMAX_TARGET:.0f} kJ/mol/nm")
    gro = (out / "em.unverified.gro").read_text().splitlines()
    atoms = gro[2:2 + int(gro[1])]
    xyz = np.array([[float(line[20:28]), float(line[28:36]), float(line[36:44])] for line in atoms])
    if [line[10:15].strip() for line in atoms] != topology["names"] or not np.isfinite(xyz).all():
        raise SystemExit("minimized atoms do not match topol.top")
    if index["System"] != len(atoms) or index["Protein_LIG"] + index["MEMB"] + index["SOL_ION"] != len(atoms):
        raise SystemExit(f"index groups do not partition the system: {index}")
    if not index["SOL_ION"] or abs(topology["charge"]) > 1e-3:
        raise SystemExit("EM system is not solvated and neutral")
    if not np.allclose([float(v) for v in gro[-1].split()][:2], read_gro(out / "boxed.gro")[1][:2], atol=1e-3):
        raise SystemExit("EM changed the membrane XY box")
    steps = re.findall(r"converged to Fmax < \S+ in (\d+) steps", text)
    return {"summary": f"Fmax {float(fmax[-1]):.1f}, Epot {float(energy[-1]):.4e}, {steps[-1] if steps else '?'} steps",
            "fmax_kj_mol_nm": float(fmax[-1]), "potential_kj_mol": float(energy[-1]),
            "steps": int(steps[-1]) if steps else None, "fmax_target": FMAX_TARGET, "atoms": len(atoms)}
