#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Energy minimization and its validation.
#//=============================================================
"""Run steepest-descent minimization into em.unverified.gro and validate the result before it may be published.

When the lipid topologies carry CHARMM-GUI dihedral restraints, a restrained minimization (emres) runs first and the
validated, unrestrained EM starts from its result."""
import math
import re
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .config import DIHRES_EM_FC, EM_MDP, FMAX_TARGET
from .runtools import run_command
from .structio import read_gro, wrap, xyz_nm
from .validation import check_inputs_unchanged

__all__ = ['has_dihedral_restraints', 'run_em', 'closest_contact', 'force_culprit', 'validate_em']

# On an Fmax failure, the residues within this distance of the atom with the largest force are named.
FORCE_NEIGHBOUR_NM = 0.3

def has_dihedral_restraints(out: Path) -> bool:
    """Whether any molecule topology in toppar/ carries [ dihedral_restraints ] rows."""
    for itp in (out / "toppar").glob("*.itp"):
        section = ""
        for raw in itp.read_text(errors="replace").splitlines():
            s = raw.split(";", 1)[0].strip()
            if s.startswith("["):
                section = s.strip("[] ").lower()
            elif section == "dihedral_restraints" and s and s[0].isdigit():
                return True
    return False


def run_em(out: Path, topology: dict, gmx: str, ntomp: int, nsteps: int) -> None:
    """Run steepest-descent minimization (PME when neutral), first with the lipid dihedral restraints when present."""
    coulomb = "PME" if abs(topology["charge"]) < 1e-3 else "Cut-off"
    start = "solv_ions.gro"
    if has_dihedral_restraints(out):
        _minimize(out, gmx, "emres", start, coulomb, ntomp, nsteps, f"-DDIHRES -DDIHRES_FC={DIHRES_EM_FC:g}")
        start = "emres.gro"
    _minimize(out, gmx, "em", start, coulomb, ntomp, nsteps)


def _minimize(out: Path, gmx: str, name: str, structure: str, coulomb: str, ntomp: int, nsteps: int, define: str = "") -> None:
    """grompp + mdrun one steepest-descent minimization from structure; em writes em.unverified.gro, others {name}.gro."""
    result = "em.unverified.gro" if name == "em" else f"{name}.gro"
    (out / f"{name}.mdp").write_text(EM_MDP.format(emtol=FMAX_TARGET, nsteps=nsteps, coulombtype=coulomb,
                                                   define=f"define          = {define}\n" if define else ""))
    for stale in (f"{name}.gro", result, f"{name}.log", f"{name}.edr", f"{name}.trr"):
        (out / stale).unlink(missing_ok=True)
    run_command(out, gmx.split() + ["grompp", "-f", f"{name}.mdp", "-c", structure, "-p", "topol.top",
                                    "-n", "index_ini.ndx", "-o", f"{name}.tpr", "-po", "mdout.mdp", "-maxwarn", "0"],
                produces=(f"{name}.tpr",))
    (out / "mdout.mdp").unlink(missing_ok=True)
    run_command(out, gmx.split() + ["mdrun", "-deffnm", name, "-c", result, "-ntmpi", "1", "-ntomp", str(ntomp)],
                produces=(result, f"{name}.log"))


def closest_contact(out: Path, number: int) -> str:
    """Name the atom closest to the one EM blew up on, to point at the overlap."""
    atoms, box = read_gro(out / "solv_ions.gro")
    label = lambda a: f"{a['resname']}{a['resid']}:{a['atom']}"
    xyz, target = wrap(xyz_nm(atoms), box), atoms[number - 1]
    others = [i for i, a in enumerate(atoms) if (a["resid"], a["resname"]) != (target["resid"], target["resname"])]
    d, k = cKDTree(xyz[others], boxsize=box[:3]).query(xyz[number - 1])
    return f"atom {number} {label(target)} starts {10 * d:.2f} A from {label(atoms[others[k]])} before EM"


def force_culprit(out: Path, number: int, structure: str = "em.unverified.gro", radius_nm: float = FORCE_NEIGHBOUR_NM) -> str:
    """Describe the atom EM ended with the largest force on: its residue, and every other residue within radius_nm of it.

    A finite Fmax that stays above the target comes from one strained spot (two molecules interlocked, a bond threaded
    through a ring, an atom pushed into another molecule); the residues around that atom name it. The residue's group
    (protein/ligand, membrane, solvent) comes from its position in the topology order: solute, then lipids, then water.
    """
    atoms, box = read_gro(out / structure)
    label = lambda a: f"{a['resname']}{a['resid']}:{a['atom']}"
    if not 1 <= number <= len(atoms):
        return f"atom {number} is not in {structure}"
    xyz, target = wrap(xyz_nm(atoms), box), atoms[number - 1]
    same = lambda a: (a["resid"], a["resname"]) == (target["resid"], target["resname"])
    near = cKDTree(xyz, boxsize=box[:3]).query_ball_point(xyz[number - 1], radius_nm)
    closest = {}
    for i in near:
        if not same(atoms[i]):
            d = float(np.linalg.norm((xyz[i] - xyz[number - 1] + 0.5 * np.array(box[:3])) % np.array(box[:3]) - 0.5 * np.array(box[:3])))
            key = f"{atoms[i]['resname']}{atoms[i]['resid']}"
            if key not in closest or d < closest[key][0]:
                closest[key] = (d, label(atoms[i]))
    neighbours = ", ".join(f"{atom} {10 * d:.2f} A" for d, atom in sorted(closest.values())[:8])
    return (f"largest force on atom {number} {label(target)} in {structure}; other residues within {10 * radius_nm:.1f} A: "
            + (neighbours or "none"))


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
        culprit = re.findall(r"Maximum force\s*=\s*\S+\s+on atom\s+(\d+)", text)
        where = (f"; {force_culprit(out, int(culprit[-1]))}; before EM: {closest_contact(out, int(culprit[-1]))}"
                 if culprit else "; em.log names no atom")
        raise SystemExit(f"EM final Fmax {float(fmax[-1]):.1f} >= {FMAX_TARGET:.0f} kJ/mol/nm{where}")
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
