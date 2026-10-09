#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Energy minimization and its validation.
#//=============================================================
"""Run steepest-descent minimization into em.unverified.gro and validate the result before it may be published.

When the lipid topologies carry CHARMM-GUI dihedral restraints, a restrained minimization (emres) runs first and the
validated, unrestrained EM starts from its result."""
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .config import AMINO, DIHRES_EM_FC, EM_MDP, FMAX_TARGET
from .runtools import run_command, run_path
from .structio import element, read_gro, read_pdb, wrap, xyz_nm
from .validation import check_inputs_unchanged

__all__ = ['has_dihedral_restraints', 'run_em', 'closest_contact', 'force_culprit', 'stage_files', 'read_structure',
           'closest_partner', 'slice_number', 'trace_clash', 'validate_em']

# On an Fmax failure, the residues within this distance of the atom with the largest force are named.
FORCE_NEIGHBOUR_NM = 0.3
# trace_clash: a heavy-atom pair closer than CLASH_HEAVY_A, or any pair closer than CLASH_ANY_A, between two molecules is
# a clash (heavy-atom contacts in a relaxed CHARMM membrane are about 3 A, H...H about 2 A); pairs are searched within
# CLASH_SEARCH_NM.
CLASH_HEAVY_A, CLASH_ANY_A, CLASH_SEARCH_NM = 2.0, 1.0, 0.5

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
                                    "-n", "index_ini.ndx", "-o", f"{name}.tpr", "-po", f"{name}_mdout.mdp", "-maxwarn", "0"],
                produces=(f"{name}.tpr",))  # grompp's processed .mdp is kept with the other intermediates
    run_command(out, gmx.split() + ["mdrun", "-deffnm", name, "-c", result, "-ntmpi", "1", "-ntomp", str(ntomp)],
                produces=(result, f"{name}.log"))


def closest_contact(out: Path, number: int) -> str:
    """Name the atom closest to the one EM blew up on, to point at the overlap."""
    atoms, box = read_gro(run_path(out, "solv_ions.gro"))
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


def stage_files(out: Path) -> list[tuple[str, Path]]:
    """The structures a build leaves behind, in the order the stages wrote them (those present)."""
    stages = [("backmapped membrane (membrane.pdb)", run_path(out, "membrane.pdb")),
              ("after the lipid scan and topology (prot-memb.pdb)", run_path(out, "prot-memb.pdb")),
              ("boxed (boxed.gro)", run_path(out, "boxed.gro")), ("solvated (solv.gro)", run_path(out, "solv.gro")),
              ("with ions (solv_ions.gro)", run_path(out, "solv_ions.gro")), ("restrained EM (emres.gro)", run_path(out, "emres.gro")),
              ("final EM (em.unverified.gro)", run_path(out, "em.unverified.gro"))]
    return [(label, path) for label, path in stages if path.is_file()]


def read_structure(path: Path) -> tuple[list[dict], np.ndarray]:
    """Atoms with coordinates in nm and the orthorhombic cell (nm) of a .gro or a .pdb with CRYST1."""
    if path.suffix == ".gro":
        atoms, box = read_gro(path)
        return atoms, np.array(box[:3], dtype=float)
    atoms, cryst = read_pdb(path)
    cell = np.array([float(cryst[6:15]), float(cryst[15:24]), float(cryst[24:33])]) / 10.0 if cryst else np.array([1e4] * 3)
    return [dict(a, x=a["x"] / 10.0, y=a["y"] / 10.0, z=a["z"] / 10.0) for a in atoms], cell


def closest_partner(atoms: list[dict], cell: np.ndarray, resname: str, resid: int, atom: str | None = None) -> dict | None:
    """The closest atom of any other residue to one residue (or one of its atoms): pairs over all atoms and heavy atoms."""
    is_target = lambda a: a["resname"] == resname and a["resid"] == resid and a["resname"] not in AMINO
    target = [i for i, a in enumerate(atoms) if is_target(a)]
    if not target:
        return None
    xyz = np.mod(xyz_nm(atoms), cell)
    tree = cKDTree(xyz, boxsize=cell)
    label = lambda a: f"{a['resname']}{a['resid']}:{a['atom']}"
    best = {}
    for kind, mine in (("any", target), ("heavy", [i for i in target if element(atoms[i]["atom"]) != "H"]),
                       ("atom", [i for i in target if atoms[i]["atom"] == atom] if atom else [])):
        found = None
        for i in mine:
            for j in tree.query_ball_point(xyz[i], CLASH_SEARCH_NM):
                if is_target(atoms[j]) or (kind == "heavy" and element(atoms[j]["atom"]) == "H"):
                    continue
                d = float(np.linalg.norm((xyz[j] - xyz[i] + 0.5 * cell) % cell - 0.5 * cell))
                if found is None or d < found[0]:
                    found = (d, label(atoms[i]), label(atoms[j]))
        if found:
            best[kind] = {"distance_A": round(10 * found[0], 3), "mine": found[1], "partner": found[2]}
    return best


def slice_number(table: Path, resname: str, resid: int) -> int | None:
    """The coarse-grained slice number (work/cg_membrane.gro) of membrane.pdb lipid `resid`, from the backmapped table.

    backmapping.assemble_membrane_pdb numbers the backmapped lipids 1, 2, 3, ... in the table's order (residue n gets
    n % 9999 + 1); the table keeps each lipid's own slice number in its resid column. None when the lipid at that place
    is not `resname` (the names do not line up) or the table is shorter.
    """
    groups, last = [], None
    for line in table.read_text().splitlines()[1:]:
        chain, number, name = line.split("\t")[:3]
        if (chain, number, name) != last:
            groups.append((int(number), name))
            last = (chain, number, name)
    candidates = [g for n, g in enumerate(groups) if n % 9999 + 1 == resid and g[1] == resname]
    return candidates[0][0] if len(candidates) == 1 else None


def trace_clash(out: Path, resname: str, resid: int, atom: str | None = None) -> dict:
    """Follow one membrane residue through every saved stage and find where it first clashes with another molecule.

    membrane.pdb numbers the lipids 1, 2, 3, ... in the order backmapping wrote them, and every .gro keeps that number, so
    the same molecule is compared at each stage; its coarse-grained slice number is read from the backmapped table
    (slice_number). A clash is a heavy-atom pair closer than
    CLASH_HEAVY_A or any atom pair closer than CLASH_ANY_A. At the coarse-grained stage its beads are measured against
    the protein that backmapping holds rigid (work/rock.pdb).
    """
    rows = []
    cg, rock, table = run_path(out, "work") / "cg_membrane.gro", run_path(out, "work") / "rock.pdb", run_path(out, "work") / "membrane_aa.tsv"
    source = slice_number(table, resname, resid) if table.is_file() else None
    if source is not None and cg.is_file() and rock.is_file():
        beads, box = read_gro(cg)
        mine = np.array([[b["x"], b["y"], b["z"]] for b in beads if b["resid"] == source])
        protein = xyz_nm([a for a in read_pdb(rock)[0] if element(a["atom"]) != "H"]) / 10.0
        if len(mine) and len(protein):
            cell = np.array(box[:3], dtype=float)
            d = cKDTree(np.mod(protein, cell), boxsize=cell).query(np.mod(mine, cell))[0].min()
            rows.append({"stage": "coarse-grained slice (work/cg_membrane.gro)", "slice_residue": source,
                         "closest_protein_heavy_atom_to_a_bead_A": round(10 * d, 3)})
    for label, path in stage_files(out):
        atoms, cell = read_structure(path)
        contact = closest_partner(atoms, cell, resname, resid, atom)
        if contact is None:
            rows.append({"stage": label, "absent": True})
            continue
        clash = (contact.get("heavy", {}).get("distance_A", 99.0) < CLASH_HEAVY_A
                 or contact.get("any", {}).get("distance_A", 99.0) < CLASH_ANY_A)
        rows.append({"stage": label, **contact, "clash": bool(clash)})
    first = next((r for r in rows if r.get("clash")), None)
    return {"residue": f"{resname}{resid}", "atom": atom, "stages": rows,
            "first_clash": {"stage": first["stage"], "partner": (first.get("heavy") or first["any"])["partner"]} if first else None}


def validate_em(system: dict, topology: dict, index: dict) -> dict:
    """Check EM converged and that the minimized structure, topology, index and box are all consistent."""
    out = system["out"]
    check_inputs_unchanged(system["inputs"])
    text = (run_path(out, "em.log")).read_text(errors="replace")
    fmax, energy = re.findall(r"Maximum force\s*=\s*(\S+)", text), re.findall(r"Potential Energy\s*=\s*(\S+)", text)
    if not fmax or not energy or not math.isfinite(float(energy[-1])) or not math.isfinite(float(fmax[-1])):
        culprit = re.findall(r"on atom (\d+)", text)
        raise SystemExit(f"non-finite EM force: {closest_contact(out, int(culprit[-1])) if culprit else 'no force report in em.log'}")
    if not float(fmax[-1]) < FMAX_TARGET:
        culprit = re.findall(r"Maximum force\s*=\s*\S+\s+on atom\s+(\d+)", text)
        where = (f"; {force_culprit(out, int(culprit[-1]))}; before EM: {closest_contact(out, int(culprit[-1]))}"
                 if culprit else "; em.log names no atom")
        if culprit:
            atoms = read_gro(run_path(out, "em.unverified.gro"))[0]
            number = int(culprit[-1])
            if 1 <= number <= len(atoms) and atoms[number - 1]["resname"] not in AMINO:
                a = atoms[number - 1]
                trace = trace_clash(out, a["resname"], a["resid"], a["atom"])
                (out / "em_clash_trace.json").write_text(json.dumps(trace, indent=2) + "\n")
                first = trace["first_clash"]
                where += (f"; {trace['residue']} first clashes in the {first['stage']} stage, with {first['partner']} "
                          "(every stage: em_clash_trace.json)" if first else
                          f"; {trace['residue']} clashes with no molecule at any saved stage (em_clash_trace.json)")
        raise SystemExit(f"EM final Fmax {float(fmax[-1]):.1f} >= {FMAX_TARGET:.0f} kJ/mol/nm{where}")
    gro = (run_path(out, "em.unverified.gro")).read_text().splitlines()
    atoms = gro[2:2 + int(gro[1])]
    xyz = np.array([[float(line[20:28]), float(line[28:36]), float(line[36:44])] for line in atoms])
    if [line[10:15].strip() for line in atoms] != topology["names"] or not np.isfinite(xyz).all():
        raise SystemExit("minimized atoms do not match topol.top")
    if index["System"] != len(atoms) or index["Protein_LIG"] + index["MEMB"] + index["SOL_ION"] != len(atoms):
        raise SystemExit(f"index groups do not partition the system: {index}")
    if not index["SOL_ION"] or abs(topology["charge"]) > 1e-3:
        raise SystemExit("EM system is not solvated and neutral")
    if not np.allclose([float(v) for v in gro[-1].split()][:2], read_gro(run_path(out, "boxed.gro"))[1][:2], atol=1e-3):
        raise SystemExit("EM changed the membrane XY box")
    steps = re.findall(r"converged to Fmax < \S+ in (\d+) steps", text)
    return {"summary": f"Fmax {float(fmax[-1]):.1f}, Epot {float(energy[-1]):.4e}, {steps[-1] if steps else '?'} steps",
            "fmax_kj_mol_nm": float(fmax[-1]), "potential_kj_mol": float(energy[-1]),
            "steps": int(steps[-1]) if steps else None, "fmax_target": FMAX_TARGET, "atoms": len(atoms)}
