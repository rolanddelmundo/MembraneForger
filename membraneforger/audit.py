#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Independent audit of a finished build directory.
#//=============================================================
"""Re-derive what a build claims from its files alone, with parsers that share no code with the builder."""
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.spatial import cKDTree

from .config import Settings
from .stereo import check_gm3_stereo

__all__ = ['gro_table', 'topology_includes', 'topology_molecules', 'itp_atoms', 'gromacs', 'check_topology',
           'check_grompp', 'check_energy', 'same_image', 'superposed_rmsd', 'check_structure', 'closest_contact_between', 'itp_bonds',
           'itp_dihedral_restraints', 'dihedrals_pbc', 'check_dihedral_restraints', 'molecule_rings',
           'ring_piercings', 'check_ring_piercing', 'audit_run']

PROTEIN = {"ALA", "ARG", "ASN", "ASP", "CYS", "CYS2", "CYSG", "CYSP", "GLN", "GLU", "GLY", "HIS", "HSD", "HSE", "HSP",
           "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "AIB", "LEM", "KTZ", "KRT", "KSM"}
WATER_IONS = {"TIP3", "SOD", "CLA", "POT"}


def gro_table(path: Path) -> tuple[list, list, np.ndarray, np.ndarray]:
    """Read residue names, atom names, coordinates (nm) and box from a .gro file."""
    lines = path.read_text().splitlines()
    n = int(lines[1])
    body = lines[2:2 + n]
    xyz = np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])] for l in body])
    return [l[5:10].strip() for l in body], [l[10:15].strip() for l in body], xyz, np.array(lines[2 + n].split(), float)


def topology_includes(top: Path, defines: set) -> list[Path]:
    """Follow #include lines the way grompp does (relative to the including file) for a set of defines."""
    seen, pending = [], [top]
    while pending:
        path = pending.pop(0)
        if path in seen:
            continue
        if not path.is_file():
            raise SystemExit(f"topology include {path} does not exist")
        seen.append(path)
        active, local = [], set(defines)
        for raw in path.read_text(errors="replace").splitlines():
            s = raw.split(";", 1)[0].split()
            if not s:
                continue
            if s[0] in ("#ifdef", "#ifndef"):
                active.append((s[1] in local) == (s[0] == "#ifdef"))
            elif s[0] == "#else":
                active[-1] = not active[-1]
            elif s[0] == "#endif":
                active.pop()
            elif all(active) and s[0] == "#define":
                local.add(s[1])
            elif all(active) and s[0] == "#include":
                name = s[1].strip('"')
                pending.append(next((p for p in (path.parent / name, top.parent / name) if p.is_file()), path.parent / name))
    return seen


def topology_molecules(top: Path) -> list[tuple[str, int]]:
    """Read the [ molecules ] table of a topology."""
    section, rows = "", []
    for raw in top.read_text().splitlines():
        s = raw.split(";", 1)[0].strip()
        if s.startswith("["):
            section = s.strip("[] ").lower()
        elif section == "molecules" and s and not s.startswith("#"):
            rows.append((s.split()[0], int(s.split()[1])))
    return rows


def itp_atoms(paths: list[Path]) -> dict:
    """Map each moleculetype defined in the given files to its (atom name, charge) list."""
    molecules, current, section = {}, None, ""
    for path in paths:
        for raw in path.read_text(errors="replace").splitlines():
            s = raw.split(";", 1)[0].strip()
            if not s or s.startswith("#"):
                continue
            if s.startswith("["):
                section, fresh = s.strip("[] \t").lower(), True
                continue
            f = s.split()
            if section == "moleculetype" and fresh:
                current, fresh = f[0], False
                molecules[current] = []
            elif section == "atoms" and current and len(f) >= 7 and f[0].isdigit():
                molecules[current].append((f[4], float(f[6])))
    return molecules


def gromacs(gmx: str, args: list, cwd: Path, stdin: str | None = None) -> tuple[int, str]:
    """Run one GROMACS tool and return its exit code and combined output."""
    result = subprocess.run(gmx.split() + args, cwd=cwd, text=True, capture_output=True, input=stdin)
    return result.returncode, result.stdout + result.stderr


def check_topology(run: Path, names: list) -> tuple[dict, list]:
    """Compare topol.top (includes, molecules, charge) with the coordinate atom names."""
    fails = []
    top = run / "topol.top"
    plain, restrained = topology_includes(top, set()), topology_includes(top, {"POSRES"})
    reached = {p.resolve() for p in plain + restrained}
    strays = sorted(str(p.relative_to(run)) for p in (run / "toppar").rglob("*") if p.is_file() and p.resolve() not in reached)
    if strays:
        fails.append(f"toppar files never included: {strays}")
    atoms = itp_atoms(plain)
    molecules = topology_molecules(top)
    unknown = [mol for mol, _ in molecules if mol not in atoms]
    if unknown:
        return {"molecules": molecules}, fails + [f"[ molecules ] entries without a moleculetype: {unknown}"]
    expected = [name for mol, n in molecules for _ in range(n) for name, _ in atoms[mol]]
    if expected != names:
        fails.append(f"coordinate atom names differ from topol.top ({len(names)} vs {len(expected)} atoms)")
    charge = sum(n * sum(q for _, q in atoms[mol]) for mol, n in molecules)
    if abs(charge) > 1e-3:
        fails.append(f"net charge {charge:+.4f} e")
    return {"molecules": molecules, "atoms": len(expected), "net_charge": round(charge, 4),
            "include_files": len(reached), "stray_toppar_files": strays}, fails


def check_grompp(run: Path, gmx: str, structure: str) -> tuple[dict, list]:
    """Rebuild the run input with -maxwarn 0, as built and with position restraints defined."""
    fails, report = [], {}
    scratch = Path(tempfile.mkdtemp(prefix="membraneforger_audit_"))
    try:
        for label, mdp_extra, extra in (("as_built", "", []), ("posres", "define = -DPOSRES\n", ["-r", str(run / structure)])):
            mdp = scratch / f"{label}.mdp"
            mdp.write_text(mdp_extra + (run / "em.mdp").read_text())
            code, output = gromacs(gmx, ["grompp", "-f", str(mdp), "-c", str(run / structure), "-p", str(run / "topol.top"),
                                         "-n", str(run / "index_ini.ndx"), "-o", str(scratch / f"{label}.tpr"),
                                         "-po", str(scratch / f"{label}_out.mdp"), "-maxwarn", "0"] + extra, run)
            warnings, notes = len(re.findall(r"(?m)^WARNING \d+", output)), len(re.findall(r"(?m)^NOTE \d+", output))
            report[label] = {"exit": code, "warnings": warnings, "notes": notes}
            if code or warnings or not (scratch / f"{label}.tpr").is_file():
                fails.append(f"grompp -maxwarn 0 ({label}) exit {code}, {warnings} warnings")
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    return report, fails


def check_energy(run: Path, gmx: str) -> tuple[dict, list]:
    """Check the energy and trajectory files with gmx check and compare gmx energy with em.log."""
    fails, report = [], {}
    for flag, name in (("-e", "em.edr"), ("-f", "em.trr")):
        code, output = gromacs(gmx, ["check", flag, name], run)
        frames = re.findall(r"(?m)^(?:Last energy frame read|Last frame)\s+(\d+)", output)
        report[f"gmx_check_{name}"] = {"exit": code, "last_frame": int(frames[-1]) if frames else None}
        if code or not frames:  # a truncated file can still exit 0, so a complete read must be reported
            fails.append(f"gmx check {flag} {name}: exit {code}, {'no' if not frames else 'a'} complete read reported")
    scratch = Path(tempfile.mkdtemp(prefix="membraneforger_audit_"))
    try:
        code, output = gromacs(gmx, ["energy", "-f", str(run / "em.edr"), "-o", str(scratch / "e.xvg")], run, "Potential\n\n")
        series = [float(l.split()[1]) for l in (scratch / "e.xvg").read_text().splitlines()
                  if l and l[0] not in "#@"] if code == 0 and (scratch / "e.xvg").is_file() else []
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    log_text = (run / "em.log").read_text(errors="replace")
    logged = re.findall(r"Potential Energy\s*=\s*(\S+)", log_text)
    fmax = re.findall(r"Maximum force\s*=\s*(\S+)", log_text)
    steps = re.findall(r"converged to Fmax < \S+ in (\d+) steps", log_text)
    if not series or not logged or not fmax:
        return report, fails + ["gmx energy or em.log gave no potential energy / maximum force"]
    report.update(potential_first=series[0], potential_last=series[-1], potential_em_log=float(logged[-1]),
                  fmax=float(fmax[-1]), em_steps=int(steps[-1]) if steps else None, energy_frames=len(series))
    if not steps:
        fails.append("em.log does not report convergence")
    if not np.isfinite(series).all() or series[-1] >= series[0]:
        fails.append(f"potential energy did not decrease ({series[0]:.4e} -> {series[-1]:.4e})")
    if abs(series[-1] - float(logged[-1])) > 1e-4 * abs(float(logged[-1])):
        fails.append(f"gmx energy ({series[-1]:.6e}) disagrees with em.log ({float(logged[-1]):.6e})")
    return report, fails


def same_image(reference: np.ndarray, xyz: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Shift every atom by whole box vectors to the periodic image closest to its reference position."""
    delta = xyz - reference
    return reference + delta - cell * np.round(delta / cell)


def superposed_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """RMSD between two coordinate sets after optimal rigid superposition."""
    P, Q = P - P.mean(axis=0), Q - Q.mean(axis=0)
    U, _, Vt = np.linalg.svd(P.T @ Q)
    R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
    return float(np.sqrt((((P @ R.T) - Q) ** 2).sum(axis=1).mean()))


def check_structure(run: Path, structure: str, settings: Settings) -> tuple[dict, list]:
    """Measure what EM did to the system: solute drift, lipid contacts, leaflets, box, solvent and ions."""
    fails = []
    resn, names, xyz, box = gro_table(run / structure)
    resn0, names0, xyz0, box0 = gro_table(run / "solv_ions.gro")
    boxed = gro_table(run / "boxed.gro")[3]
    if not np.isfinite(xyz).all():
        return {}, ["non-finite coordinates"]
    if names != names0:
        fails.append("minimized structure and solv_ions.gro list different atoms")
    cell = box[:3]
    protein = np.array([r in PROTEIN for r in resn])
    bulk = np.array([r in WATER_IONS for r in resn])
    heavy = np.array([not n.startswith("H") for n in names])
    index_text = (run / "index_ini.ndx").read_text()
    index_groups = re.findall(r"(?m)^\[ (\S+) \]", index_text)
    members = [np.array(block.split(), dtype=int) for block in re.split(r"(?m)^\[ \S+ \]", index_text)[1:]]
    sizes = dict(zip(index_groups, (len(m) for m in members)))
    solute = np.zeros(len(resn), bool)
    solute[:sizes.get("Protein_LIG", 0)] = True
    membrane = np.zeros(len(resn), bool)
    membrane[sizes.get("Protein_LIG", 0):sizes.get("Protein_LIG", 0) + sizes.get("MEMB", 0)] = True
    if index_groups != ["System", "Protein_LIG", "MEMB", "SOL_ION"]:
        return {"index_groups": sizes}, fails + [f"index groups are {index_groups}, not System/Protein_LIG/MEMB/SOL_ION"]
    everything = np.arange(1, len(resn) + 1)
    if not np.array_equal(members[0], everything) or not np.array_equal(np.concatenate(members[1:]), everything):
        fails.append(f"index groups are not consecutive blocks that partition atoms 1..{len(resn)}: {sizes}")
    if protein[~solute].any() or bulk[~(solute | membrane)].sum() != bulk.sum() or (~bulk[~(solute | membrane)]).any():
        fails.append("index groups hold the wrong kind of molecule")
    # protein CA RMSD across EM after optimal superposition
    ca = protein & np.array([n == "CA" for n in names])
    if ca.sum() < 3:
        return {"index_groups": sizes}, fails + ["no protein CA atoms found to measure drift across EM"]
    moved = same_image(xyz0, xyz, cell)  # mdrun may write a molecule in another periodic image
    ca_rmsd = superposed_rmsd(xyz0[ca], moved[ca])
    if ca_rmsd > settings.audit_max_ca_rmsd_nm:
        fails.append(f"protein CA RMSD across EM {ca_rmsd:.3f} nm > {settings.audit_max_ca_rmsd_nm} nm")
    # closest lipid-solute heavy-atom contact under periodic boundaries
    wrapped = np.mod(xyz, cell)
    wrapped[wrapped >= cell] = 0.0
    contact = float(cKDTree(wrapped[solute & heavy], boxsize=cell).query(wrapped[membrane & heavy])[0].min())
    if contact < settings.audit_min_contact_nm:
        fails.append(f"lipid heavy atom {contact:.3f} nm from the solute after EM (< {settings.audit_min_contact_nm} nm)")
    # leaflets from lipid phosphorus atoms
    phosphorus = membrane & np.array([n == "P" for n in names])
    z = wrapped[phosphorus][:, 2]
    middle = 0.5 * float(z.min() + z.max())  # leaflets may hold different numbers of lipids, so do not split at the median
    lower, upper = z[z < middle], z[z >= middle]
    if len(lower) == 0 or len(upper) == 0:
        fails.append("lipid phosphates do not form two leaflets")
        planes = [None, None]
    else:
        planes = [float(np.median(lower)), float(np.median(upper))]
        if max(lower.std(), upper.std()) > 0.6:
            fails.append(f"phosphate planes are not flat (z spread {lower.std():.2f}/{upper.std():.2f} nm)")
    water = np.array([n == "OH2" for n in names]) & np.array([r == "TIP3" for r in resn])
    core_water = 0
    if planes[0] is not None:
        zw = wrapped[water][:, 2]
        inner = (zw > planes[0] + 0.5) & (zw < planes[1] - 0.5)
        # the builder keeps cavity waters within 0.45 nm of the solute and EM moves atoms by a few tenths of a nm,
        # so only a core water more than 1.0 nm from any solute heavy atom is a water that should not be there
        near = np.array([bool(h) for h in cKDTree(wrapped[solute & heavy], boxsize=cell).query_ball_point(wrapped[water][inner], 1.0)])
        core_water = int(inner.sum() - near.sum())
        if core_water:
            fails.append(f"{core_water} waters inside the bilayer core away from the protein")
    if not np.allclose(box[:2], boxed[:2], atol=1e-3):
        fails.append("EM changed the membrane XY cell")
    solute_z = wrapped[solute][:, 2]
    report = {"atoms": len(resn), "box_nm": [float(v) for v in cell], "protein_ca_atoms": int(ca.sum()),
              "protein_ca_rmsd_nm": round(ca_rmsd, 4),
              "solute_max_displacement_nm": round(float(np.linalg.norm(moved[solute] - xyz0[solute], axis=1).max()), 4),
              "closest_lipid_solute_heavy_nm": round(contact, 4), "phosphate_planes_nm": planes,
              "leaflet_phosphates": [len(lower), len(upper)],
              "bilayer_thickness_nm": round(planes[1] - planes[0], 3) if planes[0] is not None else None,
              "water_in_bilayer_core": core_water, "waters": int(water.sum()),
              "sodium": int(sum(r == "SOD" for r in resn)), "chloride": int(sum(r == "CLA" for r in resn)),
              "water_padding_nm": [round(float(solute_z.min()), 3), round(float(cell[2] - solute_z.max()), 3)],
              "index_groups": sizes}
    return report, fails


def closest_contact_between(run: Path, structure: str, first_group_atoms: int) -> dict:
    """Find the closest pair between the first N atoms (the solute) and the rest, hydrogens included, with periodic images."""
    resn, names, xyz, box = gro_table(run / structure)
    cell = box[:3]
    wrapped = np.mod(xyz, cell)
    wrapped[wrapped >= cell] = 0.0
    distance, partner = cKDTree(wrapped[:first_group_atoms], boxsize=cell).query(wrapped[first_group_atoms:])
    k = int(np.argmin(distance))
    i, j = first_group_atoms + k, int(partner[k])
    return {"distance_nm": round(float(distance[k]), 4), "atoms": [f"{resn[i]}:{names[i]}(atom {i + 1})", f"{resn[j]}:{names[j]}(atom {j + 1})"]}


def itp_bonds(paths: list[Path]) -> dict:
    """Map each moleculetype defined in the given files to its bond list (zero-based atom pairs)."""
    molecules, current, section, fresh = {}, None, "", False
    for path in paths:
        for raw in path.read_text(errors="replace").splitlines():
            s = raw.split(";", 1)[0].strip()
            if not s or s.startswith("#"):
                continue
            if s.startswith("["):
                section, fresh = s.strip("[] \t").lower(), True
                continue
            f = s.split()
            if section == "moleculetype" and fresh:
                current, fresh = f[0], False
                molecules[current] = []
            elif section == "bonds" and current and len(f) >= 2 and f[0].isdigit() and f[1].isdigit():
                molecules[current].append((int(f[0]) - 1, int(f[1]) - 1))
    return molecules


def itp_dihedral_restraints(paths: list[Path]) -> dict:
    """Map each moleculetype to its [ dihedral_restraints ] rows: (zero-based atom quadruple, target angle in degrees)."""
    molecules, current, section, fresh = {}, None, "", False
    for path in paths:
        for raw in path.read_text(errors="replace").splitlines():
            s = raw.split(";", 1)[0].strip()
            if not s or s.startswith("#"):  # rows inside #ifdef DIHRES are read: the check needs them either way
                continue
            if s.startswith("["):
                section, fresh = s.strip("[] \t").lower(), True
                continue
            f = s.split()
            if section == "moleculetype" and fresh:
                current, fresh = f[0], False
                molecules[current] = []
            elif section == "dihedral_restraints" and current and len(f) >= 6 and all(v.isdigit() for v in f[:5]):
                molecules[current].append((tuple(int(v) - 1 for v in f[:4]), float(f[5])))
    return molecules


def dihedrals_pbc(xyz: np.ndarray, cell: np.ndarray, quads: np.ndarray) -> np.ndarray:
    """Dihedral angles (degrees) of atom quadruples, with minimum-image bond vectors in a rectangular cell."""
    image = lambda v: v - cell * np.round(v / cell)
    b0 = image(xyz[quads[:, 0]] - xyz[quads[:, 1]])
    b1 = image(xyz[quads[:, 2]] - xyz[quads[:, 1]])
    b2 = image(xyz[quads[:, 3]] - xyz[quads[:, 2]])
    b1 = b1 / np.linalg.norm(b1, axis=1)[:, None]
    v = b0 - (b0 * b1).sum(axis=1)[:, None] * b1
    w = b2 - (b2 * b1).sum(axis=1)[:, None] * b1
    return np.degrees(np.arctan2((np.cross(b1, v) * w).sum(axis=1), (v * w).sum(axis=1)))


def check_dihedral_restraints(run: Path, structure: str) -> tuple[dict, list]:
    """Evaluate every CHARMM-GUI dihedral restraint of the lipid topologies on the structure.

    The +-120 rows fix stereocentres and the 0/180 rows fix double bonds: a structure on the wrong side has the wrong
    molecule, which no MD run corrects, so it fails. The +-60 rows hold sugar and inositol rings in their chair;
    a ring outside it is a conformation MD can repair, so it is reported, not failed."""
    _, _, xyz, box = gro_table(run / structure)
    cell = box[:3]
    top = run / "topol.top"
    includes = topology_includes(top, set())
    atoms, rows = itp_atoms(includes), itp_dihedral_restraints(includes)
    offset, report, fails = 0, {}, []
    for mol, count in topology_molecules(top):
        n, restraints = len(atoms[mol]), rows.get(mol, [])
        if restraints and count:
            quads = np.array([q for q, _ in restraints])
            target = np.array([t for _, t in restraints])
            starts = offset + n * np.arange(count)
            phi = dihedrals_pbc(xyz, cell, (starts[:, None, None] + quads[None]).reshape(-1, 4)).reshape(count, -1)
            kind = np.where(np.abs(np.abs(target) - 120) < 30, "configuration",
                            np.where(np.abs(np.abs(target) - 60) < 30, "ring", "double_bond"))
            wrong = np.where(kind == "double_bond", (np.abs(phi) < 90) != (np.abs(target) < 90),
                             np.sign(phi) != np.sign(target))
            names = [a for a, _ in atoms[mol]]
            entry = {"molecules": int(count), "restraints": int(len(restraints))}
            for label in ("configuration", "double_bond"):
                bad = wrong[:, kind == label]
                entry[f"wrong_{label}s"] = int(bad.sum())
                if bad.any():
                    k, r = np.argwhere(bad)[0]
                    q = quads[kind == label][r]
                    fails.append(f"{mol}: {int(bad.sum())} {label.replace('_', ' ')} restraint(s) on the wrong side, e.g. molecule "
                                 f"{k + 1} {'-'.join(names[i] for i in q)} = {phi[k, kind == label][r]:+.0f} (target "
                                 f"{target[kind == label][r]:+.0f})")
            ring_rows = np.flatnonzero(kind == "ring")
            if len(ring_rows):
                graph = nx.Graph()
                for r in ring_rows:  # the rows of one ring share its atoms
                    graph.add_edges_from((r, ("atom", int(a))) for a in quads[r])
                rings = [sorted(x for x in c if not isinstance(x, tuple)) for c in nx.connected_components(graph)]
                out_of_chair = sum(int(wrong[:, ring].any(axis=1).sum()) for ring in rings)
                entry.update({"rings": len(rings) * int(count), "rings_out_of_chair": out_of_chair})
            report[mol] = entry
        offset += n * count
    return report, fails


def molecule_rings(bonds: list) -> list[tuple]:
    """Find the 5- and 6-membered rings of one molecule from its bond graph."""
    graph = nx.Graph(bonds)
    rings = []
    for component in nx.biconnected_components(graph):
        if len(component) < 5:
            continue
        sub = graph.subgraph(component)
        cycles = nx.minimum_cycle_basis(sub) if len(component) <= 60 else nx.cycle_basis(sub)
        for cycle in cycles:
            if len(cycle) in (5, 6):
                order = nx.find_cycle(sub.subgraph(cycle))
                rings.append(tuple(a for a, _ in order))
    return rings


def ring_piercings(xyz: np.ndarray, cell: np.ndarray, rings: np.ndarray, sizes: np.ndarray, bonds: np.ndarray) -> list:
    """Return (ring row, bond row) for every bond that passes through the inside of a ring under periodic boundaries."""
    # rings: (n, 6) atom indices padded with -1; a bond pierces a ring when its ends lie on opposite sides of the
    # ring plane and it crosses that plane inside the ring's inscribed circle.
    if not len(rings) or not len(bonds):
        return []
    mask = rings >= 0
    first = xyz[rings[:, 0]]
    members = np.where(mask[..., None], xyz[np.where(mask, rings, 0)] - first[:, None, :], 0.0)
    members -= cell * np.round(members / cell) * mask[..., None]
    centre_rel = members.sum(axis=1) / sizes[:, None]
    centres = first + centre_rel
    spokes = np.where(mask[..., None], members - centre_rel[:, None, :], 0.0)
    radius = np.linalg.norm(spokes, axis=2).sum(axis=1) / sizes
    normals = np.cross(spokes[:, 0], spokes[:, 1]) + np.cross(spokes[:, 1], spokes[:, 2])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    a = xyz[bonds[:, 0]]
    half = xyz[bonds[:, 1]] - a
    half -= cell * np.round(half / cell)
    middle = np.mod(a + 0.5 * half, cell)
    middle[middle >= cell] = 0.0
    wrapped = np.mod(centres, cell)
    wrapped[wrapped >= cell] = 0.0
    found = []
    # a bond's midpoint lies within half its length of any point it crosses, so search that far (at least 0.25 nm);
    # stretched pre-EM bonds are therefore not missed
    reach = max(0.25, 0.5 * float(np.linalg.norm(half, axis=1).max()) + 0.02)
    for r, near in enumerate(cKDTree(middle, boxsize=cell).query_ball_point(wrapped, min(reach, 0.49 * float(cell.min())))):
        ring_atoms = set(rings[r][mask[r]].tolist())
        for b in near:
            if bonds[b, 0] in ring_atoms or bonds[b, 1] in ring_atoms:
                continue
            start = a[b] - centres[r]
            start -= cell * np.round(start / cell)
            end = start + half[b]
            da, db = float(start @ normals[r]), float(end @ normals[r])
            if da * db >= 0:
                continue
            crossing = start + (da / (da - db)) * half[b]
            if np.linalg.norm(crossing) < radius[r] * np.cos(np.pi / sizes[r]):
                found.append((r, b))
    return found


def check_ring_piercing(run: Path, structure: str) -> tuple[dict, list]:
    """Fail if any covalent bond threads a 5- or 6-membered ring (lipid tails through aromatic, sterol or sugar rings)."""
    resn, names, xyz, box = gro_table(run / structure)
    top = run / "topol.top"
    files = topology_includes(top, set())
    bonds_of = itp_bonds(files)
    natoms = {mol: len(atoms) for mol, atoms in itp_atoms(files).items()}
    rings, sizes, bonds, owner, offset = [], [], [], [], 0
    for mol, count in topology_molecules(top):
        if mol in WATER_IONS or not bonds_of.get(mol):
            offset += count * natoms[mol]
            continue
        local_bonds = np.array(bonds_of[mol])
        local_rings = molecule_rings(bonds_of[mol])
        starts = offset + natoms[mol] * np.arange(count)
        bonds.append((local_bonds[None, :, :] + starts[:, None, None]).reshape(-1, 2))
        for ring in local_rings:
            padded = np.array(list(ring) + [-1] * (6 - len(ring)))
            rings.append(np.where(padded >= 0, padded[None, :] + starts[:, None], -1))
            sizes.append(np.full(count, len(ring)))
        owner.append((mol, count, len(local_rings)))
        offset += count * natoms[mol]
    rings = np.vstack(rings) if rings else np.zeros((0, 6), int)
    sizes = np.concatenate(sizes) if sizes else np.zeros(0, int)
    bonds = np.vstack(bonds) if bonds else np.zeros((0, 2), int)
    found = ring_piercings(xyz, box[:3], rings, sizes, bonds)
    label = lambda i: f"{resn[i]}:{names[i]}(atom {i + 1})"
    pierced = [{"ring": [label(i) for i in rings[r][rings[r] >= 0]], "bond": [label(bonds[b, 0]), label(bonds[b, 1])],
                "ring_atoms": [int(i) + 1 for i in rings[r][rings[r] >= 0]],
                "bond_atoms": [int(bonds[b, 0]) + 1, int(bonds[b, 1]) + 1]} for r, b in found]
    report = {"structure": structure, "rings_checked": int(len(rings)), "bonds_checked": int(len(bonds)),
              "rings_per_molecule": {mol: n for mol, _, n in owner if n}, "pierced": pierced}
    fails = [f"{len(pierced)} ring(s) pierced by a covalent bond, e.g. {pierced[0]['bond'][0]}-{pierced[0]['bond'][1]} through "
             f"the ring of {pierced[0]['ring'][0]}"] if pierced else []
    return report, fails


DEFAULT_DATA = Path(__file__).resolve().parents[1] / "backmap_data"  # map.dat: the GM3 chirality definitions the stereo check uses


def audit_run(run: Path, gmx: str, structure: str = "em.gro", settings: Settings = Settings(), data: Path | None = None) -> dict:
    """Audit one build directory and write audit.json; raises if any check fails."""
    run = run.resolve()
    names = gro_table(run / structure)[1]
    report, failures = {"structure": structure}, []
    data = Path(data) if data else DEFAULT_DATA
    for key, (part, fails) in (("topology", check_topology(run, names)), ("grompp", check_grompp(run, gmx, structure)),
                               ("energy", check_energy(run, gmx)), ("geometry", check_structure(run, structure, settings)),
                               ("ring_piercing", check_ring_piercing(run, structure)),
                               ("dihedral_restraints", check_dihedral_restraints(run, structure)),
                               ("gm3_stereochemistry", check_gm3_stereo(run / structure, data))):
        report[key] = part
        failures += fails
    report["checks"] = ["atom names and count match topol.top", "net charge is zero", "no stray toppar files",
                        "grompp -maxwarn 0 as built", "grompp -maxwarn 0 with -DPOSRES", "gmx check em.edr/em.trr",
                        "gmx energy agrees with em.log", "EM converged and lowered the potential",
                        "index groups partition the system", "protein CA RMSD across EM", "lipid-solute contacts",
                        "two flat phosphate leaflets", "no water in the bilayer core", "membrane XY cell preserved",
                        "no covalent bond threads a 5- or 6-membered ring",
                        "lipid stereocentres and double bonds satisfy their CHARMM-GUI dihedral restraints",
                        "every GM3 has its 16 sugar stereocentres and 2 ceramide trans bonds in the configuration of map.dat"]
    report["failures"] = failures
    report["result"] = "FAIL" if failures else "PASS"
    (run / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    if failures:
        raise SystemExit("independent audit failed: " + "; ".join(failures))
    return report


if __name__ == "__main__":
    result = audit_run(Path(sys.argv[1]).resolve(), sys.argv[2] if len(sys.argv) > 2 else "gmx")
    print(f"{sys.argv[1]}\t{result['result']}")
