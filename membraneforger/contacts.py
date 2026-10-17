#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Force-field-aware contact validator and bond-crossing detector.
#//=============================================================
"""Find atom pairs that overlap by the force field's own measure, and bonds that pass through each other.

Steepest descent cannot undo two molecules threaded through each other, and it stops at a finite Fmax when two atoms
are held on top of each other by their bonds. Both are geometry defects that exist before minimization and can be found
there. This module reads the topology grompp would read (topol.top and its includes) and tests a structure the way
GROMACS sees it: every pair of atoms that is not excluded (different molecules, or the same molecule more than nrexcl
bonds apart) is a non-bonded pair, and its Lennard-Jones minimum Rmin_ij comes from the pair's own parameters
([ nonbond_params ], the NBFIX table) or the combination rule over the two atom types. A pair closer than
contact_rmin_fraction * Rmin_ij is a contact violation: at 0.6 Rmin the LJ repulsion alone is several thousand
kJ/mol/nm for heavy atoms, far above the EM target, while no equilibrium contact (hydrogen bonds, packed tails) comes
below about 0.7 Rmin. Atom types without LJ repulsion (epsilon 0, such as TIP3 hydrogens in some water models) are left
out, as GROMACS leaves them out. Pairs of two solvent molecules (water-water, ion-water, ion-ion) are not examined:
gmx genion puts an ion on a water's site, so its neighbours sit at water-water distances, which for a chloride and an
oxygen is 0.6-0.7 of their LJ Rmin (a 2.4-2.8 A neighbour against a 4.0 A Rmin), and a free water or ion relaxes in
the first steps of EM; what EM cannot undo is a solvent molecule caught in the solute or the membrane, which is tested.

A bond crossing is a pair of bonds of two different molecules whose closest approach is below bond_crossing_nm with
the closest points inside both bonds (not at an atom): one chain passes through the other. Endpoint-to-endpoint
distances can miss this; the segment test does not.

The gates in the build: boxed.gro right after backmapping and topology (backmap_verdict retries another seed),
solv_ions.gro before EM (water and ions were placed by the same measure), and em.unverified.gro after EM.
"""
import json
import math
from collections import deque
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .audit import WATER_IONS, topology_includes, topology_molecules
from .config import Settings
from .runtools import run_path
from .structio import element, read_gro, read_pdb, xyz_nm

__all__ = ['lj_table', 'itp_molecules', 'system_tables', 'rmin_matrix', 'overlapping', 'pair_limits', 'violating_pairs', 'bond_crossings',
           'check_contacts', 'describe_contacts', 'record_contacts', 'check_fraction', 'CONTACTS_JSON']

CONTACTS_JSON = "contacts.json"  # one file: the report of every gate, by stage


def lj_table(paths: list[Path]) -> dict:
    """Lennard-Jones parameters of a topology: {"sigma": {type: nm}, "epsilon": {type: kJ/mol}, "pairs": {(t1, t2): sigma}}.

    Reads every [ atomtypes ] and [ nonbond_params ] block of the given files (as grompp would), and the combination
    rule from [ defaults ]. Only rules 2 and 3 (sigma/epsilon tables) are understood; CHARMM36 uses rule 2.
    """
    table, section = {"sigma": {}, "epsilon": {}, "pairs": {}, "comb_rule": 2}, ""
    for path in paths:
        for raw in path.read_text(errors="replace").splitlines():
            s = raw.split(";", 1)[0].strip()
            if not s or s.startswith("#"):
                continue
            if s.startswith("["):
                section = s.strip("[] \t").lower()
                continue
            f = s.split()
            if section == "defaults" and len(f) >= 2 and f[0].isdigit():
                table["comb_rule"] = int(f[1])
            elif section == "atomtypes" and len(f) >= 6:
                # name [bonded type] [at.num] mass charge ptype sigma epsilon: the last two numbers are sigma, epsilon
                try:
                    table["sigma"][f[0]], table["epsilon"][f[0]] = float(f[-2]), float(f[-1])
                except ValueError:
                    continue
            elif section == "nonbond_params" and len(f) >= 5:
                try:
                    table["pairs"][(f[0], f[1])] = table["pairs"][(f[1], f[0])] = (float(f[3]), float(f[4]))
                except ValueError:
                    continue
    if table["comb_rule"] not in (2, 3):
        raise SystemExit(f"combination rule {table['comb_rule']} (C6/C12 tables) is not supported by the contact check")
    return table


def itp_molecules(paths: list[Path]) -> dict:
    """Every moleculetype of the given files: {name: {"names": [...], "types": [...], "bonds": [(i, j)], "nrexcl": n}}."""
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
                molecules[current] = {"names": [], "types": [], "bonds": [], "nrexcl": int(f[1]) if len(f) > 1 and f[1].isdigit() else 3}
            elif section == "atoms" and current and len(f) >= 5 and f[0].isdigit():
                molecules[current]["names"].append(f[4])
                molecules[current]["types"].append(f[1])
            elif section in ("bonds", "settles") and current and len(f) >= 2 and f[0].isdigit() and f[1].isdigit():
                if section == "bonds":
                    molecules[current]["bonds"].append((int(f[0]) - 1, int(f[1]) - 1))
                else:  # SETTLE water: O-H1, O-H2 are the bonds, O is the first atom
                    o = int(f[0]) - 1
                    molecules[current]["bonds"] += [(o, o + 1), (o, o + 2)]
    return molecules


def excluded_keys(natoms: int, bonds: list, nrexcl: int) -> np.ndarray:
    """Keys i * natoms + j (i < j) of the atom pairs within nrexcl bonds of each other: GROMACS excludes them."""
    neighbours = [[] for _ in range(natoms)]
    for a, b in bonds:
        neighbours[a].append(b)
        neighbours[b].append(a)
    keys = []
    for start in range(natoms):
        seen, queue = {start: 0}, deque([start])
        while queue:
            atom = queue.popleft()
            if seen[atom] == nrexcl:
                continue
            for other in neighbours[atom]:
                if other not in seen:
                    seen[other] = seen[atom] + 1
                    queue.append(other)
        keys += [start * natoms + other for other in seen if other > start]
    return np.array(sorted(keys), dtype=np.int64)


def system_tables(top: Path) -> dict:
    """Per-atom tables of a whole system in topology order: types, molecule ids, bonds and exclusions, plus LJ data."""
    files = topology_includes(top, set())
    molecules, lj, rows = itp_molecules(files), lj_table(files), topology_molecules(top)
    types, mol_of, moltype_of, bonds, offset, molecule, first_atom, mol_n, excluded_by_type = [], [], [], [], 0, 0, [], [], {}
    for name, count in rows:
        if name not in molecules:
            raise SystemExit(f"moleculetype {name} of {top.name} is not defined in its includes")
        m = molecules[name]
        n = len(m["names"])
        local_bonds = np.array(m["bonds"], dtype=np.int64).reshape(-1, 2)
        if name not in excluded_by_type:
            excluded_by_type[name] = excluded_keys(n, m["bonds"], m["nrexcl"]) if n > 1 else np.zeros(0, dtype=np.int64)
        for _ in range(count):
            types += m["types"]
            mol_of.append(np.full(n, molecule, dtype=np.int64))
            moltype_of.append(name)
            first_atom.append(offset)
            mol_n.append(n)
            if len(local_bonds):
                bonds.append(local_bonds + offset)
            offset += n
            molecule += 1
    type_names = sorted(set(types))
    index, rmin = rmin_matrix(lj, type_names)
    return {"types": np.array([index[t] for t in types], dtype=np.int64), "type_names": type_names, "rmin": rmin,
            "molecule": np.concatenate(mol_of) if mol_of else np.zeros(0, dtype=np.int64), "moltype": moltype_of,
            "first_atom": np.array(first_atom, dtype=np.int64), "mol_n": np.array(mol_n, dtype=np.int64),
            "bonds": np.vstack(bonds) if bonds else np.zeros((0, 2), dtype=np.int64), "excluded_by_type": excluded_by_type,
            "natoms": offset}


def rmin_matrix(lj: dict, type_names: list[str]) -> tuple[dict, np.ndarray]:
    """Index of the given atom types and their LJ Rmin_ij matrix (nm; 0 where a pair has no LJ repulsion)."""
    index = {t: i for i, t in enumerate(type_names)}
    missing = [t for t in type_names if t not in lj["sigma"]]
    if missing:
        raise SystemExit(f"atom type(s) {missing} have no [ atomtypes ] entry in the topology includes")
    sigma = np.array([lj["sigma"][t] for t in type_names])
    epsilon = np.array([lj["epsilon"][t] for t in type_names])
    sigma_ij = 0.5 * (sigma[:, None] + sigma[None, :]) if lj["comb_rule"] == 2 else np.sqrt(sigma[:, None] * sigma[None, :])
    repulsive = (epsilon[:, None] > 0) & (epsilon[None, :] > 0)
    for (a, b), (s, e) in lj["pairs"].items():
        if a in index and b in index:
            sigma_ij[index[a], index[b]] = s
            repulsive[index[a], index[b]] = e > 0
    return index, np.where(repulsive, 2.0 ** (1.0 / 6.0) * sigma_ij, 0.0)


def overlapping(probe_xyz: np.ndarray, probe_types: list[str], ref_xyz: np.ndarray, ref_types: list[str], cell: np.ndarray,
                lj: dict, fraction: float) -> tuple[np.ndarray, tuple | None]:
    """Which probe atoms sit inside fraction * Rmin_ij of any reference atom (periodic), and the worst pair found.

    Used while placing solvent: every water atom (probe) against every solute atom (reference), and every ion against
    everything else. Returns a boolean per probe atom and (probe index, reference index, distance_nm, limit_nm) of the
    pair deepest inside its limit, or None.
    """
    names = sorted(set(probe_types) | set(ref_types))
    index, rmin = rmin_matrix(lj, names)
    pt = np.array([index[t] for t in probe_types], dtype=np.int64)
    rt = np.array([index[t] for t in ref_types], dtype=np.int64)
    cell = np.asarray(cell, dtype=float)[:3]
    probe, ref = np.mod(probe_xyz, cell), np.mod(ref_xyz, cell)
    probe[probe >= cell], ref[ref >= cell] = 0.0, 0.0
    hits = cKDTree(ref, boxsize=cell).query_ball_point(probe, float(fraction * rmin.max()))
    ia = np.repeat(np.arange(len(probe)), [len(h) for h in hits])
    ib = np.concatenate([np.asarray(h, dtype=np.int64) for h in hits]) if len(ia) else np.zeros(0, dtype=np.int64)
    bad = np.zeros(len(probe), dtype=bool)
    if not len(ia):
        return bad, None
    d = np.linalg.norm(minimum_image(ref[ib] - probe[ia], cell), axis=1)
    limit = fraction * rmin[pt[ia], rt[ib]]
    inside = (limit > 0) & (d < limit)
    if not inside.any():
        return bad, None
    bad[ia[inside]] = True
    k = int(np.argmin(np.where(inside, d / np.where(limit > 0, limit, 1.0), np.inf)))
    return bad, (int(ia[k]), int(ib[k]), float(d[k]), float(limit[k]))


def pair_limits(tables: dict, i: np.ndarray, j: np.ndarray, fraction: float) -> np.ndarray:
    """The contact limit (nm) of atom pairs (i, j): fraction of their LJ Rmin; 0 where the pair has no LJ repulsion."""
    return fraction * tables["rmin"][tables["types"][i], tables["types"][j]]


def minimum_image(delta: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Shortest periodic displacement vectors of an orthorhombic cell."""
    return delta - cell * np.round(delta / cell)


def is_excluded(tables: dict, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """Whether pairs (i, j) are topology exclusions: the same molecule, within nrexcl bonds of each other."""
    out = np.zeros(len(i), dtype=bool)
    same = np.flatnonzero(tables["molecule"][i] == tables["molecule"][j])
    if not len(same):
        return out
    mol = tables["molecule"][i[same]]
    start, n = tables["first_atom"][mol], tables["mol_n"][mol]
    a, b = np.minimum(i[same], j[same]) - start, np.maximum(i[same], j[same]) - start
    key = a * n + b
    names = np.array([tables["moltype"][m] for m in mol])
    for name in set(names):
        here = names == name
        out[same[here]] = np.isin(key[here], tables["excluded_by_type"][name])
    return out


def violating_pairs(xyz: np.ndarray, cell: np.ndarray, tables: dict, fraction: float) -> list[tuple]:
    """Non-excluded atom pairs closer than fraction * Rmin_ij: (i, j, distance_nm, limit_nm), closest ratio first."""
    search = float(fraction * tables["rmin"].max())
    wrapped = np.mod(xyz, cell)
    wrapped[wrapped >= cell] = 0.0
    pairs = cKDTree(wrapped, boxsize=cell).query_pairs(search, output_type="ndarray")
    if not len(pairs):
        return []
    i, j = pairs[:, 0], pairs[:, 1]
    d = np.linalg.norm(minimum_image(wrapped[j] - wrapped[i], cell), axis=1)
    limit = pair_limits(tables, i, j, fraction)
    hit = (limit > 0) & (d < limit)
    i, j, d, limit = i[hit], j[hit], d[hit], limit[hit]
    keep = ~is_excluded(tables, i, j)
    order = np.argsort(d[keep] / limit[keep])
    return [(int(a), int(b), float(x), float(y)) for a, b, x, y in zip(i[keep][order], j[keep][order], d[keep][order], limit[keep][order])]


def segment_distances(p1: np.ndarray, p2: np.ndarray, q1: np.ndarray, q2: np.ndarray, cell: np.ndarray) -> tuple:
    """Closest approach of segments p1-p2 and q1-q2 (minimum image on the segment offset): distance, s, t in [0, 1]."""
    u, v = p2 - p1, q2 - q1
    w0 = minimum_image(p1 - q1, cell)
    a, b, c = (u * u).sum(1), (u * v).sum(1), (v * v).sum(1)
    d, e = (u * w0).sum(1), (v * w0).sum(1)
    den = a * c - b * b
    s = np.where(den > 1e-12, (b * e - c * d) / np.where(den > 1e-12, den, 1.0), 0.0).clip(0.0, 1.0)
    t = np.where(c > 1e-12, (b * s + e) / np.where(c > 1e-12, c, 1.0), 0.0).clip(0.0, 1.0)
    s = np.where(a > 1e-12, (b * t - d) / np.where(a > 1e-12, a, 1.0), 0.0).clip(0.0, 1.0)
    gap = w0 + s[:, None] * u - t[:, None] * v
    return np.linalg.norm(gap, axis=1), s, t


def bond_crossings(xyz: np.ndarray, cell: np.ndarray, tables: dict, threshold_nm: float, interior: float = 0.1) -> list[tuple]:
    """Bond pairs of different molecules whose closest approach is below threshold_nm inside both bonds.

    Returns (bond_a, bond_b, distance_nm) with bonds as (i, j) atom index pairs, closest first. `interior` keeps the
    closest points at least that fraction of a bond length away from either atom: an approach at an atom is a contact
    (reported by violating_pairs), an approach between the atoms is a crossing. Bonds of water and ions are left out.
    """
    bonds = tables["bonds"]
    if len(bonds):  # water and ions cannot be threaded; a water bond across a lipid bond puts its O inside a contact limit anyway
        solvent = np.array([tables["moltype"][m] in WATER_IONS for m in range(len(tables["moltype"]))])
        bonds = bonds[~solvent[tables["molecule"][bonds[:, 0]]]]
    if not len(bonds):
        return []
    wrapped = np.mod(xyz, cell)
    wrapped[wrapped >= cell] = 0.0
    p1, p2 = wrapped[bonds[:, 0]], wrapped[bonds[:, 0]] + minimum_image(wrapped[bonds[:, 1]] - wrapped[bonds[:, 0]], cell)
    mid, half = np.mod(0.5 * (p1 + p2), cell), 0.5 * np.linalg.norm(p2 - p1, axis=1)
    mid[mid >= cell] = 0.0
    pairs = cKDTree(mid, boxsize=cell).query_pairs(float(2.0 * half.max() + threshold_nm), output_type="ndarray")
    if not len(pairs):
        return []
    a, b = pairs[:, 0], pairs[:, 1]
    other = tables["molecule"][bonds[a, 0]] != tables["molecule"][bonds[b, 0]]
    a, b = a[other], b[other]
    if not len(a):
        return []
    d, s, t = segment_distances(p1[a], p2[a], p1[b], p2[b], cell)
    hit = (d < threshold_nm) & (s > interior) & (s < 1.0 - interior) & (t > interior) & (t < 1.0 - interior)
    order = np.argsort(d[hit])
    return [((int(bonds[x, 0]), int(bonds[x, 1])), (int(bonds[y, 0]), int(bonds[y, 1])), float(z))
            for x, y, z in zip(a[hit][order], b[hit][order], d[hit][order])]


def read_coordinates(path: Path) -> tuple[list[dict], np.ndarray]:
    """Atoms (nm) and the orthorhombic cell (nm) of a .gro, or of a .pdb with CRYST1."""
    if path.suffix == ".gro":
        atoms, box = read_gro(path)
        return atoms, np.array(box[:3], dtype=float)
    atoms, cryst = read_pdb(path)
    if not cryst:
        raise SystemExit(f"{path.name} has no CRYST1 box")
    cell = np.array([float(cryst[6:15]), float(cryst[15:24]), float(cryst[24:33])]) / 10.0
    return [dict(a, x=a["x"] / 10.0, y=a["y"] / 10.0, z=a["z"] / 10.0) for a in atoms], cell


def check_contacts(out: Path, structure: str, solute_atoms: int, settings: Settings = Settings(), top: Path | None = None) -> dict:
    """Test one structure of a build against its topology: contact violations and bond crossings, with who is involved.

    `solute_atoms` is the number of leading atoms that are protein/ligand (topology order puts them first). Each finding
    is classed by the molecules it involves: "solvent" (a water or ion against the solute or a lipid: it can be replaced),
    "membrane" (a lipid is involved: another backmapping seed can change it), "solute" (both atoms in the all-atom
    input). Pairs of two solvent molecules are not findings (see the module docstring).
    """
    check_fraction(settings.contact_rmin_fraction)
    path = run_path(out, structure)
    atoms, cell = read_coordinates(path)
    tables = system_tables(top or run_path(out, "topol.top"))
    if len(atoms) != tables["natoms"]:
        raise SystemExit(f"{path.name} has {len(atoms)} atoms, topol.top {tables['natoms']}")
    xyz = xyz_nm(atoms)
    label = lambda i: f"{atoms[i]['resname']}{atoms[i]['resid']}:{atoms[i]['atom']}({element(atoms[i]['atom'])}, atom {i + 1})"
    solvent = lambda i: tables["moltype"][tables["molecule"][i]] in WATER_IONS

    def classify(*idx) -> str:  # a water or ion can be replaced, a lipid re-backmapped, the solute is the input's own
        if any(solvent(i) for i in idx):
            return "solvent"
        return "membrane" if any(i >= solute_atoms for i in idx) else "solute"

    contacts = [{"atoms": [label(i), label(j)], "distance_A": round(10 * d, 3), "limit_A": round(10 * lim, 3),
                 "rmin_A": round(10 * lim / settings.contact_rmin_fraction, 3), "class": classify(i, j),
                 "molecules": [int(tables["molecule"][i]) + 1, int(tables["molecule"][j]) + 1]}
                for i, j, d, lim in violating_pairs(xyz, cell, tables, settings.contact_rmin_fraction)
                if not (solvent(i) and solvent(j))]
    crossings = [{"bonds": [[label(a[0]), label(a[1])], [label(b[0]), label(b[1])]], "distance_A": round(10 * d, 3),
                  "class": classify(a[0], b[0]), "molecules": [int(tables["molecule"][a[0]]) + 1, int(tables["molecule"][b[0]]) + 1]}
                 for a, b, d in bond_crossings(xyz, cell, tables, settings.bond_crossing_nm)]
    by_class = lambda rows: {c: sum(r["class"] == c for r in rows) for c in ("solute", "membrane", "solvent")}
    return {"structure": structure, "atoms": len(atoms), "molecules": len(tables["moltype"]),
            "limits": {"contact_rmin_fraction": settings.contact_rmin_fraction, "bond_crossing_A": round(10 * settings.bond_crossing_nm, 2),
                       "nonbonded_rule": "pairs not excluded by the topology (different molecules, or more than nrexcl bonds apart)"},
            "contacts": contacts, "crossings": crossings,
            "contacts_by_class": by_class(contacts), "crossings_by_class": by_class(crossings)}


def describe_contacts(report: dict, classes: tuple = ("solute", "membrane", "solvent"), limit: int = 3) -> str:
    """One sentence naming the worst findings of the given classes, for a log line or a failure message."""
    contacts = [c for c in report["contacts"] if c["class"] in classes]
    crossings = [c for c in report["crossings"] if c["class"] in classes]
    parts = []
    if contacts:
        worst = "; ".join(f"{c['atoms'][0]} and {c['atoms'][1]} are {c['distance_A']:.2f} A apart (limit {c['limit_A']:.2f} A, "
                          f"{report['limits']['contact_rmin_fraction']:g} of their LJ Rmin {c['rmin_A']:.2f} A)" for c in contacts[:limit])
        parts.append(f"{len(contacts)} non-bonded pair(s) inside their LJ contact limit: {worst}")
    if crossings:
        worst = "; ".join(f"bond {c['bonds'][0][0]}-{c['bonds'][0][1]} passes {c['distance_A']:.2f} A from bond "
                          f"{c['bonds'][1][0]}-{c['bonds'][1][1]}" for c in crossings[:limit])
        parts.append(f"{len(crossings)} bond crossing(s) (chains threaded through each other; minimization cannot undo this): {worst}")
    return "; ".join(parts)


def record_contacts(out: Path, stage: str, report: dict) -> Path:
    """Add one gate's report to <out>/contacts.json under its stage name (the file holds every gate of the build)."""
    path = out / CONTACTS_JSON
    data = json.loads(path.read_text()) if path.is_file() else {}
    data[stage] = report
    path.write_text(json.dumps(data, indent=2) + "\n")
    return path


def check_fraction(fraction: float) -> None:
    """Refuse a contact fraction that would flag equilibrium contacts (above about 0.75 Rmin) or nothing (at or below 0)."""
    if not 0.0 < fraction <= 0.75 or not math.isfinite(fraction):
        raise SystemExit(f"Settings.contact_rmin_fraction must be in (0, 0.75], not {fraction}")
