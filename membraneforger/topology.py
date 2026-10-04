#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Structure preparation and CHARMM36 topology construction.
#//=============================================================
"""Split an assembled PDB into molecules, repair it, and build toppar/, topol.top and topology-ordered coordinates."""
import itertools
import math
import os
import pty
import re
import select
import shutil
import subprocess
import time
from collections import OrderedDict, defaultdict, deque
from pathlib import Path

import networkx as nx
import numpy as np
from networkx.algorithms import isomorphism as iso
from scipy.spatial import cKDTree

from .config import (AMINO, C_CAP, CYSG_HDB, DISULFIDE_MAX_A, DISULFIDE_OK_A, GLPA_MAX_BOND_A, GM3_XML_TO_GLPA, H_BOND_RANGE,
                     LIPID_DELETE_A, LIPID_SCAN_A, LIPIDATED, LYS_BACKBONE, N_CAP, RENAME_ATOM, RENAME_MOLECULE, RENAME_RESIDUE,
                     SOLVENT, TERMINUS_MENU)
from .runtools import log, sha256
from .structio import dist, element, read_itp, read_pdb, residues_in_order, source_element, write_pdb, xyz_nm

__all__ = ['prepare_structure', 'find_disulfides', 'repair_structure', 'set_disulfide_geometry', 'lipid_clashes',
           'resolve_lipid_clashes', 'setup_forcefield', 'add_lipidated_residue', 'charmm_patches', 'chain_termini',
           'split_caps', 'run_pdb2gmx', 'restore_source_coordinates', 'moleculetype_block', 'bonded_terms',
           'patch_caps', 'write_posre', 'check_disulfides', 'build_protein', 'place_h', 'h_parents', 'check_h_bonds',
           'match_by_name', 'match_by_graph', 'match_glpa', 'molecule_itp', 'forcefield_includes', 'build_topology']

def prepare_structure(source: Path, out: Path, name: str, forcefield: Path, ligands: set) -> dict:
    """Split the input into protein chains and other residues (ligands, lipids), keeping altloc A."""
    atoms, cryst1 = read_pdb(source)
    if not cryst1:
        raise SystemExit(f"{source.name} has no CRYST1 box")
    chains: OrderedDict[str, list[dict]] = OrderedDict()
    others: OrderedDict[tuple, list[dict]] = OrderedDict()
    for a in atoms:
        if a["altloc"] not in " A":
            continue
        if a["resname"] in AMINO:
            chains.setdefault(a["chain"] or a["segid"] or "A", []).append(a)
        else:
            others.setdefault((RENAME_MOLECULE.get(a["resname"], a["resname"]), a["chain"], a["segid"], a["resid"]), []).append(a)
    unknown = sorted({k[0] for k in others} - set(ligands) - {p.stem for p in (forcefield / "toppar").glob("*.itp")})
    if unknown:
        raise SystemExit(f"no topology source for residue(s) {unknown}")
    if not chains:
        raise SystemExit(f"{source.name} has no protein residues")
    dropped = [a for a in atoms if a["altloc"] not in " A"]
    return {"source": source, "out": out, "name": name, "sha": sha256(source), "cryst1": cryst1, "chains": chains,
            "others": others, "natoms": len(atoms) - len(dropped),
            "alternates": sorted({f"{a['chain']}:{a['resname']}{a['resid']}" for a in dropped}),
            "disulfides": find_disulfides(chains), "forcefield": forcefield, "ligands": set(ligands)}


def find_disulfides(chains: dict) -> dict[str, list[tuple[int, int]]]:
    """Disulfides are read from the structure itself: any two cysteine SG atoms within DISULFIDE_MAX_A are bonded."""
    sg = [(label, a) for label, chain in chains.items() for a in chain if a["resname"] == "CYS" and a["atom"] == "SG"]
    found: dict[str, list[tuple[int, int]]] = {label: [] for label in chains}
    partner: dict[tuple, tuple] = {}
    for (la, a), (lb, b) in itertools.combinations(sg, 2):
        if dist(a, b) >= DISULFIDE_MAX_A:
            continue
        ka, kb = (la, a["resid"]), (lb, b["resid"])
        if la != lb:
            raise SystemExit(f"inter-chain disulfide {la}:CYS{a['resid']}-{lb}:CYS{b['resid']} ({dist(a, b):.2f} A) "
                             "is not supported: each chain gets its own topology")
        if ka in partner or kb in partner:
            raise SystemExit(f"{la}:CYS{a['resid'] if ka in partner else b['resid']} is within {DISULFIDE_MAX_A} A "
                             "of two cysteine SG atoms; the disulfide pattern is ambiguous")
        partner[ka], partner[kb] = kb, ka
        found[la].append(tuple(sorted((a["resid"], b["resid"]))))
    return found


def repair_structure(system: dict) -> list[str]:
    """Fix stretched disulfides and turn the acyl-grafted Lys26 into KSM with an amide NZ."""
    repairs = []
    for label, pairs in system["disulfides"].items():
        by = {(a["resid"], a["atom"]): a for a in system["chains"][label] if a["resname"] == "CYS"}
        for left, right in pairs:
            if not all((r, n) in by for r in (left, right) for n in ("CB", "SG")):
                raise SystemExit(f"{label} CYS{left}-CYS{right}: disulfide cysteine lacks CB or SG")
            sg_sg = dist(by[(left, "SG")], by[(right, "SG")])
            if not DISULFIDE_OK_A[0] <= sg_sg <= DISULFIDE_OK_A[1]:
                set_disulfide_geometry(by, left, right)
                repairs.append(f"{label} CYS{left}-CYS{right} SG-SG {sg_sg:.2f} A -> 2.03 A")
    for label, chain in system["chains"].items():
        for rid in sorted({a["resid"] for a in chain if a["resname"] == "LYS" and a["atom"] == "C11"}):
            lys = {a["atom"]: a for a in chain if a["resid"] == rid and a["resname"] == "LYS"}
            hz = [n for n in ("HZ1", "HZ2", "HZ3") if n in lys]
            keep = max(hz, key=lambda n: dist(lys[n], lys["C11"]))
            drop = {id(lys[n]) for n in hz if n != keep}
            for a in lys.values():
                a["resname"] = "KSM"
            lys[keep]["atom"] = "HZ1"
            system["chains"][label] = chain = [a for a in chain if id(a) not in drop]
            repairs.append(f"{label} LYS{rid}: acyl graft -> KSM, NZ amide keeps {keep} as HZ1, dropped {len(drop)} H")
    return repairs


def set_disulfide_geometry(by: dict, left: int, right: int) -> None:
    """Rebuild both SG atoms of a disulfide at 2.03 A, keeping them on the side they already point to."""
    xyz = lambda a: np.array([a["x"], a["y"], a["z"]])
    cb1, sg1, cb2, sg2 = by[(left, "CB")], by[(left, "SG")], by[(right, "CB")], by[(right, "SG")]
    axis = xyz(cb2) - xyz(cb1)
    cb_cb = float(np.linalg.norm(axis))
    axis /= cb_cb
    half_gap = (cb_cb - 2.03) / 2.0
    if not 0 < half_gap < 1.81:
        raise SystemExit(f"cannot repair disulfide C{left}-C{right}: CB-CB {cb_cb:.2f} A")
    side = (xyz(sg1) - xyz(cb1)) + (xyz(sg2) - xyz(cb2))
    side -= axis * side.dot(axis)
    side = side if np.linalg.norm(side) >= 1e-6 else np.cross(axis, [0.0, 0.0, 1.0])
    side /= np.linalg.norm(side)
    lift = math.sqrt(1.81 ** 2 - half_gap ** 2)
    for atom, origin, sign in ((sg1, cb1, 1.0), (sg2, cb2, -1.0)):
        atom["x"], atom["y"], atom["z"] = map(float, xyz(origin) + sign * axis * half_gap + side * lift)


def lipid_clashes(system: dict) -> dict[tuple, tuple[float, bool]]:
    """Find lipids near protein/ligand heavy atoms, including across the periodic XY boundary."""
    cryst = system["cryst1"].ljust(80)
    box = np.array([float(cryst[6:15]), float(cryst[15:24]), 0.0])
    heavy = lambda atoms: xyz_nm([a for a in atoms if element(a["atom"]) != "H"])
    solute = heavy([a for c in system["chains"].values() for a in c]
                   + [a for k, atoms in system["others"].items() if k[0] in system["ligands"] for a in atoms])
    in_cell = cKDTree(solute)
    images = cKDTree(np.vstack([solute + box * (i, j, 0) for i in (-1, 0, 1) for j in (-1, 0, 1)]))
    found = {}
    for key, atoms in system["others"].items():
        if key[0] not in system["ligands"]:
            closest = float(images.query(heavy(atoms))[0].min())
            if closest < LIPID_SCAN_A:
                found[key] = (closest, float(in_cell.query(heavy(atoms))[0].min()) >= LIPID_SCAN_A)
    return found


def resolve_lipid_clashes(system: dict, clashes: dict[tuple, tuple[float, bool]]) -> list[str]:
    """Drop the one lipid that sits on top of the solute; fail if more than one does."""
    fatal = sorted((d, key) for key, (d, _image) in clashes.items() if d < LIPID_DELETE_A)
    if len(fatal) > 1:
        raise SystemExit(f"{len(fatal)} lipids within {LIPID_DELETE_A} A of the solute; only one may be removed: "
                         + ", ".join(f"{k[0]}{k[3]} {d:.2f} A" for d, k in fatal))
    notes = [f"lipid scan: {len(clashes)} lipids within {LIPID_SCAN_A} A of the solute (XY minimum image)"]
    for d, key in fatal:
        del system["others"][key]
        notes.append(f"removed {key[0]}{key[3]}: heavy atom {d:.2f} A from the solute "
                     + ("via a periodic XY image" if clashes[key][1] else "in cell"))
    return notes


def setup_forcefield(work: Path, resnames: set[str], forcefield: Path) -> None:
    """Make a private copy of charmm36.ff for pdb2gmx with CYSG hydrogens and extra residue types."""
    ff = work / "charmm36.ff"
    ff.mkdir()
    for f in (forcefield / "charmm36.ff").iterdir():
        if f.name != "merged.hdb":
            (ff / f.name).symlink_to(f)
    hdb = (forcefield / "charmm36.ff" / "merged.hdb").read_text().rstrip("\n") + "\n"
    (ff / "merged.hdb").write_text(hdb + (CYSG_HDB if "CYSG" in resnames and "\nCYSG" not in "\n" + hdb else ""))
    types = (forcefield / "charmm36.ff" / "residuetypes.dat").read_text().rstrip("\n") + "\n"
    known = {line.split()[0] for line in types.splitlines() if line.split()}
    (work / "residuetypes.dat").write_text(types + "".join(f"{n}    Protein\n" for n in sorted(resnames - known)))


def add_lipidated_residue(work: Path, name: str, body: list[dict], forcefield: Path) -> str:
    """Turn a CGenFF lipidated lysine into an RTP residue (CHARMM backbone, CGenFF side chain) plus its extra parameters."""
    spec, package = LIPIDATED[name], forcefield / LIPIDATED[name]["package"]
    donor = [a for a in body if a["resname"] == name]
    names = [a["atom"] for a in donor]
    if len(names) != spec["count"] - len(spec["unmapped"]):
        raise SystemExit(f"{name}: {len(names)} atoms, expected {spec['count'] - len(spec['unmapped'])}")
    atoms, bonds, section = {}, [], ""
    for line in (package / spec["top"]).read_text().splitlines():
        if line.lstrip().startswith("["):
            section = line.strip(" []\t")
            continue
        f = line.split(";", 1)[0].split()
        if section == "atoms" and len(f) >= 8:
            atoms[int(f[0])] = (f[4], f[1], float(f[6]))
        elif section == "bonds" and len(f) >= 2:
            bonds.append((int(f[0]), int(f[1])))
    if len(atoms) != spec["count"] or abs(sum(row[2] for row in atoms.values()) + 2) > 1e-5:
        raise SystemExit(f"{name}: CGenFF free-residue atom count or charge changed")
    mapping = dict(spec["mapping"])
    if spec["tail"] is None:
        mapping.update({i: atoms[i][0] for i in atoms if i not in mapping and i not in spec["unmapped"]})
    else:
        mapping.update({i: names[i + spec["offset"]] for i in spec["tail"]})
    if set(mapping.values()) != set(names):
        raise SystemExit(f"{name}: CGenFF atom order does not map onto the construct")
    coord = {a["atom"]: a for a in donor}
    kept_bonds = [(mapping[a], mapping[b]) for a, b in bonds if a in mapping and b in mapping]
    lengths = [dist(coord[a], coord[b]) for a, b in kept_bonds]
    if min(lengths) < 0.85 or max(lengths) > 1.85 or {"NZ", spec["junction"]} not in [set(b) for b in kept_bonds]:
        raise SystemExit(f"{name}: CGenFF bond graph disagrees with construct coordinates")
    side_charge = sum(atoms[i][2] for i, n in mapping.items() if n not in LYS_BACKBONE and n != "CA")
    ca_charge = -2 - side_charge - sum(q for _, q in LYS_BACKBONE.values())
    if abs(ca_charge - spec["ca_charge"]) > 1e-5:
        raise SystemExit(f"{name}: backbone charge closure changed (CA={ca_charge:.5f})")
    backbone = dict(LYS_BACKBONE, CA=("CT1", ca_charge))
    block = [f"\n[ {name} ]", "  [ atoms ]"]
    block += [f"{mapping[i]:>6s} {backbone.get(mapping[i], atoms[i][1:3])[0]:>10s} "
              f"{backbone.get(mapping[i], atoms[i][1:3])[1]:10.6f} {i:4d}" for i in sorted(mapping)]
    block += ["  [ bonds ]"] + [f"{a:>6s} {b:>6s}" for a, b in kept_bonds]
    block += ["     C     +N", "  [ impropers ]", "     N     -C     CA     HN", "     C     CA     +N     O",
              "  [ cmap ]", "    -C      N     CA      C    +N", ""]
    ff = work / "charmm36.ff"
    (ff / "merged.rtp").unlink()
    (ff / "merged.rtp").write_text((forcefield / "charmm36.ff" / "merged.rtp").read_text() + "\n".join(block))
    angle = ("CG2O1", "CG321", "OG301")
    dihedrals = {("NG2S1", "CG2O1", "CG321", "OG301"), ("OG2D1", "CG2O1", "CG321", "OG301"),
                 ("CG2O1", "CG321", "OG301", "CG321"), ("NG2S1", "CG321", "CG321", "OG301")}
    extra_angles, extra_dihedrals = [], []
    for line in (package / "charmm36.ff" / "ffbonded.itp").read_text().splitlines(True):
        f = line.split(";", 1)[0].split()
        if len(f) > 3 and f[3] == "5" and tuple(f[:3]) in (angle, angle[::-1]):
            extra_angles.append(line)
        if len(f) > 4 and f[4] == "9" and (tuple(f[:4]) in dihedrals or tuple(f[3::-1]) in dihedrals):
            extra_dihedrals.append(line)
    if len(extra_angles) != 1 or len(extra_dihedrals) != 8:
        raise SystemExit(f"{name}: supplemental CGenFF bonded parameters changed")
    (ff / "ffbonded.itp").unlink()
    (ff / "ffbonded.itp").write_text((forcefield / "charmm36.ff" / "ffbonded.itp").read_text() + "\n[ angletypes ]\n"
                                     + "".join(extra_angles) + "\n[ dihedraltypes ]\n" + "".join(extra_dihedrals))
    return f"{name} polymer RTP: {len(names)} atoms, {len(kept_bonds)} bonds, NZ-{spec['junction']}, charge -2"


def charmm_patches(forcefield: Path) -> dict:
    """Pull the ACE and CT3 patch atoms (type, charge, mass) and bonds from the CHARMM protein RTF."""
    masses, patches, current = {}, {}, None
    for raw in (forcefield / "protein_parameters" / "top_all36_prot.rtf").read_text().splitlines():
        p = raw.split("!", 1)[0].split()
        if not p:
            continue
        key = p[0].upper()
        if key == "MASS" and len(p) >= 4:
            masses[p[2].upper()] = float(p[3])
        elif key in {"RESI", "PRES"}:
            current = p[1].upper() if p[1].upper() in {"ACE", "CT3"} else None
            if current:
                patches[current] = {"atoms": {}, "bonds": []}
        elif current and key == "ATOM":
            patches[current]["atoms"][p[1].upper()] = (p[2].upper(), float(p[3]), masses[p[2].upper()])
        elif current and key in {"BOND", "DOUBLE"}:
            vals = [x.upper().strip("+-") for x in p[1:]]
            patches[current]["bonds"] += list(zip(vals[::2], vals[1::2]))
    return patches


def chain_termini(residues: list) -> tuple[str, str]:
    """Work out each chain's termini from the atoms on its first and last residues."""
    first, last = {a["atom"] for a in residues[0][2]}, {a["atom"] for a in residues[-1][2]}
    nterm = "ACE" if set(N_CAP) <= first else "NTER"
    if first & set(N_CAP) and nterm != "ACE":
        raise SystemExit(f"incomplete ACE cap on {residues[0][1]}{residues[0][0]}")
    if set(C_CAP) <= last:
        return nterm, "CT3"
    if {"NT", "HT1", "HT2"} <= last and "CAT" not in last:
        return nterm, "CT2"
    if last & {"OXT", "OT1", "OT2"}:
        return nterm, "CTER"
    raise SystemExit(f"unrecognized C terminus on {residues[-1][1]}{residues[-1][0]}: {sorted(last)}")


def split_caps(residues: list) -> tuple[list[dict], dict, dict]:
    """Separate cap atoms from the chain body and normalize residue and atom names for pdb2gmx."""
    body, ncap, ccap = [], {}, {}
    for _rid, rn, atoms in residues:
        rn = RENAME_RESIDUE.get(rn, rn)
        for a in atoms:
            a = dict(a, resname=rn)
            if a["atom"] in N_CAP:
                ncap[a["atom"]] = a
            elif a["atom"] in C_CAP:
                ccap[a["atom"]] = a
            else:
                a["atom"] = RENAME_ATOM.get(rn, {}).get(a["atom"], "HN" if a["atom"] == "H" else a["atom"])
                body.append(a)
    return body, ncap, ccap


def run_pdb2gmx(cmd: list[str], cwd: Path, env: dict, selections: list[str], timeout: int = 900) -> tuple[int, str]:
    """Drive pdb2gmx through a pty, picking each terminus menu entry by its label."""
    menu = re.compile(r"^\s*(\d+)\s*:\s*(.+?)\s*$", re.M)
    master, slave = pty.openpty()
    proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdin=slave, stdout=slave, stderr=slave, close_fds=True)
    os.close(slave)
    out, sent, menu_start = "", 0, 0
    deadline = time.time() + timeout
    try:
        while True:
            if time.time() > deadline:
                proc.kill()
                raise SystemExit("pdb2gmx timed out waiting for a terminus menu")
            if not select.select([master], [], [], 0.2)[0]:
                if proc.poll() is not None:
                    break
                continue
            try:
                chunk = os.read(master, 8192).decode("utf-8", "replace")
            except OSError:
                break
            if not chunk:
                break
            out += chunk
            if sent < len(selections) and "terminus type" in out[menu_start:].lower():
                choice = [n for n, label in menu.findall(out[menu_start:][-4000:]) if label.lower() == selections[sent].lower()]
                if choice:
                    os.write(master, (choice[-1] + "\n").encode())
                    sent, menu_start = sent + 1, len(out)
    finally:
        os.close(master)
        proc.wait()
    return proc.returncode, out


def restore_source_coordinates(processed: Path, body: list[dict]) -> None:
    """Put source coordinates back on residues pdb2gmx moved, and on CYSG hydrogens."""
    atoms = read_pdb(processed)[0]
    source = {(a["resid"], a["resname"], a["atom"]): a for a in body}
    shifted = {(a["resid"], a["resname"]) for a in atoms if element(a["atom"]) != "H"
               and (a["resid"], a["resname"], a["atom"]) in source
               and dist(a, source[(a["resid"], a["resname"], a["atom"])]) > 0.25}
    for a in atoms:
        src = source.get((a["resid"], a["resname"], a["atom"]))
        if src and ((a["resid"], a["resname"]) in shifted or a["resname"] == "CYSG"):
            if a["resname"] == "CYSG" and element(a["atom"]) != "H" and dist(a, src) > 0.02:
                raise SystemExit(f"CYSG heavy atom {a['atom']} moved {dist(a, src):.3f} A in pdb2gmx")
            a["x"], a["y"], a["z"] = src["x"], src["y"], src["z"]
    write_pdb(atoms, processed)


def moleculetype_block(top: Path, mol: str) -> str:
    """Cut the [ moleculetype ] block out of a .top and rename the molecule."""
    body = re.search(r"(?s)\[ moleculetype \].*?(?=\n; Include Position restraint|\n#ifdef POSRES|\n\[ system \])",
                     top.read_text())
    if not body:
        raise SystemExit(f"{mol}: no [ moleculetype ] in {top.name}")
    return re.sub(r"(\[ moleculetype \]\n(?:;.*\n)*)\S+", rf"\g<1>{mol}", body.group(0), count=1) + "\n"


def bonded_terms(natoms: int, bonds: set) -> tuple[list, list, list]:
    """Derive 1-4 pairs, angles and proper dihedrals from a bond list."""
    adj = defaultdict(set)
    for a, b in bonds:
        adj[a].add(b)
        adj[b].add(a)
    angles = sorted({(a, c, b) for c in range(1, natoms + 1) for a in adj[c] for b in adj[c] if a < b})
    near, pairs = bonds | {(a, b) for a, _, b in angles}, set()
    for start in range(1, natoms + 1):
        queue, depth = deque([start]), {start: 0}
        while queue:
            x = queue.popleft()
            if depth[x] == 3:
                if start < x and (start, x) not in near:
                    pairs.add((start, x))
                continue
            for y in adj[x] - depth.keys():
                depth[y] = depth[x] + 1
                queue.append(y)
    dihedrals = {min((a, b, c, d), (d, c, b, a)) for b, c in bonds for a in adj[b] - {c} for d in adj[c] - {b, a}}
    return sorted(pairs), angles, sorted(dihedrals)


def patch_caps(mol: str, processed: Path, body_itp: Path, ncap: dict, ccap: dict, patches: dict, out_itp: Path) -> None:
    """Write a chain .itp with CHARMM ACE/CT3 patch atoms folded into the terminal residues."""
    pdb_atoms, itp = read_pdb(processed)[0], read_itp(body_itp)
    if len(pdb_atoms) != len(itp["atoms"]):
        raise SystemExit(f"{mol}: pdb2gmx coordinate/topology atom count mismatch")
    first, last = pdb_atoms[0], pdb_atoms[-1]
    coords, rows, old_to_new = [], [], {}

    def add(atom: dict, typ: str, charge: float, mass: float) -> int:
        """Append one atom and its topology row; returns its new index."""
        coords.append(atom)
        rows.append({"type": typ, "resnr": atom["resid"], "residue": atom["resname"], "atom": atom["atom"],
                     "charge": charge, "mass": mass})
        return len(coords)

    for name in N_CAP if ncap else []:
        add(dict(ncap[name], resid=first["resid"], resname=first["resname"]), *patches["ACE"]["atoms"][name])
    has_hn = any(a["resid"] == first["resid"] and a["atom"] == "HN" for a in pdb_atoms)
    synthetic_hn = None
    for old, (atom, row) in enumerate(zip(pdb_atoms, itp["atoms"]), 1):
        old_to_new[old] = add(atom, row["type"], row["charge"], row["mass"])
        if ncap and not has_hn and atom["resid"] == first["resid"] and atom["atom"] == "N":
            cy = ncap["CY"]
            direction = np.array([atom["x"] - cy["x"], atom["y"] - cy["y"], atom["z"] - cy["z"]])
            xyz = np.array([atom["x"], atom["y"], atom["z"]]) + 1.01 * direction / np.linalg.norm(direction)
            hn = dict(atom, atom="HN", elem="H", x=float(xyz[0]), y=float(xyz[1]), z=float(xyz[2]))
            synthetic_hn = add(hn, *patches["CT3"]["atoms"]["HNT"])
    for name in C_CAP if ccap else []:
        add(dict(ccap[name], resid=last["resid"], resname=last["resname"]), *patches["CT3"]["atoms"][name])

    idx = {(a["resid"], a["atom"]): i for i, a in enumerate(coords, 1)}
    bonds = {tuple(sorted((old_to_new[a], old_to_new[b]))) for a, b in itp["bonds"]}
    if synthetic_hn:
        bonds.add((idx[(first["resid"], "N")], synthetic_hn))
    side = lambda n: first["resid"] if n in N_CAP or n in {"N", "CA", "HN", "CD"} else last["resid"]
    for a, b in patches["ACE"]["bonds"] + patches["CT3"]["bonds"]:
        if (side(a), a) in idx and (side(b), b) in idx:
            bonds.add(tuple(sorted((idx[(side(a), a)], idx[(side(b), b)]))))
    impropers = [tuple(old_to_new[i] for i in row) for row in itp["impropers"]]
    for names, rid in ((("CY", "CAY", "N", "OY"), first["resid"]), (("N", "CY", "CA", "HN"), first["resid"]),
                       (("NT", "C", "CAT", "HNT"), last["resid"]), (("C", "CA", "NT", "O"), last["resid"])):
        if all((rid, n) in idx for n in names):
            impropers.append(tuple(idx[(rid, n)] for n in names))
    order = list(OrderedDict.fromkeys(a["resid"] for a in coords))
    cmap = []
    for i, rid in enumerate(order):
        prev = (first["resid"], "CY") if i == 0 else (order[i - 1], "C")
        nxt = (last["resid"], "NT") if i == len(order) - 1 else (order[i + 1], "N")
        keys = [prev, (rid, "N"), (rid, "CA"), (rid, "C"), nxt]
        if all(k in idx for k in keys):
            cmap.append(tuple(idx[k] for k in keys))

    pairs, angles, dihedrals = bonded_terms(len(rows), bonds)
    with out_itp.open("w") as fh:
        fh.write(f"[ moleculetype ]\n; Name            nrexcl\n{mol:<20s} 3\n\n[ atoms ]\n")
        for i, r in enumerate(rows, 1):
            fh.write(f"{i:6d} {r['type']:>10s} {r['resnr']:6d} {r['residue']:>6s} {r['atom']:>6s} {i:6d} "
                     f"{r['charge']:10.5f} {r['mass']:10.5f}\n")
        for name, entries, funct in [("bonds", sorted(bonds), 1), ("pairs", pairs, 1), ("angles", angles, 5),
                                     ("dihedrals", dihedrals, 9), ("dihedrals", impropers, 2), ("cmap", cmap, 1)]:
            fh.write(f"\n[ {name} ]\n" + "".join(" ".join(f"{i:5d}" for i in e) + f" {funct:5d}\n" for e in entries))
    write_pdb(coords, processed)


def write_posre(itp: Path, mol: str) -> None:
    """Write backbone/side-chain position restraints for a chain and include them from its .itp."""
    backbone = {"N", "CA", "C", "O", "OT1", "OT2", "CY", "OY", "NT"}
    rows = [f"{a['nr']:6d}     1    {m:<12s} {m:<12s} {m}\n" for a in read_itp(itp)["atoms"] if element(a["atom"]) != "H"
            for m in ["POSRES_FC_BB" if a["atom"] in backbone else "POSRES_FC_SC"]]
    (itp.parent / f"posre_{mol}.itp").write_text("[ position_restraints ]\n" + "".join(rows))
    with itp.open("a") as fh:
        fh.write(f'\n#ifdef POSRES\n#include "posre_{mol}.itp"\n#endif\n')


def check_disulfides(itp: Path, mol: str, expected: list[tuple[int, int]]) -> None:
    """Make sure pdb2gmx bonded exactly the disulfides found in the input structure."""
    data = read_itp(itp)
    sg = {a["nr"]: a["resnr"] for a in data["atoms"] if a["atom"] == "SG" and a["residue"] in {"CYS", "CYS2"}}
    actual = {tuple(sorted((sg[a], sg[b]))) for a, b in data["bonds"] if a in sg and b in sg}
    if actual != set(expected):
        raise SystemExit(f"{mol} disulfides {sorted(actual)} != expected {sorted(expected)}")


def build_protein(system: dict, label: str, work: Path, toppar: Path, gmx: str, patches: dict) -> list[str]:
    """Build one protein chain's .itp, posre and coordinates with pdb2gmx (plus cap patching)."""
    out, mol = system["out"], f"PRO{label}"
    residues = residues_in_order(system["chains"][label])
    nterm, cterm = chain_termini(residues)
    body, ncap, ccap = split_caps(residues)
    lipidated = {a["resname"] for a in body} & set(LIPIDATED)
    notes = [f"{mol}: {residues[0][1]}{residues[0][0]}-{residues[-1][1]}{residues[-1][0]} {nterm}/{cterm}"]
    notes += [add_lipidated_residue(work, name, body, system["forcefield"]) for name in sorted(lipidated)]
    pdb_in, pdb_out, top = work / f"{mol}_in.pdb", work / f"{mol}.pdb", work / f"{mol}.top"
    write_pdb(body, pdb_in)
    cmd = gmx.split() + ["pdb2gmx", "-f", str(pdb_in), "-o", str(pdb_out), "-p", str(top), "-i", str(work / f"{mol}_pr.itp"),
                         "-ff", "charmm36", "-water", "tip3p", "-ter"] + ([] if lipidated else ["-ignh", "-missing"])
    rc, transcript = run_pdb2gmx(cmd, work, dict(os.environ, GMXLIB=str(work), GMX_MAXBACKUP="-1"),
                                 [TERMINUS_MENU[nterm], TERMINUS_MENU[cterm]])
    log(out, f"$ {' '.join(cmd)}\n{transcript}", "DEBUG")
    if rc != 0:
        errors = re.findall(r"(?s)Fatal error:\s*\n(.*?)\n\s*\n", transcript)
        raise SystemExit(f"pdb2gmx failed for {mol}: {' '.join(errors[-1].split()) if errors else transcript[-600:]}")
    restore_source_coordinates(pdb_out, body)
    itp = toppar / f"{mol}.itp"
    if nterm == "ACE" or cterm == "CT3":
        (work / f"{mol}_body.itp").write_text(moleculetype_block(top, mol))
        patch_caps(mol, pdb_out, work / f"{mol}_body.itp", ncap if nterm == "ACE" else {},
                   ccap if cterm == "CT3" else {}, patches, itp)
    else:
        itp.write_text(moleculetype_block(top, mol))
    write_posre(itp, mol)
    check_disulfides(itp, mol, system["disulfides"][label])
    if system["disulfides"][label]:
        notes.append(f"{mol}: disulfides " + ", ".join(f"CYS{a}-CYS{b}" for a, b in sorted(system["disulfides"][label])))
    topology_atoms = read_itp(itp)["atoms"]
    heavy_source = sum(element(a["atom"]) != "H" for a in system["chains"][label])
    heavy_topology = sum(element(a["atom"]) != "H" for a in topology_atoms)
    if heavy_source != heavy_topology or len(read_pdb(pdb_out)[0]) != len(topology_atoms):
        raise SystemExit(f"{mol}: heavy atoms source={heavy_source} topology={heavy_topology}")
    return notes


def place_h(parent: dict, name: str, n: int) -> dict:
    """Place a missing hydrogen on a fixed direction from its parent; EM sorts out the geometry."""
    dx, dy, dz = [(1.0, 0.0, 0.0), (-0.333, 0.943, 0.0), (-0.333, -0.471, 0.816), (-0.333, -0.471, -0.816)][n % 4]
    bond = 0.98 if source_element(parent) in {"N", "O"} else 1.09
    return dict(parent, atom=name, elem="H", x=parent["x"] + bond * dx, y=parent["y"] + bond * dy, z=parent["z"] + bond * dz)


def h_parents(itp: dict) -> dict[int, int]:
    """Map each hydrogen in an .itp to its bonded heavy atom."""
    names = {a["nr"]: a["atom"] for a in itp["atoms"]}
    return {h: p for i, j in itp["bonds"] for h, p in ((i, j), (j, i)) if element(names[h]) == "H" and element(names[p]) != "H"}


def check_h_bonds(out: list[dict], itp: dict, mol: str) -> None:
    """Fail if any X-H bond is outside a sane length range."""
    for h, p in h_parents(itp).items():
        if not H_BOND_RANGE[0] <= dist(out[h - 1], out[p - 1]) <= H_BOND_RANGE[1]:
            raise SystemExit(f"{mol}: {out[h - 1]['atom']}-{out[p - 1]['atom']} bond {dist(out[h - 1], out[p - 1]):.2f} A")


def match_by_name(res: list[dict], itp: dict, mol: str) -> tuple[list[dict], int]:
    """Order a residue's atoms like its .itp by atom name, adding hydrogens that are missing."""
    by = defaultdict(deque)
    for a in res:
        by[a["atom"]].append(a)
    parent, out, synthetic = h_parents(itp), [], 0
    for spec in itp["atoms"]:
        if by[spec["atom"]]:
            out.append(by[spec["atom"]].popleft())
        elif spec["nr"] in parent and parent[spec["nr"]] < spec["nr"]:
            out.append(place_h(out[parent[spec["nr"]] - 1], spec["atom"], synthetic))
            synthetic += 1
        else:
            raise SystemExit(f"{mol} {res[0]['chain']}{res[0]['resid']}: missing atom {spec['atom']}")
    leftover = [a["atom"] for q in by.values() for a in q]
    if leftover:
        raise SystemExit(f"{mol} {res[0]['chain']}{res[0]['resid']}: atoms not in topology {leftover[:10]}")
    check_h_bonds(out, itp, mol)
    return out, synthetic


def match_by_graph(res: list[dict], itp: dict, mol: str) -> tuple[list[dict], int]:
    """Order a ligand's atoms like its .itp by matching bond graphs, since ligand atom names differ between sources."""
    radii = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07, "F": 0.57, "CL": 1.02}
    source = nx.Graph()
    heavy = [(i, a) for i, a in enumerate(res, 1) if source_element(a) != "H"]
    source.add_nodes_from((i, {"elem": source_element(a)}) for i, a in heavy)
    for k, (i, a) in enumerate(heavy):
        for j, b in heavy[k + 1:]:
            if dist(a, b) <= 1.25 * (radii.get(source_element(a), 0.7) + radii.get(source_element(b), 0.7)):
                source.add_edge(i, j)
    names = {a["nr"]: a["atom"] for a in itp["atoms"]}
    target = nx.Graph()
    target.add_nodes_from((nr, {"elem": element(n)}) for nr, n in names.items() if element(n) != "H")
    target.add_edges_from((i, j) for i, j in itp["bonds"] if i in target and j in target)
    matcher = iso.GraphMatcher(target, source, node_match=iso.categorical_node_match("elem", None))
    candidates = list(itertools.islice(matcher.isomorphisms_iter(), 2000))
    if not candidates:
        raise SystemExit(f"{mol} {res[0]['chain']}{res[0]['resid']}: heavy-atom graph is not isomorphic to its ITP")
    mapping = max(candidates, key=lambda m: sum(names[t] == res[s - 1]["atom"] for t, s in m.items()))
    parent, free_h, out, synthetic = h_parents(itp), [a for a in res if source_element(a) == "H"], [], 0
    for spec in itp["atoms"]:
        if spec["nr"] in mapping:
            out.append(dict(res[mapping[spec["nr"]] - 1], atom=spec["atom"]))
            continue
        anchor = res[mapping[parent[spec["nr"]]] - 1]
        near = sorted((h for h in free_h if H_BOND_RANGE[0] <= dist(h, anchor) <= H_BOND_RANGE[1]),
                      key=lambda h: (h["atom"] != spec["atom"], dist(h, anchor)))
        if near:
            free_h.remove(near[0])
            out.append(dict(near[0], atom=spec["atom"]))
        else:
            out.append(place_h(anchor, spec["atom"], synthetic))
            synthetic += 1
    if free_h:
        raise SystemExit(f"{mol}: {len(free_h)} source hydrogens have no topology counterpart")
    check_h_bonds(out, itp, mol)
    return out, synthetic


def match_glpa(res: list[dict], itp: dict) -> list[dict]:
    """Check a GLPA (GM3) residue is complete, in GLPA.itp order and unbroken, and tag its sub-residue names."""
    label = f"GLPA {res[0]['chain']}{res[0]['resid']}"
    order = [spec["atom"] for spec in itp["atoms"]]
    if {a["resname"] for a in res} == {"GM3"}:  # backmapped GM3 still carries GM3.xml names and order
        by = {a["atom"]: a for a in res}
        if len(by) != len(res) or set(by) != {src for src, _ in GM3_XML_TO_GLPA}:
            raise SystemExit(f"{label}: GM3 atom names do not match the GM3.xml -> GLPA plan")
        if [dst for _, dst in GM3_XML_TO_GLPA] != order:
            raise SystemExit("GM3_XML_TO_GLPA does not follow GLPA.itp atom order")
        res = [dict(by[src], atom=dst) for src, dst in GM3_XML_TO_GLPA]
    if [a["atom"] for a in res] != order:
        raise SystemExit(f"{label}: atoms are not in GLPA.itp order")
    long = [(res[i - 1]["atom"], res[j - 1]["atom"], dist(res[i - 1], res[j - 1])) for i, j in itp["bonds"]
            if dist(res[i - 1], res[j - 1]) > GLPA_MAX_BOND_A]
    if long:
        raise SystemExit(f"{label}: {len(long)} bonds over {GLPA_MAX_BOND_A} A, e.g. {long[0][0]}-{long[0][1]} {long[0][2]:.2f} A")
    check_h_bonds(res, itp, "GLPA")
    return [dict(a, topology_resname=spec["residue"]) for a, spec in zip(res, itp["atoms"])]


def molecule_itp(name: str, toppar: Path, forcefield: Path, ligands: set) -> Path:
    """Copy a molecule's .itp into toppar/ and give it heavy-atom position restraints if it has none."""
    path = toppar / f"{name}.itp"
    if name == "V6G":
        path.write_text(moleculetype_block(forcefield / "V6G_gromacs" / "V6G_gmx.top", name))
    else:
        shutil.copy(forcefield / "toppar" / f"{name}.itp", path)
    data = read_itp(path)
    if not data["posres"] and name not in SOLVENT:
        macro = "POSRES_FC_LIG" if name in ligands else "POSRES_FC_LIPID"
        rows = [f"{a['nr']:6d}  1  {macro}  {macro}  {macro}\n" for a in data["atoms"] if element(a["atom"]) != "H"]
        (toppar / f"posre_{name}.itp").write_text("[ position_restraints ]\n" + "".join(rows))
        with path.open("a") as fh:
            fh.write(f'\n#ifdef POSRES\n#include "posre_{name}.itp"\n#endif\n')
    return path


def forcefield_includes(ff: Path) -> list[str]:
    """List forcefield.itp and every file it includes by default, so only those get copied."""
    needed, pending, defines = [], ["forcefield.itp"], set()
    while pending:
        name = pending.pop(0)
        if name in needed:
            continue
        needed.append(name)
        active = []
        for line in (ff / name).read_text(errors="replace").splitlines():
            s = line.split(";", 1)[0].split()
            if not s:
                continue
            if s[0] in ("#ifdef", "#ifndef"):
                active.append((s[1] in defines) == (s[0] == "#ifdef"))
            elif s[0] == "#else":
                active[-1] = not active[-1]
            elif s[0] == "#endif":
                active.pop()
            elif all(active) and s[0] == "#define":
                defines.add(s[1])
            elif all(active) and s[0] == "#include":
                pending.append(Path(s[1].strip('"')).name)
    return needed


def build_topology(system: dict, work: Path, gmx: str) -> dict:
    """Write toppar/, topol.top and prot-memb.pdb, with coordinates in exact topology order."""
    out, forcefield, ligands = system["out"], system["forcefield"], system["ligands"]
    toppar = out / "toppar"
    shutil.rmtree(toppar, ignore_errors=True)
    (toppar / "charmm36.ff").mkdir(parents=True)
    protein_names = {RENAME_RESIDUE.get(a["resname"], a["resname"]) for c in system["chains"].values() for a in c}
    setup_forcefield(work, protein_names, forcefield)
    patches = charmm_patches(forcefield)
    molecules, coords, notes = [], [], []
    for label in system["chains"]:
        notes += build_protein(system, label, work, toppar, gmx, patches)
        molecules.append((f"PRO{label}", 1, read_itp(toppar / f"PRO{label}.itp")["atoms"], "Protein_LIG"))
        coords += [dict(a, chain=label, segid=f"PRO{label}", group="Protein_LIG") for a in read_pdb(work / f"PRO{label}.pdb")[0]]
    residues: OrderedDict[str, list] = OrderedDict()
    for key, atoms in system["others"].items():
        residues.setdefault(key[0], []).append(atoms)
    for name in ("TIP3", "SOD", "CLA"):
        residues.setdefault(name, [])
    solvent_itps = {}
    for name, blocks in residues.items():
        itp = read_itp(molecule_itp(name, toppar, forcefield, ligands))
        group = "Protein_LIG" if name in ligands else "SOL_ION" if name in SOLVENT else "MEMB"
        if not blocks:
            solvent_itps[name] = itp
            continue
        synthetic = 0
        for res in blocks:
            if len(itp["atoms"]) == 1:
                if len(res) != 1:
                    raise SystemExit(f"{name} {res[0]['resid']}: {len(res)} atoms for a one-atom topology")
                mapped, n = [dict(res[0], atom=itp["atoms"][0]["atom"])], 0
            elif name == "GLPA":
                mapped, n = match_glpa(res, itp), 0
            else:
                mapped, n = (match_by_graph if name in ligands else match_by_name)(res, itp, name)
            synthetic += n
            coords += [dict(a, resname=a.get("topology_resname", name), segid=name[:4], group=group) for a in mapped]
        molecules.append((itp["mol"], len(blocks), itp["atoms"], group))
        if synthetic:
            notes.append(f"{name}: placed {synthetic} H absent from the input")
    active_ff = work / "charmm36.ff" if protein_names & set(LIPIDATED) else forcefield / "charmm36.ff"
    for name in forcefield_includes(active_ff):
        shutil.copy(active_ff / name, toppar / "charmm36.ff" / name)
    v6g_params = [forcefield / "V6G_gromacs/charmm36.ff/V6G_ffbonded.itp", forcefield / "V6G_gromacs/V6G_base_ff_params.itp"]
    v6g_params = v6g_params if "V6G" in residues else []
    for f in v6g_params:
        shutil.copy(f, toppar / f.name)
    lines = [f"#ifndef {m}\n#define {m} 0\n#endif" for m in ("POSRES_FC_SC", "POSRES_FC_BB", "POSRES_FC_LIG", "POSRES_FC_LIPID")]
    lines += ['#include "toppar/charmm36.ff/forcefield.itp"'] + [f'#include "toppar/{f.name}"' for f in v6g_params]
    lines += [f'#include "toppar/{name}.itp"' for name in [f"PRO{c}" for c in system["chains"]] + list(residues)]
    lines += ["", "[ system ]", system["name"], "", "[ molecules ]"] + [f"{mol:<12s} {n}" for mol, n, _, _ in molecules]
    (out / "topol.top").write_text("\n".join(lines) + "\n")
    expected = [a["atom"] for _, n, atoms, _ in molecules for _ in range(n) for a in atoms]
    if [a["atom"] for a in coords] != expected:
        raise SystemExit(f"coordinate/topology atom order mismatch ({len(coords)} vs {len(expected)} atoms)")
    if len({(round(a["x"], 3), round(a["y"], 3), round(a["z"], 3)) for a in coords}) != len(coords):
        raise SystemExit("duplicate atom coordinates in the topology-ordered structure")
    write_pdb(coords, out / "prot-memb.pdb", system["cryst1"])
    charge = sum(n * sum(a["charge"] for a in atoms) for _, n, atoms, _ in molecules)
    notes.append(f"topol.top: {len(coords)} atoms, net charge {charge:+.3f}; " + ", ".join(f"{m} {n}" for m, n, _, _ in molecules))
    return {"molecules": molecules, "names": expected, "charge": charge, "notes": notes, "coords": coords,
            "solvent_itps": solvent_itps}
