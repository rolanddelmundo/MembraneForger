#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// PDB / GRO / ITP readers and writers.
#//=============================================================
"""Minimal structure and topology file readers and writers used by the build stages."""
import math
from collections import OrderedDict
from pathlib import Path

import numpy as np

__all__ = ['read_pdb', 'write_pdb', 'read_gro', 'write_gro', 'xyz_nm', 'wrap', 'read_itp', 'element',
           'source_element', 'dist', 'residues_in_order']

def read_pdb(path: Path) -> tuple[list[dict], str | None]:
    """Read ATOM/HETATM records and the CRYST1 line from a PDB."""
    atoms, cryst1 = [], None
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("CRYST1") and cryst1 is None:
            cryst1 = line.rstrip()
        if line[:6] in ("ATOM  ", "HETATM"):
            line = line.ljust(80)
            atoms.append({"atom": line[12:16].strip(), "altloc": line[16], "resname": line[17:21].strip(),
                          "chain": line[21].strip(), "resid": int(line[22:26]), "x": float(line[30:38]),
                          "y": float(line[38:46]), "z": float(line[46:54]), "segid": line[72:76].strip(),
                          "elem": line[76:78].strip().upper()})
    return atoms, cryst1


def write_pdb(atoms: list[dict], path: Path, cryst1: str | None = None) -> None:
    """Write atoms as a PDB, optionally with a box."""
    with path.open("w") as fh:
        fh.write(cryst1 + "\n" if cryst1 else "")
        for i, a in enumerate(atoms, 1):
            name = f" {a['atom']:<3s}" if len(a["atom"]) < 4 else a["atom"][:4]
            fh.write(f"ATOM  {i % 100000:5d} {name} {a['resname'][:4]:<4s}{(a.get('chain') or ' ')[:1]}"
                     f"{a['resid'] % 10000:4d}    {a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}  1.00  0.00      "
                     f"{(a.get('segid') or '')[:4]:<4s}{(a.get('elem') or element(a['atom']))[:2]:>2s}\n")
        fh.write("END\n")


def read_gro(path: Path) -> tuple[list[dict], list[float]]:
    """Read a .gro file into atoms (nm) and its box."""
    lines = path.read_text().splitlines()
    n = int(lines[1])
    atoms = [{"resid": int(l[0:5]), "resname": l[5:10].strip(), "atom": l[10:15].strip(),
              "x": float(l[20:28]), "y": float(l[28:36]), "z": float(l[36:44])} for l in lines[2:2 + n]]
    return atoms, [float(v) for v in lines[2 + n].split()]


def write_gro(atoms: list[dict], box: list[float], path: Path, title: str) -> None:
    """Write atoms (nm) and a box as a .gro file."""
    with path.open("w") as fh:
        fh.write(f"{title}\n{len(atoms)}\n")
        for i, a in enumerate(atoms, 1):
            fh.write(f"{a['resid'] % 100000:5d}{a['resname'][:5]:<5s}{a['atom'][:5]:>5s}{i % 100000:5d}"
                     f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}\n")
        fh.write("".join(f"{v:10.5f}" for v in box) + "\n")


def xyz_nm(atoms: list[dict]) -> np.ndarray:
    """Stack atom coordinates into an array."""
    return np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=float)


def wrap(xyz: np.ndarray, box: list[float]) -> np.ndarray:
    """Wrap coordinates into the box."""
    out = np.mod(xyz, box[:3])
    out[out >= np.asarray(box[:3])] = 0.0
    return out


def read_itp(path: Path) -> dict:
    """Parse the parts of an .itp we use: molecule name, atoms, bonds, impropers, posres flag."""
    data = {"mol": None, "atoms": [], "bonds": [], "impropers": [], "posres": False}
    state, dihedral_block = None, 0
    for raw in path.read_text(errors="replace").splitlines():
        s = raw.split(";", 1)[0].strip()
        data["posres"] |= s.startswith("#ifdef POSRES")
        if not s or s.startswith("#"):
            continue
        if s.startswith("["):
            state = s.strip("[] \t").lower()
            dihedral_block += state == "dihedrals"
            continue
        p = s.split()
        if state == "moleculetype" and data["mol"] is None:
            data["mol"] = p[0]
        elif state == "atoms" and len(p) >= 8 and p[0].isdigit():
            data["atoms"].append({"nr": int(p[0]), "type": p[1], "resnr": int(p[2]), "residue": p[3],
                                  "atom": p[4], "charge": float(p[6]), "mass": float(p[7])})
        elif state == "bonds" and len(p) >= 2:
            data["bonds"].append((int(p[0]), int(p[1])))
        elif state == "dihedrals" and len(p) >= 4 and ((len(p) >= 5 and p[4] == "2") or dihedral_block >= 2):
            data["impropers"].append(tuple(int(x) for x in p[:4]))
    return data


def element(name: str) -> str:
    """Guess an element from an atom name."""
    n = name.strip().upper().lstrip("0123456789")
    return "CL" if n.startswith("CL") else "MG" if n.startswith("MG") else n[:1]


def source_element(atom: dict) -> str:
    """Element from the PDB element column, falling back to the atom name."""
    return "H" if element(atom["atom"]) == "H" else (atom.get("elem") or element(atom["atom"]))


def dist(a: dict, b: dict) -> float:
    """Distance between two atoms in their own units."""
    return math.dist((a["x"], a["y"], a["z"]), (b["x"], b["y"], b["z"]))


def residues_in_order(atoms: list[dict]) -> list[tuple[int, str, list[dict]]]:
    """Group atoms into residues, keeping file order."""
    seen: OrderedDict[tuple, list[dict]] = OrderedDict()
    for a in atoms:
        seen.setdefault((a["resid"], a["resname"]), []).append(a)
    return [(rid, rn, vals) for (rid, rn), vals in seen.items()]
