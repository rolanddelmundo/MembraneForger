#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// EM-only protein backbone phi/psi restraints.
#//=============================================================
"""EM-only flat-bottom phi/psi restraints for every protein chain (#ifdef BBDIHRES in the chain's .itp).

Minimization relieves the overlaps left by backmapping and embedding, and on its own it moves backbone dihedrals freely: on 140 builds
without restraints the median build had 10 residues drift by more than 60 degrees (backbone_check). The restraints hold each phi and
psi near the input model's value, with zero force inside a window and DIHRES_EM_FC beyond it. Production topologies never define
BBDIHRES, so MD is unrestrained.

Modes (--backbone-restraints):
  off      no restraints (the default)
  fixed    every window is the input angle +-BB_DIHRES_DPHI
  clipped  every window is the input angle +-BB_DIHRES_DPHI cut back to the allowed Ramachandran region of the backbone check
           (backbone_check.rama_outlier), shrunk inward by --backbone-restraint-inset degrees. Without the clip, a residue that starts
           close to a region boundary can be pushed across it while staying inside its +-30 degree window, which neither the drift test
           nor the restraint catches; with it, the zero-force region of each residue lies inside the allowed region.

Clipping works per residue on the box phi x psi: the largest box around the input point, inside the +-BB_DIHRES_DPHI square, in which
every grid point is allowed, chosen by the larger smaller half-width and then by area. Angles are taken unwrapped around the input
value, so a window may run through +-180; GROMACS restraints are symmetric, so each side of the box becomes a centre (its midpoint,
wrapped to [-180, 180)) and a half-width. The window always contains the input angle.

A residue keeps the +-BB_DIHRES_DPHI window, and the record says why, when: the backbone check does not judge it (no standard parent
residue, a Gly parent, or no phi or no psi because it ends a chain or a chain break); its input is already a Ramachandran outlier
(input_outlier); its input lies inside the inset band at a boundary (flag_margin); or the best box has a half-width below
BB_DIHRES_MIN_HALF_DEG (flag_narrow). flag_* rows mean the model is too close to a boundary for a safe window and should be refined.

Known gap: in a chain with ACE/CT3 caps the cap atoms (CY, NT) are folded into the first and last residue, so those residues have
no neighbouring residue for phi (first) or psi (last) and get no restraint.

The restraint is finite, so a window bounds where the restraint starts to push, not where the angle ends: check the result with the
backbone check (--strict-segments for residues that were modelled rather than observed).
"""
from pathlib import Path

import numpy as np

from .backbone_check import PARENT_RESIDUE
from .config import BB_DIHRES_DPHI, BB_DIHRES_MIN_HALF_DEG
from .structio import read_itp, read_pdb

__all__ = ['RESTRAINT_MODES', 'RECORD_NAME', 'rama_allowed', 'backbone_angles', 'restraint_windows', 'write_backbone_dihres',
           'write_restraint_record', 'summarize_record']

RESTRAINT_MODES = ("off", "fixed", "clipped")
RECORD_NAME = "backbone_restraints.tsv"
RECORD_COLUMNS = ("molecule", "index", "resid", "resname", "angle", "input", "centre", "half_width", "low", "high", "status")
PEPTIDE_BOND_A = 2.0    # C(i)-N(i+1) shorter than this (A) is a peptide bond, as in backbone_check (0.20 nm)
STEP_DEG = 0.5          # grid of the box search, in degrees from the input angle


def rama_allowed(parent: str, phi, psi, inset: float = 0.0) -> np.ndarray:
    """Vectorized `not backbone_check.rama_outlier(parent, phi, psi, margin=-inset)` for arrays of angles in [-180, 180).

    inset shrinks every allowed region by that many degrees (rama_outlier's margin widens them). The test suite checks this function
    against rama_outlier point by point, so the two cannot drift apart."""
    phi, psi = np.broadcast_arrays(np.asarray(phi, float), np.asarray(psi, float))
    m = -inset
    if parent == "PRO":
        return (-110 - m <= phi) & (phi <= -35 + m)
    left = phi <= 0
    left_out = (phi > -20 + m) | ((phi < -150 - m) & (-60 + m < psi) & (psi < 30 - m))
    beta_wrap = (phi >= 150 - m) & ((psi >= 90 - m) | (psi <= -150 + m))
    right_in = (25 - m <= phi) & (phi <= 115 + m) & (-70 - m <= psi) & (psi <= 120 + m)
    return np.where(left, ~left_out, beta_wrap | right_in)


def _wrap(x):
    return (np.asarray(x, float) + 180.0) % 360.0 - 180.0


def _dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - b0 @ b1 * b1, b2 - b2 @ b1 * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def backbone_angles(atoms: list, xyz: list) -> list:
    """Residues of one chain in file order, each {"resid", "resname", "phi": (atom numbers, degrees) or None, "psi": ...}.

    atoms are .itp atoms (nr, resnr, residue, atom) and xyz the matching coordinates (A), in the same order. A residue starts wherever
    (resnr, residue) changes or an atom name repeats, so repeated or restarting residue numbers (insertion codes) stay separate. phi
    needs a peptide bond to the previous residue and psi one to the next, judged by the C-N distance: chain breaks give no angle."""
    residues, key, names = [], None, set()
    for a, b in zip(atoms, xyz):
        if key != (a["resnr"], a["residue"]) or a["atom"] in names or not residues:
            residues.append({"resid": a["resnr"], "resname": a["residue"], "atoms": {}})
            key, names = (a["resnr"], a["residue"]), set()
        names.add(a["atom"])
        if a["atom"] in ("N", "CA", "C"):
            residues[-1]["atoms"][a["atom"]] = (a["nr"], np.array([b["x"], b["y"], b["z"]], float))
    linked = lambda c, n: np.linalg.norm(c[1] - n[1]) < PEPTIDE_BOND_A
    out = []
    for i, r in enumerate(residues):
        at, phi, psi = r["atoms"], None, None
        if all(k in at for k in ("N", "CA", "C")):
            prev = residues[i - 1]["atoms"] if i else {}
            nxt = residues[i + 1]["atoms"] if i + 1 < len(residues) else {}
            if "C" in prev and linked(prev["C"], at["N"]):
                q = (prev["C"], at["N"], at["CA"], at["C"])
                phi = (tuple(x[0] for x in q), _dihedral(*(x[1] for x in q)))
            if "N" in nxt and linked(at["C"], nxt["N"]):
                q = (at["N"], at["CA"], at["C"], nxt["N"])
                psi = (tuple(x[0] for x in q), _dihedral(*(x[1] for x in q)))
        out.append({"resid": r["resid"], "resname": r["resname"], "phi": phi, "psi": psi})
    return out


def _run(mask: np.ndarray, centre: int) -> tuple:
    """(first, last) index of the run of True in a 1-D mask that contains `centre` (mask[centre] is True)."""
    false = np.flatnonzero(~mask)
    first = false[false < centre].max() + 1 if np.any(false < centre) else 0
    last = false[false > centre].min() - 1 if np.any(false > centre) else len(mask) - 1
    return int(first), int(last)


def _best_box(allowed: np.ndarray, c: int) -> tuple:
    """Largest all-True box of a square mask (phi rows, psi columns) containing the centre cell (c, c): (phi0, phi1, psi0, psi1).

    Every phi range [lo, hi] around c is tried; its psi range is the run around c of the rows ANDed together. Scored by the smaller
    half-width, then by area."""
    best, best_score = (c, c, c, c), (-1.0, -1.0)
    n = allowed.shape[0]
    for lo in range(c, -1, -1):
        if not allowed[lo, c]:
            break
        column = np.logical_and.reduce(allowed[lo:c + 1], axis=0)
        for hi in range(c, n):
            if hi > c:
                column = column & allowed[hi]
            if not column[c]:
                break
            s0, s1 = _run(column, c)
            score = (min(hi - lo, s1 - s0), (hi - lo + 1) * (s1 - s0 + 1))
            if score > best_score:
                best, best_score = (lo, hi, s0, s1), score
    return best


def restraint_windows(residues: list, mode: str = "clipped", inset: float = 5.0, half: float = BB_DIHRES_DPHI,
                      min_half: float = BB_DIHRES_MIN_HALF_DEG) -> list:
    """One row per phi/psi of `residues` (backbone_angles): {"index", "resid", "resname", "angle", "atoms", "input", "centre",
    "half_width", "low", "high", "status"}, low/high unwrapped around the input. See the module docstring for the statuses."""
    if mode not in ("fixed", "clipped"):
        raise ValueError(f"backbone restraint mode {mode!r} writes no windows")
    if not 0 <= inset < half:
        raise ValueError(f"inset {inset} must be at least 0 and below the window half-width {half}")
    steps = int(round(half / STEP_DEG))
    offsets = np.arange(-steps, steps + 1) * STEP_DEG
    rows = []
    for index, r in enumerate(residues):
        kinds = [k for k in ("phi", "psi") if r[k] is not None]
        if not kinds:
            continue
        box, status = None, "fixed"
        if mode == "clipped":
            parent = PARENT_RESIDUE.get(r["resname"])
            if len(kinds) < 2 or parent in (None, "GLY"):
                status = "not_checked"
            else:
                phi0, psi0 = r["phi"][1], r["psi"][1]
                if not rama_allowed(parent, _wrap(phi0), _wrap(psi0)):
                    status = "input_outlier"
                elif not rama_allowed(parent, _wrap(phi0), _wrap(psi0), inset):
                    status = "flag_margin"
                else:
                    allowed = rama_allowed(parent, _wrap(phi0 + offsets)[:, None], _wrap(psi0 + offsets)[None, :], inset)
                    if allowed.all():
                        status = "unchanged"
                    else:
                        p0, p1, s0, s1 = _best_box(allowed, steps)
                        if not allowed[p0:p1 + 1, s0:s1 + 1].all():   # cannot happen: the search only returns all-True boxes
                            raise RuntimeError(f"backbone restraint window of {r['resname']}{r['resid']} is not inside the allowed region")
                        lohi = {"phi": (phi0 + offsets[p0], phi0 + offsets[p1]), "psi": (psi0 + offsets[s0], psi0 + offsets[s1])}
                        if min(hi - lo for lo, hi in lohi.values()) / 2 < min_half:
                            status = "flag_narrow"
                        else:
                            box = lohi
        for k in kinds:
            x0 = r[k][1]
            lo, hi = box[k] if box else (x0 - half, x0 + half)
            st = status
            if box:
                st = "clipped" if hi - lo < 2 * half - 1e-9 else "unchanged"
            rows.append({"index": index, "resid": r["resid"], "resname": r["resname"], "angle": k, "atoms": r[k][0], "input": x0,
                         "centre": float(_wrap((lo + hi) / 2)), "half_width": (hi - lo) / 2, "low": lo, "high": hi, "status": st})
    return rows


def write_backbone_dihres(itp: Path, coords: Path, mode: str = "clipped", inset: float = 5.0) -> list:
    """Append the EM-only restraints (#ifdef BBDIHRES) to a chain's .itp and return its rows (restraint_windows), each with "molecule".

    coords is the pdb2gmx output with the source coordinates restored; its atoms are the .itp atoms in the same order."""
    atoms, xyz = read_itp(itp)["atoms"], read_pdb(coords)[0]
    if len(atoms) != len(xyz) or any(a["atom"] != b["atom"] for a, b in zip(atoms, xyz)):
        raise SystemExit(f"{itp.name}: atom order differs from {coords.name}; cannot place backbone dihedral restraints")
    rows = restraint_windows(backbone_angles(atoms, xyz), mode, inset)
    molecule = read_itp(itp)["mol"] or itp.stem
    for row in rows:
        row["molecule"] = molecule
    with itp.open("a") as fh:
        fh.write(f"\n#ifdef BBDIHRES\n[ dihedral_restraints ]\n; ai aj ak al type phi dphi kfac   (EM only, {mode} windows on input phi/psi)\n")
        fh.writelines(f"{a:6d}{b:6d}{c:6d}{d:6d}     1 {row['centre']:8.2f} {row['half_width']:6.2f}  DIHRES_FC\n"
                      for row in rows for a, b, c, d in [row["atoms"]])
        fh.write("#endif\n")
    return rows


def write_restraint_record(out: Path, rows: list) -> Path:
    """Write every restraint window (all chains) to <out>/backbone_restraints.tsv, outside toppar/ (the audit fails on toppar files
    the topology does not include)."""
    path = Path(out) / RECORD_NAME
    fmt = lambda v: f"{v:.2f}" if isinstance(v, float) else str(v)
    with path.open("w") as fh:
        fh.write("\t".join(RECORD_COLUMNS) + "\n")
        fh.writelines("\t".join(fmt(row[c]) for c in RECORD_COLUMNS) + "\n" for row in rows)
    return path


def summarize_record(rows: list) -> tuple:
    """(status counts as text, residues flagged for refinement as 'MOL:RESNAMEresid')."""
    counts = {}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    flagged = sorted({f"{row['molecule']}:{row['resname']}{row['resid']}" for row in rows if row["status"].startswith("flag_")})
    return ", ".join(f"{k} {counts[k]}" for k in sorted(counts)), flagged
