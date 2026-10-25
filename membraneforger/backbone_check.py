#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Protein backbone check after minimization.
#//=============================================================
"""Compare the protein backbone after minimization with the structure minimization started from.

The audit looks at C-alpha positions and contacts; nothing there looks at peptide planes, alpha-carbon handedness or phi/psi,
so a build could PASS with a twisted peptide bond or an inverted centre. This check measures what EM did, per residue, against
solv_ions.gro (the same atoms in the same order, before EM). Nothing here needs to know the protein, the lipids or the force field.

FAIL (the build stops; --backbone-check yes):
  omega        a peptide bond that was trans (or, for Xaa-Pro, cis or trans) and is twisted or has changed isomer after EM
  chirality    a C-alpha (GLY and residues without CB excepted) or a Thr/Ile C-beta whose handedness differs from the start; the
               start decides, so D residues in the input are not an error
  drift        a phi or psi that moved by more than BB_DIHRES_DPHI + --backbone-drift-margin degrees; applied only when the EM-only
               backbone restraints (#ifdef BBDIHRES) were active, since without them any drift is possible
  strict       a residue of a --strict-segments range (segments the caller knows were modelled, not observed) that is a
               Ramachandran outlier after EM, whatever it was before
WARN (reported, never stops a build): new Ramachandran outliers elsewhere, new omega between 120 and 150 degrees, bonds that were
already twisted in the start.

What it does not do: tell a modelled residue from an observed one (give --strict-segments for that), judge the input, or catch a
residue that drifts inside the restraint width into a borderline region (--backbone-restraints clipped keeps the restraint windows
inside the allowed region, backbone_restraints.py). Measured on five restrained builds (GIPR/Gs, 1214 protein
residues each) the largest drift was 35.5 degrees and none moved more than 45; on 140 older builds without restraints the median build
had 10 residues drift by more than 60 degrees, and 65 of them gained at least one Ramachandran outlier.
"""
import json
from pathlib import Path

import numpy as np

from .config import AMINO, BB_DIHRES_DPHI
from .structio import read_gro

__all__ = ['PARENT_RESIDUE', 'parse_segments', 'rama_outlier', 'backbone_residues', 'chain_labels', 'check_backbone',
           'write_backbone_report']

# Parent amino acid of every residue name the build accepts; names with no standard parent (AIB, LEM, KSM...) keep their own
# geometry checks but are exempt from the Ramachandran regions, which are those of the 20 standard residues.
PARENT_RESIDUE = {n: n for n in "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()}
PARENT_RESIDUE.update({"HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "CYSG": "CYS", "CYSM": "CYS", "CYSP": "CYS", "GLYP": "GLY",
                       "KSM": "LYS", "KTZ": "LYS", "KRT": "LYS"})
NEEDED = {"N", "CA", "C", "CB", "OG1", "CG1", "CG2"}
BOND_NM = 0.20      # C(i)-N(i+1) shorter than this is a peptide bond (the genuine one is 0.133)
TWISTED_DEG, WARN_DEG = 120.0, 150.0
CIS_DEG = 45.0      # Xaa-Pro: |omega| below this is cis, above TWISTED_DEG is trans, between is twisted


def parse_segments(text: str) -> tuple:
    """'R:338-345,A:5-9,B:12' -> (('R', 338, 345), ('A', 5, 9), ('B', 12, 12)); ValueError when a part is malformed."""
    segments = []
    for part in (text or "").replace(" ", "").split(","):
        if not part:
            continue
        chain, _, span = part.partition(":")
        lo, _, hi = span.partition("-")
        if not chain or not lo.lstrip("-").isdigit() or (hi and not hi.lstrip("-").isdigit()):
            raise ValueError(f"segment {part!r} is not CHAIN:FIRST-LAST")
        first, last = int(lo), int(hi or lo)
        if last < first:
            raise ValueError(f"segment {part!r} ends before it starts")
        segments.append((chain, first, last))
    return tuple(segments)


def rama_outlier(resname: str, phi: float, psi: float, margin: float = 0.0) -> bool:
    """Generous Ramachandran regions (the ones the md_cognate geometry gate uses); `margin` widens every allowed region by that many degrees.

    GLY is the caller's to exempt. PRO has its own phi window; every other residue is outside only in the phi > -20 band of the
    left half, the phi ~ -180 / psi ~ -60..30 corner, and the right half outside the alpha-L / bridge window."""
    if resname == "PRO":
        return not (-110 - margin <= phi <= -35 + margin)
    if phi <= 0:
        return phi > -20 + margin or (phi < -150 - margin and -60 + margin < psi < 30 - margin)
    if phi >= 150 - margin and (psi >= 90 - margin or psi <= -150 + margin):
        return False
    return not (25 - margin <= phi <= 115 + margin and -70 - margin <= psi <= 120 + margin)


def _dih(p0, p1, p2, p3, box):
    """Dihedral (degrees) of four points, every difference taken as the shortest periodic image when a box is given."""
    image = (lambda v: v - box * np.round(v / box)) if box is not None else (lambda v: v)
    b0, b1, b2 = image(p0 - p1), image(p2 - p1), image(p3 - p2)
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - b0 @ b1 * b1, b2 - b2 @ b1 * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def _volume(centre, a, b, c, box):
    """Signed volume of the three bonds a, b, c around `centre` (shortest images)."""
    image = (lambda v: v - box * np.round(v / box)) if box is not None else (lambda v: v)
    return float(image(a - centre) @ np.cross(image(b - centre), image(c - centre)))


def backbone_residues(path: Path) -> dict:
    """Protein residues of a .gro in file order: {"box": [x, y, z] nm, "residues": [{"name", "resid", "atoms": {N, CA, C, CB...}}]}."""
    atoms, box = read_gro(Path(path))
    residues, key = [], None
    for a in atoms:
        if a["resname"] not in AMINO:
            continue
        if key != (a["resid"], a["resname"]) or not residues:
            residues.append({"name": a["resname"], "resid": a["resid"], "atoms": {}})
            key = (a["resid"], a["resname"])
        if a["atom"] in NEEDED:
            residues[-1]["atoms"][a["atom"]] = np.array([a["x"], a["y"], a["z"]])
    return {"box": np.array(box[:3], float), "residues": residues}


def chain_labels(system: dict) -> list:
    """[(chain, resid, resname)] for the protein residues of `system["chains"]` (input order), the labels of the .gro residues."""
    labels = []
    for chain, atoms in system["chains"].items():
        last = None
        for a in atoms:
            if a["resname"] in AMINO or a["resname"] in PARENT_RESIDUE:
                k = (a["resid"], a["resname"])
                if k != last:
                    labels.append((chain, a["resid"], a["resname"]))
                    last = k
    return labels


def _angles(structure: dict) -> list:
    """Per residue: phi, psi, omega (the bond to the next residue) in degrees or None, the C-alpha and C-beta handedness (+1 L, -1 D, 0 n/a)."""
    box, res = structure["box"], structure["residues"]
    image = lambda v: v - box * np.round(v / box)
    out = []
    for i, r in enumerate(res):
        at = r["atoms"]
        phi = psi = omega = None
        has = all(k in at for k in ("N", "CA", "C"))
        prev = res[i - 1]["atoms"] if i else None
        nxt = res[i + 1]["atoms"] if i + 1 < len(res) else None
        if has and prev and "C" in prev and np.linalg.norm(image(prev["C"] - at["N"])) < BOND_NM:
            phi = _dih(prev["C"], at["N"], at["CA"], at["C"], box)
        if has and nxt and all(k in nxt for k in ("N", "CA")) and np.linalg.norm(image(at["C"] - nxt["N"])) < BOND_NM:
            psi = _dih(at["N"], at["CA"], at["C"], nxt["N"], box)
            omega = _dih(at["CA"], at["C"], nxt["N"], nxt["CA"], box)
        ca = cb = 0
        parent = PARENT_RESIDUE.get(r["name"])
        if has and "CB" in at and parent != "GLY":
            ca = 1 if _volume(at["CA"], at["N"], at["C"], at["CB"], box) > 0 else -1
        sub = {"THR": "OG1", "ILE": "CG1"}.get(parent)
        if sub and all(k in at for k in ("CB", "CA", sub, "CG2")):
            cb = 1 if _volume(at["CB"], at["CA"], at[sub], at["CG2"], box) > 0 else -1
        out.append((phi, psi, omega, ca, cb))
    return out


def _omega_class(omega: float, proline_next: bool) -> str:
    """trans / cis (Xaa-Pro only) / twisted (anything else short of trans)."""
    a = abs(omega)
    if a >= TWISTED_DEG:
        return "trans"
    if proline_next:
        return "cis" if a <= CIS_DEG else "twisted"
    return "twisted"


def _circ(a: float, b: float) -> float:
    """Angular difference in degrees, 0..180."""
    return abs((a - b + 180.0) % 360.0 - 180.0)


def _pp(phi: float, psi: float | None) -> str:
    """'phi/psi' for a message; a missing psi (chain end) prints as '-'."""
    return f"{phi:+.0f}/{psi:+.0f}" if psi is not None else f"{phi:+.0f}/-"


def check_backbone(before: dict, after: dict, labels: list | None = None, restrained: bool = True, margin_deg: float = 30.0,
                   strict: tuple = (), mode: str = "yes") -> dict:
    """Judge the backbone after EM against the structure before it; returns the report (see the module docstring).

    before/after come from backbone_residues; labels (chain, resid, resname) name the residues when they match one for one.
    mode "yes": failures stop the build; "report": failures are listed as warnings; "no" is handled by the caller."""
    rb, ra = before["residues"], after["residues"]
    report = {"status": "PASS", "mode": mode, "residues": len(ra), "restraints_active": bool(restrained),
              "drift_limit_deg": BB_DIHRES_DPHI + margin_deg, "failures": [], "warnings": []}
    if len(rb) != len(ra) or [r["name"] for r in rb] != [r["name"] for r in ra]:
        report.update(status="skipped", reason=f"the protein residues before ({len(rb)}) and after ({len(ra)}) EM do not correspond")
        report["warnings"].append(report["reason"])
        return report
    if not ra:
        report.update(status="skipped", reason="no protein residues")
        return report
    named = labels is not None and len(labels) == len(ra)
    label = (lambda i: f"{labels[i][0]}:{labels[i][2]}{labels[i][1]}") if named else (lambda i: f"{ra[i]['name']}{ra[i]['resid']}")
    where = (lambda i: (labels[i][0], labels[i][1])) if named else (lambda i: ("?", ra[i]["resid"]))
    report["labels"] = "input chain and residue number" if named else "residue numbers of the .gro"
    A0, A1 = _angles(before), _angles(after)
    made_twisted, inherited, new_warn, flipped, drift, new_outliers, strict_out = [], 0, [], [], [], [], []
    compared = 0
    for i, ((p0, s0, o0, c0, b0), (p1, s1, o1, c1, b1)) in enumerate(zip(A0, A1)):
        parent = PARENT_RESIDUE.get(ra[i]["name"])
        if o0 is not None and o1 is not None:
            compared += 1
            pro = PARENT_RESIDUE.get(ra[i + 1]["name"]) == "PRO"
            k0, k1 = _omega_class(o0, pro), _omega_class(o1, pro)
            bond = f"{label(i)}-{label(i + 1)} omega {o0:+.0f} -> {o1:+.0f}"
            if k0 != k1 and k0 != "twisted":      # EM made a trans or cis bond twisted, or changed its isomer; a bond that was
                made_twisted.append(f"{bond} ({k0} -> {k1})")       # already twisted may relax, which is no failure
            elif k0 == k1 == "twisted":
                inherited += 1
            elif k0 == k1 == "trans" and abs(o1) < WARN_DEG <= abs(o0):
                new_warn.append(bond)
        if c0 and c1 and c0 != c1:
            flipped.append(f"{label(i)} C-alpha {'L -> D' if c0 > 0 else 'D -> L'}")
        if b0 and b1 and b0 != b1:
            flipped.append(f"{label(i)} C-beta inverted")
        if p0 is not None and p1 is not None:
            drift.append((max(_circ(p0, p1), _circ(s0, s1)) if s0 is not None and s1 is not None else _circ(p0, p1), i, p0, s0, p1, s1))
            if parent and parent != "GLY" and s0 is not None and s1 is not None:
                out0, out1 = rama_outlier(parent, p0, s0), rama_outlier(parent, p1, s1)
                if out1 and not out0:
                    new_outliers.append(f"{label(i)} phi/psi {_pp(p0, s0)} -> {_pp(p1, s1)}")
                chain, number = where(i)
                if out1 and any(chain == c and lo <= number <= hi for c, lo, hi in strict):
                    strict_out.append(f"{label(i)} phi/psi {_pp(p1, s1)} (before {_pp(p0, s0)}) in a strict segment")
    drift.sort(reverse=True)
    report["bonds_compared"] = compared
    report["max_drift_deg"] = round(drift[0][0], 1) if drift else 0.0
    report["worst_drift"] = [f"{label(i)} phi/psi {_pp(p0, s0)} -> {_pp(p1, s1)} ({d:.0f} deg)" for d, i, p0, s0, p1, s1 in drift[:5]]
    over = [x for x in drift if x[0] > report["drift_limit_deg"]]
    report["omega_made_twisted"], report["omega_twisted_before_and_after"] = made_twisted, inherited
    report["chirality_changed"], report["strict_segments"] = flipped, [f"{c}:{lo}-{hi}" for c, lo, hi in strict]
    report["strict_segment_outliers"], report["new_ramachandran_outliers"] = strict_out, new_outliers
    failures = [f"peptide bond {m}" for m in made_twisted] + [f"stereochemistry changed: {m}" for m in flipped] + \
               [f"strict segment: {m}" for m in strict_out]
    if restrained and over:
        failures += [f"phi/psi moved more than {report['drift_limit_deg']:g} deg with the backbone restraints on: {label(i)} "
                     f"{_pp(p0, s0)} -> {_pp(p1, s1)} ({d:.0f} deg)" for d, i, p0, s0, p1, s1 in over[:5]]
        failures += [f"... and {len(over) - 5} more residues moved that far"] if len(over) > 5 else []
    elif not restrained and over:
        report["warnings"].append(f"{len(over)} residue(s) moved more than {report['drift_limit_deg']:g} deg, with no backbone restraints in "
                                  f"this build (largest {drift[0][0]:.0f} deg)")
    if new_outliers:
        report["warnings"].append(f"{len(new_outliers)} residue(s) became Ramachandran outliers, e.g. {new_outliers[0]}")
    if new_warn:
        report["warnings"].append(f"{len(new_warn)} peptide bond(s) moved below {WARN_DEG:g} deg from trans, e.g. {new_warn[0]}")
    if inherited:
        report["warnings"].append(f"{inherited} peptide bond(s) were already twisted before EM and still are")
    if mode == "report":
        report["warnings"] += [f"(would fail) {m}" for m in failures]
        failures = []
    report["failures"] = failures
    report["status"] = "FAIL" if failures else "PASS"
    return report


def write_backbone_report(out: Path, report: dict) -> None:
    """Save the report as backbone_check.json in the output directory."""
    Path(out, "backbone_check.json").write_text(json.dumps(report, indent=1, default=str) + "\n")
