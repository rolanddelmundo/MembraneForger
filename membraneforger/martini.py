#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// The Martini 3 classification layer: the only place that decides what a CG residue is.
#//=============================================================
"""Martini 3 residue tables, classification of every coarse-grained residue, and name aliases for backmapping."""
import math
from collections import Counter

import numpy as np

from .config import RENAME_MOLECULE
from .structio import xyz_nm

__all__ = ['MARTINI3_PROTEIN', 'MARTINI2_PROTEIN', 'MARTINI3_PROTEIN_ALIASES', 'MARTINI3_SOLVENT',
           'MARTINI3_ION_BEADS', 'MARTINI2_ONLY', 'MARTINI3_MEMBRANE', 'MARTINI3_GLYCOLIPIDS', 'gro_safe_aliases',
           'structure_resname', 'cg_residues', 'membrane_entry', 'classify_cg_residue', 'classify_glycolipid',
           'classify_cg', 'make_membrane_whole']

# Side-chain beads per residue in Martini 3. Martini 2 differs (ALA 0, TYR 3, TRP 4) and is rejected.
MARTINI3_PROTEIN = {"GLY": 0, "ALA": 1, "CYS": 1, "VAL": 1, "LEU": 1, "ILE": 1, "MET": 1, "PRO": 1, "SER": 1,
                    "THR": 1, "ASN": 1, "GLN": 1, "ASP": 1, "GLU": 1, "LYS": 2, "ARG": 2, "HIS": 3, "PHE": 3,
                    "TYR": 4, "TRP": 5}


MARTINI2_PROTEIN = dict(MARTINI3_PROTEIN, ALA=0, TYR=3, TRP=4)


MARTINI3_PROTEIN_ALIASES = {"HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "HIH": "HIS", "HID": "HIS", "HIE": "HIS",
                            "HIP": "HIS", "ASH": "ASP", "ASPP": "ASP", "GLH": "GLU", "GLUP": "GLU", "LSN": "LYS",
                            "LYN": "LYS"}


MARTINI3_SOLVENT = {"W": ("W",), "SW": ("SW",), "TW": ("TW",)}


MARTINI3_ION_BEADS = {"NA", "CL", "CA", "K", "MG", "BR", "IOD", "TMA"}  # residue ION, or named after the bead


MARTINI2_ONLY = {"PW": "Martini 2 polarizable water", "WF": "Martini 2 antifreeze water",
                 "NA+": "Martini 2 ion naming", "CL-": "Martini 2 ion naming", "CA+": "Martini 2 ion naming"}


# Membrane residues whose Martini 3 name or bead set differs from the all-atom mapping. "unmapped" beads carry no
# atoms in the mapping and are dropped (and reported) before backmapping. Any other residue is a lipid only if the
# installed mapping has a residue of the same name with exactly the same beads.
MARTINI3_MEMBRANE = {
    "CHOL": dict(cls="sterol", aa="CHL1", unmapped=("R6",),
                 beads=("ROH", "R1", "R2", "R3", "R4", "R5", "R6", "C1", "C2")),
    "SAP6": dict(cls="phospholipid", aa="SAPI25", unmapped=("C4",),
                 beads=("C1", "C2", "C3", "C4", "PO4", "P4", "P5", "GL1", "GL2", "C1A", "D2A", "D3A", "D4A", "C5A",
                        "C1B", "C2B", "C3B", "C4B")),
}


# Martini 3 glycolipids come as consecutive sugar/ceramide residues with a virtual ring-centre bead V.
# Each is merged into one mapped residue; None marks a bead that is dropped.
MARTINI3_GLYCOLIPIDS = {
    "GM3": (("GLC", {"A": "G1", "B": "G2", "C": "G3", "V": None}),
            ("GAL", {"A": "A1", "B": "A2", "C": "A3", "V": None}),
            ("NMC", {"A": "S1", "B": "S2", "C": "S3", "D": "S4", "E": "S5", "V": None}),
            ("CER", {n: n for n in ("AM1", "AM2", "T1A", "C2A", "C3A", "C1B", "C2B", "C3B", "C4B")})),
}


def gro_safe_aliases(names: list[str]) -> dict[str, str]:
    """Give every residue name longer than the five-character .gro field a short alias that maps back uniquely."""
    long = sorted({n for n in names if len(n) > 5})
    alias = {name: f"M{k:03d}" for k, name in enumerate(long)}
    if set(alias.values()) & set(names):
        raise SystemExit(f"residue names {sorted(set(alias.values()) & set(names))} collide with backmapping aliases")
    return alias


def structure_resname(name: str) -> str:
    """Return the four-character PDB spelling of an all-atom membrane residue and require it to round-trip."""
    short = name[:4]
    if RENAME_MOLECULE.get(short, short) not in (name, RENAME_MOLECULE.get(name, name)):
        raise SystemExit(f"residue {name} does not survive the 4-character PDB field (reads back as "
                         f"{RENAME_MOLECULE.get(short, short)}); add it to RENAME_MOLECULE")
    return short


def cg_residues(atoms: list[dict]) -> list[list[dict]]:
    """Group consecutive beads that share chain, residue number and residue name."""
    residues, last = [], None
    for a in atoms:
        key = (a["chain"], a["resid"], a["resname"])
        if key != last:
            residues.append([])
            last = key
        residues[-1].append(a)
    return residues


def membrane_entry(cg: str, cls: str, aa: str, beads: list[dict], dropped: list[str], mapping: dict) -> tuple[str, object]:
    """Build one membrane molecule record after checking the mapping covers exactly the beads that are kept."""
    kept = [b["atom"] for b in beads]
    if aa not in mapping or sorted(mapping[aa]["beads"]) != sorted(kept):
        return "unsupported", f"{cg}: the installed mapping for {aa} does not cover Martini 3 beads {kept}"
    return cls, {"cg": cg, "cls": cls, "aa": aa, "beads": [dict(b, resname=aa) for b in beads], "dropped": dropped}


def classify_cg_residue(beads: list[dict], mapping: dict) -> tuple[str, object]:
    """Classify one Martini 3 residue from its residue name and bead inventory."""
    # Returns (class, payload): payload is the reason for "unsupported", a molecule record for membrane classes.
    resname, names = beads[0]["resname"], [b["atom"] for b in beads]
    label = f"{resname} [{' '.join(names)}]"
    base = MARTINI3_PROTEIN_ALIASES.get(resname, resname)
    old = next((MARTINI2_ONLY[n] for n in [resname] + names if n in MARTINI2_ONLY), None)
    if old:
        return "unsupported", f"{label}: {old}"
    if base in MARTINI3_PROTEIN and "BB" in names:
        side = lambda table: ["BB"] + [f"SC{k}" for k in range(1, table[base] + 1)]
        if names == side(MARTINI3_PROTEIN):
            return "protein", None
        if names == side(MARTINI2_PROTEIN):
            return "unsupported", f"{label}: Martini 2 protein bead inventory"
        return "unsupported", f"{label}: not the Martini 3 beads {' '.join(side(MARTINI3_PROTEIN))}"
    if resname in MARTINI3_SOLVENT and tuple(names) == MARTINI3_SOLVENT[resname]:
        return "solvent", None
    if (resname == "ION" or resname in MARTINI3_ION_BEADS) and len(names) == 1 and names[0] in MARTINI3_ION_BEADS:
        return "ion", None
    if resname in MARTINI3_MEMBRANE:
        spec = MARTINI3_MEMBRANE[resname]
        if sorted(names) != sorted(spec["beads"]):
            legacy = sorted(names) == sorted(mapping.get(resname, {}).get("beads", ["-"]))
            return "unsupported", (f"{label}: not the Martini 3 bead set {' '.join(spec['beads'])}"
                                   + (" (Martini 2 inventory)" if legacy else ""))
        return membrane_entry(resname, spec["cls"], spec["aa"], [b for b in beads if b["atom"] not in spec["unmapped"]],
                              list(spec["unmapped"]), mapping)
    if resname in mapping and base not in MARTINI3_PROTEIN and sorted(names) == sorted(mapping[resname]["beads"]):
        return membrane_entry(resname, "phospholipid" if "PO4" in names else "lipid", resname, beads, [], mapping)
    return "unsupported", f"{label}: unknown residue or bead set, no all-atom mapping"


def classify_glycolipid(group: list[list[dict]], name: str, mapping: dict) -> tuple[str, object]:
    """Merge the consecutive sugar/ceramide residues of one Martini 3 glycolipid into a single mapped molecule."""
    parts = MARTINI3_GLYCOLIPIDS[name]
    complete = len(group) == len(parts) and all(
        g[0]["resname"] == part and sorted(b["atom"] for b in g) == sorted(rename) for g, (part, rename) in zip(group, parts))
    if not complete:
        first = group[0]
        return "unsupported", (f"{first[0]['resname']} [{' '.join(b['atom'] for b in first)}]: not a complete Martini 3 "
                               f"{name} ({'+'.join(part for part, _ in parts)})")
    merged = [dict(b, atom=rename[b["atom"]]) for g, (_, rename) in zip(group, parts) for b in g if rename[b["atom"]]]
    dropped = [f"{part}:{bead}" for part, rename in parts for bead, new in rename.items() if new is None]
    return membrane_entry(name, "glycolipid", name, merged, dropped, mapping)


def classify_cg(atoms: list[dict], mapping: dict) -> dict:
    """Sort every coarse-grained residue into a Martini 3 class and refuse anything unsupported."""
    # `mapping` is {all-atom residue: {"beads": [...], "atoms": [...]}} as read from the installed mapping files.
    residues = cg_residues(atoms)
    protein, membrane, counts, problems = [], [], Counter(), Counter()
    i = 0
    while i < len(residues):
        glyco = next((g for g, parts in MARTINI3_GLYCOLIPIDS.items() if parts[0][0] == residues[i][0]["resname"]), None)
        if glyco:
            group = residues[i:i + len(MARTINI3_GLYCOLIPIDS[glyco])]
            cls, payload = classify_glycolipid(group, glyco, mapping)
            i += len(group) if cls != "unsupported" else 1
        else:
            cls, payload = classify_cg_residue(residues[i], mapping)
            if cls == "protein":
                protein.append(residues[i])
            i += 1
        if cls == "unsupported":
            problems[payload] += 1
            continue
        counts[cls] += 1
        if isinstance(payload, dict):
            membrane.append(payload)
    if problems:
        raise SystemExit("unsupported coarse-grained content (Martini 3 is required): "
                         + "; ".join(f"{n} x {what}" for what, n in sorted(problems.items())))
    if not protein:
        raise SystemExit("the coarse-grained input has no Martini 3 protein (no BB beads)")
    if not membrane:
        raise SystemExit("the coarse-grained input has no membrane lipids")
    return {"protein": protein, "membrane": membrane, "counts": counts,
            "composition": Counter(m["cg"] for m in membrane),
            "dropped": {m["cg"]: m["dropped"] for m in membrane if m["dropped"]}}


def make_membrane_whole(membrane: list[dict], box: list[float]) -> tuple[float, float]:
    """Make every membrane molecule whole and bring the bilayer into one periodic image in z."""
    # Returns its z range (nm).
    cell = np.array(box)
    centres = []
    for mol in membrane:
        xyz = xyz_nm(mol["beads"])
        xyz = xyz[0] + (xyz - xyz[0]) - cell * np.round((xyz - xyz[0]) / cell)
        centres.append(xyz.mean(axis=0)[2])
        mol["xyz"] = xyz
    angle = 2.0 * math.pi * np.array(centres) / cell[2]
    middle = (math.atan2(np.sin(angle).mean(), np.cos(angle).mean()) / (2.0 * math.pi) * cell[2]) % cell[2]
    for mol, z in zip(membrane, centres):
        mol["xyz"][:, 2] -= cell[2] * round((z - middle) / cell[2])
        for bead, (x, y, zz) in zip(mol["beads"], mol["xyz"]):
            bead["x"], bead["y"], bead["z"] = float(x), float(y), float(zz)
    low = min(float(mol["xyz"][:, 2].min()) for mol in membrane)
    high = max(float(mol["xyz"][:, 2].max()) for mol in membrane)
    if high - low > 0.5 * cell[2]:
        raise SystemExit(f"membrane beads span {high - low:.1f} nm of a {cell[2]:.1f} nm box in z; "
                         "expected one planar bilayer normal to z")
    return low, high
