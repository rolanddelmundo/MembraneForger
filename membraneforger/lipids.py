#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Centralized lipid metadata: the one reference site of every lipid at both resolutions.
#//=============================================================
"""Which bead or atom stands for a lipid when it is sliced, assigned to a leaflet, tessellated or correlated.

Slicing, the area-per-lipid analysis and the lateral radial distribution functions all need ONE representative
in-plane position per lipid. It is the headgroup site, which is what "where does this lipid sit in its leaflet"
means chemically; tail beads wander by up to a lipid length and must not decide anything. The Martini 3 bead names
are those of the installed mapping (backmap_data/map.dat, membraneforger/martini.py); the atomistic names are those
of the CHARMM36 topologies under forcefield/toppar/.

    phospholipids (POPC, DOPC, POPE, DOPE, POPS, DOPS, PSM)   PO4                     <->  P   (the phosphate)
    PIP2 (Martini SAP6, atomistic SAPI25)                     PO4 (the diester, not P4/P5) <->  P
    cholesterol (CHOL, atomistic CHL1)                        ROH                     <->  O3  (the 3-beta hydroxyl)
    GM3 (Martini DPG3, atomistic GLPA)                        centroid of the three sugar residues' beads
                                                              (GLC G1-G3, GAL A1-A3, NMC S1-S5) <-> the amide N (NF)
    anything else                                             centre of geometry of the whole residue

An anchor given as several site names means the centroid of those sites. A lipid that has no entry falls back to the
centre of geometry of its (whole) residue, so an unlisted species is analysed, not refused. Names here are the
classifier's spellings (martini.MARTINI3_MEMBRANE, embedding.LIPID_ALIASES): SAP6 is PIP2 and GM3 is DPG3; the
all-atom residues CHL1, SAPI25 and GLPA map back to them through AA_TO_CG_NAME.
"""
import numpy as np

__all__ = ['LIPID_SLICE_ANCHORS', 'LIPID_AA_ANCHORS', 'AA_TO_CG_NAME', 'STEROLS', 'LIPID_DISPLAY_NAMES', 'cg_name',
           'display_name', 'anchor_bead', 'anchor_indices', 'anchor_xyz', 'is_sterol', 'leaflet_of']

# Martini 3 site(s) that represent each lipid in the membrane plane (slicing anchor, Voronoi generator, RDF site).
LIPID_SLICE_ANCHORS = {
    "POPC": "PO4", "DOPC": "PO4", "POPE": "PO4", "DOPE": "PO4", "POPS": "PO4", "DOPS": "PO4", "PSM": "PO4",
    "SAP6": "PO4",   # PIP2: the diester phosphate, not the inositol ring phosphates P4/P5
    "CHOL": "ROH",
    "GM3": ("G1", "G2", "G3", "A1", "A2", "A3", "S1", "S2", "S3", "S4", "S5"),  # DPG3: centroid of the three sugars
}

# The atomistic site(s) that correspond to each anchor, by CHARMM36 residue name (forcefield/toppar/*.itp).
LIPID_AA_ANCHORS = {
    "POPC": "P", "DOPC": "P", "POPE": "P", "DOPE": "P", "POPS": "P", "DOPS": "P", "PSM": "P",
    "SAPI25": "P", "SAPI": "P",
    "CHL1": "O3",
    "GLPA": "NF",    # the ceramide amide nitrogen, the atom the sugars are attached to (sugar atom names repeat in GLPA)
}

# All-atom residue name -> the classifier's lipid name, so that every stage reports the same species labels.
AA_TO_CG_NAME = {"CHL1": "CHOL", "SAPI25": "SAP6", "SAPI": "SAP6", "GLPA": "GM3", "GM3": "GM3"}

# Lipids without a headgroup in the interfacial plane: left out of the headgroup-plane tessellation, where their
# anchor would sit about 0.5 nm below the phosphates.
STEROLS = frozenset({"CHOL"})

# How the species are written in tables and figures (the INSANE/user-facing names of embedding.LIPID_NAMES).
LIPID_DISPLAY_NAMES = {"SAP6": "PIP2 (SAP6)", "GM3": "GM3 (DPG3)"}


def cg_name(resname: str) -> str:
    """The classifier's lipid name for an all-atom or coarse-grained residue name."""
    return AA_TO_CG_NAME.get(resname, resname)


def display_name(cg: str) -> str:
    """The user-facing name of a lipid species."""
    return LIPID_DISPLAY_NAMES.get(cg, cg)


def anchor_bead(cg: str, atomistic: bool = False) -> str | tuple | None:
    """The anchor site name(s) of a lipid species at the given resolution, or None (centre-of-geometry fallback)."""
    return (LIPID_AA_ANCHORS if atomistic else LIPID_SLICE_ANCHORS).get(cg)


def anchor_indices(names: list[str], anchor: str | tuple | None) -> list[int]:
    """Indices of the anchor site(s) among the bead/atom names; empty when they are absent (centre of geometry)."""
    wanted = (anchor,) if isinstance(anchor, str) else tuple(anchor or ())
    found = [i for i, n in enumerate(names) if n in wanted]
    return found if len(found) == len(wanted) else []


def anchor_xyz(xyz: np.ndarray, names: list[str], anchor: str | tuple | None) -> np.ndarray:
    """Representative position of one WHOLE residue: its anchor site(s), or its centre of geometry when there is none."""
    i = anchor_indices(names, anchor)
    return xyz[i].mean(axis=0).astype(float) if i else xyz.mean(axis=0).astype(float)


def is_sterol(cg: str) -> bool:
    """True for lipids whose anchor lies below the headgroup plane (left out of the headgroup-only tessellation)."""
    return cg in STEROLS


def leaflet_of(anchor_z: float, midplane: float) -> str:
    """'lower' or 'upper' from the z of a lipid's anchor relative to the bilayer midplane."""
    return "lower" if anchor_z < midplane else "upper"
