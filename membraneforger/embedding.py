#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Embedding an all-atom complex into a membrane simulated with another protein.
#//=============================================================
"""Place an all-atom complex into a pre-equilibrated Martini 3 membrane, trim the patch, and edit its lipids."""
import math
from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

from .config import AMINO, MIN_Z_PAD_TOTAL_NM, SLAB_Z_PAD_NM
from .lipids import anchor_bead, anchor_xyz, leaflet_of
from .slicing import relax_seam
from .structio import element, wrap, xyz_nm

__all__ = ['LIPID_NAMES', 'LIPID_ALIASES', 'CONVERTIBLE', 'BELT_NM', 'OVERLAP_NM', 'HYDROPHOBIC', 'CHARGED',
           'EXPOSURE_RADIUS_NM', 'BURIED_ABOVE', 'HARD_CORE_NM', 'MAX_VOID_NM3', 'lipid_name', 'hydrophobic_belt',
           'periodic_mean', 'phosphate_planes', 'membrane_voids', 'EMBED_SITES', 'FRAME_PROTEIN_CLEARANCE_NM',
           'image_onto_slab', 'frame_protein_beads', 'free_site', 'rectangle_distance', 'slice_window',
           'frame_protein_in_slice',
           'embed_complex', 'trim_membrane', 'edit_lipids', 'check_box_z']

# User-facing Martini 3 lipid names (INSANE spelling) and how they are spelled inside the classifier.
LIPID_NAMES = ("CHOL", "POPC", "DOPC", "POPE", "DOPE", "POPS", "DOPS", "PSM", "DPG3", "SAP6")
LIPID_ALIASES = {"DPG3": "GM3", "GM3": "GM3", "PIP2": "SAP6", "SAPI25": "SAP6", "CHL1": "CHOL", "DPSM": "PSM"}
# Lipids with the same twelve-bead Martini 3 layout: one can be turned into another by renaming its beads.
CONVERTIBLE = ("POPC", "DOPC", "POPE", "DOPE", "POPS", "DOPS", "PSM")
BELT_NM = 3.0  # hydrophobic thickness of the bilayer: the window that locates the transmembrane region
# Making room for the complex (embed_complex). Lipids closer than OVERLAP_NM to a heavy atom of the complex are pushed
# away at the coarse-grained level (repulsion range OVERLAP_NM + PUSH_MARGIN_NM, so restrained beads settle at about
# OVERLAP_NM) for at most RELAX_STEPS steps, without creating inter-lipid contacts closer than SEAM_CONTACT_NM
# (config.Settings.seam_min_bead_nm). Only lipids still closer than HARD_CORE_NM afterwards are removed: the complex
# occupies their place. On the 6WHC receptor-Gs complex the earlier rule (remove every lipid with any bead within
# OVERLAP_NM) deleted about 99 lipids per bundled frame, 39 of them only for touching the Ggamma geranylgeranyl chain
# or the Galpha N-terminus, and left an empty column of about 4 nm^3 in the inner leaflet. Removing lipids whose
# headgroup anchor the complex occupies was also tried and rejected: next to a thin lipid anchor it still deletes
# lipids whose volume nothing replaces (3 nm^3 pocket on GPR1). Over the 18 bundled frames this embedding removes
# 24-41 lipids instead of 88-110. RELAX_STEPS: the largest pocket stops shrinking after about 1200 steps on GPR1
# (600 steps: 1.5 nm^3, 1200: 0.44, 4000: 0.44); one embedding takes about a minute.
OVERLAP_NM = 0.40
PUSH_MARGIN_NM = 0.05
HARD_CORE_NM = 0.30
RELAX_STEPS = 1200
SEAM_CONTACT_NM = 0.30
# Empty pockets in the acyl region (membrane_voids): grid points farther than VOID_EMPTY_NM from every lipid bead and
# solute atom, on a VOID_GRID_NM grid between VOID_MARGIN_NM off the midplane and VOID_MARGIN_NM inside the PO4 plane.
# The embedding stops when a pocket exceeds MAX_VOID_NM3, less than the volume of one phospholipid. With the 6WHC
# complex in the 18 bundled frames: the uncut frames' largest pocket is 0.00-0.20 nm^3, this embedding leaves
# 0.02-0.76 nm^3, and the earlier remove-on-contact rule left 1.2-6.8 nm^3 in the inner leaflet.
# Where the complex goes in the frame (embed_complex site): "hole", where the frame's own protein was, or "free", the
# point farthest from that protein, searched on a FREE_SITE_GRID_NM grid.
EMBED_SITES = ("hole", "free")
FREE_SITE_GRID_NM = 0.25
# With the free site, the slice window must keep this far from every bead of the frame's receptor. In the bundled
# frames no lipid bead comes closer than 0.35-0.37 nm to a receptor bead (GPR1, KOR1, GPR5), so the hole reaches that
# far past the bead centres; 1.0 nm keeps the cut edge a further ~0.6 nm (about one lipid bead) off it.
FRAME_PROTEIN_CLEARANCE_NM = 1.0
VOID_GRID_NM = 0.1
VOID_EMPTY_NM = 0.60
VOID_MARGIN_NM = 0.30
MAX_VOID_NM3 = 1.0
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "CYS", "GLY", "PRO"}
CHARGED = {"ASP", "GLU", "LYS", "ARG"}
# A residue is on the surface when at most this many protein heavy atoms lie within this radius of its side chain.
# On the example receptor-G protein complex surface side chains have 40-60 such neighbours and buried ones 80-110.
EXPOSURE_RADIUS_NM, BURIED_ABOVE = 0.8, 60


def lipid_name(name: str) -> str:
    """Return the classifier's spelling of a user-given lipid name, refusing names that are not supported."""
    upper = name.strip().upper()
    canonical = LIPID_ALIASES.get(upper, upper)
    if canonical not in {LIPID_ALIASES.get(n, n) for n in LIPID_NAMES}:
        raise SystemExit(f"unknown lipid {name}; choose from {', '.join(LIPID_NAMES)}")
    return canonical


def hydrophobic_belt(aa_atoms: list[dict], width_nm: float = BELT_NM, step_nm: float = 0.05,
                     exposure_radius_nm: float = EXPOSURE_RADIUS_NM, buried_above: int = BURIED_ABOVE) -> dict:
    """Find the z (A) of the slab of width `width_nm` whose exposed residues are the most hydrophobic."""
    # The all-atom input must already have its membrane normal along z; this only decides where along z the
    # bilayer centre lies, which is where the membrane's midplane is put. Only surface residues count (a side chain
    # with at most `buried_above` protein heavy atoms within `exposure_radius_nm` of its centroid): the surface of a
    # transmembrane segment is hydrophobic, the surface of a soluble domain is not, whatever their cores hold.
    heavy = [a for a in aa_atoms if a["resname"] in AMINO and element(a["atom"]) != "H"]
    residues = {}
    for a in heavy:
        residues.setdefault((a["chain"], a.get("segid", ""), a["resid"], a["resname"]), []).append(a)
    if len(residues) < 10:
        raise SystemExit("the all-atom input has fewer than 10 protein residues; cannot locate its transmembrane region")
    xyz = xyz_nm(heavy) / 10.0
    tree = cKDTree(xyz)
    centroids, kind = [], []
    for key, atoms in residues.items():
        side = [a for a in atoms if a["atom"] not in ("N", "CA", "C", "O", "OXT")] or [a for a in atoms if a["atom"] == "CA"] or atoms
        centroids.append(xyz_nm(side).mean(axis=0) / 10.0)
        kind.append(1 if key[3] in HYDROPHOBIC else -1 if key[3] in CHARGED else 0)
    centroids, kind = np.array(centroids), np.array(kind)
    neighbours = np.array([len(hits) for hits in tree.query_ball_point(centroids, exposure_radius_nm)])
    exposed = neighbours <= buried_above
    z = centroids[:, 2]
    order = np.argsort(z)
    z, weight = z[order], np.where(exposed[order], kind[order], 0)
    cumulative = np.concatenate([[0], np.cumsum(weight)])
    centres = np.arange(z.min() + width_nm / 2.0, z.max() - width_nm / 2.0 + step_nm, step_nm)
    if not len(centres):
        centres = np.array([0.5 * (z.min() + z.max())])
    low = np.searchsorted(z, centres - width_nm / 2.0, side="left")
    high = np.searchsorted(z, centres + width_nm / 2.0, side="right")
    scores = cumulative[high] - cumulative[low]
    best = int(np.argmax(scores))
    inside = weight[low[best]:high[best]]
    return {"z_a": round(float(centres[best]) * 10.0, 2), "width_nm": width_nm, "score": int(scores[best]),
            "hydrophobic": int((inside == 1).sum()), "charged": int((inside == -1).sum()),
            "residues": int(high[best] - low[best]), "exposed_residues": int(exposed.sum()),
            "exposure": {"radius_nm": exposure_radius_nm, "buried_above": buried_above}}


def periodic_mean(values: np.ndarray, length: float) -> float:
    """Mean of coordinates on a periodic axis of the given length, inside [0, length)."""
    angle = 2.0 * math.pi * values / length
    return (math.atan2(np.sin(angle).mean(), np.cos(angle).mean()) / (2.0 * math.pi) * length) % length


def phosphate_planes(membrane: list[dict], midplane_nm: float) -> dict:
    """Median z (nm) of the PO4 beads of each leaflet; 2 nm from the midplane for a leaflet without phosphate."""
    z = np.array([p[2] for mol in membrane for bead, p in zip(mol["beads"], mol["xyz"]) if bead["atom"] == "PO4"])
    lower, upper = z[z < midplane_nm], z[z >= midplane_nm]
    return {"lower": float(np.median(lower)) if len(lower) else midplane_nm - 2.0,
            "upper": float(np.median(upper)) if len(upper) else midplane_nm + 2.0}


def membrane_voids(membrane: list[dict], obstacles_nm: np.ndarray, box: list[float], midplane_nm: float,
                   grid_nm: float = VOID_GRID_NM, empty_nm: float = VOID_EMPTY_NM) -> dict:
    """Empty pockets in the acyl region of each leaflet: grid points with no lipid bead and no solute within empty_nm.

    The acyl region of a leaflet runs from VOID_MARGIN_NM off the midplane to VOID_MARGIN_NM inside the median plane
    of its phosphate (PO4) beads. Empty points that touch on the grid form one pocket; the largest pocket volume
    (nm^3) and the number of empty points per leaflet are returned. In the uncut bundled frames (equilibrated, with
    their own protein as the obstacle) the largest pocket is at most 0.2 nm^3, so a larger one is a hole the build made.
    """
    cell = np.asarray(box[:3], float)
    beads = np.vstack([np.asarray(mol["xyz"], float) for mol in membrane])
    planes = phosphate_planes(membrane, midplane_nm)
    solids = np.vstack([beads, obstacles_nm]) if len(obstacles_nm) else beads
    tree = cKDTree(wrap(solids, cell), boxsize=cell)
    nx, ny = max(1, int(round(cell[0] / grid_nm))), max(1, int(round(cell[1] / grid_nm)))
    gx, gy = (np.arange(nx) + 0.5) * cell[0] / nx, (np.arange(ny) + 0.5) * cell[1] / ny
    out = {"grid_nm": grid_nm, "empty_nm": empty_nm}
    for leaflet, sign in (("lower", -1.0), ("upper", 1.0)):
        z0, z1 = sorted((midplane_nm + sign * VOID_MARGIN_NM, planes[leaflet] - sign * VOID_MARGIN_NM))
        gz = np.arange(z0, z1 + 1e-9, grid_nm) if z1 > z0 else np.array([0.5 * (z0 + z1)])
        X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")
        points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
        distance = tree.query(wrap(points, cell), distance_upper_bound=empty_nm)[0]
        empty = np.isinf(distance).reshape(X.shape)
        pockets, sizes = label_periodic(empty)
        voxel = (cell[0] / nx) * (cell[1] / ny) * grid_nm
        largest = int(max(sizes)) if sizes else 0
        where = None
        if largest:
            k = int(np.argmax(sizes)) + 1
            centre = points[(pockets == k).ravel()]
            where = [round(periodic_mean(centre[:, 0], cell[0]), 2), round(periodic_mean(centre[:, 1], cell[1]), 2),
                     round(float(centre[:, 2].mean()), 2)]
        out[leaflet] = {"acyl_z_nm": [round(z0, 3), round(z1, 3)], "empty_points": int(empty.sum()),
                        "pockets": len(sizes), "largest_pocket_nm3": round(largest * voxel, 4), "largest_pocket_centre_nm": where}
    out["largest_pocket_nm3"] = max(out["lower"]["largest_pocket_nm3"], out["upper"]["largest_pocket_nm3"])
    return out


def label_periodic(mask: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """Connected components of a boolean grid (6-neighbour), periodic in the first two axes; labels and sizes."""
    labels = np.zeros(mask.shape, dtype=int)
    sizes, shape = [], mask.shape
    for start in zip(*np.nonzero(mask)):
        if labels[start]:
            continue
        sizes.append(0)
        labels[start] = len(sizes)
        stack = [start]
        while stack:
            i, j, k = stack.pop()
            sizes[-1] += 1
            for di, dj, dk in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                n = ((i + di) % shape[0], (j + dj) % shape[1], k + dk)
                if 0 <= n[2] < shape[2] and mask[n] and not labels[n]:
                    labels[n] = len(sizes)
                    stack.append(n)
    return labels, sizes


def image_onto_slab(xyz: np.ndarray, slab: tuple[float, float], cell_z: float) -> np.ndarray:
    """Move coordinates by whole box lengths in z to the image nearest the bilayer's middle.

    martini.make_membrane_whole brings the lipids into one periodic image in z but leaves the frame's protein where
    the file had it, so a bilayer that crossed the z boundary would otherwise be compared with a protein one box
    length away from it.
    """
    out = np.array(xyz, dtype=float).reshape(-1, 3)
    middle = 0.5 * (slab[0] + slab[1])
    out[:, 2] -= cell_z * np.round((out[:, 2] - middle) / cell_z)
    return out


def frame_protein_beads(cg_protein: list[list[dict]], slab: tuple[float, float], box: list[float],
                        names: tuple | None = None) -> np.ndarray:
    """Coordinates (nm) of the frame's own protein beads (all, or those named) within the bilayer's z range, in its image."""
    beads = np.array([[b["x"], b["y"], b["z"]] for res in cg_protein for b in res if names is None or b["atom"] in names],
                     dtype=float).reshape(-1, 3)
    beads = image_onto_slab(beads, slab, float(box[2]))
    return beads[(beads[:, 2] >= slab[0]) & (beads[:, 2] <= slab[1])]


def free_site(cg_protein: list[list[dict]], slab: tuple[float, float], box: list[float], window: tuple | None = None,
              grid_nm: float = FREE_SITE_GRID_NM) -> tuple[np.ndarray, float | None]:
    """Where to put the complex so that the membrane cut around it stays farthest from the frame's own protein.

    A bundled frame was equilibrated around a receptor; where the receptor was, the frame has a hole of the receptor's
    shape. A complex much smaller than that receptor (a single transmembrane helix) cannot fill it, so it is placed
    on unbroken membrane instead. `window` is the rectangle the slice will keep, as (lower, upper) x/y offsets (nm)
    from the placement point; every grid point of the cell is scored by the distance from the window placed there to
    the nearest protein bead (periodic in x and y; an axis the window spans entirely contributes nothing), and the
    best point and its distance are returned. Without a window the score is the distance from the point itself.
    Without protein in the bilayer the centre of the patch is returned with no distance.
    """
    cell = np.asarray(box[:2], float)
    beads = frame_protein_beads(cg_protein, slab, box)
    if not len(beads):
        return 0.5 * cell, None
    nx, ny = max(1, int(round(cell[0] / grid_nm))), max(1, int(round(cell[1] / grid_nm)))
    X, Y = np.meshgrid((np.arange(nx) + 0.5) * cell[0] / nx, (np.arange(ny) + 0.5) * cell[1] / ny, indexing="ij")
    points = np.column_stack([X.ravel(), Y.ravel()])
    lower, upper = (np.zeros(2), np.zeros(2)) if window is None else (np.asarray(window[0], float), np.asarray(window[1], float))
    clearance = rectangle_distance(beads[:, :2], points + lower, upper - lower, cell).min(axis=1)
    best = int(np.argmax(clearance))
    return points[best], float(clearance[best])


def rectangle_distance(points_xy: np.ndarray, lowers: np.ndarray, size: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Periodic Euclidean distance (nm) from each point to each axis-aligned rectangle [lower, lower + size); 0 inside.

    lowers holds one lower corner per rectangle (shape (R, 2)); the result has shape (R, number of points). An axis
    the rectangle spans entirely contributes no gap. free_site and frame_protein_in_slice share it, so the placement
    check and the check of the cut that follows measure the same distance.
    """
    lowers, size, cell = np.atleast_2d(lowers).astype(float), np.asarray(size, float), np.asarray(cell, float)
    gaps = []
    for axis in range(2):
        if size[axis] >= cell[axis]:
            gaps.append(np.zeros((len(lowers), len(points_xy))))
            continue
        relative = np.mod(points_xy[None, :, axis] - lowers[:, axis, None], cell[axis])
        gaps.append(np.where(relative <= size[axis], 0.0, np.minimum(relative - size[axis], cell[axis] - relative)))
    return np.hypot(gaps[0], gaps[1])


def slice_window(xyz_nm: np.ndarray, anchor_xy: np.ndarray, cell: np.ndarray, buffer_nm: float, slack_nm: float,
                 requested_xy_nm: tuple | None = None) -> tuple[np.ndarray, np.ndarray]:
    """The rectangle a slice will keep around a complex, as (lower, upper) x/y offsets (nm) from anchor_xy.

    Mirrors slicing.crop_windows: the complex's x/y extent plus buffer_nm on each side, or the requested width centred
    on the extent, widened by slack_nm on each side for the crop offset the slicer may choose. An axis at least as
    wide as the cell keeps the whole cell.
    """
    low, high = xyz_nm[:, :2].min(axis=0), xyz_nm[:, :2].max(axis=0)
    centre = 0.5 * (low + high)
    width = (high - low) + 2.0 * buffer_nm if requested_xy_nm is None else np.asarray(requested_xy_nm, float)
    width = np.where(width < cell - 0.01, width + 2.0 * slack_nm, np.maximum(width, cell))
    return centre - 0.5 * width - anchor_xy, centre + 0.5 * width - anchor_xy


def frame_protein_in_slice(cg_protein: list[list[dict]], slab: tuple[float, float], box: list[float], placed: list[dict],
                           cut: dict, clearance_nm: float = FRAME_PROTEIN_CLEARANCE_NM) -> int:
    """How many beads of the frame's own protein (in the bilayer) lie inside, or within clearance_nm of, the slice window.

    The slicer moves the complex by the window's lower corner; the same shift maps the frame's protein into the new
    cell. The receptor's hole is wider than its bead centres (a bead radius plus the lipid-free gap around it), so a
    bead within clearance_nm of the window counts too. An axis that was not cut keeps the whole cell, so any bead
    counts there.
    """
    beads = frame_protein_beads(cg_protein, slab, box)
    if not len(beads):
        return 0
    cell, size = np.asarray(box[:2], float), np.asarray(cut["box"][:2], float)
    size = np.where(size < cell - 1e-9, size, cell)  # an uncut axis keeps the whole cell
    lower = np.array([placed[0]["x"] - cut["placed"][0]["x"], placed[0]["y"] - cut["placed"][0]["y"]]) / 10.0
    distance = rectangle_distance(beads[:, :2], lower, size, cell)[0]
    return int(((distance < clearance_nm) | (distance == 0.0)).sum())


def embed_complex(aa_atoms: list[dict], cg_protein: list[list[dict]], membrane: list[dict], box: list[float],
                  slab: tuple[float, float], bilayer_z_a: float | None = None, midplane_nm: float | None = None,
                  overlap_nm: float = OVERLAP_NM, hard_core_nm: float = HARD_CORE_NM, relax_steps: int = RELAX_STEPS,
                  max_void_nm3: float = MAX_VOID_NM3, site: str = "hole", window_buffer_nm: float = 1.0,
                  window_slack_nm: float = 0.5, requested_xy_nm: tuple | None = None) -> dict:
    """Put the complex into the frame's membrane with its hydrophobic belt on the midplane, and make room for it.

    site "hole" (default) puts the complex where the frame's own protein was; site "free" puts it on the stretch of
    membrane farthest from that protein (free_site), where the frame has no hole, for a complex much smaller than
    the frame's protein. With "free" the frame's protein still holds its own place: the site is chosen so that the
    slice window (the complex's extent plus window_buffer_nm, or requested_xy_nm, plus window_slack_nm for the crop
    offset; slice_window) stays farthest from it, and must keep FRAME_PROTEIN_CLEARANCE_NM from it; its beads count as
    solid in the pocket search; and the pipeline refuses a slice whose window still comes that close
    (frame_protein_in_slice).

    Room is made by moving lipids, not by deleting every lipid the complex touches: the lipids closer than overlap_nm
    to a heavy atom are pushed out of the way at the coarse-grained level (slicing.relax_seam: soft repulsion from the
    heavy atoms until no bead is closer than overlap_nm, weak position restraints, intramolecular shape restraints,
    and repulsion between lipid beads that the push brings together). Only a lipid still closer than hard_core_nm to a
    heavy atom afterwards is removed: the complex occupies its place in the leaflet. A lipid anchor of the complex (a
    prenyl chain, a palmitoylated N-terminus) therefore displaces the lipids around it instead of emptying a column of
    the leaflet. The acyl region is then searched for empty pockets (membrane_voids); a pocket larger than
    max_void_nm3 stops the build.
    """
    # Returns the identity rotation and the translation (A) applied to the complex, the placed atoms, the membrane
    # molecules that remain (at the positions the push gave them), the indices of the removed ones in `membrane`, and
    # notes/metrics for the log and manifest. `membrane` itself is not modified.
    cell = np.array(box)
    # midplane_nm: the bilayer midplane to put the bilayer centre on (default: the middle of the membrane's z extent).
    midplane = 0.5 * (slab[0] + slab[1]) if midplane_nm is None else float(midplane_nm)
    belt = hydrophobic_belt(aa_atoms) if bilayer_z_a is None else {"z_a": float(bilayer_z_a), "given": True}
    centre_z = belt["z_a"] / 10.0
    inside = frame_protein_beads(cg_protein, slab, box, names=("BB",))
    if site not in EMBED_SITES:
        raise SystemExit(f"unknown embedding site {site}; choose from {', '.join(EMBED_SITES)}")
    away = clear = None
    if site == "free":
        target_xy, where = None, ""  # chosen below, once the complex's transmembrane centre and extent are known
    elif len(inside):
        target_xy = np.array([periodic_mean(inside[:, 0], cell[0]), periodic_mean(inside[:, 1], cell[1])])
        where = f"the centre of the frame's own transmembrane protein ({len(inside)} BB beads in the bilayer)"
    else:
        target_xy = 0.5 * cell[:2]
        where = "the centre of the membrane patch (the frame has no protein in the bilayer)"
    heavy = np.array([element(a["atom"]) != "H" for a in aa_atoms])
    xyz = xyz_nm(aa_atoms) / 10.0
    tm = heavy & (np.abs(xyz[:, 2] - centre_z) <= BELT_NM / 2.0)
    if tm.sum() < 10:
        raise SystemExit(f"fewer than 10 heavy atoms within {BELT_NM / 2.0:.1f} nm of the bilayer centre z = "
                         f"{belt['z_a']:.1f} A; orient the input with the membrane normal along z or pass --bilayer-z")
    if site == "free":
        tm_xy = xyz[tm, :2].mean(axis=0)
        window = slice_window(xyz, tm_xy, cell[:2], window_buffer_nm, window_slack_nm, requested_xy_nm)
        target_xy, clear = free_site(cg_protein, slab, box, window)
        if clear is not None:
            offsets = frame_protein_beads(cg_protein, slab, box)[:, :2] - target_xy
            away = float(np.linalg.norm(offsets - cell[:2] * np.round(offsets / cell[:2]), axis=1).min())
            if clear < FRAME_PROTEIN_CLEARANCE_NM:
                size = window[1] - window[0]
                raise SystemExit(f"no place in this {cell[0]:.1f} x {cell[1]:.1f} nm frame keeps the membrane cut around this complex "
                                 f"({min(size[0], cell[0]):.1f} x {min(size[1], cell[1]):.1f} nm with the crop slack) "
                                 f"{FRAME_PROTEIN_CLEARANCE_NM} nm from the frame's own protein (best {clear:.2f} nm): use a smaller "
                                 "--xy-buffer or --box, or --embed-site hole")
        where = ("the membrane farthest from the frame's own protein for the cut around the complex "
                 + (f"(cut edge {clear:.2f} nm, transmembrane centre {away:.2f} nm from its nearest bead in the bilayer)"
                    if clear is not None else "(the frame has no protein in the bilayer)"))
    shift_nm = np.array([target_xy[0] - xyz[tm, 0].mean(), target_xy[1] - xyz[tm, 1].mean(), midplane - centre_z])
    t = shift_nm * 10.0
    placed = [dict(a, x=a["x"] + t[0], y=a["y"] + t[1], z=a["z"] + t[2]) for a in aa_atoms]
    protein_nm = xyz[heavy] + shift_nm
    whole = [np.asarray(mol["xyz"], float) if "xyz" in mol else xyz_nm(mol["beads"]) for mol in membrane]
    tree = cKDTree(wrap(protein_nm, box), boxsize=box)
    closest_start = np.array([tree.query(wrap(x, box))[0].min() for x in whole])
    relaxed = relax_seam(whole, cell[:2], float(cell[2]), [False, False], protein_nm, SEAM_CONTACT_NM, relax_steps,
                         protein_repulsion_nm=overlap_nm + PUSH_MARGIN_NM, protein_clearance_nm=overlap_nm)
    closest = np.array([tree.query(wrap(x, box))[0].min() for x in relaxed["xyz"]])
    kept, removed, removed_leaflets, removed_indices, displacement, shape = [], Counter(), Counter(), [], [], []
    for index, (mol, new, gap) in enumerate(zip(membrane, relaxed["xyz"], closest)):
        if gap < hard_core_nm:
            removed[mol["cg"]] += 1
            removed_indices.append(index)
            anchor = anchor_xyz(new, [b["atom"] for b in mol["beads"]], anchor_bead(mol["cg"]))
            removed_leaflets[leaflet_of(float(anchor[2]), midplane)] += 1
            continue
        beads = [dict(b, x=float(p[0]), y=float(p[1]), z=float(p[2])) for b, p in zip(mol["beads"], new)]
        kept.append(dict(mol, beads=beads, xyz=new.copy()))
        old = whole[index]
        i, j = np.triu_indices(len(new), 1)
        displacement.append(float(np.linalg.norm(new - old, axis=1).max()))
        shape.append(float(np.abs(np.linalg.norm(new[i] - new[j], axis=1) - np.linalg.norm(old[i] - old[j], axis=1)).max())
                     if len(i) else 0.0)
    if not kept:
        raise SystemExit("every membrane molecule overlaps the placed complex; the input is not oriented along z")
    clearance = float(closest[closest >= hard_core_nm].min())
    frame_protein = frame_protein_beads(cg_protein, slab, box) if site == "free" else np.zeros((0, 3))
    voids = membrane_voids(kept, np.vstack([protein_nm, frame_protein]), box, midplane)
    touching = int((closest_start < overlap_nm).sum())
    notes = [f"bilayer centre of the all-atom input at z = {belt['z_a']:.1f} A "
             + ("(given)" if belt.get("given") else f"(hydrophobic belt: {belt['hydrophobic']} hydrophobic and "
                                                    f"{belt['charged']} charged residues in {belt['width_nm']:.1f} nm)"),
             f"complex moved by ({t[0]:+.1f}, {t[1]:+.1f}, {t[2]:+.1f}) A onto {where}, midplane z = {midplane:.2f} nm; "
             "it is not rotated: orient the input with the membrane normal along z before the build",
             f"{touching} membrane molecule(s) within {overlap_nm:.2f} nm of the complex were pushed aside "
             f"({relaxed['steps']} steps, largest bead displacement of a kept lipid {max(displacement):.2f} nm); "
             f"{sum(removed.values())} still closer than {hard_core_nm:.2f} nm removed: "
             + (", ".join(f"{name} {n}" for name, n in sorted(removed.items())) or "none")
             + f"; closest lipid bead to a heavy atom {clearance:.3f} nm",
             f"largest empty pocket in the acyl region {voids['largest_pocket_nm3']:.3f} nm^3 (limit {max_void_nm3} nm^3)"]
    metrics = {"mode": "embed", "site": site, "distance_to_frame_protein_nm": round(away, 3) if away is not None else None,
               "slice_window_clearance_nm": round(clear, 3) if clear is not None else None,
               "bilayer_centre_input_A": belt["z_a"], "hydrophobic_belt": belt,
               "translation_A": [round(float(v), 3) for v in t], "target_xy_nm": [round(float(v), 3) for v in target_xy],
               "midplane_nm": round(midplane, 3), "overlap_cutoff_nm": overlap_nm, "hard_core_nm": hard_core_nm,
               "lipids_touching_complex": touching, "removed_lipids": dict(removed),
               "removed_by_leaflet": dict(removed_leaflets),
               # displacement and shape change of the lipids that are kept
               "relaxation": {"steps": relaxed["steps"], "converged": relaxed["converged"],
                              "lipids_moved_over_0.1_nm": int(sum(d > 0.1 for d in displacement)),
                              "max_bead_displacement_nm": round(max(displacement), 4),
                              "max_intramolecular_distance_change_nm": round(max(shape), 4),
                              "closest_created_pair_nm": relaxed["closest_created_pair_nm"]},
               "closest_bead_to_complex_nm": round(clearance, 4), "voids": voids, "membrane_molecules_kept": len(kept)}
    if voids["largest_pocket_nm3"] > max_void_nm3:
        side = max(("lower", "upper"), key=lambda k: voids[k]["largest_pocket_nm3"])
        raise SystemExit(f"embedding left an empty pocket of {voids['largest_pocket_nm3']:.2f} nm^3 in the acyl region of the "
                         f"{side} leaflet at {voids[side]['largest_pocket_centre_nm']} nm (limit {max_void_nm3} nm^3): the "
                         "frame's lipids cannot be moved around this complex without a hole; "
                         + ("for a complex much smaller than the frame's own protein (e.g. one transmembrane helix) use "
                            "--embed-site free, which places it away from the frame's protein hole; otherwise "
                            if site == "hole" else "")
                         + "use a coarse-grained frame of this complex (--cg FILE) or another frame")
    return {"R": np.eye(3), "t": t, "placed": placed, "membrane": kept, "removed": removed, "removed_indices": removed_indices,
            "notes": notes, "metrics": metrics}


def trim_membrane(membrane: list[dict], placed: list[dict], box: list[float], xy_nm: tuple[float, float]) -> dict:
    """Cut the membrane patch down to xy_nm around the placed complex and move everything into the new cell."""
    # Molecules whose centre falls outside the new cell are dropped; the complex ends up centred in x and y.
    (x, y), cell = xy_nm, np.array(box)
    if x > cell[0] + 1e-6 or y > cell[1] + 1e-6:
        raise SystemExit(f"requested box {x * 10:.0f} x {y * 10:.0f} A exceeds the membrane patch "
                         f"{cell[0] * 10:.0f} x {cell[1] * 10:.0f} A; the box can only be trimmed, not enlarged")
    if x <= 0 or y <= 0:
        raise SystemExit("box edges must be positive")
    heavy = [a for a in placed if element(a["atom"]) != "H"]
    centre = xyz_nm(heavy).mean(axis=0)[:2] / 10.0
    shift = np.array([0.5 * x - centre[0], 0.5 * y - centre[1], 0.0])
    kept, removed = [], Counter()
    for mol in membrane:
        xyz = np.array([[b["x"], b["y"], b["z"]] for b in mol["beads"]]) + shift
        middle = xyz.mean(axis=0)
        image = np.array([-cell[0] * math.floor(middle[0] / cell[0]), -cell[1] * math.floor(middle[1] / cell[1]), 0.0])
        xyz, middle = xyz + image, middle + image
        if middle[0] < x and middle[1] < y:
            for bead, (bx, by, bz) in zip(mol["beads"], xyz):
                bead["x"], bead["y"], bead["z"] = float(bx), float(by), float(bz)
            kept.append(mol)
        else:
            removed[mol["cg"]] += 1
    if not kept:
        raise SystemExit("no membrane molecule is left inside the requested box")
    moved = [dict(a, x=a["x"] + shift[0] * 10.0, y=a["y"] + shift[1] * 10.0) for a in placed]
    note = (f"membrane trimmed to {x:.2f} x {y:.2f} nm around the complex: {sum(removed.values())} molecule(s) outside "
            f"removed ({', '.join(f'{n} {k}' for n, k in sorted(removed.items())) or 'none'}); lipids across the new "
            "periodic seam were not equilibrated together, so equilibrate the system before production")
    return {"membrane": kept, "placed": moved, "box": [x, y, float(cell[2])], "removed": removed, "note": note,
            "shift_A": [round(float(v) * 10.0, 3) for v in shift]}


def edit_lipids(membrane: list[dict], delete: list[str], add: str | None, mapping: dict) -> tuple[list[dict], dict]:
    """Remove every molecule of the named lipids or, with `add`, turn them into that lipid by renaming their beads."""
    names = {lipid_name(n) for n in delete}
    target = lipid_name(add) if add else None
    if target and not names:
        raise SystemExit("--addlipid replaces the lipids removed by --dellipid; say which lipid to replace with --dellipid")
    if target and (not {*names, target} <= set(CONVERTIBLE)):
        raise SystemExit(f"cannot turn {', '.join(sorted(names))} into {target}: only {', '.join(CONVERTIBLE)} share a "
                         "Martini 3 bead layout; cholesterol, DPG3 and SAP6 can be removed but not replaced")
    if target and target in names:
        raise SystemExit(f"--addlipid {target} is also named by --dellipid")
    if target and target not in mapping:
        raise SystemExit(f"the installed mapping has no residue {target}")
    kept, removed, converted = [], Counter(), Counter()
    for mol in membrane:
        if mol["cg"] not in names:
            kept.append(mol)
        elif target:
            position = {name: i for i, name in enumerate(mapping[mol["aa"]]["beads"])}
            beads = [dict(b, atom=mapping[target]["beads"][position[b["atom"]]], resname=target) for b in mol["beads"]]
            kept.append(dict(mol, cg=target, aa=target, beads=beads, dropped=[]))
            converted[mol["cg"]] += 1
        else:
            removed[mol["cg"]] += 1
    if not kept:
        raise SystemExit("every membrane molecule was removed by --dellipid")
    missing = sorted(n for n in names if not removed[n] and not converted[n])
    if missing:
        raise SystemExit(f"the membrane has no {', '.join(missing)} to remove")
    return kept, {"deleted": dict(removed), "converted_to_" + target if target else "converted": dict(converted)}


def check_box_z(placed: list[dict], slab: tuple[float, float], requested_z_nm: float) -> float:
    """Refuse a requested box z that cannot hold the placed complex, before the slow backmapping is paid for."""
    # Same rule as rebox_system (water padding above and below the complex, measured from the bilayer centre),
    # estimated from the coarse-grained bilayer; rebox_system re-checks on the all-atom membrane. Returns the need (nm).
    z = np.array([a["z"] for a in placed if element(a["atom"]) != "H"]) / 10.0
    centre = 0.5 * (slab[0] + slab[1])
    need = max(2.0 * (max(z.max() - centre, centre - z.min()) + SLAB_Z_PAD_NM),
               max(z.max(), slab[1]) - min(z.min(), slab[0]) + MIN_Z_PAD_TOTAL_NM)
    if requested_z_nm < need - 1e-6:
        raise SystemExit(f"requested box z {requested_z_nm * 10:.0f} A is below the {need * 10:.0f} A this complex needs "
                         f"({SLAB_Z_PAD_NM} nm of water above and below it); checked before backmapping")
    return need
