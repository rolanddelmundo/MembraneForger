#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// All-atom <-> coarse-grained correspondence and rigid placement.
#//=============================================================
"""Match all-atom chains to CG protein segments by sequence and fit the complex onto the CG protein."""
import itertools
import math
from collections import OrderedDict, defaultdict

import numpy as np

from .config import AMINO, ONE_LETTER, Settings
from .martini import MARTINI3_PROTEIN_ALIASES
from .structio import residues_in_order, xyz_nm

__all__ = ['align_sequences', 'kabsch', 'kabsch_about_z', 'rmsd', 'all_atom_chains', 'cg_protein_segments',
           'candidate_matches', 'assign_matches', 'fit_assignment', 'resolve_interchangeable_chains',
           'choose_periodic_image', 'map_all_atom_to_cg']

def align_sequences(a: str, b: str) -> list[tuple[int, int]]:
    """Semi-global alignment (end gaps are free): returns the aligned index pairs of two one-letter sequences."""
    n, m = len(a), len(b)
    score = np.zeros((n + 1, m + 1), dtype=np.int32)
    for i in range(1, n + 1):
        row, prev, ai = score[i], score[i - 1], a[i - 1]
        for j in range(1, m + 1):
            diag = prev[j - 1] + (2 if ai == b[j - 1] and ai != "X" else -1)
            row[j] = max(diag, prev[j] - 2, row[j - 1] - 2)
    i, j = max([(n, j) for j in range(m + 1)] + [(i, m) for i in range(n + 1)], key=lambda p: score[p])
    pairs = []
    while i > 0 and j > 0:
        if score[i, j] == score[i - 1, j - 1] + (2 if a[i - 1] == b[j - 1] and a[i - 1] != "X" else -1):
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif score[i, j] == score[i - 1, j] - 2:
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


def kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Least-squares rotation and translation taking points P onto Q."""
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    R = Vt.T @ np.diag([1.0, 1.0, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
    return R, qc - R @ pc


def kabsch_about_z(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Least-squares rotation about z plus translation taking P onto Q; the z axis (membrane normal) is kept."""
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc)[:, :2].T @ (Q - qc)[:, :2]
    theta = math.atan2(H[0, 1] - H[1, 0], H[0, 0] + H[1, 1])
    c, s = math.cos(theta), math.sin(theta)
    R = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    return R, qc - R @ pc


def rmsd(error: np.ndarray) -> float:
    """Root-mean-square of a vector of per-pair deviations."""
    return float(np.sqrt((error ** 2).mean()))


def all_atom_chains(aa_atoms: list[dict]) -> OrderedDict:
    """Collect each all-atom protein chain's one-letter sequence, residue labels and backbone centres (A)."""
    chains: OrderedDict[str, list[dict]] = OrderedDict()
    for a in aa_atoms:
        if a["resname"] in AMINO:
            chains.setdefault(a["chain"] or a["segid"] or "A", []).append(a)
    described = OrderedDict()
    for label, atoms in chains.items():
        residues = residues_in_order(atoms)
        centres = []
        for rid, rn, members in residues:
            backbone = [x for x in members if x["atom"] in ("N", "CA", "C", "O")]
            if not any(x["atom"] == "CA" for x in backbone):
                raise SystemExit(f"all-atom {label}:{rn}{rid} has no CA atom")
            centres.append(xyz_nm(backbone).mean(axis=0))
        described[label] = {"seq": "".join(ONE_LETTER.get(rn, "X") for _, rn, _ in residues),
                            "ids": [f"{rn}{rid}" for rid, rn, _ in residues], "xyz": np.array(centres)}
    return described


def cg_protein_segments(cg_protein: list[list[dict]], cell: np.ndarray, chain_break_nm: float) -> list[dict]:
    """Split the CG protein into continuous BB traces, unwrapping each across the periodic boundary (A)."""
    segments, last = [], None
    for beads in cg_protein:
        bb = next(b for b in beads if b["atom"] == "BB")
        point = np.array([bb["x"], bb["y"], bb["z"]]) * 10.0
        step = None if last is None else point - last["raw"] - cell * np.round((point - last["raw"]) / cell)
        if (last is None or bb["chain"] != last["chain"] or bb["resid"] <= last["resid"]
                or np.linalg.norm(step) > chain_break_nm * 10.0):
            segments.append({"seq": "", "ids": [], "xyz": []})
            whole = point
        else:
            whole = last["whole"] + step
        base = MARTINI3_PROTEIN_ALIASES.get(bb["resname"], bb["resname"])
        segments[-1]["seq"] += ONE_LETTER[base]
        segments[-1]["ids"].append(f"{bb['resname']}{bb['resid']}")
        segments[-1]["xyz"].append(whole)
        last = {"raw": point, "whole": whole, "chain": bb["chain"], "resid": bb["resid"]}
    for seg in segments:
        seg["xyz"] = np.array(seg["xyz"])
    return segments


def candidate_matches(aa: OrderedDict, segments: list[dict], settings: Settings) -> dict:
    """Align every all-atom chain with every CG segment and keep the pairs that pass the identity criteria."""
    candidates = {}
    for label, chain in aa.items():
        for k, seg in enumerate(segments):
            pairs = align_sequences(chain["seq"], seg["seq"])
            same = sum(chain["seq"][i] == seg["seq"][j] for i, j in pairs)
            if (len(pairs) >= settings.fit_min_pairs and same >= settings.fit_min_identity * len(pairs)
                    and len(pairs) >= settings.fit_min_coverage * min(len(chain["seq"]), len(seg["seq"]))):
                candidates[(label, k)] = {"pairs": pairs, "same": same}
    return candidates


def assign_matches(candidates: dict) -> tuple[list, dict]:
    """Greedily accept the best sequence matches whose all-atom and CG residues do not overlap earlier ones."""
    assigned, covered = [], defaultdict(set)
    for label, k in sorted(candidates, key=lambda key: -candidates[key]["same"]):
        used_aa = {i for i, _ in candidates[(label, k)]["pairs"]}
        used_cg = {j for _, j in candidates[(label, k)]["pairs"]}
        if not covered[label] & used_aa and not covered[k] & used_cg:
            assigned.append((label, k))
            covered[label] |= used_aa
            covered[k] |= used_cg
    return assigned, covered


def fit_assignment(assignment: list, aa: OrderedDict, segments: list[dict], candidates: dict, cell: np.ndarray,
                   settings: Settings, about_z: bool = False) -> dict:
    """Image every assigned CG segment next to the largest one, then find one trimmed rigid transform for all pairs."""
    # about_z restricts the transform to a rotation about z plus a translation: used when the all-atom complex is
    # already membrane-oriented, so that the CG registration can never tilt it away from the membrane normal.
    fit = kabsch_about_z if about_z else kabsch
    blocks = [(k, label, np.array([aa[label]["xyz"][i] for i, _ in candidates[(label, k)]["pairs"]]),
               np.array([segments[k]["xyz"][j] for _, j in candidates[(label, k)]["pairs"]]))
              for label, k in assignment]
    anchor = max(blocks, key=lambda b: len(b[2]))
    R, t = fit(anchor[2], anchor[3])
    shifts = [cell * np.round(((P @ R.T + t).mean(axis=0) - Q.mean(axis=0)) / cell) for _, _, P, Q in blocks]
    P = np.vstack([b[2] for b in blocks])
    Q = np.vstack([b[3] + shift for b, shift in zip(blocks, shifts)])
    keep = np.ones(len(P), dtype=bool)
    for _ in range(8):
        R, t = fit(P[keep], Q[keep])
        error = np.linalg.norm(P @ R.T + t - Q, axis=1)
        trimmed = error <= max(settings.fit_trim_factor * float(np.median(error)), settings.fit_trim_floor_a)
        if trimmed.sum() < 3 or (trimmed == keep).all():
            break
        keep = trimmed
    R, t = fit(P[keep], Q[keep])  # the returned transform always belongs to the returned core
    error = np.linalg.norm(P @ R.T + t - Q, axis=1)
    return {"R": R, "t": t, "error": error, "keep": keep, "Q": Q, "blocks": blocks, "shifts": shifts,
            "rmsd": rmsd(error), "core": rmsd(error[keep])}


def resolve_interchangeable_chains(assigned: list, candidates: dict, fit, settings: Settings) -> tuple[list, list[str]]:
    """Let geometry decide between chains that share a sequence, and refuse a dead heat."""
    rivals = defaultdict(set)
    for label, k in assigned:
        for (other, k2), cand in candidates.items():
            if k2 == k and other != label and cand["same"] >= settings.fit_rival_identity * candidates[(label, k)]["same"]:
                rivals[label].add(other)
                rivals[other].add(label)
    notes, done = [], set()
    for label in list(rivals):
        if label in done:
            continue
        group, stack = set(), [label]
        while stack:
            x = stack.pop()
            if x not in group:
                group.add(x)
                stack += rivals[x]
        done |= group
        slots = sorted(k for owner, k in assigned if owner in group)
        if len(set(slots)) != len(slots):
            raise SystemExit(f"ambiguous all-atom to coarse-grained chain mapping: chains {sorted(group)} share a "
                             "sequence and one coarse-grained segment matches more than one of them")
        if math.perm(len(group), len(slots)) > 720:
            raise SystemExit(f"ambiguous all-atom to coarse-grained chain mapping: {len(group)} interchangeable chains "
                             f"{sorted(group)} are too many to resolve by geometry")
        options = []
        for order in itertools.permutations(sorted(group), len(slots)):
            if all((owner, k) in candidates for owner, k in zip(order, slots)):
                trial = [pair for pair in assigned if pair[0] not in group] + list(zip(order, slots))
                options.append((fit(trial)["rmsd"], trial))
        options.sort(key=lambda o: o[0])
        if len(options) > 1 and options[1][0] - options[0][0] < settings.fit_dead_heat_a:
            raise SystemExit(f"ambiguous all-atom to coarse-grained chain mapping: chains {sorted(group)} are "
                             f"interchangeable by sequence and fit equally well ({options[0][0]:.2f} vs {options[1][0]:.2f} A)")
        assigned = options[0][1]
        notes.append(f"chains {','.join(sorted(group))} share a sequence; assigned by geometry "
                     f"(RMSD {options[0][0]:.2f} A, next best {options[1][0]:.2f} A)" if len(options) > 1 else
                     f"chains {','.join(sorted(group))} share a sequence; only one assignment is possible")
    return assigned, notes


def choose_periodic_image(points: np.ndarray, cell: np.ndarray, slab: tuple[float, float]) -> tuple[np.ndarray, int]:
    """Pick the z image that puts the most BB beads inside the bilayer and keep the complex in the cell in XY."""
    low, high = slab[0] * 10.0, slab[1] * 10.0
    inside = lambda n: int(((points[:, 2] + n * cell[2] >= low) & (points[:, 2] + n * cell[2] <= high)).sum())
    image = max(range(-2, 3), key=lambda n: (inside(n), -abs(points[:, 2].mean() + n * cell[2] - 0.5 * (low + high))))
    centre = points.mean(axis=0)
    shift = np.array([-cell[0] * math.floor(centre[0] / cell[0]), -cell[1] * math.floor(centre[1] / cell[1]),
                      image * cell[2]])
    return shift, inside(image)


def map_all_atom_to_cg(aa_atoms: list[dict], cg_protein: list[list[dict]], box: list[float], slab: tuple[float, float],
                       settings: Settings = Settings(), midplane_nm: float | None = None, anchors: tuple = (),
                       embedded_half_a: float | None = None) -> dict:
    """Match all-atom chains to CG protein segments by sequence and find the rigid transform that places the complex."""
    # Chain IDs and residue numbers are never compared. Backbone centres (N, CA, C, O) are fitted onto BB beads,
    # which is where Martini 3 places BB. One transform moves the whole complex so its interfaces stay intact.
    # With midplane_nm (the CG bilayer midplane) and anchors, the all-atom complex is taken to be membrane-oriented
    # already (normal +z, midplane z = 0): the registration is then a rotation about z plus a translation fitted on
    # the anchor chain(s), with z = 0 put on the CG midplane. The orientation is kept; the CG frame is used laterally.
    cell = np.array(box) * 10.0
    about_z = midplane_nm is not None
    aa = all_atom_chains(aa_atoms)
    segments = cg_protein_segments(cg_protein, cell, settings.cg_chain_break_nm)
    candidates = candidate_matches(aa, segments, settings)
    assigned, covered = assign_matches(candidates)
    orphans = [k for k in range(len(segments)) if k not in {k for _, k in assigned}]
    if orphans:
        raise SystemExit("coarse-grained protein segment(s) with no all-atom counterpart: " + ", ".join(
            f"{segments[k]['ids'][0]}-{segments[k]['ids'][-1]} ({len(segments[k]['ids'])} residues)" for k in orphans))
    fit = lambda assignment, z=False: fit_assignment(assignment, aa, segments, candidates, cell, settings, z)
    assigned, notes = resolve_interchangeable_chains(assigned, candidates, fit, settings)
    ambiguity = "; ".join(notes) or "none: every chain has a unique sequence match"
    # The free (unconstrained) fit of the whole complex establishes that the all-atom complex IS the complex of the
    # CG frame: same assembly, every chain where the CG frame has it. In the plain workflow it is also the placement.
    best = fit(assigned)
    core_fraction = float(best["keep"].mean())
    if core_fraction < settings.fit_min_core_fraction or best["core"] > settings.fit_max_core_rmsd_a:
        raise SystemExit(f"all-atom and coarse-grained protein do not superpose: core RMSD {best['core']:.2f} A over "
                         f"{int(best['keep'].sum())}/{len(best['keep'])} backbone pairs (limits "
                         f"{settings.fit_max_core_rmsd_a} A, {settings.fit_min_core_fraction:.0%})")
    offset, displaced = 0, []
    for k, label, P, _ in best["blocks"]:  # no matched chain may sit somewhere else than its CG counterpart
        median = float(np.median(best["error"][offset:offset + len(P)]))
        offset += len(P)
        if median > settings.fit_max_chain_median_a:
            raise SystemExit(f"all-atom chain {label} does not sit where the coarse-grained frame has it: median backbone "
                             f"deviation {median:.1f} A over {len(P)} pairs in the complex fit (limit "
                             f"{settings.fit_max_chain_median_a} A)")
        if median > settings.fit_max_core_rmsd_a:
            displaced.append(f"chain {label} deviates from its CG counterpart by a median {median:.1f} A in the complex fit; "
                             "it keeps its all-atom pose relative to the rest of the complex")
    registration = None
    if about_z:
        # Orientation mode: the complex is already in the membrane frame (normal +z, midplane z = 0) and that frame
        # is authoritative. The CG frame only supplies the lateral registration, read from the anchor chain(s) with
        # a rotation about z plus an xy translation; the z translation puts z = 0 on the CG bilayer midplane. The
        # tilt and depth by which the equilibrated CG pose differs are measured and reported, never adopted.
        anchor_assignment = [(label, k) for label, k in assigned if label in anchors]
        if not anchor_assignment:
            raise SystemExit(f"anchor chain(s) {list(anchors)} have no coarse-grained counterpart, so the oriented complex cannot "
                             "be registered laterally to the CG membrane")
        anchor_free = fit(anchor_assignment)
        # The lateral registration must put the membrane-embedded part of the anchor into the lipid cavity; soluble
        # domains above or below the membrane are a lever arm that only amplifies the (reported) tilt discrepancy.
        # So the constrained fit uses the anchor residues inside the hydrophobic slab of the orientation when there
        # are enough of them, and every anchor residue otherwise (a peripheral anchor).
        embedded = {key: dict(candidates[key], pairs=[(i, j) for i, j in candidates[key]["pairs"]
                                                      if embedded_half_a is None or abs(aa[key[0]]["xyz"][i][2]) <= embedded_half_a])
                    for key in anchor_assignment}
        n_embedded = sum(len(c["pairs"]) for c in embedded.values())
        use_embedded = embedded_half_a is not None and n_embedded >= settings.register_min_embedded_pairs
        reg = fit_assignment(anchor_assignment, aa, segments, embedded if use_embedded else candidates, cell, settings, True)
        tilt = float(np.degrees(np.arccos(np.clip(anchor_free["R"][2, 2], -1.0, 1.0))))
        complex_tilt = float(np.degrees(np.arccos(np.clip(best["R"][2, 2], -1.0, 1.0))))
        midplane_a = midplane_nm * 10.0
        depth_offset_a = float(reg["t"][2] - midplane_a)  # the CG pose puts the anchor this much higher than the midplane does
        registration = {"constrained_to": "rotation about z + translation, fitted on the anchor chain(s) only; z translation "
                                          "from the CG bilayer midplane",
                        "anchor_chains": list(anchors), "anchor_backbone_pairs": int(len(reg["error"])),
                        "fitted_on": (f"anchor residues with backbone centre inside +-{embedded_half_a} A of the midplane"
                                      if use_embedded else "every matched anchor residue"),
                        "embedded_pairs_available": int(n_embedded),
                        "rotation_about_z_deg": round(float(np.degrees(np.arctan2(reg["R"][1, 0], reg["R"][0, 0]))), 3),
                        "cg_midplane_nm": round(midplane_nm, 4),
                        "depth_offset_cg_pose_minus_midplane_A": round(depth_offset_a, 3),
                        "tilt_between_cg_anchor_pose_and_orientation_deg": round(tilt, 3),
                        "tilt_between_cg_complex_pose_and_orientation_deg": round(complex_tilt, 3),
                        "anchor_free_fit_core_rmsd_A": round(anchor_free["core"], 3),
                        "anchor_constrained_core_rmsd_A": round(reg["core"], 3), "anchor_constrained_rmsd_A": round(reg["rmsd"], 3),
                        "anchor_constrained_core_fraction": round(float(reg["keep"].mean()), 4),
                        "complex_free_fit_core_rmsd_A": round(best["core"], 3),
                        "limits": {"register_max_tilt_deg": settings.register_max_tilt_deg,
                                   "register_max_depth_offset_nm": settings.register_max_depth_offset_nm,
                                   "fit_max_core_rmsd_a": settings.fit_max_core_rmsd_a}}
        if tilt > settings.register_max_tilt_deg:
            raise SystemExit(f"the coarse-grained anchor pose is tilted {tilt:.1f} deg from the membrane orientation of the "
                             f"all-atom complex (limit {settings.register_max_tilt_deg} deg); the CG membrane cavity does not "
                             "fit the oriented complex, so no scientifically defensible registration exists")
        if abs(depth_offset_a) > settings.register_max_depth_offset_nm * 10.0:
            raise SystemExit(f"the coarse-grained anchor sits {depth_offset_a / 10:.2f} nm from where the membrane orientation "
                             f"puts it (limit {settings.register_max_depth_offset_nm} nm)")
        if reg["core"] > settings.fit_max_core_rmsd_a or float(reg["keep"].mean()) < settings.fit_min_core_fraction:
            raise SystemExit(f"the oriented anchor cannot be registered laterally to its coarse-grained counterpart: core RMSD "
                             f"{reg['core']:.2f} A over {int(reg['keep'].sum())}/{len(reg['keep'])} {'membrane-embedded ' if use_embedded else ''}"
                             f"backbone pairs with a rotation about z only (limits {settings.fit_max_core_rmsd_a} A, "
                             f"{settings.fit_min_core_fraction:.0%}); the CG pose is tilted {tilt:.1f} deg from the orientation")
        best["R"], best["t"] = reg["R"], np.array([reg["t"][0], reg["t"][1], midplane_a])
        # Deviations reported from here on are those of the actual placement; the core flags keep describing the
        # correspondence (the free complex fit), which is what the rows' "fit" column has always meant.
        best["error"] = np.linalg.norm(np.vstack([b[2] for b in best["blocks"]]) @ best["R"].T + best["t"] - best["Q"], axis=1)
        registration["complex_rmsd_after_registration_A"] = round(rmsd(best["error"]), 3)
    shift, inside = choose_periodic_image(best["Q"], cell, slab)
    R, t = best["R"], best["t"] + shift

    rows, chains, offset = [], [], 0
    for k, label, P, Q in best["blocks"]:
        match = candidates[(label, k)]
        error, keep = best["error"][offset:offset + len(P)], best["keep"][offset:offset + len(P)]
        offset += len(P)
        own_R, own_t = kabsch(P, Q)
        own = rmsd(np.linalg.norm(P @ own_R.T + own_t - Q, axis=1))
        rows += [(label, aa[label]["ids"][i], k + 1, segments[k]["ids"][j], f"{e:.2f}", "core" if c else "trimmed")
                 for (i, j), e, c in zip(match["pairs"], error, keep)]
        chains.append({"aa_chain": label, "cg_segment": k + 1, "aa_residues": len(aa[label]["seq"]),
                       "cg_residues": len(segments[k]["seq"]), "matched_residues": len(P),
                       "identity": round(match["same"] / len(P), 4),
                       "aa_coverage": round(len(P) / len(aa[label]["seq"]), 4),
                       "cg_coverage": round(len(P) / len(segments[k]["seq"]), 4),
                       "rmsd_complex_fit_A": round(rmsd(error), 3), "rmsd_core_pairs_A": round(rmsd(error[keep]), 3) if keep.any() else None,
                       "rmsd_own_fit_A": round(own, 3), "core_pairs": int(keep.sum()),
                       "median_deviation_A": round(float(np.median(error)), 3)})
        notes.append(f"chain {label} ({len(aa[label]['seq'])} residues) -> CG segment {k + 1} "
                     f"{segments[k]['ids'][0]}-{segments[k]['ids'][-1]} ({len(segments[k]['seq'])} residues): "
                     f"{len(P)} pairs, identity {match['same'] / len(P):.0%}, RMSD {rmsd(error):.2f} A "
                     f"{'after registration' if registration else 'in the complex fit'}, {own:.2f} A fitted alone")
    for k, seg in enumerate(segments):
        if len(seg["seq"]) > len(covered[k]):
            notes.append(f"CG segment {k + 1}: {len(seg['seq']) - len(covered[k])} of {len(seg['seq'])} residues have no "
                         "all-atom partner (the all-atom input decides what is built)")
    unmatched = [label for label in aa if label not in {label for label, _ in assigned}]
    for label in unmatched:
        notes.append(f"chain {label} ({len(aa[label]['seq'])} residues) has no coarse-grained counterpart; "
                     "it is carried rigidly with the complex")
    notes += displaced
    notes.append(f"{'correspondence (free complex fit)' if registration else 'rigid fit'}: {len(best['error'])} backbone pairs, "
                 f"RMSD {best['rmsd']:.2f} A, core RMSD {best['core']:.2f} A over {int(best['keep'].sum())} pairs; "
                 f"{inside} BB beads inside the bilayer")
    if registration:
        notes.append(f"registration: rotation about z ({registration['rotation_about_z_deg']} deg) and translation fitted on "
                     f"anchor {','.join(anchors)} ({registration['anchor_backbone_pairs']} pairs, {registration['fitted_on']}, core RMSD "
                     f"{registration['anchor_constrained_core_rmsd_A']} A); the CG anchor pose is tilted "
                     f"{registration['tilt_between_cg_anchor_pose_and_orientation_deg']} deg from the membrane orientation and sits "
                     f"{registration['depth_offset_cg_pose_minus_midplane_A']:+.1f} A from the midplane placement (measured, not "
                     f"adopted); whole-complex RMSD after registration {registration['complex_rmsd_after_registration_A']} A")
    metrics = {"chains": chains, "unmatched_aa_chains": unmatched, "backbone_pairs": len(best["error"]),
               "raw_rmsd_A": round(best["rmsd"], 3), "core_rmsd_A": round(best["core"], 3),
               "core_fraction": round(core_fraction, 4), "rotation": [[round(float(v), 6) for v in row] for row in R],
               "translation_A": [round(float(v), 3) for v in t], "bb_beads_in_bilayer": inside, "ambiguity": ambiguity,
               "limits": {"fit_max_core_rmsd_a": settings.fit_max_core_rmsd_a,
                          "fit_min_core_fraction": settings.fit_min_core_fraction,
                          "fit_min_identity": settings.fit_min_identity,
                          "fit_max_chain_median_a": settings.fit_max_chain_median_a, "fit_min_coverage": settings.fit_min_coverage,
                          "fit_trim_factor": settings.fit_trim_factor, "fit_trim_floor_a": settings.fit_trim_floor_a},
               "registration": registration}
    return {"R": R, "t": t, "rows": rows, "notes": notes, "assigned": sorted(assigned), "metrics": metrics}
