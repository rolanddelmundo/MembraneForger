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

__all__ = ['align_sequences', 'kabsch', 'rmsd', 'all_atom_chains', 'cg_protein_segments', 'candidate_matches',
           'assign_matches', 'fit_assignment', 'resolve_interchangeable_chains', 'choose_periodic_image',
           'map_all_atom_to_cg']

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
                   settings: Settings) -> dict:
    """Image every assigned CG segment next to the largest one, then find one trimmed rigid transform for all pairs."""
    blocks = [(k, label, np.array([aa[label]["xyz"][i] for i, _ in candidates[(label, k)]["pairs"]]),
               np.array([segments[k]["xyz"][j] for _, j in candidates[(label, k)]["pairs"]]))
              for label, k in assignment]
    anchor = max(blocks, key=lambda b: len(b[2]))
    R, t = kabsch(anchor[2], anchor[3])
    shifts = [cell * np.round(((P @ R.T + t).mean(axis=0) - Q.mean(axis=0)) / cell) for _, _, P, Q in blocks]
    P = np.vstack([b[2] for b in blocks])
    Q = np.vstack([b[3] + shift for b, shift in zip(blocks, shifts)])
    keep = np.ones(len(P), dtype=bool)
    for _ in range(8):
        R, t = kabsch(P[keep], Q[keep])
        error = np.linalg.norm(P @ R.T + t - Q, axis=1)
        trimmed = error <= max(settings.fit_trim_factor * float(np.median(error)), settings.fit_trim_floor_a)
        if trimmed.sum() < 3 or (trimmed == keep).all():
            break
        keep = trimmed
    R, t = kabsch(P[keep], Q[keep])  # the returned transform always belongs to the returned core
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
                       settings: Settings = Settings()) -> dict:
    """Match all-atom chains to CG protein segments by sequence and find the rigid transform that places the complex."""
    # Chain IDs and residue numbers are never compared. Backbone centres (N, CA, C, O) are fitted onto BB beads,
    # which is where Martini 3 places BB. One transform moves the whole complex so its interfaces stay intact.
    cell = np.array(box) * 10.0
    aa = all_atom_chains(aa_atoms)
    segments = cg_protein_segments(cg_protein, cell, settings.cg_chain_break_nm)
    candidates = candidate_matches(aa, segments, settings)
    assigned, covered = assign_matches(candidates)
    orphans = [k for k in range(len(segments)) if k not in {k for _, k in assigned}]
    if orphans:
        raise SystemExit("coarse-grained protein segment(s) with no all-atom counterpart: " + ", ".join(
            f"{segments[k]['ids'][0]}-{segments[k]['ids'][-1]} ({len(segments[k]['ids'])} residues)" for k in orphans))
    fit = lambda assignment: fit_assignment(assignment, aa, segments, candidates, cell, settings)
    assigned, notes = resolve_interchangeable_chains(assigned, candidates, fit, settings)
    ambiguity = "; ".join(notes) or "none: every chain has a unique sequence match"
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
                     f"{len(P)} pairs, identity {match['same'] / len(P):.0%}, RMSD {rmsd(error):.2f} A in the complex fit, "
                     f"{own:.2f} A fitted alone")
    for k, seg in enumerate(segments):
        if len(seg["seq"]) > len(covered[k]):
            notes.append(f"CG segment {k + 1}: {len(seg['seq']) - len(covered[k])} of {len(seg['seq'])} residues have no "
                         "all-atom partner (the all-atom input decides what is built)")
    unmatched = [label for label in aa if label not in {label for label, _ in assigned}]
    for label in unmatched:
        notes.append(f"chain {label} ({len(aa[label]['seq'])} residues) has no coarse-grained counterpart; "
                     "it is carried rigidly with the complex")
    notes += displaced
    notes.append(f"rigid fit: {len(best['error'])} backbone pairs, RMSD {best['rmsd']:.2f} A, core RMSD {best['core']:.2f} A "
                 f"over {int(best['keep'].sum())} pairs; {inside} BB beads inside the bilayer")
    metrics = {"chains": chains, "unmatched_aa_chains": unmatched, "backbone_pairs": len(best["error"]),
               "raw_rmsd_A": round(best["rmsd"], 3), "core_rmsd_A": round(best["core"], 3),
               "core_fraction": round(core_fraction, 4), "rotation": [[round(float(v), 6) for v in row] for row in R],
               "translation_A": [round(float(v), 3) for v in t], "bb_beads_in_bilayer": inside, "ambiguity": ambiguity,
               "limits": {"fit_max_core_rmsd_a": settings.fit_max_core_rmsd_a,
                          "fit_min_core_fraction": settings.fit_min_core_fraction,
                          "fit_min_identity": settings.fit_min_identity,
                          "fit_max_chain_median_a": settings.fit_max_chain_median_a, "fit_min_coverage": settings.fit_min_coverage,
                          "fit_trim_factor": settings.fit_trim_factor, "fit_trim_floor_a": settings.fit_trim_floor_a}}
    return {"R": R, "t": t, "rows": rows, "notes": notes, "assigned": sorted(assigned), "metrics": metrics}
