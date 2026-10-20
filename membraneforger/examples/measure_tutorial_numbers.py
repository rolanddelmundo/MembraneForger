"""Recompute the frame-dependent numbers of the MembraneForger tutorial on the bundled membranes.

    python membraneforger/examples/measure_tutorial_numbers.py --repo . --out tutorial50_out [--gmx env/bin/gmx]
        [--frames KOR_Membrane1 GPR139_Membrane5]

Only the coarse-grained stages run (embedding, slicing, the slice gate): nothing is backmapped, so a frame takes
minutes. For every bundled frame it records

  * `frame`        the frame as read: rotation found about z, overlapping pairs, largest empty pocket of the uncut frame;
  * `6whc`         the oriented 6WHC complex embedded in the frame (the tutorial's KOR_Membrane1 table), with the slice window
                   centred on the complex (the default);
  * `6whc_search`  the same with --slice-offset-search (the tutorial's offset-search comparison);
  * `helix_free`   TM6 of the GCGR receptor alone on unbroken membrane (--embed-site free, --xy-buffer 3);
  * `helix_hole`   the same helix in the frame's receptor hole, which the empty-pocket gate is expected to stop.

Every case keeps its folder under <out>/cases/ (log, membrane_validation.md/json); a finished case is skipped when the
script runs again, so a job that ran out of time continues where it stopped. The numbers go to <out>/summary.md and
<out>/summary.json.
"""
import argparse
import hashlib
import json
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np

CASES = ("6whc", "6whc_search", "helix_free", "helix_hole")


def jsonable(value):
    """Make a record JSON-safe (numpy numbers, paths, tuples)."""
    return json.loads(json.dumps(value, default=lambda o: o.item() if hasattr(o, "item") else str(o)))


def get(tree, *path):
    """tree[path[0]][path[1]]... or None when any step is missing."""
    for key in path:
        if not isinstance(tree, dict) or key not in tree:
            return None
        tree = tree[key]
    return tree


def frame_numbers(mf, frame, mapping):
    """The uncut frame: alignment, overlapping pairs and the largest empty pocket with its own receptor as obstacle."""
    raw, box = mf.read_cg(frame)
    raw, alignment = mf.align_frame_to_box(raw, box)
    out = {"box_nm": [round(float(b), 4) for b in box], "beads": len(raw),
           "rotation_about_z_deg": alignment.get("rotation_about_z_deg"),
           "overlapping_pairs_as_read": alignment.get("overlapping_pairs_as_read")}
    try:
        cg = mf.classify_cg(raw, mapping)
        slab = mf.make_membrane_whole(cg["membrane"], box)
        midplane = mf.bilayer_midplane(cg["membrane"], slab)[0]
        obstacles = mf.frame_protein_beads(cg["protein"], slab, box)
        voids = mf.membrane_voids(cg["membrane"], obstacles, box, midplane)
        out["largest_pocket_nm3"] = voids["largest_pocket_nm3"]
        out["lipids"] = len(cg["membrane"])
    except Exception:
        out["pocket_error"] = traceback.format_exc(limit=3)
    return out


def summarise(session, prepared):
    """The numbers the tutorial quotes, from the session record of one prepare_inputs call."""
    rec = session.record
    emb, cut, chk = rec.get("embedding") or {}, rec.get("slice") or {}, rec.get("slice_check") or {}
    mv = rec.get("membrane_validation") or {}
    leaflets = {name: {k: get(chk, "leaflets", name, k)
                       for k in ("lipids", "apl_A2", "region_reference_apl_A2", "construction_delta_percent",
                                 "representativeness_delta_percent", "construction_pass", "composition_status")}
                for name in ("upper", "lower")}
    rdf = {name: get(chk, "rdf", "leaflets", name, "all") for name in ("upper", "lower")}
    return {
        "embedding": {"site": emb.get("site"), "distance_to_frame_protein_nm": emb.get("distance_to_frame_protein_nm"),
                      "lipids_touching_complex": emb.get("lipids_touching_complex"),
                      "removed_lipids": emb.get("removed_lipids"), "removed_total": sum((emb.get("removed_lipids") or {}).values()),
                      "largest_pocket_nm3": get(emb, "voids", "largest_pocket_nm3"),
                      "membrane_molecules_kept": emb.get("membrane_molecules_kept"),
                      "relaxation": emb.get("relaxation")},
        "slice": {k: cut.get(k) for k in ("mode", "cropped", "new_cell_nm", "area_fraction_kept", "lipids_before", "lipids_after",
                                          "leaflets_after", "membrane_margin_around_complex_nm", "frame_protein_beads_in_window")},
        "crop_offset": cut.get("crop_offset"),
        "selection": cut.get("selection"),
        "seam": cut.get("seam"),
        "slice_check": {"pass": chk.get("pass"), "status": chk.get("status"), "leaflets": leaflets,
                        "composition": chk.get("composition"), "rdf": rdf},
        "old_rule_estimate": mv.get("old_rule_estimate"),
        "stages": [{"label": s.get("label"), "upper": get(s, "leaflets", "upper", "lipids"), "lower": get(s, "leaflets", "lower", "lipids"),
                    "upper_apl_A2": get(s, "leaflets", "upper", "apl_A2"), "lower_apl_A2": get(s, "leaflets", "lower", "apl_A2")}
                   for s in mv.get("stages", [])],
    }


def run_case(mf, args, frame, case, folder):
    """Run prepare_inputs for one frame and case; returns the summary or the failure."""
    repo = Path(args.repo).resolve()
    forcefield = repo / "membraneforger/forcefield"
    aa_pdb = repo / "membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_prot-lig.pdb"
    opm = repo / "membraneforger/tests/data/6whc_opm_chainR.pdb"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "work").mkdir(exist_ok=True)
    (folder / mf.LOG_NAME).write_text("")
    kwargs, orientation, aa = {"embed": True}, mf.OrientationRequest(mode="none"), aa_pdb
    if case.startswith("6whc"):
        orientation = mf.OrientationRequest(mode="opm", chains=("R",), pdb_id="6WHC", opm_file=opm)
        if case == "6whc_search":
            kwargs["settings"] = mf.Settings(slice_optimize_offset=True)
    else:
        atoms = mf.read_all_atom(aa_pdb, forcefield)[0]
        aa = folder / "tm6.pdb"
        mf.write_pdb([a for a in atoms if a["chain"] == "R" and 349 <= a["resid"] <= 376], aa)
        kwargs["settings"] = mf.Settings(box_xy_buffer_nm=3.0)       # one helix plus 1 nm is a ~4 nm cell: the cut edges alone fail the gate
        if case == "helix_free":
            kwargs["embed_site"] = "free"
    session = mf.Session(out=folder, name=f"{frame.stem}_{case}", gmx=args.gmx, forcefield=forcefield, data=mf.locate_data(None),
                         python=sys.executable, ntomp=1, nsteps=1, orientation=orientation, work=folder / "work", **kwargs)
    started = time.time()
    try:
        prepared = mf.prepare_inputs(session, aa, frame)
        result = {"ok": True, **summarise(session, prepared)}
    except BaseException as exc:  # a gate that stops the build (StageFailure, SystemExit) is a result, not a crash
        if isinstance(exc, KeyboardInterrupt):
            raise
        result = {"ok": False, "stopped_by": type(exc).__name__, "stage": getattr(exc, "stage", None),
                  "message": getattr(exc, "what", None) or str(exc.code if isinstance(exc, SystemExit) else exc),
                  **summarise(session, {})}
    result["seconds"] = round(time.time() - started, 1)
    (folder / "session_record.json").write_text(json.dumps(jsonable(session.record), indent=1) + "\n")
    if not args.keep_work:
        import shutil
        shutil.rmtree(folder / "work", ignore_errors=True)
    return result


def fingerprint(repo, frame):
    """What a cached result depends on: the frame's bytes and the repository's commit (and uncommitted changes)."""
    def git(*a):
        return subprocess.run(["git", "-C", str(repo), *a], capture_output=True, text=True).stdout.strip()

    return {"frame_md5": hashlib.md5(Path(frame).read_bytes()).hexdigest(), "commit": git("rev-parse", "HEAD"),
            "dirty": bool(git("status", "--porcelain", "--untracked-files=no"))}


def load_cached(saved, stamp):
    """A cached result for the same fingerprint, else None (a changed frame or commit is measured again)."""
    if not saved.exists():
        return None
    cached = json.loads(saved.read_text())
    return cached if cached.get("fingerprint") == stamp and not stamp["dirty"] else None


def percentile(values, q):
    """q-th percentile of a list, or None when it is empty."""
    return round(float(np.percentile(values, q)), 2) if values else None


def span(values):
    """[minimum, maximum] of a list, or None when it is empty."""
    return [min(values), max(values)] if values else None


def aggregate(results):
    """Across-frame statistics for the tutorial's sentences."""
    agg = {}
    for case in ("6whc", "6whc_search"):
        rows = [r[case] for r in results.values() if r.get(case, {}).get("ok")]
        construction = [abs(r["slice_check"]["leaflets"][l]["construction_delta_percent"]) for r in rows for l in ("upper", "lower")
                        if r["slice_check"]["leaflets"][l]["construction_delta_percent"] is not None]
        whole = {l: [r["slice_check"]["leaflets"][l]["representativeness_delta_percent"] for r in rows
                     if r["slice_check"]["leaflets"][l]["representativeness_delta_percent"] is not None] for l in ("upper", "lower")}
        agg[case] = {
            "frames": len(rows), "passed_slice_gate": sum(bool(r["slice_check"]["pass"]) for r in rows),
            "lipids_removed_by_embedding": span([r["embedding"]["removed_total"] for r in rows]),
            "embedded_largest_pocket_nm3": span([r["embedding"]["largest_pocket_nm3"] for r in rows]),
            "construction_abs_percent": {"median": percentile(construction, 50), "p95": percentile(construction, 95),
                                         "max": percentile(construction, 100), "leaflets": len(construction)},
            "whole_cell_percent_sd": {l: round(float(np.std(v)), 1) if v else None for l, v in whole.items()},
            "whole_cell_percent_worst_abs": round(max((abs(x) for v in whole.values() for x in v), default=float("nan")), 1),
            "seam_lipids_removed": [get(r, "seam", "lipids_removed") for r in rows],
        }
    pockets = [r["frame"]["largest_pocket_nm3"] for r in results.values() if r.get("frame", {}).get("largest_pocket_nm3") is not None]
    agg["uncut_frames_largest_pocket_nm3"] = [min(pockets), max(pockets)] if pockets else None
    agg["frame_rotation_deg"] = {n: r["frame"]["rotation_about_z_deg"] for n, r in results.items() if "frame" in r}
    return agg


def write_summary(out, results, aggregates):
    """summary.json and a Markdown summary meant to be pasted back."""
    (out / "summary.json").write_text(json.dumps(jsonable({"results": results, "aggregate": aggregates}), indent=1) + "\n")
    lines = ["# Tutorial numbers on the bundled frames", "", "## Across frames", "", "```json", json.dumps(jsonable(aggregates), indent=1), "```", ""]
    lines += ["## Per frame, oriented GCGR embedded (slice window centred)", "",
              "| Frame | Removed by embedding | Embedded pocket (nm3) | Uncut pocket (nm3) | Cell (nm) | Lipids kept | "
              "dAPL upper / lower vs uncut (%) | dAPL upper / lower vs embed (%) | Old rule lipids (up / lo) | Slice gate |", "|" + "---|" * 10]
    for name, r in sorted(results.items()):
        c, f = r.get("6whc", {}), r.get("frame", {})
        if not c.get("ok"):
            lines.append(f"| {name} | stopped: {c.get('message')} | | | | | | | | |")
            continue
        lv, old = c["slice_check"]["leaflets"], c.get("old_rule_estimate") or {}
        lines.append(f"| {name} | {c['embedding']['removed_total']} | {c['embedding']['largest_pocket_nm3']} | {f.get('largest_pocket_nm3')} | "
                     f"{c['slice']['new_cell_nm']} | {c['slice']['lipids_after']} / {c['slice']['lipids_before']} | "
                     f"{lv['upper']['construction_delta_percent']} / {lv['lower']['construction_delta_percent']} | "
                     f"{lv['upper']['representativeness_delta_percent']} / {lv['lower']['representativeness_delta_percent']} | "
                     f"{get(old, 'leaflets', 'upper', 'lipids')} / {get(old, 'leaflets', 'lower', 'lipids')} | {c['slice_check']['status']} |")
    lines += ["", "## Single helix (TM6 of GCGR)", "", "| Frame | free site: distance to frame receptor (nm) | removed lipids | pocket (nm3) | "
              "window clear of receptor | hole site stopped? | hole message |", "|" + "---|" * 7]
    for name, r in sorted(results.items()):
        free, hole = r.get("helix_free", {}), r.get("helix_hole", {})
        lines.append(f"| {name} | {get(free, 'embedding', 'distance_to_frame_protein_nm')} | {get(free, 'embedding', 'removed_total')} | "
                     f"{get(free, 'embedding', 'largest_pocket_nm3')} | {get(free, 'slice', 'membrane_margin_around_complex_nm')} | "
                     f"{'yes' if hole and not hole.get('ok') else 'no'} | {str(hole.get('message', ''))[:160]} |")
    lines += ["", "## Seam and offset search (GCGR embedded)", ""]
    for name, r in sorted(results.items()):
        c = r.get("6whc", {})
        lines.append(f"- {name}: seam {json.dumps(jsonable(c.get('seam')))}; offset {json.dumps(jsonable(c.get('crop_offset')))}")
    (out / "summary.md").write_text("\n".join(lines) + "\n")


def main(argv=None):
    """Run every case for every bundled frame and write the summary."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repo", default=".", help="the MembraneForger repository (default: the current folder)")
    parser.add_argument("--out", required=True, help="output folder (created; finished cases in it are skipped)")
    parser.add_argument("--gmx", default="gmx", help="GROMACS command (not used by these stages; recorded in the session)")
    parser.add_argument("--frames", nargs="*", help="bundled frame names (default: all, see --list-membranes)")
    parser.add_argument("--keep-work", action="store_true", help="keep each case's work folder")
    args = parser.parse_args(argv)
    sys.path.insert(0, str(Path(args.repo).resolve()))
    import membraneforger as mf
    out = Path(args.out).resolve()
    (out / "cases").mkdir(parents=True, exist_ok=True)
    mapping = json.loads((Path(args.repo).resolve() / "membraneforger/tests/data/mapping_residues.json").read_text())
    frames = {n: p for n, p in mf.bundled_frames().items() if not args.frames or n in {f.upper() for f in args.frames}}
    if not frames:
        raise SystemExit("no bundled frame matches --frames")
    results = {}
    for name, path in frames.items():
        results[name] = {}
        done = out / "cases" / f"{name}__frame.json"
        stamp = fingerprint(args.repo, path)
        results[name]["frame"] = load_cached(done, stamp)
        if results[name]["frame"] is None:
            results[name]["frame"] = {**frame_numbers(mf, path, mapping), "fingerprint": stamp}
            done.write_text(json.dumps(jsonable(results[name]["frame"]), indent=1) + "\n")
        for case in CASES:
            saved = out / "cases" / f"{name}__{case}.json"
            results[name][case] = load_cached(saved, stamp)
            if results[name][case] is not None:
                print(f"{name} {case}: done earlier", flush=True)
                continue
            print(f"{name} {case}: running", flush=True)
            results[name][case] = {**run_case(mf, args, path, case, out / "cases" / f"{name}__{case}"), "fingerprint": stamp}
            saved.write_text(json.dumps(jsonable(results[name][case]), indent=1) + "\n")
            print(f"{name} {case}: {'ok' if results[name][case]['ok'] else 'stopped'} in {results[name][case]['seconds']} s", flush=True)
    aggregates = aggregate(results)
    write_summary(out, results, aggregates)
    print(f"wrote {out / 'summary.md'}")
    # only the single helix at the hole site is expected to be stopped by a gate; any other stopped case would
    # leave the cross-frame statistics partial
    failed = [f"{n} {c}" for n, r in results.items() for c in CASES if c != "helix_hole" and not r[c].get("ok")]
    if failed:
        print("ERROR: these cases stopped and are missing from the statistics: " + ", ".join(failed))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
