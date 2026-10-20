"""Compare the existing headgroup g(r) with the accessible-shell-area normalization on the embedded and sliced bundled membranes.

    python membraneforger/examples/validate_rdf_normalization.py --repo . --out rdf_shell_out
        [--frames KOR_Membrane1 KOR_Membrane5] [--centred|--search]

For every frame the oriented GCGR complex is embedded and the membrane sliced (the CG stages only, as in
measure_tutorial_numbers.py). Per leaflet it reports

  * the existing g(r) (area = the Voronoi lipid-accessible area) and the accessible-shell g(r) of the embedded membrane and
    of the slice, and their mean over 1.0-2.5 nm (1 = uncorrelated beyond the first shell);
  * the RMS difference between the slice and the embedded curves with each normalization;
  * the same two curves for random anchors in the geometry of the slice (the same cell, the protein cross-section at headgroup
    depth and the number of anchors): an ideal gas and hard discs of 0.45 nm. A correct normalization gives 1 for the
    ideal gas, so the departure of the existing curve from 1 there is the geometric part of the offset and what remains
    on the real leaflets is lipid structure.

Nothing here changes a gate: the existing g(r) keeps deciding the slice check. Output: <out>/rdf_normalization.json and .md.
"""
import argparse
import json
import sys
import traceback
from pathlib import Path

import numpy as np

BAND_NM = (1.0, 2.5)


def jsonable(value):
    """Make a record JSON-safe (numpy numbers and arrays)."""
    return json.loads(json.dumps(value, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))


def band_mean(r, g, band=BAND_NM):
    """Mean of g over the band."""
    inside = (r >= band[0]) & (r <= band[1])
    return float(np.mean(g[inside]))


def rms_difference(a, b, band=(0.3, 2.0)):
    """RMS difference of two curves on their common r grid within band."""
    n = min(len(a["g"]), len(b["g"]))
    r = a["r_nm"][:n]
    inside = (r >= band[0]) & (r <= band[1])
    return float(np.sqrt(np.mean((a["g"][:n][inside] - b["g"][:n][inside]) ** 2)))


def stages_for(mf, repo, frame, out, search, gmx):
    """Embed GCGR in the frame and slice it; returns the MembraneValidation holding the embed and slice stages."""
    forcefield = repo / "membraneforger/forcefield"
    aa = repo / "membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_prot-lig.pdb"
    opm = repo / "membraneforger/tests/data/6whc_opm_chainR.pdb"
    out.mkdir(parents=True, exist_ok=True)
    (out / "work").mkdir(exist_ok=True)
    (out / mf.LOG_NAME).write_text("")
    session = mf.Session(out=out, name=f"{frame.stem}_rdf", gmx=gmx, forcefield=forcefield, data=mf.locate_data(None), python=sys.executable,
                         ntomp=1, nsteps=1, embed=True, work=out / "work", settings=mf.Settings(slice_optimize_offset=search),
                         orientation=mf.OrientationRequest(mode="opm", chains=("R",), pdb_id="6WHC", opm_file=opm))
    mf.prepare_inputs(session, aa, frame)
    return session.analysis


def analyse(mf, analysis, replicates):
    """The comparison for one frame, from the embed and slice stages of its MembraneValidation."""
    result = {}
    embed, cut = analysis.stages["embed"], analysis.stages["slice"]
    shell = {name: mf.leaflet_rdfs_shell(entry["stage"]) for name, entry in (("embed", embed), ("slice", cut))}
    for leaflet in ("lower", "upper"):
        raw = {name: entry["rdfs"][leaflet]["all"] for name, entry in (("embed", embed), ("slice", cut))}
        sh = {name: shell[name][leaflet]["all"] for name in shell}
        stage = cut["stage"]
        rows = [r for r in stage["lipids"] if r["leaflet"] == leaflet]
        plane = float(np.median([r["z_nm"] for r in rows]))
        protein = np.asarray(stage["protein_nm"], dtype=float).reshape(-1, 3)
        section = protein[np.abs(protein[:, 2] - plane) <= 0.6][:, :2]
        grid = mf.open_area_grid(section, np.array(stage["box_nm"][:2]))
        r_max = float(raw["slice"]["r_nm"][-1] + 0.5 * raw["slice"]["dr_nm"])
        refs = {kind: mf.random_reference_rdfs(grid, len(rows), replicates=replicates, hard_core_nm=core, seed=11, r_max=r_max)
                for kind, core in (("ideal gas", 0.0), ("hard discs 0.45 nm", 0.45))}
        r = raw["slice"]["r_nm"]
        entry = {"lipids_slice": len(rows), "lipids_embed": embed["measure"][leaflet]["lipids"],
                 "accessible_nm2_voronoi_slice": cut["measure"][leaflet]["accessible_area_nm2"],
                 "accessible_nm2_grid_slice": sh["slice"]["accessible_nm2"],
                 "mean_g_1_to_2.5_nm": {name: {"existing": band_mean(raw[name]["r_nm"], raw[name]["g"]),
                                               "shell": band_mean(sh[name]["r_nm"], sh[name]["g"])} for name in raw},
                 "rms_slice_vs_embed_0.3_to_2.0_nm": {"existing": rms_difference(raw["embed"], raw["slice"]),
                                                      "shell": rms_difference(sh["embed"], sh["slice"])},
                 "random_reference_of_the_slice_geometry": {
                     kind: {"existing": band_mean(ref["r_nm"], ref["raw_mean"]), "shell": band_mean(ref["r_nm"], ref["shell_mean"]),
                            "replicates": ref["replicates"], "shell_sd_over_replicates": band_mean(ref["r_nm"], ref["shell_sd"])}
                     for kind, ref in refs.items()},
                 "curves": {"r_nm": r, "embed_existing": raw["embed"]["g"], "slice_existing": raw["slice"]["g"],
                            "embed_shell": sh["embed"]["g"], "slice_shell": sh["slice"]["g"],
                            "ideal_existing": refs["ideal gas"]["raw_mean"], "ideal_shell": refs["ideal gas"]["shell_mean"]}}
        result[leaflet] = entry
    return result


def write_report(out, results):
    """rdf_normalization.json and a Markdown summary."""
    (out / "rdf_normalization.json").write_text(json.dumps(jsonable(results), indent=1) + "\n")
    lines = ["# Accessible-shell normalization of the headgroup g(r)", "",
             "Mean g(r) over 1.0-2.5 nm. 1 = no structure beyond the first shell. 'Ideal gas' = random anchors in the slice's open area "
             "(same cell, protein cross-section and number of anchors): any departure from 1 there is geometry, not lipid structure.", "",
             "| Frame | Leaflet | Lipids | Embed existing | Embed shell | Slice existing | Slice shell | Ideal gas existing | Ideal gas shell | "
             "Hard discs existing | Hard discs shell | RMS slice-embed existing | RMS slice-embed shell |", "|" + "---|" * 13]
    for frame, leaflets in sorted(results.items()):
        if "error" in leaflets:
            lines.append(f"| {frame} | stopped: {leaflets['error']} | | | | | | | | | | | |")
            continue
        for leaflet, e in leaflets.items():
            m, ref, rms = e["mean_g_1_to_2.5_nm"], e["random_reference_of_the_slice_geometry"], e["rms_slice_vs_embed_0.3_to_2.0_nm"]
            lines.append(f"| {frame} | {leaflet} | {e['lipids_slice']} | {m['embed']['existing']:.3f} | {m['embed']['shell']:.3f} | "
                         f"{m['slice']['existing']:.3f} | {m['slice']['shell']:.3f} | {ref['ideal gas']['existing']:.3f} | "
                         f"{ref['ideal gas']['shell']:.3f} | {ref['hard discs 0.45 nm']['existing']:.3f} | "
                         f"{ref['hard discs 0.45 nm']['shell']:.3f} | {rms['existing']:.3f} | {rms['shell']:.3f} |")
    (out / "rdf_normalization.md").write_text("\n".join(lines) + "\n")


def main(argv=None):
    """Run the comparison for the chosen frames and write the report."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repo", default=".", help="the MembraneForger repository (default: the current folder)")
    parser.add_argument("--out", required=True, help="output folder (created)")
    parser.add_argument("--frames", nargs="*", help="bundled frame names (default: KOR_Membrane1)")
    parser.add_argument("--search", action="store_true", help="slice with the offset search instead of centred")
    parser.add_argument("--replicates", type=int, default=100, help="random placements per reference (default: %(default)s)")
    parser.add_argument("--gmx", default="gmx")
    args = parser.parse_args(argv)
    repo = Path(args.repo).resolve()
    sys.path.insert(0, str(repo))
    import membraneforger as mf
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    frames = {n: p for n, p in mf.bundled_frames().items() if n.upper() in {f.upper() for f in (args.frames or ["KOR_Membrane1"])}}
    if not frames:
        raise SystemExit("no bundled frame matches --frames")
    results = {}
    for name, frame in frames.items():
        print(f"{name}: embedding, slicing, comparing", flush=True)
        try:
            analysis = stages_for(mf, repo, frame, out / "cases" / name, args.search, args.gmx)
            results[name] = analyse(mf, analysis, args.replicates)
        except BaseException as exc:  # a gate that stops the build is a result for that frame
            if isinstance(exc, KeyboardInterrupt):
                raise
            results[name] = {"error": str(getattr(exc, "what", None) or (exc.code if isinstance(exc, SystemExit) else exc)),
                             "traceback": traceback.format_exc(limit=3)}
    write_report(out, results)
    print(f"wrote {out / 'rdf_normalization.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
