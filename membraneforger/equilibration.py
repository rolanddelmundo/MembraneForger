#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Gates 3 and 4: equilibration convergence and equilibrium physical properties of a finished trajectory.
#//=============================================================
"""Grade an atomistic equilibration / production run of a built system and append the result to its validation report.

    python -m membraneforger.equilibration --out BUILD_DIR --energy ener.xvg [--frames f1.gro f2.gro ...]
        [--production-fraction 0.5] [--temperature 310] [--stress profile.dat] [--reference-ka 250] [--reference-d 8]

The build directory is the one MembraneForger wrote (em.gro, membrane_validation.json). `--energy` is a `gmx energy`
output with Box-X, Box-Y and Box-Z (area and volume series); `--frames` are evenly spaced frames of the trajectory
written as .gro files (`gmx trjconv -sep`, whole molecules), from which the per-frame APL, thickness, core hydration,
protein tilt and depth, tail order parameters and the lipid MSD are computed. Convergence (gate 3) is judged on the
final analysis window of every series; the equilibrium properties (gate 4) on the production fraction only. A
lateral pressure profile (`--stress`, z and P_N - P_L in bar, from a local-stress tool) gives the leaflet tensions.

Everything is reported with qc.metric records; what cannot be measured from the data given is NOT RUN, what has no
matched reference is REFERENCE NEEDED, and what is too short is INSUFFICIENT SAMPLING. Nothing is turned into a PASS
or a FAIL that the data do not support.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

from .config import Settings
from .membrane_report import AA_LIPIDS, CORE_HALF_NM, matrix_markdown, stage_from_gro, stage_matrix
from .packing import measure_packing
from .qc import classify, metric, overall_status, records_table
from .structure_metrics import bilayer_thickness, core_hydration, layered_order, protein_orientation, tail_interdigitation, z_density_profiles
from .trajectory import (
    CHARMM36_CHAINS,
    area_compressibility,
    box_series_from_xvg,
    convergence,
    convergence_record,
    diffusion_fit,
    diffusion_vs_reference,
    ka_vs_reference,
    leaflet_tension,
    msd_xy,
    read_stress_profile,
    scd_profile_from_residues,
    scd_profile_rmse,
    unwrap_xy,
)

__all__ = ['frame_series', 'convergence_gate', 'equilibrium_gate', 'write_equilibration_report', 'main']


def frame_series(paths: list[Path], settings: Settings) -> dict:
    """Per-frame APL, thickness, core hydration and protein orientation from .gro frames, plus anchor and box series."""
    series = {"frame": [], "apl_upper_A2": [], "apl_lower_A2": [], "thickness_A": [], "core_water_percent": [], "tilt_deg": [],
              "depth_A": [], "area_nm2": [], "water_path": [], "interdigitation": []}
    anchors, boxes = [], []
    for k, path in enumerate(paths):
        stage = stage_from_gro(Path(path), name="equilibrated")
        measure = measure_packing(stage)
        thickness = bilayer_thickness(stage)
        orientation = protein_orientation(stage["protein_nm"], stage["midplane_nm"]) if len(stage["protein_nm"]) else None
        hydration = core_hydration(stage["water_nm"], stage["protein_nm"], stage["midplane_nm"], CORE_HALF_NM, stage["box_nm"],
                                   thickness_nm=thickness["thickness_nm"], stage="equilibrated") if len(stage["water_nm"]) else None
        series["frame"].append(k)
        series["apl_upper_A2"].append(measure["upper"]["apl_A2"])
        series["apl_lower_A2"].append(measure["lower"]["apl_A2"])
        series["thickness_A"].append(thickness["thickness_A"])
        series["core_water_percent"].append(hydration["ratio_percent"] if hydration else np.nan)
        series["tilt_deg"].append(orientation["tilt_deg"] if orientation else np.nan)
        series["depth_A"].append(orientation["depth_A"] if orientation else np.nan)
        series["area_nm2"].append(stage["box_nm"][0] * stage["box_nm"][1])
        series["water_path"].append(bool(hydration["transmembrane_water_path"]) if hydration else None)
        series["interdigitation"].append(tail_interdigitation(stage["tails"][0], stage["tails"][1], stage["midplane_nm"])["overlap"]
                                         if "tails" in stage else np.nan)
        if k == 0:
            series["first_profiles"] = z_density_profiles(stage["z_groups"], stage["box_nm"]) if "z_groups" in stage else None
        anchors.append(np.array([[r["x_nm"], r["y_nm"]] for r in stage["lipids"]]))
        boxes.append(stage["box_nm"][:2])
        if k == len(paths) - 1:
            series["last_stage"] = stage
            series["last_profiles"] = z_density_profiles(stage["z_groups"], stage["box_nm"]) if "z_groups" in stage else None
    series["anchors"] = anchors
    series["boxes"] = np.array(boxes)
    series["species"] = [r["resname"] for r in series["last_stage"]["lipids"]] if paths else []
    return series


SERIES_NAMES = {"apl_upper_A2": "APL", "apl_lower_A2": "APL", "thickness_A": "bilayer thickness", "core_water_percent": "core water",
                "tilt_deg": "protein tilt", "depth_A": "protein insertion depth"}
# Protein orientation drifts are judged in their own units, not in percent of a mean angle: (PASS, WARNING) ceilings.
ORIENTATION_DRIFT = {"tilt_deg": (3.0, 7.0), "depth_A": (1.0, 2.0)}


def absolute_drift_record(name: str, stage: str, result: dict, units: str, ceilings: tuple) -> dict:
    """A convergence record graded on the absolute drift over the final window (slope x window length) in the series' units."""
    drift = result["slope_per_ps"] * result["window_ps"] if result.get("slope_per_ps") is not None else None
    status = "INSUFFICIENT SAMPLING" if result["status"] == "INSUFFICIENT SAMPLING" else classify(drift, *ceilings)
    if status != "INSUFFICIENT SAMPLING" and result.get("monotonic") and drift is not None and abs(drift) > ceilings[0]:
        status = "FAIL"
    reason = (f"drift {drift:+.2f} {units} over the final window (PASS <= {ceilings[0]}, WARNING <= {ceilings[1]})"
              + ("; monotonic trend" if result.get("monotonic") else "")) if drift is not None else result["reason"]
    return metric(name, stage, result["mean"], units, uncertainty=result["sem"], n=result["n"], deviation=drift, status=status, reason=reason)


def water_path_persistence(flags: list, stage: str = "equilibrated") -> dict:
    """Fraction of frames with a water column through the core: none PASS, transient (< half) WARNING, persistent FAIL."""
    seen = [f for f in flags if f is not None]
    if not seen:
        return metric("transmembrane water path persistence", stage, None, "fraction of frames", status="NOT RUN", reason="no water in the frames")
    fraction = sum(bool(f) for f in seen) / len(seen)
    status = "PASS" if fraction == 0 else "WARNING" if fraction < 0.5 else "FAIL"
    return metric("transmembrane water path persistence", stage, fraction, "fraction of frames", n=len(seen), status=status,
                  reason=f"a spanning water column in {sum(bool(f) for f in seen)} of {len(seen)} frames"
                  + (" (transient)" if 0 < fraction < 0.5 else " (persistent)" if fraction >= 0.5 else ""))


def convergence_gate(energy: dict | None, frames: dict | None, frame_time_ps: np.ndarray | None, settings: Settings) -> list[dict]:
    """Gate 3 records: drift and half-window tests of area, volume, APL, thickness, core hydration, tilt and depth."""
    records = []
    if energy is not None:
        for key, units in (("area_nm2", "nm^2"), ("volume_nm3", "nm^3")):
            records.append(convergence_record(f"box {key.split('_')[0]}", "equilibrated", convergence(energy["time_ps"], energy[key]), units))
    else:
        records.append(metric("box area", "equilibrated", None, "nm^2", status="NOT RUN", reason="no --energy xvg given"))
    if frames is not None and frame_time_ps is not None and len(frame_time_ps) == len(frames["frame"]):
        for key, units, leaflet in (("apl_upper_A2", "A^2", "upper"), ("apl_lower_A2", "A^2", "lower"), ("thickness_A", "A", None),
                                    ("core_water_percent", "% of bulk", None), ("tilt_deg", "deg", None), ("depth_A", "A", None)):
            values = np.array(frames[key], dtype=float)
            if not np.isfinite(values).all():
                records.append(metric(key, "equilibrated", None, units, leaflet, status="NOT RUN", reason="not measurable in these frames"))
                continue
            result = convergence(frame_time_ps, values)
            if key in ORIENTATION_DRIFT:
                records.append(absolute_drift_record(SERIES_NAMES[key], "equilibrated", result, units, ORIENTATION_DRIFT[key]))
            else:
                records.append(convergence_record(SERIES_NAMES[key], "equilibrated", result, units, leaflet))
        records.append(water_path_persistence(frames["water_path"]))
        if frames.get("last_profiles"):
            layered = layered_order(frames["last_profiles"], frames["last_stage"]["midplane_nm"], stage="equilibrated")
            records.append(layered["record"])
            if frames.get("first_profiles"):
                first, last = frames["first_profiles"]["density_nm3"], frames["last_profiles"]["density_nm3"]
                change = {k: float(np.sqrt(np.mean((last[k] - first[k]) ** 2)) / max(float(first[k].max()), 1e-9)) for k in last if k in first}
                records.append(metric("Z-density profile change first-to-last frame", "equilibrated", max(change.values()) if change else None,
                                      "fraction of peak", status="NOT RUN", reason="largest RMS change of a group profile over its peak; reported"))
    else:
        records.append(metric("APL", "equilibrated", None, "A^2", status="NOT RUN", reason="no --frames given (or their times are unknown)"))
    return records


def equilibrium_gate(energy: dict | None, frames: dict | None, frame_time_ps: np.ndarray | None, args, settings: Settings) -> list[dict]:
    """Gate 4 records on the production fraction: lateral diffusion, K_A, leaflet tension, final hydration and order."""
    records = []
    if energy is not None:
        n = len(energy["time_ps"])
        production = slice(int(n * (1.0 - args.production_fraction)), n)
        ka = area_compressibility(energy["area_nm2"][production], args.temperature)
        versus = ka_vs_reference(ka["ka_mn_per_m"], args.reference_ka)
        status = versus["status"] if ka["status"] == "PASS" else ka["status"]
        records.append(metric("area compressibility K_A", "final", ka["ka_mn_per_m"], "mN/m", reference=args.reference_ka,
                              uncertainty=ka["uncertainty_mn_per_m"], n=ka["n"], deviation=versus.get("deviation_percent"), status=status,
                              reason=f"{ka['reason']}; {versus['reason']}"))
    else:
        records.append(metric("area compressibility K_A", "final", None, "mN/m", status="NOT RUN", reason="no --energy xvg given"))
    if frames is not None and frame_time_ps is not None and len(frames["anchors"]) >= 10:
        n = len(frames["anchors"])
        start = int(n * (1.0 - args.production_fraction))
        positions = unwrap_xy(np.array(frames["anchors"][start:]), frames["boxes"][start:])
        time = frame_time_ps[start:] - frame_time_ps[start]
        species = np.array(frames["species"])
        groups = [("all lipids", np.ones(len(species), bool))] + [(s, species == s) for s in sorted(set(species)) if (species == s).sum() >= 20]
        for label, mask in groups:
            msd = msd_xy(positions[:, mask], time)
            fit = diffusion_fit(msd["lag_ps"], msd["msd_nm2"])
            versus = (diffusion_vs_reference(fit["d_um2_per_s"], args.reference_d) if fit["status"] in ("PASS", "WARNING")
                      else {"status": fit["status"], "reason": fit["reason"]})
            records.append(metric(f"lateral diffusion D ({label})", "final", fit["d_um2_per_s"], "um^2/s", reference=args.reference_d,
                                  n=int(mask.sum()), deviation=fit["alpha"], status=versus["status"] if fit["status"] == "PASS" else fit["status"],
                                  reason=f"alpha {fit['alpha']}, R^2 {fit['r2']}; {versus['reason']}"))
        last = frames["last_stage"]
        thickness = float(np.nanmean(frames["thickness_A"][-max(1, len(frames["thickness_A"]) // 2):]))
        if args.reference_thickness:
            deviation = 100.0 * (thickness - args.reference_thickness) / args.reference_thickness
            records.append(metric("bilayer thickness vs reference", "final", thickness, "A", reference=args.reference_thickness, deviation=deviation,
                                  status=classify(deviation, 3.0, 5.0), reason=f"{deviation:+.1f} % vs the matched reference (production mean)"))
        else:
            records.append(metric("bilayer thickness vs reference", "final", thickness, "A", status="REFERENCE NEEDED",
                                  reason="production-window mean; no matched reference given (--reference-thickness)"))
        inter = float(np.nanmean(frames["interdigitation"][-max(1, len(frames["interdigitation"]) // 2):]))
        if np.isfinite(inter):
            if args.reference_interdigitation:
                deviation = 100.0 * (inter - args.reference_interdigitation) / args.reference_interdigitation
                records.append(metric("tail interdigitation vs reference", "final", inter, "overlap", reference=args.reference_interdigitation,
                                      deviation=deviation, status=classify(deviation, 10.0, 20.0),
                                      reason=f"{deviation:+.1f} % vs the matched reference"))
            else:
                records.append(metric("tail interdigitation vs reference", "final", inter, "overlap", status="REFERENCE NEEDED",
                                      reason="production-window mean; no matched reference given (--reference-interdigitation)"))
        if "core_hydration" in last:
            records.append(metric("final core water", "final", last["core_hydration"]["ratio_percent"], "% of bulk",
                                  status=last["core_hydration"]["status"], reason="last frame; see the convergence series for persistence"))
        for name in sorted({r["resname"] for r in last["lipids"]}):
            residues = [res for res in last.get("residues", []) if res and res[0]["resname"] == name]
            if len(residues) >= 20 and name in AA_LIPIDS and name not in ("CHL1",):
                profile = scd_profile_from_residues(residues, CHARMM36_CHAINS)
                reference_scd = (args.reference_scd or {}).get(name, {})
                for chain, values in profile.items():
                    if values["n"]:
                        graded = scd_profile_rmse(reference_scd.get(chain), values["scd"])
                        records.append(metric(f"S_CD {name} {chain}", "final", float(np.nanmean(values["scd"])), "-", n=int(values["n"]),
                                              deviation=graded["rmse"], status=graded["status"],
                                              reason="mean S_CD over the chain; " + graded["reason"]))
    else:
        records.append(metric("lateral diffusion D", "final", None, "um^2/s", status="INSUFFICIENT SAMPLING",
                              reason="fewer than 10 frames (or no --frames)"))
    if args.stress:
        profile = read_stress_profile(Path(args.stress))
        midplane = frames["last_stage"]["midplane_nm"] if frames else float(np.median(profile["z_nm"]))
        tension = leaflet_tension(profile["z_nm"], profile["pi_bar"], midplane, profile.get("error_bar"))
        for key, leaflet in (("gamma_upper_mn_per_m", "upper"), ("gamma_lower_mn_per_m", "lower")):
            records.append(metric("leaflet tension", "final", tension[key], "mN/m", leaflet, uncertainty=tension[f"uncertainty_{leaflet}_mn_per_m"],
                                  status="NOT RUN", reason="per-leaflet value; the difference is graded"))
        records.append(metric("leaflet differential tension", "final", tension["delta_mn_per_m"], "mN/m",
                              uncertainty=tension["uncertainty_delta_mn_per_m"], status=tension["status"], reason=tension["reason"]))
    else:
        records.append(metric("leaflet differential tension", "final", None, "mN/m", status="NOT RUN", reason="no --stress profile given"))
    return records


def write_equilibration_report(out: Path, records: list[dict], series: dict | None) -> Path:
    """Write equilibration_validation.json/.md and append gates 3 and 4 to membrane_validation.json when it exists."""
    gate3 = [r for r in records if r["stage"] == "equilibrated"]
    gate4 = [r for r in records if r["stage"] == "final"]
    keys = ("frame", "apl_upper_A2", "apl_lower_A2", "thickness_A", "core_water_percent", "tilt_deg", "depth_A", "area_nm2")
    clean = lambda v: None if (isinstance(v, float) and not np.isfinite(v)) else v
    result = {"records": records,
              "gates": {"3 equilibration convergence": overall_status(gate3), "4 equilibrium physical properties": overall_status(gate4)},
              "series": {k: [clean(v) for v in series[k]] for k in keys} if series else None}
    (out / "equilibration_validation.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    text = ["# Equilibration and equilibrium validation\n",
            f"Gate 3 (convergence): **{result['gates']['3 equilibration convergence']['status']}**; "
            f"gate 4 (equilibrium properties): **{result['gates']['4 equilibrium physical properties']['status']}**\n",
            "## 3. Equilibration convergence\n", records_table(gate3), "## 4. Equilibrium physical validation\n", records_table(gate4)]
    (out / "equilibration_validation.md").write_text("\n".join(text))
    validation = out / "membrane_validation.json"
    if validation.is_file():
        record = json.loads(validation.read_text())
        record["records"] = [r for r in record.get("records", []) if r["stage"] not in ("equilibrated", "final")] + records
        record["gates"] = {**record.get("gates", {}), **result["gates"]}
        graded = [r for r in record["records"] if r["stage"] in ("slice", "backmap", "minimize", "equilibrated", "final")]
        record["overall"] = overall_status(graded)
        stages = {s["name"]: {} for s in record.get("stages", [])}
        stages["equilibrated"] = {"measure": True, "thickness": True, "orientation": True, "core_hydration": True}
        record["stage_matrix"] = stage_matrix(stages)
        validation.write_text(json.dumps(record, indent=2, default=str) + "\n")
        md = out / "membrane_validation.md"
        if md.is_file():
            body = md.read_text()
            sections = (("## 3. Equilibration convergence (gate 3)\n", gate3), ("## 4. Equilibrium physical validation (gate 4)\n", gate4))
            for heading, table in sections:
                body = body.replace(heading, heading + "\n" + records_table(table) + "\n", 1)
            previous = matrix_markdown(stage_matrix({s["name"]: {} for s in record.get("stages", [])}))
            body = body.replace(previous, matrix_markdown(record["stage_matrix"]), 1)
            md.write_text(body)
    return out / "equilibration_validation.md"


def main(argv: list | None = None) -> int:
    """Command line: grade an equilibration run and append it to the build's validation report."""
    parser = argparse.ArgumentParser(prog="python -m membraneforger.equilibration", description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", type=Path, required=True, help="the MembraneForger build directory")
    parser.add_argument("--energy", type=Path, help="gmx energy .xvg with Box-X, Box-Y, Box-Z (time in ps)")
    parser.add_argument("--frames", type=Path, nargs="*", default=[], help="evenly spaced .gro frames of the trajectory (whole molecules)")
    parser.add_argument("--frame-dt", type=float, default=None, metavar="PS",
                        help="time between frames in ps (default: spread over the energy series)")
    parser.add_argument("--production-fraction", type=float, default=0.5, help="final fraction of the run treated as production (default 0.5)")
    parser.add_argument("--temperature", type=float, default=310.0, metavar="K", help="simulation temperature for K_A (default 310 K)")
    parser.add_argument("--stress", type=Path, help="lateral pressure profile: z (nm), P_N - P_L (bar)[, error]")
    parser.add_argument("--reference-ka", type=float, default=None, metavar="MN_PER_M", help="matched reference K_A (mN/m)")
    parser.add_argument("--reference-d", type=float, default=None, metavar="UM2_PER_S", help="matched reference lipid diffusion (um^2/s)")
    parser.add_argument("--reference-thickness", type=float, default=None, metavar="A", help="matched reference bilayer thickness (A)")
    parser.add_argument("--reference-interdigitation", type=float, default=None, metavar="OVERLAP",
                        help="matched reference tail interdigitation (overlap fraction)")
    parser.add_argument("--reference-scd", type=Path, default=None, metavar="JSON",
                        help='matched S_CD profiles: {"POPC": {"sn1": [...], "sn2": [...]}, ...} in CHARMM36 carbon order')
    args = parser.parse_args(argv)
    args.reference_scd = json.loads(args.reference_scd.read_text()) if args.reference_scd else None
    settings = Settings()
    energy = box_series_from_xvg(args.energy) if args.energy else None
    frames = frame_series(args.frames, settings) if args.frames else None
    times = None
    if frames is not None:
        if args.frame_dt:
            times = np.arange(len(args.frames)) * args.frame_dt
        elif energy is not None:
            times = np.linspace(energy["time_ps"][0], energy["time_ps"][-1], len(args.frames))
    records = convergence_gate(energy, frames, times, settings) + equilibrium_gate(energy, frames, times, args, settings)
    path = write_equilibration_report(args.out, records, frames)
    verdict = overall_status(records)
    print(f"{args.out.name}\t{verdict['status']}\t{path}")
    for line in verdict["failures"] + verdict["warnings"]:
        print("  " + line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
