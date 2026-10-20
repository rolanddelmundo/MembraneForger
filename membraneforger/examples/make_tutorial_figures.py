#!/usr/bin/env python
"""Draw the tutorial figures from a folder written by measure_tutorial_numbers.py (tutorial50_out).

    python examples/make_tutorial_figures.py TUTORIAL50_OUT [--out membraneforger/docs/figures]

Writes, per figure, a PNG and a CSV with the numbers behind it, and two tables:

    KOR_Membrane1_before_after.png / .csv   KOR_Membrane1 as bundled, embedded and sliced: Voronoi area per lipid, top-down and side views
    protein_density_by_frame.png / .csv     density of all lipids within R of the protein over the leaflet mean, every frame
    table_apl.csv                           leaflet APL by stage, same lipids uncut, head-plane APL, every frame and leaflet
    table_protein_density.csv               the 10 A first-shell ratio, embedded and slice, every frame and leaflet

Input: TUTORIAL50_OUT/summary.json and TUTORIAL50_OUT/cases/<FRAME>__6whc/membrane_validation.{json} and
membrane_validation_lipids.tsv (each case is a build stopped before backmapping). Needs numpy and matplotlib.
"""
import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 15, "axes.labelsize": 15, "axes.titlesize": 15,
    "legend.fontsize": 13.5, "xtick.labelsize": 14, "ytick.labelsize": 14, "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})
UPPER, LOWER = "#0072B2", "#E69F00"
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9", "#F0E442", "#6e6e6e"]


def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


def before_after(src, out):
    """KOR_Membrane1: the bundled frame, the membrane after embedding the GCGR complex, and the slice cut around it."""
    case = src / "cases/KOR_Membrane1__6whc"
    lipids = list(csv.DictReader(open(case / "membrane_validation_lipids.tsv"), delimiter="\t"))
    for r in lipids:
        for key in ("x_A", "y_A", "z_A", "voronoi_area_A2"):
            r[key] = float(r[key])
    stages = {s["name"]: s for s in json.load(open(case / "membrane_validation.json"))["stages"]}
    columns = [("cg_frame", "CG frame (50 µs,\nframe's own receptor)"), ("embed", "Embedded\n(GCGR complex)"), ("slice", "Slice\n(cut around GCGR)")]
    fig = plt.figure(figsize=(13, 12))
    grid = fig.add_gridspec(3, 4, width_ratios=[1, 1, 1, 0.06], height_ratios=[1, 1, 0.55], hspace=0.28, wspace=0.12)
    scatter = None
    for col, (stage, title) in enumerate(columns):
        box = stages[stage]["box_nm"]
        for row, leaflet in enumerate(("upper", "lower")):
            ax = fig.add_subplot(grid[row, col])
            pts = [r for r in lipids if r["stage"] == stage and r["leaflet"] == leaflet]
            scatter = ax.scatter([r["x_A"] / 10 for r in pts], [r["y_A"] / 10 for r in pts], c=[r["voronoi_area_A2"] for r in pts], s=26,
                                 marker="s", cmap="coolwarm", vmin=20, vmax=80, linewidths=0)
            ax.add_patch(plt.Rectangle((0, 0), box[0], box[1], fill=False, ec="#333333", lw=1.2))
            ax.set_xlim(-0.5, box[0] + 0.5)
            ax.set_ylim(-0.5, box[1] + 0.5)
            ax.set_aspect("equal")
            ax.set_xticks([0, 5, 10, 15])
            ax.set_yticks([0, 5, 10, 15])
            if row == 0:
                ax.set_title(title, fontweight="bold", fontsize=14)
            if col == 0:
                ax.set_ylabel(f"{leaflet.capitalize()} leaflet\ny (nm)")
            else:
                ax.set_yticklabels([])
            if row == 1:
                ax.set_xlabel("x (nm)")
            else:
                ax.set_xticklabels([])
            mean = np.mean([r["voronoi_area_A2"] for r in pts])
            ax.text(0.03, 0.03, f"N = {len(pts)}   mean {mean:.1f} Å$^2$", transform=ax.transAxes, fontsize=11,
                    bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))
        ax = fig.add_subplot(grid[2, col])
        for leaflet, colour in (("upper", UPPER), ("lower", LOWER)):
            pts = [r for r in lipids if r["stage"] == stage and r["leaflet"] == leaflet]
            ax.scatter([r["x_A"] / 10 for r in pts], [r["z_A"] / 10 for r in pts], s=5, color=colour, label=leaflet + " leaflet", linewidths=0)
        ax.set_xlim(-0.5, box[0] + 0.5)
        ax.set_xlabel("x (nm)")
        z = [r["z_A"] / 10 for r in lipids if r["stage"] == stage]
        ax.set_ylim(min(z) - 0.5, max(z) + 1.8)
        if col == 0:
            ax.set_ylabel("z (nm)")
            ax.legend(frameon=False, loc="upper left", fontsize=11, markerscale=3, ncol=2)
        else:
            ax.set_yticklabels([])
    colourbar = fig.colorbar(scatter, cax=fig.add_subplot(grid[0:2, 3]))
    colourbar.set_label("Voronoi area per lipid (Å$^2$)")
    fig.suptitle("KOR_Membrane1: the bundled 50 µs frame and the membrane it becomes in the build", fontsize=15, fontweight="bold", y=0.94)
    fig.savefig(out / "KOR_Membrane1_before_after.png")
    plt.close(fig)
    write_csv(out / "KOR_Membrane1_before_after.csv", ["stage", "leaflet", "resid", "resname", "x_A", "y_A", "z_A", "voronoi_area_A2"],
              [[r["stage"], r["leaflet"], r["resid"], r["resname"], r["x_A"], r["y_A"], r["z_A"], r["voronoi_area_A2"]] for r in lipids])


def protein_density_and_apl(src, out, frames):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    curve_rows, apl_rows, density_rows = [], [], []
    for number, frame in enumerate(frames):
        report = json.load(open(src / f"cases/{frame}__6whc/membrane_validation.json"))
        stages = {s["name"]: s for s in report["stages"]}
        check = report["slice_check"]["leaflets"]
        for ax, leaflet in zip(axes, ("upper", "lower")):
            curve = stages["slice"]["protein_density"][leaflet]["species"]["all"]
            ax.plot(curve["R_A"], curve["cumulative_ratio"], color=PALETTE[number % len(PALETTE)], lw=2, label=frame)
            curve_rows += [[frame, leaflet, r, v] for r, v in zip(curve["R_A"], curve["cumulative_ratio"])]
            shell = {name: stages[name]["protein_density"][leaflet]["species"]["all"] for name in ("embed", "slice")}
            density_rows.append([frame, leaflet, shell["embed"]["first_shell_ratio"], shell["embed"]["first_shell_n"], shell["embed"]["n"],
                                 shell["slice"]["first_shell_ratio"], shell["slice"]["first_shell_n"], shell["slice"]["n"]])
            apl_rows.append([frame, leaflet] + [stages[n]["leaflets"][leaflet]["apl_A2"] for n in ("cg_frame", "embed", "slice")] +
                            [check[leaflet]["region_reference_apl_A2"], check[leaflet]["construction_delta_percent"],
                             stages["embed"]["leaflets"][leaflet]["head_apl_A2"], stages["slice"]["leaflets"][leaflet]["head_apl_A2"]])
    for ax, leaflet in zip(axes, ("upper", "lower")):
        ax.axhline(1.0, color="#333333", lw=0.8, ls=":")
        ax.axvline(10, color="#999999", lw=0.8, ls=":")
        ax.set_xlabel("R, distance from the protein (Å)")
        ax.set_title(f"{leaflet.capitalize()} leaflet, all lipids", fontweight="bold")
    axes[0].set_ylabel("Lipid density within R / leaflet mean")
    axes[1].legend(frameon=False, ncol=2, fontsize=11, loc="upper right")
    fig.tight_layout()
    fig.savefig(out / "protein_density_by_frame.png")
    plt.close(fig)
    write_csv(out / "protein_density_by_frame.csv", ["frame", "leaflet", "R_A", "cumulative_ratio"], curve_rows)
    write_csv(out / "table_apl.csv", ["frame", "leaflet", "apl_cg_frame_A2", "apl_embedded_A2", "apl_slice_A2", "apl_same_lipids_uncut_A2",
                                      "construction_delta_percent", "head_plane_apl_embedded_A2", "head_plane_apl_slice_A2"], apl_rows)
    write_csv(out / "table_protein_density.csv", ["frame", "leaflet", "ratio_10A_embedded", "n_shell_embedded", "n_leaflet_embedded",
                                                  "ratio_10A_slice", "n_shell_slice", "n_leaflet_slice"], density_rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tutorial50_out", help="the folder written by measure_tutorial_numbers.py")
    parser.add_argument("--out", default="membraneforger/docs/figures", help="folder for the PNG and CSV files (created)")
    args = parser.parse_args()
    src, out = Path(args.tutorial50_out), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    results = json.load(open(src / "summary.json"))["results"]
    frames = [f for f in sorted(results) if f.startswith("KOR")] + [f for f in sorted(results) if not f.startswith("KOR")]
    before_after(src, out)
    protein_density_and_apl(src, out, frames)
    print("wrote", ", ".join(sorted(p.name for p in out.iterdir())))


if __name__ == "__main__":
    main()
