#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Equilibration convergence and equilibrium physical properties of an atomistic membrane trajectory.
#//=============================================================
"""Has the built membrane equilibrated, and does it behave like the lipid bilayer it is supposed to be?

Two families of questions, both answered from numpy arrays so that every estimator can be tested on synthetic data;
the readers (GROMACS .xvg, .gro frames, an MDAnalysis trajectory, a stress profile) are thin and separate.

    convergence             drift, half-difference, monotonic trend and block SEM of any scalar time series
    order_parameters        S_CD = <(3 cos^2 theta - 1) / 2> of C-H vectors relative to the bilayer normal (z)
    msd_xy / diffusion_fit  lateral MSD, its exponent alpha, and the 2D diffusion coefficient D = slope / 4
    area_compressibility    K_A = k_B T <A> / var(A) with a block-bootstrap uncertainty
    leaflet_tension         per-leaflet surface tension from a lateral pressure profile and its upper - lower difference

Units: nm, ps and degrees; areas nm^2; K_A and tensions in mN/m. Diffusion: an MSD slope in nm^2/ps is converted
with 1 nm^2/ps = 1e3 nm^2/ns = 1e6 um^2/s (1 nm^2 = 1e-6 um^2, 1 ns = 1e-9 s), and D is reported both in um^2/s and
in 1e-8 cm^2/s (1 um^2/s = 1e-12 m^2/s = 1e-8 cm^2/s, so the two numbers coincide). Tension: pi(z) in bar
integrated over z in nm gives bar nm = 1e5 Pa * 1e-9 m = 1e-4 N/m = 0.1 mN/m.

Every graded quantity ends in a qc.metric record; poor sampling is reported as INSUFFICIENT SAMPLING, a missing
experimental or control value as REFERENCE NEEDED, never as a PASS or a FAIL.
"""
from pathlib import Path

import numpy as np

from . import qc
from .structio import read_gro

__all__ = ['TENSION_THRESHOLDS_MN_PER_M', 'CHARMM36_CHAINS', 'KB_J_PER_K', 'NM2_PER_PS_TO_UM2_PER_S', 'BAR_NM_TO_MN_PER_M',
           'read_xvg', 'box_series_from_xvg', 'read_gro_frames', 'trajectory_frames', 'read_stress_profile',
           'block_sem', 'convergence', 'convergence_record',
           'order_parameters', 'scd_profile_from_residues', 'scd_profile_rmse',
           'unwrap_xy', 'msd_xy', 'diffusion_fit', 'diffusion_vs_reference',
           'area_compressibility', 'ka_vs_reference', 'leaflet_tension']

KB_J_PER_K = 1.380649e-23
NM2_PER_PS_TO_UM2_PER_S = 1.0e6          # 1 nm^2/ps = 1e3 nm^2/ns = 1e6 um^2/s
UM2_PER_S_TO_1E8_CM2_PER_S = 1.0         # 1 um^2/s = 1e-8 cm^2/s
BAR_NM_TO_MN_PER_M = 0.1                 # 1 bar * 1 nm = 1e5 Pa * 1e-9 m = 1e-4 N/m = 0.1 mN/m

# |gamma_upper - gamma_lower| ceilings (PASS, WARNING) and the magnitude above which one leaflet's tension is flagged.
TENSION_THRESHOLDS_MN_PER_M = {"pass": 5.0, "warning": 10.0, "leaflet_flag": 10.0}

# CHARMM36 chain atom names of a POPC-like phospholipid: sn1 carbons C3k carry H{k}X/H{k}Y (H{k}Z on the terminal
# methyl); sn2 carbons C2k carry H{k}R/H{k}S (H{k}T on the methyl, a single H{k}1 on an unsaturated carbon). Only
# the hydrogens that exist in a residue are used, so one table serves saturated and unsaturated chains.
CHARMM36_CHAINS = {
    "sn1": [(f"C3{k}", [f"H{k}X", f"H{k}Y", f"H{k}Z"]) for k in range(2, 19)],
    "sn2": [(f"C2{k}", [f"H{k}R", f"H{k}S", f"H{k}T", f"H{k}1"]) for k in range(2, 19)],
}


# ---------------------------------------------------------------------------------------------------- readers
def read_xvg(path: Path) -> dict:
    """A GROMACS .xvg: {"legends": [column names after time], "labels": {xaxis, yaxis}, "data": (n, columns) array, time first}."""
    legends, labels, rows = {}, {}, []
    for raw in Path(path).read_text(errors="replace").splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        if s.startswith("@"):
            parts = s[1:].split()
            if len(parts) >= 3 and parts[0].startswith("s") and parts[0][1:].isdigit() and parts[1] == "legend":
                legends[int(parts[0][1:])] = s.split('"')[1] if '"' in s else " ".join(parts[2:])
            elif len(parts) >= 3 and parts[0] in ("xaxis", "yaxis") and parts[1] == "label":
                labels[parts[0]] = s.split('"')[1] if '"' in s else " ".join(parts[2:])
            continue
        rows.append([float(v) for v in s.split()])
    data = np.array(rows, dtype=float) if rows else np.zeros((0, 1 + len(legends)))
    return {"legends": [legends.get(i, f"s{i}") for i in range(max(data.shape[1] - 1, len(legends)))], "labels": labels, "data": data}


def box_series_from_xvg(path: Path) -> dict:
    """Time (ps) and box_x, box_y, box_z (nm) from a `gmx energy` .xvg with Box-X/Box-Y/Box-Z, plus area (nm^2) and volume (nm^3)."""
    xvg = read_xvg(path)
    columns = {}
    for axis in ("X", "Y", "Z"):
        hits = [i for i, name in enumerate(xvg["legends"]) if name.strip().lower() == f"box-{axis.lower()}"]
        if not hits:
            raise SystemExit(f"{path}: no Box-{axis} column (legends: {xvg['legends']}); run `gmx energy` with Box-X Box-Y Box-Z")
        columns[axis] = xvg["data"][:, hits[0] + 1]
    out = {"time_ps": xvg["data"][:, 0], "box_x_nm": columns["X"], "box_y_nm": columns["Y"], "box_z_nm": columns["Z"]}
    out["area_nm2"] = out["box_x_nm"] * out["box_y_nm"]
    out["volume_nm3"] = out["area_nm2"] * out["box_z_nm"]
    return out


def read_gro_frames(paths: list) -> list[tuple[list[dict], list[float]]]:
    """(atoms, box) of each of several .gro files (nm), in the order given."""
    return [read_gro(Path(p)) for p in paths]


def trajectory_frames(topology: Path, trajectory: Path, selection: str | None = None):
    """Iterate (atoms, box) in nm over an MDAnalysis-readable trajectory; MDAnalysis is optional and imported here only."""
    try:
        import MDAnalysis as mda
    except ImportError:
        raise SystemExit("reading a trajectory needs MDAnalysis, which is not installed: `pip install MDAnalysis` "
                         "(or convert the frames to .gro files and use read_gro_frames)") from None
    universe = mda.Universe(str(topology), str(trajectory))
    group = universe.select_atoms(selection) if selection else universe.atoms

    def frames():
        for _ in universe.trajectory:
            xyz = group.positions / 10.0
            atoms = [{"resid": int(rid), "resname": str(rn), "atom": str(an), "x": float(p[0]), "y": float(p[1]), "z": float(p[2])}
                     for rid, rn, an, p in zip(group.resids, group.resnames, group.names, xyz)]
            yield atoms, [float(v) / 10.0 for v in universe.dimensions[:3]]
    return frames()


def read_stress_profile(path: Path) -> dict:
    """A lateral pressure profile text file (z nm, pi bar[, error bar]; # comments): {"z_nm", "pi_bar", "error_bar" or None}."""
    rows = [[float(v) for v in line.split("#", 1)[0].split()] for line in Path(path).read_text(errors="replace").splitlines()
            if line.split("#", 1)[0].strip()]
    if not rows or min(len(r) for r in rows) < 2:
        raise SystemExit(f"{path}: expected two or three numeric columns (z, lateral pressure[, error])")
    width = min(len(r) for r in rows)
    data = np.array([r[:width] for r in rows], dtype=float)
    return {"z_nm": data[:, 0], "pi_bar": data[:, 1], "error_bar": data[:, 2] if width >= 3 else None}


# ---------------------------------------------------------------------------------------------------- convergence
def block_sem(series: np.ndarray, blocks: int = 5) -> float:
    """Standard error of the mean from `blocks` consecutive block averages (ddof = 1); nan with fewer than 2 blocks."""
    s = np.asarray(series, dtype=float)
    n = len(s) // blocks
    if blocks < 2 or n < 1:
        return float("nan")
    means = s[:n * blocks].reshape(blocks, n).mean(axis=1)
    return float(means.std(ddof=1) / np.sqrt(blocks))


def convergence(time_ps: np.ndarray, series: np.ndarray, final_fraction: float = 0.5, blocks: int = 5) -> dict:
    """Is a scalar time series (area, volume, APL, thickness, hydration, tilt, depth) stationary over its final window?

    Drift: a linear fit over the last `final_fraction` of the run, slope * window length as percent of the window
    mean. Half difference: second-half minus first-half mean as percent of the overall mean. Monotonic: the
    `blocks` consecutive block means of the whole series rise or fall strictly. SEM: block SEM over the window.
    Status: the worse of drift (PASS <= 1 %, WARNING <= 3 %) and half difference (PASS <= 2 %, WARNING <= 5 %),
    FAIL also when the trend is monotonic and the drift exceeds 1 %.
    """
    t, s = np.asarray(time_ps, dtype=float), np.asarray(series, dtype=float)
    keep = np.isfinite(t) & np.isfinite(s)
    t, s = t[keep], s[keep]
    n = len(s)
    if n < 2 * blocks or t[-1] <= t[0]:
        return {"n": int(n), "mean": float(s.mean()) if n else None, "sem": None, "drift_percent": None, "half_difference_percent": None,
                "monotonic": False, "slope_per_ps": None, "window_ps": None, "status": "INSUFFICIENT SAMPLING",
                "reason": f"{n} samples, fewer than {2 * blocks}" if n < 2 * blocks else "time does not advance"}
    window = t >= t[-1] - final_fraction * (t[-1] - t[0])
    tw, sw = t[window], s[window]
    mean_w = float(sw.mean())
    slope = float(np.polyfit(tw, sw, 1)[0]) if len(sw) >= 2 and tw[-1] > tw[0] else 0.0
    length = float(tw[-1] - tw[0])
    drift = 100.0 * slope * length / mean_w if mean_w else float("inf")
    overall = float(s.mean())
    half = 100.0 * (s[n // 2:].mean() - s[:n // 2].mean()) / overall if overall else float("inf")
    width = n // blocks
    block_means = s[:width * blocks].reshape(blocks, width).mean(axis=1)
    steps = np.diff(block_means)
    monotonic = bool(np.all(steps > 0) or np.all(steps < 0))
    sem = block_sem(sw, blocks) if len(sw) >= 2 * blocks else block_sem(s, blocks)
    status = qc.worst([qc.classify(drift, 1.0, 3.0), qc.classify(half, 2.0, 5.0)])
    reasons = [f"drift {drift:+.2f} % of the mean over the final {int(round(100 * final_fraction))} %",
               f"second half - first half {half:+.2f} %"]
    if monotonic and abs(drift) > 1.0:
        status = "FAIL"
        reasons.append(f"block means strictly monotonic over {blocks} blocks")
    elif monotonic:
        reasons.append("block means monotonic but the drift is within 1 %")
    return {"n": int(n), "mean": mean_w, "sem": sem, "drift_percent": float(drift), "half_difference_percent": float(half),
            "monotonic": monotonic, "slope_per_ps": slope, "window_ps": length, "block_means": block_means.tolist(),
            "status": status, "reason": "; ".join(reasons)}


def convergence_record(name: str, stage: str, result: dict, units: str, leaflet: str | None = None) -> dict:
    """The qc.metric record of a convergence() result: value = window mean, uncertainty = block SEM, deviation = drift %."""
    return qc.metric(name, stage, result["mean"], units, leaflet=leaflet, uncertainty=result["sem"], n=result["n"],
                     deviation=result["drift_percent"], status=result["status"], reason=result["reason"])


# ---------------------------------------------------------------------------------------------------- order parameters
def order_parameters(ch_vectors) -> np.ndarray:
    """S_CD per carbon from C->H vectors relative to z: an (n, n_carbons, n_H, 3) array or a list of (..., 3) arrays per carbon."""
    def scd(v):
        v = np.asarray(v, dtype=float).reshape(-1, 3)
        norm = np.linalg.norm(v, axis=1)
        v = v[norm > 0]
        cos = v[:, 2] / norm[norm > 0]
        return float(np.mean(1.5 * cos ** 2 - 0.5)) if len(v) else float("nan")
    if isinstance(ch_vectors, np.ndarray) and ch_vectors.ndim == 4:
        return np.array([scd(ch_vectors[:, c]) for c in range(ch_vectors.shape[1])])
    return np.array([scd(v) for v in ch_vectors])


def scd_profile_from_residues(residues: list[list[dict]], chain_spec: dict = CHARMM36_CHAINS) -> dict:
    """Per-chain S_CD profile of one lipid species: {chain: {"carbon", "index", "scd", "sem", "n"}} over the residues (nm)."""
    # Each residue contributes one value per carbon (its hydrogens averaged); the SEM is over residues (ddof = 1).
    out = {}
    for chain, carbons in chain_spec.items():
        names, index, means, sems, counts = [], [], [], [], []
        for k, (carbon, hydrogens) in enumerate(carbons, 2):
            per_molecule = []
            for residue in residues:
                atoms = {a["atom"]: a for a in residue}
                c = atoms.get(carbon)
                hs = [atoms[h] for h in hydrogens if h in atoms]
                if c is None or not hs:
                    continue
                v = np.array([[h["x"] - c["x"], h["y"] - c["y"], h["z"] - c["z"]] for h in hs])
                per_molecule.append(order_parameters([v])[0])
            if not per_molecule:
                continue
            values = np.array(per_molecule)
            names.append(carbon)
            index.append(int(carbon[2:]) if carbon[2:].isdigit() else k)
            means.append(float(values.mean()))
            sems.append(float(values.std(ddof=1) / np.sqrt(len(values))) if len(values) > 1 else float("nan"))
            counts.append(int(len(values)))
        out[chain] = {"carbon": names, "index": index, "scd": means, "sem": sems, "n": counts}
    return out


def scd_profile_rmse(reference, sample) -> dict:
    """RMS difference of two S_CD profiles on the same carbons: PASS <= 0.03, WARNING <= 0.05, FAIL beyond; no reference -> REFERENCE NEEDED."""
    s = np.asarray(sample, dtype=float)
    if reference is None:
        return {"rmse": None, "n": int(np.isfinite(s).sum()), "status": "REFERENCE NEEDED", "reason": "no reference S_CD profile"}
    r = np.asarray(reference, dtype=float)
    n = min(len(r), len(s))
    ok = np.isfinite(r[:n]) & np.isfinite(s[:n])
    if not ok.any():
        return {"rmse": None, "n": 0, "status": "INSUFFICIENT SAMPLING", "reason": "no carbons in common"}
    rmse = float(np.sqrt(np.mean((s[:n][ok] - r[:n][ok]) ** 2)))
    return {"rmse": rmse, "n": int(ok.sum()), "status": qc.classify(rmse, 0.03, 0.05),
            "reason": f"RMS |S_CD - reference| = {rmse:.3f} over {int(ok.sum())} carbons"}


# ---------------------------------------------------------------------------------------------------- diffusion
def unwrap_xy(positions_wrapped: np.ndarray, box_xy_per_frame: np.ndarray) -> np.ndarray:
    """Undo periodic jumps of (n_frames, n, 2) wrapped x, y positions; the box is (n_frames, 2) or one (2,) cell."""
    p = np.asarray(positions_wrapped, dtype=float)
    box = np.broadcast_to(np.asarray(box_xy_per_frame, dtype=float).reshape(-1, 1, 2), (p.shape[0], 1, 2))
    steps = np.diff(p, axis=0)
    steps -= box[1:] * np.round(steps / box[1:])
    return np.concatenate([p[:1], p[:1] + np.cumsum(steps, axis=0)], axis=0)


def msd_xy(positions: np.ndarray, time_ps: np.ndarray, remove_drift: bool = True, reference: np.ndarray | None = None) -> dict:
    """Lateral MSD(tau) of unwrapped (n_frames, n_lipids, 2) positions, averaged over lipids and time origins up to half the run.

    The per-frame centre of mass of the set (or of `reference`, e.g. the whole leaflet) is subtracted first so
    that a drifting leaflet is not mistaken for diffusion. Lags are in ps from the (uniform) time spacing.
    """
    p = np.asarray(positions, dtype=float)
    t = np.asarray(time_ps, dtype=float)
    if remove_drift:
        p = p - np.asarray(reference if reference is not None else p, dtype=float).mean(axis=1, keepdims=True)
    n_frames = p.shape[0]
    max_lag = n_frames // 2
    dt = float(np.median(np.diff(t))) if n_frames > 1 else float("nan")
    # Kneller's FFT route: MSD(m) = S1(m) - 2 S2(m), with S2 the position autocorrelation summed over x, y.
    length = 1 << (2 * n_frames - 1).bit_length()
    fft = np.fft.rfft(p, n=length, axis=0)
    s2 = np.fft.irfft(fft * np.conj(fft), n=length, axis=0)[:n_frames].sum(axis=2)
    s2 /= (n_frames - np.arange(n_frames))[:, None]
    squares = (p ** 2).sum(axis=2)
    total = 2.0 * squares.sum(axis=0)
    s1 = np.empty_like(s2)
    for m in range(n_frames):
        if m:
            total = total - squares[m - 1] - squares[n_frames - m]
        s1[m] = total / (n_frames - m)
    msd = (s1 - 2.0 * s2).mean(axis=1)[1:max_lag + 1]
    return {"lag_ps": dt * np.arange(1, max_lag + 1), "msd_nm2": msd, "n_lipids": int(p.shape[1]), "n_frames": int(n_frames)}


def diffusion_fit(lag_ps: np.ndarray, msd: np.ndarray, window: tuple = (0.2, 0.5)) -> dict:
    """alpha from log MSD vs log tau and D = slope / 4 (2D) from MSD = 4 D tau + c, both over the given fraction window of the lags.

    D is in um^2/s and in 1e-8 cm^2/s (numerically equal). PASS: alpha in [0.9, 1.1] and R^2 >= 0.98; WARNING:
    alpha in [0.8, 0.9) or (1.1, 1.2]; otherwise INSUFFICIENT SAMPLING (not in a diffusive regime).
    """
    tau, m = np.asarray(lag_ps, dtype=float), np.asarray(msd, dtype=float)
    inside = (tau >= window[0] * tau.max()) & (tau <= window[1] * tau.max()) & (m > 0) & (tau > 0)
    if inside.sum() < 3:
        return {"alpha": None, "d_um2_per_s": None, "d_1e8_cm2_per_s": None, "slope_nm2_per_ps": None, "r2": None, "n_lags": int(inside.sum()),
                "window": list(window), "status": "INSUFFICIENT SAMPLING", "reason": "fewer than 3 lags in the fit window"}
    alpha = float(np.polyfit(np.log(tau[inside]), np.log(m[inside]), 1)[0])
    slope, intercept = np.polyfit(tau[inside], m[inside], 1)
    fitted = slope * tau[inside] + intercept
    ss_res, ss_tot = float(((m[inside] - fitted) ** 2).sum()), float(((m[inside] - m[inside].mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    d_um2 = float(slope) / 4.0 * NM2_PER_PS_TO_UM2_PER_S
    if 0.9 <= alpha <= 1.1 and r2 >= 0.98:
        status, reason = "PASS", f"alpha = {alpha:.2f}, R^2 = {r2:.3f}: diffusive"
    elif 0.8 <= alpha < 0.9 or 1.1 < alpha <= 1.2:
        status, reason = "WARNING", f"alpha = {alpha:.2f}, R^2 = {r2:.3f}: marginally diffusive"
    else:
        status, reason = "INSUFFICIENT SAMPLING", f"alpha = {alpha:.2f}, R^2 = {r2:.3f}: not in a diffusive regime"
    return {"alpha": alpha, "d_um2_per_s": d_um2, "d_1e8_cm2_per_s": d_um2 * UM2_PER_S_TO_1E8_CM2_PER_S, "slope_nm2_per_ps": float(slope),
            "r2": float(r2), "n_lags": int(inside.sum()), "window": list(window), "status": status, "reason": reason}


def diffusion_vs_reference(d, d_ref) -> dict:
    """Grade D against a reference by the factor max(D/D_ref, D_ref/D): <= 2 PASS, <= 5 WARNING, beyond FAIL; no reference -> REFERENCE NEEDED."""
    if d_ref is None:
        return {"factor": None, "status": "REFERENCE NEEDED", "reason": "no reference diffusion coefficient"}
    if d is None or not np.isfinite(d) or d <= 0 or d_ref <= 0:
        return {"factor": None, "status": "INSUFFICIENT SAMPLING", "reason": "no positive diffusion coefficient to compare"}
    factor = float(max(d / d_ref, d_ref / d))
    return {"factor": factor, "status": qc.classify(factor, 2.0, 5.0), "reason": f"D differs from the reference by a factor {factor:.2f}"}


# ---------------------------------------------------------------------------------------------------- compressibility
def area_compressibility(area_nm2_series: np.ndarray, temperature_k: float, blocks: int = 5, bootstrap: int = 200, seed: int = 0) -> dict:
    """K_A = k_B T <A> / var(A) in mN/m, with its SD over a block bootstrap (blocks resampled with replacement).

    Relative uncertainty <= 20 % is acceptable (PASS), 20-40 % WARNING, beyond INSUFFICIENT SAMPLING.
    """
    a = np.asarray(area_nm2_series, dtype=float)
    a = a[np.isfinite(a)]
    ka = lambda x: KB_J_PER_K * temperature_k * x.mean() * 1e-18 / (x.var(ddof=1) * 1e-36) * 1000.0
    if len(a) < 2 * blocks or a.var(ddof=1) <= 0:
        return {"ka_mn_per_m": None, "uncertainty_mn_per_m": None, "relative_uncertainty_percent": None,
                "mean_area_nm2": float(a.mean()) if len(a) else None, "n": int(len(a)), "status": "INSUFFICIENT SAMPLING",
                "reason": f"{len(a)} area samples, fewer than {2 * blocks} or no fluctuation"}
    value = float(ka(a))
    width = len(a) // blocks
    chunks = a[:width * blocks].reshape(blocks, width)
    rng = np.random.default_rng(seed)
    resampled = np.array([ka(chunks[rng.integers(0, blocks, blocks)].ravel()) for _ in range(bootstrap)])
    sd = float(resampled.std(ddof=1)) if bootstrap > 1 else float("nan")
    relative = 100.0 * sd / value
    status = qc.classify(relative, 20.0, 40.0)
    status = "INSUFFICIENT SAMPLING" if status == "FAIL" else status
    return {"ka_mn_per_m": value, "uncertainty_mn_per_m": sd, "relative_uncertainty_percent": float(relative), "mean_area_nm2": float(a.mean()),
            "n": int(len(a)), "status": status,
            "reason": f"K_A = {value:.0f} +/- {sd:.0f} mN/m ({relative:.0f} % from {bootstrap} block-bootstrap resamples)"}


def ka_vs_reference(ka, ka_ref) -> dict:
    """Grade K_A against a reference by relative deviation: <= 20 % PASS, <= 40 % WARNING, beyond FAIL; no reference -> REFERENCE NEEDED."""
    if ka_ref is None:
        return {"deviation_percent": None, "status": "REFERENCE NEEDED", "reason": "no reference area compressibility modulus"}
    if ka is None or not np.isfinite(ka) or not ka_ref:
        return {"deviation_percent": None, "status": "INSUFFICIENT SAMPLING", "reason": "no K_A to compare"}
    deviation = 100.0 * (ka - ka_ref) / ka_ref
    return {"deviation_percent": float(deviation), "status": qc.classify(deviation, 20.0, 40.0),
            "reason": f"K_A deviates {deviation:+.0f} % from the reference"}


# ---------------------------------------------------------------------------------------------------- tension
def leaflet_tension(z_nm: np.ndarray, lateral_pressure_bar: np.ndarray, midplane_nm: float, uncertainty_bar: np.ndarray | None = None) -> dict:
    """Per-leaflet tension from a lateral pressure profile pi(z) = P_N - P_L: gamma = sum pi dz over each side of the midplane, in mN/m.

    delta = gamma_upper - gamma_lower is graded with TENSION_THRESHOLDS_MN_PER_M (PASS, WARNING ceilings); beyond
    the WARNING ceiling it is a FAIL only when the propagated uncertainty is smaller than |delta|, otherwise
    INSUFFICIENT SAMPLING. A leaflet with |gamma| above the leaflet_flag ceiling is a WARNING on its own.
    """
    z, pi = np.asarray(z_nm, dtype=float), np.asarray(lateral_pressure_bar, dtype=float)
    err = None if uncertainty_bar is None else np.asarray(uncertainty_bar, dtype=float)
    dz = float(np.median(np.diff(np.sort(z)))) if len(z) > 1 else float("nan")
    upper, lower = z > midplane_nm, z < midplane_nm
    gamma = lambda side: float(BAR_NM_TO_MN_PER_M * (pi[side] * dz).sum())
    sigma = lambda side: None if err is None else float(BAR_NM_TO_MN_PER_M * dz * np.sqrt((err[side] ** 2).sum()))
    g_up, g_lo = gamma(upper), gamma(lower)
    u_up, u_lo = sigma(upper), sigma(lower)
    delta = g_up - g_lo
    u_delta = None if err is None else float(np.sqrt(u_up ** 2 + u_lo ** 2))
    thresholds = TENSION_THRESHOLDS_MN_PER_M
    if not (upper.any() and lower.any()) or not np.isfinite(dz):
        status, reason = "INSUFFICIENT SAMPLING", "the profile does not cover both sides of the midplane"
    else:
        status = qc.classify(delta, thresholds["pass"], thresholds["warning"])
        reason = f"gamma_upper - gamma_lower = {delta:+.2f} mN/m"
        if status == "FAIL" and u_delta is not None and u_delta >= abs(delta):
            status, reason = "INSUFFICIENT SAMPLING", reason + f" within its uncertainty {u_delta:.2f} mN/m"
        flagged = [name for name, g in (("upper", g_up), ("lower", g_lo)) if abs(g) > thresholds["leaflet_flag"]]
        if flagged:
            status = qc.worst([status, "WARNING"])
            reason += "; |gamma| > " + f"{thresholds['leaflet_flag']:g} mN/m in the " + " and ".join(flagged) + " leaflet"
    return {"gamma_upper_mn_per_m": g_up, "gamma_lower_mn_per_m": g_lo, "delta_mn_per_m": float(delta),
            "uncertainty_upper_mn_per_m": u_up, "uncertainty_lower_mn_per_m": u_lo, "uncertainty_delta_mn_per_m": u_delta,
            "dz_nm": dz, "n_bins": int(len(z)), "status": status, "reason": reason}
