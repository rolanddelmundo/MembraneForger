#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// One record format and one classification rule for every validation metric of the build.
#//=============================================================
"""The common currency of the validation gates: a metric record and its PASS / WARNING / FAIL classification.

Every analysis (packing, RDF, thickness, hydration, convergence, diffusion, ...) ends in records of this shape so
that the report tables, the stage matrix and the overall verdict can be assembled without knowing the analysis:

    metric(name, stage, value, units, leaflet=..., reference=..., uncertainty=..., n=..., deviation=..., status=..., reason=...)

Statuses: PASS, WARNING, FAIL, INSUFFICIENT SAMPLING (the data cannot decide), REFERENCE NEEDED (no matched control
to judge against: the value is reported, never graded), NOT RUN. Poor sampling is never turned into a PASS or a FAIL.

classify() grades an absolute deviation against a (pass, warning) ceiling pair; classify_percentile() grades a value
against the empirical distribution of a matched control (central 95 % PASS, 95-99 % WARNING, beyond FAIL).
overall_status() reduces a list of records to one of PASS, PASS WITH WARNINGS, FAIL or INSUFFICIENT SAMPLING.
"""
import numpy as np

__all__ = ['STATUSES', 'metric', 'classify', 'classify_percentile', 'overall_status', 'worst', 'records_table']

STATUSES = ("PASS", "WARNING", "FAIL", "INSUFFICIENT SAMPLING", "REFERENCE NEEDED", "NOT RUN")
RANK = {"PASS": 0, "NOT RUN": 0, "REFERENCE NEEDED": 1, "WARNING": 2, "INSUFFICIENT SAMPLING": 3, "FAIL": 4}


def metric(name: str, stage: str, value, units: str, leaflet: str | None = None, reference=None, uncertainty=None,
           n: int | None = None, deviation=None, status: str = "NOT RUN", reason: str = "") -> dict:
    """One validation metric with everything a reader needs to judge it."""
    if status not in STATUSES:
        raise ValueError(f"unknown status {status}")
    clean = lambda v: None if v is None else (float(v) if isinstance(v, (int, float, np.floating, np.integer)) and not isinstance(v, bool) else v)
    return {"metric": name, "stage": stage, "leaflet": leaflet, "value": clean(value), "units": units, "reference": clean(reference),
            "uncertainty": clean(uncertainty), "n": None if n is None else int(n), "deviation": clean(deviation),
            "status": status, "reason": reason}


def classify(deviation, pass_max: float, warning_max: float, absolute: bool = True) -> str:
    """PASS when |deviation| <= pass_max, WARNING up to warning_max, FAIL beyond; None is INSUFFICIENT SAMPLING."""
    if deviation is None or (isinstance(deviation, float) and not np.isfinite(deviation)):
        return "INSUFFICIENT SAMPLING"
    d = abs(deviation) if absolute else deviation
    return "PASS" if d <= pass_max else "WARNING" if d <= warning_max else "FAIL"


def classify_percentile(value: float, control: np.ndarray, minimum: int = 20) -> tuple[str, float | None]:
    """Grade a value against a control distribution: inside its central 95 % PASS, 95-99 % WARNING, beyond FAIL."""
    # Returns the status and the two-sided empirical percentile of the value (0.5 = median, 0 or 1 = extreme).
    control = np.asarray(control, dtype=float)
    control = control[np.isfinite(control)]
    if len(control) < minimum or value is None or not np.isfinite(value):
        return "INSUFFICIENT SAMPLING", None
    below = float((control < value).mean())
    two_sided = 2.0 * min(below, 1.0 - below)  # 1 at the median, 0 beyond the extremes
    return ("PASS" if two_sided > 0.05 else "WARNING" if two_sided > 0.01 else "FAIL"), round(1.0 - two_sided, 4)


def worst(statuses: list[str]) -> str:
    """The most severe of several statuses."""
    return max(statuses, key=lambda s: RANK[s]) if statuses else "NOT RUN"


def overall_status(records: list[dict], required: tuple = ()) -> dict:
    """Reduce metric records to one verdict: PASS, PASS WITH WARNINGS, FAIL or INSUFFICIENT SAMPLING."""
    # `required` names metrics that must be present and graded for a plain PASS (an absent or ungraded required
    # metric gives INSUFFICIENT SAMPLING at best). FAIL anywhere is FAIL; otherwise poor sampling of any graded
    # metric is INSUFFICIENT SAMPLING; otherwise any WARNING or REFERENCE NEEDED is PASS WITH WARNINGS.
    statuses = [r["status"] for r in records]
    present = {r["metric"] for r in records if r["status"] not in ("NOT RUN",)}
    missing = [name for name in required if name not in present]
    if "FAIL" in statuses:
        verdict = "FAIL"
    elif "INSUFFICIENT SAMPLING" in statuses or missing:
        verdict = "INSUFFICIENT SAMPLING"
    elif "WARNING" in statuses or "REFERENCE NEEDED" in statuses:
        verdict = "PASS WITH WARNINGS"
    else:
        verdict = "PASS"
    counts = {s: statuses.count(s) for s in STATUSES if statuses.count(s)}
    where = lambda r: f"{r['metric']} ({r['stage']}{', ' + r['leaflet'] if r['leaflet'] else ''}): {r['reason']}"
    return {"status": verdict, "counts": counts, "missing_required": missing,
            "failures": [where(r) for r in records if r["status"] == "FAIL"],
            "warnings": [where(r) for r in records if r["status"] in ("WARNING", "REFERENCE NEEDED", "INSUFFICIENT SAMPLING")]}


def records_table(records: list[dict]) -> str:
    """Metric records as a Markdown table."""
    cell = lambda v, d=3: "-" if v is None else (f"{v:.{d}g}" if isinstance(v, float) else str(v))
    head = ("| Metric | Stage | Leaflet | Value | Units | Reference | Uncertainty | N | Deviation | Status | Reason |\n"
            "|---|---|---|---:|---|---:|---:|---:|---:|---|---|\n")
    rows = [f"| {r['metric']} | {r['stage']} | {r['leaflet'] or '-'} | {cell(r['value'], 4)} | {r['units']} | {cell(r['reference'], 4)} | "
            f"{cell(r['uncertainty'])} | {cell(r['n'])} | {cell(r['deviation'])} | {r['status']} | {r['reason']} |" for r in records]
    return head + "\n".join(rows) + "\n"
