#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from membraneforger.cli import cmd_stage


def configure_local_resources() -> None:
    threads = os.environ.get("MEMBRANEFORGER_LOCAL_THREADS") or str(os.cpu_count() or 1)
    for key in (
        "OMP_NUM_THREADS",
        "OPENMM_CPU_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[key] = threads
    requested_platform = os.environ.get("MEMBRANEFORGER_OPENMM_PLATFORM")
    if requested_platform:
        os.environ["OPENMM_DEFAULT_PLATFORM"] = requested_platform
        return
    try:
        from openmm import Platform

        names = {Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())}
        if "OpenCL" in names:
            os.environ["OPENMM_DEFAULT_PLATFORM"] = "OpenCL"
    except Exception:
        pass


if __name__ == "__main__":
    configure_local_resources()
    import argparse

    parser = argparse.ArgumentParser(description="Run MembraneForger Stage 3 backmapping")
    parser.add_argument("--config")
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--mode")
    parser.add_argument("--accept-review", action="store_true")
    args = parser.parse_args(sys.argv[1:])
    args.stage = "stage3"
    raise SystemExit(cmd_stage(args))
