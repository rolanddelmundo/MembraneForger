#!/bin/bash
# One-time MembraneForger setup on a cluster login node (needs internet; no Anaconda, modules or root needed):
#
#   bash slurm/setup.sh                              # Python, mstool, OpenMM and GROMACS in <repo>/env
#   PPM_SRC=/path/to/ppm3_code bash slurm/setup.sh   # also compile PPM 3.0 (default source: <repo>/ppm3_code)
#
# Everything goes into the repository folder: micromamba (a single binary) in .micromamba/, the environment in env/.
# Nothing is activated or added to your shell start-up files; run_membraneforger.sbatch calls env/bin/python and
# env/bin/gmx directly. Running the script again is safe: it only adds what is missing.

set -euo pipefail

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ENV_DIR="${MFORGER_ENV:-$REPO/env}"
PPM_SRC="${PPM_SRC:-$REPO/ppm3_code}"
TOOLS="$REPO/.micromamba"
MICROMAMBA_URL=https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64

say() { echo "== $*"; }

# 1. micromamba: a stand-alone conda-forge installer
if [ ! -x "$TOOLS/micromamba" ]; then
    say "downloading micromamba to $TOOLS"
    mkdir -p "$TOOLS"
    curl -fsSL "$MICROMAMBA_URL" -o "$TOOLS/micromamba"
    chmod +x "$TOOLS/micromamba"
fi
export MAMBA_ROOT_PREFIX="$TOOLS/root"

# 2. the environment: Python, the scientific stack, OpenMM, GROMACS (CPU build) and gfortran for PPM
if [ ! -x "$ENV_DIR/bin/python" ]; then
    say "creating $ENV_DIR (a few minutes)"
    "$TOOLS/micromamba" create -y -p "$ENV_DIR" -c conda-forge --override-channels \
        "python=3.11" numpy scipy networkx "pandas<3" openmm "gromacs=*=nompi_*" gfortran make pip
fi
say "installing mstool"
"$ENV_DIR/bin/python" -m pip install --disable-pip-version-check -q "mstool==0.3.9" "pandas<3" matplotlib

# 3. check what the build needs
"$ENV_DIR/bin/python" -c "import numpy, scipy, networkx, pandas, openmm, mstool" \
    || { echo "the Python environment is incomplete; delete $ENV_DIR and run this script again"; exit 1; }
"$ENV_DIR/bin/gmx" --version 2>/dev/null | grep -E "GROMACS version|SIMD instructions"

# 4. PPM 3.0 (not redistributable, so not bundled): compile it when its source is here
if [ -x "$PPM_SRC/immers" ]; then
    say "PPM 3.0 already compiled: $PPM_SRC/immers"
elif [ -d "$PPM_SRC" ]; then
    say "compiling PPM 3.0 in $PPM_SRC"
    (cd "$PPM_SRC" && PATH="$ENV_DIR/bin:$PATH" make)
    [ -x "$PPM_SRC/immers" ] || { echo "make did not produce $PPM_SRC/immers"; exit 1; }
else
    say "PPM 3.0 source not found at $PPM_SRC: skipped (see slurm/README.md)"
fi
if [ -x "$PPM_SRC/immers" ] && [ ! -f "$PPM_SRC/res.lib" ]; then
    echo "WARNING: $PPM_SRC/res.lib is missing; PPM needs it next to immers"
fi

say "done. In run_membraneforger.sbatch set REPO=$REPO"
[ -x "$PPM_SRC/immers" ] && say "and PPM=$PPM_SRC/immers"
exit 0
