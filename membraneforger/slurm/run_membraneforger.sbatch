#!/bin/bash
# MembraneForger on Slurm: one array task per all-atom complex -> bundled Martini 3 membrane -> validated,
# energy-minimized all-atom CHARMM36 system in <name>/em.gro, next to the input PDB.
#
#   1. once, on a login node:   bash /path/to/MembraneForger/membraneforger/slurm/setup.sh
#   2. copy this file next to your PDB files and fill in blocks 1 to 3 below
#   3. from that folder:         sbatch run_membraneforger.sbatch
#
# Submitted as is, it builds the bundled example complex (GCGR) as a check of the installation.
# The last line of each task's log, and of <name>/membranebuilder.log, says PASS or FAIL.
# Any #SBATCH line can also be overridden at submission time: sbatch -p PARTITION -A ACCOUNT run_membraneforger.sbatch
#
# ============================ 1. SCHEDULER: your cluster's queue ==============================================
#SBATCH --job-name=mforger
#SBATCH --output=%x_%A_%a.log         # one log per task, written in the folder you submit from
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8             # cores for backmapping (OpenMM) and GROMACS; about 15-30 min per build on 8
#SBATCH --mem=16G                     # peak RSS of a 270k-atom GPCR build is about 1.5 GB (was 70G); 16G leaves room for larger systems
#SBATCH --time=06:00:00               # per task; add about 30 minutes per extra run when RUNS > 1 (block 3)
#SBATCH --array=0-0                   # 0-(number of NAMES - 1): one task per entry of NAMES in block 3
#SBATCH --mail-type=NONE
##SBATCH --partition=PARTITION        # uncomment and fill in if your cluster has no usable default partition
##SBATCH --account=ACCOUNT            # uncomment and fill in if your cluster bills jobs to an account
##SBATCH --qos=QOS                    # likewise for a quality-of-service name
# ==============================================================================================================

set -uo pipefail

# ============================ 2. SOFTWARE: where things are on your cluster ===================================
REPO="${MFORGER_REPO:-$HOME/MembraneForger}"   # the cloned MembraneForger repository, where setup.sh ran
MODULES=""                                     # modules to load first, space-separated; "" if none.
                                               # e.g. "Gromacs/2026.3-Container" on Gemini
GMX=""                                         # GROMACS command. "" = $REPO/env/bin/gmx from setup.sh if present,
                                               # otherwise gmx from a module. A container also works:
                                               # "singularity exec /path/image.sif gmx"
PPM="${MFORGER_PPM:-$REPO/ppm3_code/immers}"   # PPM 3.0 executable, res.lib next to it; only used for orientation
PY=""                                          # Python with the requirements installed. "" = $REPO/env/bin/python
                                               # from setup.sh; otherwise e.g. "python3" after activating your own
                                               # environment in MODULES or below
# ==============================================================================================================

# ============================ 3. THIS RUN: your structures and build options ==================================
WORK="${SLURM_SUBMIT_DIR:-$PWD}"               # folder with the input PDBs; results go to $WORK/<name>/
NAMES=(GCGR)                                   # $WORK/<name>.pdb for each entry, e.g. (receptor_A receptor_B);
                                               # keep --array in block 1 in step. GCGR is the bundled example
                                               # complex and is copied here if the file is absent
ORIENT="--orientation none"                    # already oriented (OPM, CHARMM-GUI, the bundled GCGR): as is.
                                               # a model or any structure not yet in a membrane: name the chain
                                               # that crosses the membrane and the side of its N terminus, e.g.
                                               # "--orient-chain R --nterm-side out" (GPCRs: out); uses PPM
MEMBRANE=1                                     # 1 = a random bundled GPR139 frame, 2 = a random KOR frame, or
                                               # one frame by name for every model (python -m membraneforger --list-membranes)
RUNS=1                                         # > 1 with MEMBRANE=2 (GPR139 has one frame): build each model many times, once
                                               # in each of the first RUNS frames ($WORK/<name>/<frame>/em.gro);
                                               # the runs are sequential, about 15 to 30 minutes each: raise --time
EXTRA=""                                       # other options, e.g. "--xy-buffer 2.0" for production systems
# ==============================================================================================================

task="${SLURM_ARRAY_TASK_ID:-0}"
die() { echo "ERROR: $*"; exit 1; }

if [ -n "$MODULES" ]; then
    command -v module >/dev/null 2>&1 || die "'module' command not available but MODULES=\"$MODULES\" is set"
    module load $MODULES || die "module load $MODULES failed"
fi

# Resolve the "" defaults of block 2.
[ -n "$PY" ] || PY="$REPO/env/bin/python"
if [ -z "$GMX" ]; then
    if [ -x "$REPO/env/bin/gmx" ]; then GMX="$REPO/env/bin/gmx"; else GMX=gmx; fi
fi

# Fail in seconds, not after an hour of backmapping, when something is missing.
[ -d "$REPO/membraneforger" ] || die "REPO=$REPO does not contain the membraneforger package (block 2)"
command -v "${PY%% *}" >/dev/null 2>&1 || [ -x "$PY" ] \
    || die "Python '$PY' not found: run 'bash $REPO/membraneforger/slurm/setup.sh' on a login node first, or set PY (block 2)"
"$PY" -c "import numpy, scipy, networkx, pandas, openmm, mstool" 2>/dev/null \
    || die "Python '$PY' lacks the requirements (numpy scipy networkx pandas openmm mstool): see membraneforger/slurm/README.md"
if [ "$task" -ge "${#NAMES[@]}" ]; then
    echo "task $task: NAMES has only ${#NAMES[@]} entries, nothing to do"
    exit 0
fi
if [ -n "${SLURM_ARRAY_TASK_MAX:-}" ] && [ "$SLURM_ARRAY_TASK_MAX" -lt $(( ${#NAMES[@]} - 1 )) ] \
        && [ "$task" -eq "${SLURM_ARRAY_TASK_MIN:-0}" ]; then
    echo "WARNING: NAMES has ${#NAMES[@]} entries but --array stops at $SLURM_ARRAY_TASK_MAX;" \
         "the rest will not run (use --array=0-$(( ${#NAMES[@]} - 1 )))"
fi
name="${NAMES[$task]%.pdb}"
pdb="$WORK/$name.pdb"
if [ ! -f "$pdb" ] && [ "$name" = GCGR ]; then
    cp "$REPO/membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_prot-lig.pdb" "$pdb" && echo "copied the bundled example to $pdb"
fi
[ -f "$pdb" ] || die "input $pdb not found (NAMES, block 3)"
version=$($GMX --version 2>/dev/null | grep "GROMACS version") || die "GROMACS command '$GMX' does not run (block 2)"

read -ra orient <<<"$ORIENT"
read -ra extra <<<"$EXTRA"
args=(--aa "$pdb" --cg "$MEMBRANE" --out "$WORK/$name" --gmx "$GMX" --ntomp "${SLURM_CPUS_PER_TASK:-8}")
[ "${RUNS:-1}" -gt 1 ] && args+=(--runs "$RUNS")
if [[ " $ORIENT $EXTRA " != *" --orientation none "* && " $ORIENT $EXTRA " != *" --orientation opm "* ]]; then
    [ -x "$PPM" ] || die "PPM 3.0 executable $PPM not found (block 2; see $REPO/membraneforger/slurm/README.md)"
    [ -f "$(dirname "$PPM")/res.lib" ] || die "res.lib must sit next to $PPM"
    args+=(--ppm-exe "$PPM")
fi

echo "task $task $name host $(hostname) $(date '+%F %T')"
echo "$version"
cd "$REPO" || exit 1
"$PY" -m membraneforger "${args[@]}" ${orient[@]+"${orient[@]}"} ${extra[@]+"${extra[@]}"}
status=$?

if [ "${RUNS:-1}" -gt 1 ]; then
    for log in "$WORK/$name"/*/membranebuilder.log; do [ -f "$log" ] && echo "$(basename "$(dirname "$log")"): $(tail -n 1 "$log")"; done
else
    tail -n 1 "$WORK/$name/membranebuilder.log" 2>/dev/null
fi
echo "exit $status $(date '+%F %T')"
exit $status
