#!/bin/bash
#SBATCH --job-name=mforger
#SBATCH --mail-type=none
#SBATCH -n 1
#SBATCH -c 8                          # cores for backmapping (OpenMM) and gmx
#SBATCH -N 1
#SBATCH -p fast,all                   # partition(s): change to yours
#SBATCH --mem=70G
#SBATCH --time=06:00:00
#SBATCH --array=0-1                   # one task per entry of NAMES below: 0-(number of names - 1)
#SBATCH --output=%x_%A_%a.log         # written in the folder you submit from
#
# MembraneForger: each all-atom complex in NAMES -> bundled Martini 3 membrane -> validated, energy-minimized
# all-atom CHARMM36 system in $WORK/<name>/em.gro.
#
#   1. once, on a login node:   bash /path/to/MembraneForger/slurm/setup.sh
#   2. copy this file next to your PDB files and edit the USER-DEFINED block
#   3. from that folder:         sbatch run_membraneforger.sbatch
#
# The last line of each task's log, and of <name>/membranebuilder.log, says PASS or FAIL.

set -uo pipefail

# ================================== USER-DEFINED ==================================
WORK="${SLURM_SUBMIT_DIR:-$PWD}"           # folder with the input PDBs; results go to $WORK/<name>/
REPO="$HOME/MembraneForger"                 # the cloned MembraneForger repository (where setup.sh ran)
GMX="$REPO/env/bin/gmx"                     # GROMACS command (setup.sh installs one here)
PPM="$REPO/ppm3_code/immers"                # PPM 3.0 executable with res.lib next to it
NAMES=(best_model_clean best_model_core_engaged_clean)   # $WORK/<name>.pdb; keep --array above in step

ORIENT="--orient-chain R --nterm-side out"  # chain that crosses the membrane, and the side of its N terminus
                                            # (GPCRs: out). Already oriented (OPM/CHARMM-GUI): "--orientation none"
MEMBRANE=1                                  # 1 = a random bundled GPR139 frame, 2 = a random KOR frame. Same frame
                                            # for every model: "$REPO/examples/preeq_cg_cellmem/GPR1_cg_cellmem.gro"
                                            # and add --embed to EXTRA
EXTRA=""                                    # other options, e.g. "--xy-buffer 2.0" for production systems
# ==================================================================================

PY="$REPO/env/bin/python"
task="${SLURM_ARRAY_TASK_ID:-0}"
die() { echo "ERROR: $*"; exit 1; }

# Fail in seconds, not after an hour of backmapping, when something is missing.
[ -x "$PY" ] || die "$PY not found: run 'bash $REPO/slurm/setup.sh' on a login node first (or fix REPO)"
[ -d "$REPO/membraneforger" ] || die "REPO=$REPO does not contain the membraneforger package"
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
[ -f "$pdb" ] || die "input $pdb not found"
version=$($GMX --version 2>/dev/null | grep "GROMACS version") || die "GROMACS command '$GMX' does not run"

read -ra orient <<<"$ORIENT"
read -ra extra <<<"$EXTRA"
args=(--aa "$pdb" --cg "$MEMBRANE" --out "$WORK/$name" --gmx "$GMX" --ntomp "${SLURM_CPUS_PER_TASK:-8}")
if [[ " $ORIENT $EXTRA " != *" --orientation none "* && " $ORIENT $EXTRA " != *" --orientation opm "* ]]; then
    [ -x "$PPM" ] || die "PPM 3.0 executable $PPM not found (see $REPO/slurm/README.md)"
    [ -f "$(dirname "$PPM")/res.lib" ] || die "res.lib must sit next to $PPM"
    args+=(--ppm-exe "$PPM")
fi

echo "task $task $name host $(hostname) $(date '+%F %T')"
echo "$version"
cd "$REPO" || exit 1
"$PY" -m membraneforger "${args[@]}" ${orient[@]+"${orient[@]}"} ${extra[@]+"${extra[@]}"}
status=$?

tail -n 1 "$WORK/$name/membranebuilder.log" 2>/dev/null
echo "exit $status $(date '+%F %T')"
exit $status
