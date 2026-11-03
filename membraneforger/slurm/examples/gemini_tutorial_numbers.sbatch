#!/bin/bash
# Re-measure the frame-dependent numbers of the tutorial on every bundled membrane, on Gemini (City of Hope).
# membraneforger/examples/measure_tutorial_numbers.py runs the coarse-grained stages only (embedding, slicing, the slice
# gate; nothing is backmapped), so it needs no GROMACS, but embedding runs through mstool. Tasks 0-17 each measure one bundled
# frame into the shared OUT folder; task 18 reads every finished case back (nothing is re-run) and writes OUT/summary.md and summary.json.
#
#   checkcpu                                       # pick the shortest CPU queue, never overflow; set -q below to it
#   git clone -b feat/bundle-18-cg-frames https://github.com/rolanddelmundo/MembraneForger.git /scratch/rdelmundo/mforger_tutorial18/MembraneForger
#   cd /scratch/rdelmundo/mforger_tutorial18
#   cp MembraneForger/membraneforger/slurm/examples/gemini_tutorial_numbers.sbatch .
#   jid=$(sbatch --parsable --array=0-17 gemini_tutorial_numbers.sbatch)
#   sbatch --array=18 --dependency=afterany:${jid} gemini_tutorial_numbers.sbatch
#
# Expected: each tutorial18_<job>_<n>.log ends with "exit 0" (helix_hole is the one case a gate is expected to stop), and
# tutorial18_<job>_18.log ends with "wrote .../summary.md" and "exit 0". A task that ran out of time can be resubmitted
# with the same --array index: finished cases are skipped.
#
#SBATCH --job-name=tutorial18
#SBATCH --output=%x_%A_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00               # four coarse-grained cases on one frame take minutes; 4 h leaves a wide margin
#SBATCH -p compute
#SBATCH -q compute-short              # from checkcpu on 2026-11-03 (8 h limit, nothing competing); not overflow
#SBATCH --mail-type=NONE

set -uo pipefail

REPO="${MFORGER_REPO:-/scratch/rdelmundo/mforger_tutorial18/MembraneForger}"   # the clone on the PR branch (above)
OUT="${MFORGER_OUT:-/scratch/rdelmundo/mforger_tutorial18/tutorial18_out}"
CONDA_ENV="mstool"                     # mstool plus numpy, scipy, networkx and pandas
NFRAMES=18

die() { echo "ERROR: $*"; exit 1; }
task="${SLURM_ARRAY_TASK_ID:-$NFRAMES}"

source "$HOME/anaconda3/etc/profile.d/conda.sh" && conda activate "$CONDA_ENV" || die "conda activate $CONDA_ENV failed"
[ -d "$REPO/membraneforger" ] || die "REPO=$REPO does not contain the membraneforger package: clone it first (see the header)"
cd "$REPO" || exit 1
mapfile -t frames < <(python -c "import membraneforger as mf; print('\n'.join(mf.bundled_frames()))")
[ "${#frames[@]}" -eq "$NFRAMES" ] || die "$REPO bundles ${#frames[@]} frames, not $NFRAMES: is it on feat/bundle-18-cg-frames (or a main that has it)?"

echo "task $task host $(hostname) $(date '+%F %T') commit $(git rev-parse --short HEAD)"
if [ "$task" -lt "$NFRAMES" ]; then
    echo "frame ${frames[$task]}"
    python membraneforger/examples/measure_tutorial_numbers.py --repo . --out "$OUT" --frames "${frames[$task]}"
else
    python membraneforger/examples/measure_tutorial_numbers.py --repo . --out "$OUT"   # every case cached: aggregates only
fi
status=$?
echo "exit $status $(date '+%F %T')"
exit $status
