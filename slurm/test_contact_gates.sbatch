#!/bin/bash
#SBATCH --job-name=mforger-gates
#SBATCH --mail-type=none
#SBATCH -n 1
#SBATCH -c 8                          # cores for backmapping (OpenMM) and gmx
#SBATCH -N 1
#SBATCH -p fast,all                   # partition(s): change to yours
#SBATCH --mem=70G
#SBATCH --time=06:00:00
#SBATCH --output=%x_%j.log            # written in the folder you submit from
#
# Test a MembraneForger branch on a build that failed before, without touching your main checkout:
# the branch is checked out as a git worktree next to the repository, the unit tests run there, then the failed
# build is repeated from its own run_manifest.json command into a new output folder.
#
#   1. edit the USER-DEFINED block (REPO, BRANCH, PREVIOUS)
#   2. from the folder that holds the failed output folder:   sbatch test_contact_gates.sbatch
#
# The log ends with the new build's PASS/FAIL line and the contact-gate summary (contacts.json by stage).

set -uo pipefail
# ---------------------------------------------------------------- USER-DEFINED
REPO="$HOME/MembraneForger"                 # the cloned repository (where slurm/setup.sh built env/)
BRANCH="claude/project-thread-oyatv0"       # the branch to test
PREVIOUS="${SLURM_SUBMIT_DIR:-$PWD}/test_out"   # the failed build: its run_manifest.json gives the command to repeat
OUT="${PREVIOUS}_gates"                     # where the repeated build goes (removed first if it exists)
GMX_MODULE="Gromacs/2026.3-Container"       # GROMACS on Gemini; "" to use $REPO/env/bin/gmx instead
RUN_UNIT_TESTS=1                            # 1 = run the unit tests first (about 5 minutes), 0 = skip
# ----------------------------------------------------------------------------
die() { echo "ERROR: $*"; exit 1; }
PY="$REPO/env/bin/python"
WORKTREE="$REPO-$(echo "$BRANCH" | tr '/' '-')"
[ -x "$PY" ] || die "$PY not found: run 'bash $REPO/slurm/setup.sh' on a login node first (or fix REPO)"
[ -f "$PREVIOUS/run_manifest.json" ] || die "$PREVIOUS/run_manifest.json not found: PREVIOUS must be a finished build folder"

echo "host $(hostname) $(date '+%F %T')"
# --- the branch, as a worktree beside the repository (your checkout and its env/ are left as they are)
git -C "$REPO" fetch origin "$BRANCH" || die "git fetch of $BRANCH failed"
if [ -d "$WORKTREE/.git" ] || [ -f "$WORKTREE/.git" ]; then
    git -C "$WORKTREE" checkout -q "$BRANCH" && git -C "$WORKTREE" reset -q --hard "origin/$BRANCH" || die "cannot update $WORKTREE"
else
    git -C "$REPO" worktree add "$WORKTREE" "$BRANCH" || die "cannot create the worktree $WORKTREE"
fi
echo "testing $(git -C "$WORKTREE" log --oneline -1) in $WORKTREE"
for link in env backmap_data ppm3_code; do   # untracked resources of the main checkout, shared by link
    [ -e "$REPO/$link" ] && [ ! -e "$WORKTREE/$link" ] && ln -s "$REPO/$link" "$WORKTREE/$link"
done

# --- GROMACS
if [ -n "$GMX_MODULE" ]; then
    module load "$GMX_MODULE" || die "module load $GMX_MODULE failed"
    GMX="gmx"
else
    GMX="$REPO/env/bin/gmx"
fi
$GMX --version 2>/dev/null | grep "GROMACS version" || die "GROMACS command '$GMX' does not run"

cd "$WORKTREE" || exit 1
# --- unit tests of the branch (the contact validator has its own: membraneforger/tests/test_contacts.py)
if [ "$RUN_UNIT_TESTS" = 1 ]; then
    echo "--- unit tests $(date '+%F %T')"
    "$PY" -m unittest discover -s membraneforger/tests -t . 2>&1 | grep -E "^(FAIL|ERROR):|^(Ran|OK|FAILED)" | tail -n 40
fi

# --- repeat the failed build with the branch, from the previous run's own command line
mapfile -t args < <("$PY" - "$PREVIOUS/run_manifest.json" "$OUT" "$GMX" <<'EOF'
import json, sys
manifest, out, gmx = json.loads(open(sys.argv[1]).read()), sys.argv[2], sys.argv[3]
argv = manifest["command"]
argv = argv[next(i for i, a in enumerate(argv) if a.startswith("--")):]  # drop the interpreter and the module path
for flag in ("--out", "--gmx"):
    while flag in argv:
        i = argv.index(flag)
        del argv[i:i + 2]
print("\n".join(argv + ["--out", out, "--gmx", gmx]))
EOF
)
[ "${#args[@]}" -gt 0 ] || die "could not read the previous command from $PREVIOUS/run_manifest.json"
rm -rf "$OUT"
echo "--- build $(date '+%F %T'): python -m membraneforger ${args[*]}"
"$PY" -m membraneforger "${args[@]}"
status=$?

echo "--- result"
tail -n 1 "$OUT/membranebuilder.log" 2>/dev/null
grep -E "contacts|solvate:|genion attempt|ions:|EM validated|first clashes|force-field overlaps" "$OUT/membranebuilder.log" 2>/dev/null | tail -n 12
if [ -f "$OUT/contacts.json" ]; then
    "$PY" - "$OUT/contacts.json" <<'EOF'
import json, sys
for stage, report in json.loads(open(sys.argv[1]).read()).items():
    print(f"contacts.json {stage}: {len(report['contacts'])} contact violation(s) {report['contacts_by_class']}, "
          f"{len(report['crossings'])} bond crossing(s) {report['crossings_by_class']}")
    for c in report["contacts"][:3]:
        print(f"    {c['atoms'][0]} - {c['atoms'][1]} {c['distance_A']} A (limit {c['limit_A']} A) [{c['class']}]")
    for c in report["crossings"][:3]:
        print(f"    bond {c['bonds'][0][0]}-{c['bonds'][0][1]} x bond {c['bonds'][1][0]}-{c['bonds'][1][1]} {c['distance_A']} A [{c['class']}]")
EOF
fi
echo "exit $status $(date '+%F %T')"
exit $status
