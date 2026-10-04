#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH --job-name=btop
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=36:00:01
#SBATCH --output=logs/build_topol/%x_%A_%a.log
#SBATCH --open-mode=truncate

set -euo pipefail

PYTHON=${PYTHON:-python3}

PROJECT_ROOT=${PROJECT_ROOT:-$PWD}
REPAIR_TOOL=${REPAIR_TOOL:-${PROJECT_ROOT}/md_scripts/fix_dyna_termini.py}
VMD_COMMAND=${VMD_COMMAND:-vmd}
GMX_COMMAND=${GMX_COMMAND:-${GMX_BIN:-gmx}}
TASK_ID=${SLURM_ARRAY_TASK_ID:-}
TASK_MIN=${TASK_MIN:-1}
TASK_MAX=${TASK_MAX:-130}
PREFLIGHT=0

usage() {
    cat <<'EOF'
usage: build_topol.sh --task-id N [options]

Wrapper for a site-provided topology repair tool. This script does not provide
CHARMM36 membrane parameters by itself.

Options:
  --task-id N             Task/index to pass to the repair tool.
  --project-root PATH     Directory containing step1_pdbreader.pdb/.psf.
  --repair-tool PATH      Python repair tool to execute.
  --vmd-command COMMAND   VMD command used by the repair tool.
  --gmx-command COMMAND   GROMACS command used by the repair tool.
  --task-min N            Minimum allowed task id (default: 1).
  --task-max N            Maximum allowed task id (default: 130).
  --preflight-only        Ask the repair tool to run its preflight mode.
  -h, --help              Show this help.
EOF
}

while (($#)); do
    case $1 in
        -h|--help)
            usage
            exit 0
            ;;
        --task-id)
            [[ $# -ge 2 ]] || { echo "ERROR: --task-id requires a value" >&2; exit 2; }
            TASK_ID=$2
            shift 2
            ;;
        --project-root)
            [[ $# -ge 2 ]] || { echo "ERROR: --project-root requires a value" >&2; exit 2; }
            PROJECT_ROOT=$2
            shift 2
            ;;
        --repair-tool)
            [[ $# -ge 2 ]] || { echo "ERROR: --repair-tool requires a value" >&2; exit 2; }
            REPAIR_TOOL=$2
            shift 2
            ;;
        --vmd-command)
            [[ $# -ge 2 ]] || { echo "ERROR: --vmd-command requires a value" >&2; exit 2; }
            VMD_COMMAND=$2
            shift 2
            ;;
        --gmx-command)
            [[ $# -ge 2 ]] || { echo "ERROR: --gmx-command requires a value" >&2; exit 2; }
            GMX_COMMAND=$2
            shift 2
            ;;
        --task-min)
            [[ $# -ge 2 ]] || { echo "ERROR: --task-min requires a value" >&2; exit 2; }
            TASK_MIN=$2
            shift 2
            ;;
        --task-max)
            [[ $# -ge 2 ]] || { echo "ERROR: --task-max requires a value" >&2; exit 2; }
            TASK_MAX=$2
            shift 2
            ;;
        --preflight-only)
            PREFLIGHT=1
            shift
            ;;
        *)
            echo "ERROR: unsupported argument: $1" >&2
            exit 2
            ;;
    esac
done

[[ ${TASK_MIN} =~ ^[0-9]+$ && ${TASK_MAX} =~ ^[0-9]+$ && ${TASK_ID} =~ ^[0-9]+$ ]] && ((TASK_ID >= TASK_MIN && TASK_ID <= TASK_MAX)) || {
    echo "ERROR: task ID must be ${TASK_MIN}..${TASK_MAX}" >&2
    exit 2
}
[[ -f ${REPAIR_TOOL} ]] || { echo "ERROR: missing repair tool: ${REPAIR_TOOL}" >&2; exit 2; }
[[ -f ${PROJECT_ROOT}/step1_pdbreader.pdb && -f ${PROJECT_ROOT}/step1_pdbreader.psf ]] || {
    echo "ERROR: missing root CHARMM-GUI convention reference" >&2
    exit 2
}

command=(
    "$PYTHON" "$REPAIR_TOOL"
    --task-id "$TASK_ID"
    --vmd-command "$VMD_COMMAND"
    --gmx-command "$GMX_COMMAND"
)
((PREFLIGHT == 0)) || command+=(--preflight-only)
"${command[@]}"
