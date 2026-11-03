#!/bin/bash
# Equilibrate a MembraneForger build with the bundled step6.1-6.6 schedule and prepare the production tpr.
#
#   bash membraneforger/mdp/run_equilibration.sh BUILD_DIR [GMX] [NTOMP]
#
# BUILD_DIR is the --out of a build that passed (em.gro, topol.top, index_ini.ndx). GMX is the GROMACS command (default
# gmx) and NTOMP the OpenMP threads (default 8). Every step restrains to the positions of em.gro (-r em.gro), as the
# CHARMM-GUI schedule does, and starts from the structure the previous step wrote. The six steps are 1.875 ns in all.
# The run leaves step6.*.gro/.cpt/.xtc/.edr/.log in BUILD_DIR; a later rebuild into the same folder is refused by
# MembraneForger unless --overwrite --overwrite-restart is given, and never deletes these files.
# See membraneforger/docs/equilibration.md for what each step does and how to grade the result.

set -euo pipefail

BUILD="${1:?usage: run_equilibration.sh BUILD_DIR [GMX] [NTOMP]}"
GMX="${2:-gmx}"
NTOMP="${3:-8}"
MDP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd "$BUILD"
for f in em.gro topol.top index_ini.ndx; do
    [ -f "$f" ] || { echo "ERROR: $BUILD has no $f: equilibrate a build that passed (its last log line says PASS)"; exit 1; }
done

previous=em.gro
for step in 6.1 6.2 6.3 6.4 6.5 6.6; do
    name="step${step}_equilibration"
    echo "== $name ($(date '+%F %T'))"
    $GMX grompp -f "$MDP/$name.mdp" -c "$previous" -r em.gro -p topol.top -n index_ini.ndx \
        -o "$name.tpr" -po "${name}_mdout.mdp" -maxwarn 0
    $GMX mdrun -deffnm "$name" -ntomp "$NTOMP"
    previous="$name.gro"
done

echo "== step7_production: grompp only ($(date '+%F %T'))"
$GMX grompp -f "$MDP/step7_production.mdp" -c step6.6_equilibration.gro -t step6.6_equilibration.cpt \
    -p topol.top -n index_ini.ndx -o step7_production.tpr -po step7_production_mdout.mdp -maxwarn 0
echo "equilibration done; start production with:"
echo "    $GMX mdrun -deffnm step7_production -ntomp $NTOMP"
