# Running MembraneForger on a Slurm cluster

Two files: `setup.sh` (once per user, on a login node) and `run_membraneforger.sbatch` (one array task per structure).
A third, `test_contact_gates.sbatch`, repeats a failed build with a branch checked out as a worktree beside the
repository (edit `REPO`, `BRANCH`, `PREVIOUS`; on Gemini it loads `Gromacs/2026.3-Container`), runs the unit tests
first and prints the new build's result and its `contacts.json` summary.

```bash
git clone https://github.com/rolanddelmundo/MembraneForger.git ~/MembraneForger
bash ~/MembraneForger/slurm/setup.sh          # about 5 minutes; needs internet, no Anaconda or modules

mkdir ~/my_project && cd ~/my_project         # put your PDB files here
cp ~/MembraneForger/slurm/run_membraneforger.sbatch .
# edit the USER-DEFINED block: NAMES, ORIENT, and --array to match the number of NAMES
sbatch run_membraneforger.sbatch
```

`setup.sh` downloads micromamba (a single binary) into `.micromamba/` and builds `env/` in the repository from
conda-forge: Python 3.11, numpy, scipy, networkx, pandas < 3, OpenMM, GROMACS (CPU build, picks AVX2/AVX-512 at run
time), gfortran and make, then mstool 0.3.9 from PyPI. The sbatch script calls `env/bin/python` and `env/bin/gmx`
directly, so nothing has to be activated and your shell start-up files are not touched. To use another GROMACS
(a module, or a container such as `singularity exec image.sif gmx`), set `GMX` in the sbatch script.

Each task writes `<name>/` next to the PDB files and a log `mforger_<job>_<task>.log` in the submit folder. The
script checks the environment, the input file, GROMACS and PPM before it starts, so a mistake fails in seconds.

## Several membranes for one structure

`RUNS=5` in the sbatch script builds each structure five times, once in each of the first five frames of the chosen
membrane series (`--runs 5`): `<name>/GPR1/em.gro` to `<name>/GPR5/em.gro`, five replicate systems that differ in
their lipid arrangement. The runs are sequential within the task, so raise `--time` accordingly (about 1 hour per run
on 8 cores). To run them in parallel instead, give each array task its own frame by name
(`MEMBRANE=GPR1`, `GPR2`, ...; `python -m membraneforger --list-membranes` lists the 18).

On a cluster whose GROMACS is a module (Gemini: `Gromacs/2026.3-Container`), set `MODULES="Gromacs/2026.3-Container"`
and `GMX=gmx` in the USER-DEFINED block; the Python environment from `setup.sh` is still used for everything else.

## PPM 3.0

A structure that is not yet oriented in a membrane (no OPM entry, e.g. a model) is oriented with PPM 3.0, the
OPM team's program. It is not redistributed with MembraneForger.

1. Get the PPM 3.0 source folder `ppm3_code/` (it contains a Makefile and `res.lib`) from the OPM team, or from a
   lab member who already has it.
2. Copy it to `~/MembraneForger/ppm3_code/` and run `bash slurm/setup.sh` again: it compiles `immers` with the
   environment's gfortran. For a source folder elsewhere: `PPM_SRC=/path/to/ppm3_code bash slurm/setup.sh`.
3. The default `PPM="$REPO/ppm3_code/immers"` in the sbatch script then works. A lab can share one compiled copy:
   point `PPM` at it (`res.lib` must stay next to `immers`).

PPM is not needed when the structure is already oriented (membrane normal along z, as from OPM or CHARMM-GUI): set
`ORIENT="--orientation none"`, or use `--pdb-id ID` for a deposited structure with an OPM entry. Another option is to
orient the complex on the PPM 3.0 web server (https://opm.phar.umich.edu/ppm_server3), remove the `DUM` atoms from
the result, check that all chains and hydrogens are still there, and run with `--orientation none`.
