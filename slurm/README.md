# Running MembraneForger on a Slurm cluster

| File | What it is for | Run it |
|---|---|---|
| `setup.sh` | Builds a self-contained Python environment (plus a CPU GROMACS and gfortran) inside the repository, and compiles PPM 3.0 if its source is present. No Anaconda, modules or root needed. | once per user, on a login node |
| `run_membraneforger.sbatch` | Builds one all-atom complex per array task in a bundled membrane: orientation, Martini 3 slice, backmapping, validation, energy minimization. Writes `<name>/em.gro` next to the input PDB. | copy next to your PDB files, fill in the three blocks, `sbatch` it |
| `examples/gemini_6WHC.sbatch` | The same script with the blocks filled in for the Gemini cluster (City of Hope): builds the bundled 6WHC example in five membranes. A worked example and an installation check. | `sbatch` as is on Gemini; adapt on another cluster |

```bash
git clone https://github.com/rolanddelmundo/MembraneForger.git ~/MembraneForger
bash ~/MembraneForger/slurm/setup.sh          # about 5 minutes; needs internet

mkdir ~/my_project && cd ~/my_project         # put your PDB files here
cp ~/MembraneForger/slurm/run_membraneforger.sbatch .
# fill in blocks 1 to 3 (below)
sbatch run_membraneforger.sbatch
```

Submitted unchanged, the script builds the bundled example complex (6WHC), so the first submission checks the
installation. Each task writes `<name>/` next to the PDB files and a log `mforger_<job>_<task>.log` in the submit
folder; the last line of the log, and of `<name>/membranebuilder.log`, says `PASS` or `FAIL`. The script checks the
environment, the input file, GROMACS and PPM before it starts, so a mistake fails in seconds and the error names the
block to fix.

## Filling in the script

The script has three marked blocks at the top. Nothing below them needs editing.

**Block 1, scheduler** (`#SBATCH` lines): what your cluster's queue needs. The resources are set for one build of a
GPCR complex (8 cores, 70 GB, about 1 hour per run). Add a partition, account or QOS by uncommenting the last lines
of the block if your cluster requires them. `--array=0-(number of NAMES - 1)` gives one task per structure; raise
`--time` by about 1 hour per extra run when `RUNS > 1`. Any `#SBATCH` line can also be overridden at submission time
without editing the file:

```bash
sbatch -p mypartition -A mylab --array=0-3 run_membraneforger.sbatch
```

**Block 2, software**: where things are on your cluster.

| Variable | Meaning | Default |
|---|---|---|
| `REPO` | the cloned repository (where `setup.sh` ran) | `~/MembraneForger`, or `$MFORGER_REPO` |
| `MODULES` | modules to load first, space-separated | none |
| `GMX` | GROMACS command | `$REPO/env/bin/gmx` from `setup.sh` if it exists, otherwise `gmx` from a module. A container also works: `"singularity exec image.sif gmx"` |
| `PPM` | PPM 3.0 executable, `res.lib` next to it; used only when a structure must be oriented | `$REPO/ppm3_code/immers`, or `$MFORGER_PPM` |
| `PY` | Python with the requirements installed | `$REPO/env/bin/python` from `setup.sh` |

Typical setups:

- No modules, nothing installed: run `setup.sh` and leave the block as is.
- GROMACS is a module (Gemini: `Gromacs/2026.3-Container`): `MODULES="Gromacs/2026.3-Container"` and `GMX="gmx"`;
  the Python environment from `setup.sh` is still used for everything else.
- Your own Python environment: load or activate it in `MODULES` (a module) or right after the block (a
  `source .../activate` line), and set `PY="python3"`. It needs numpy, scipy, networkx, pandas < 3, OpenMM, mstool
  0.3.9 and matplotlib, as in the top-level README.

`setup.sh` itself takes two optional variables: `MFORGER_ENV=/path` builds the environment elsewhere than
`<repo>/env`, and `PPM_SRC=/path/to/ppm3_code` compiles PPM from a source folder elsewhere than `<repo>/ppm3_code`.

**Block 3, this run**: your structures and build options.

| Variable | Meaning |
|---|---|
| `NAMES` | the PDB files in the submit folder, without `.pdb`: `(receptor_A receptor_B)`. Keep `--array` in block 1 in step |
| `ORIENT` | `"--orientation none"` for a structure already placed in a membrane (OPM, CHARMM-GUI). Otherwise the chain that crosses the membrane and the side of its N terminus, `"--orient-chain R --nterm-side out"` (GPCRs: out); this uses PPM |
| `MEMBRANE` | `1` for a random GPR139 frame, `2` for a random KOR frame, or one frame by name (`python -m membraneforger --list-membranes` lists them) |
| `RUNS` | build each structure this many times, once in each of the first `RUNS` frames of the series (`MEMBRANE=2`: `<name>/KOR1/em.gro`, `KOR2`, `KOR3`, `KOR5`, `KOR6`; the GPR139 series has one frame, so `RUNS` above 1 needs `MEMBRANE=2`). Sequential, about 1 hour each |
| `EXTRA` | any other option from the [command-line reference](../docs/cli_reference.md), for example `"--xy-buffer 2.0"` for production box sizes |

`RUNS=5` gives five replicate systems that differ in their lipid arrangement. To build them in parallel instead,
give each array task its own frame by name (`MEMBRANE=KOR1`, `KOR2`, ...; `python -m membraneforger --list-membranes`
lists the 8).

## PPM 3.0

A structure that is not yet oriented in a membrane (no OPM entry, e.g. a model) is oriented with PPM 3.0, the
OPM team's program. It is not redistributed with MembraneForger.

1. Get the PPM 3.0 source folder `ppm3_code/` (it contains a Makefile and `res.lib`) from the OPM team, or from a
   lab member who already has it.
2. Copy it to `~/MembraneForger/ppm3_code/` and run `bash slurm/setup.sh` again: it compiles `immers` with the
   environment's gfortran. For a source folder elsewhere: `PPM_SRC=/path/to/ppm3_code bash slurm/setup.sh`.
3. The default `PPM` in block 2 then works. A lab can share one compiled copy: point `PPM` at it (`res.lib` must
   stay next to `immers`).

PPM is not needed when the structure is already oriented (membrane normal along z, as from OPM or CHARMM-GUI): set
`ORIENT="--orientation none"`, or use `--pdb-id ID` in `EXTRA` for a deposited structure with an OPM entry. Another
option is to orient the complex on the PPM 3.0 web server (https://opm.phar.umich.edu/ppm_server3), remove the `DUM`
atoms from the result, check that all chains and hydrogens are still there, and run with `--orientation none`.
