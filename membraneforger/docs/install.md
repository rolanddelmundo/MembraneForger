# Install

<div class="mf-chips" markdown>
<span>About 5 minutes</span>
<span>Linux or macOS</span>
<span>Python 3.10 or later</span>
<span>GROMACS</span>
</div>

## What you need

- **Linux or macOS with Python 3.10 or later.** The Python packages (numpy, scipy, networkx, pandas < 3, openmm and
  mstool 0.3.9 or 0.3.10) are installed by `pip` below. matplotlib is optional and draws the validation figures.
- **GROMACS** (tested with 2023.3 and 2025.3). If it is not on your `PATH` as `gmx`, pass it with `--gmx`. On a
  cluster, load its module first.
- **PPM 3.0**, only if your structure is not already oriented and has no matching OPM entry
  (see [Advanced options](advanced.md#ppm-30)). The default orientation mode downloads the OPM entry of a deposited
  structure, so it needs network access or `--opm-file`.
- **Time and memory.** A GPCR complex takes about 15 to 30 minutes on 8 cores and about 1.5 GB of memory. `--ntomp`
  sets the threads (default: every core of the machine).

## Install with pip

The force field, the bundled membranes and the equilibration templates are part of the package.

```bash
git clone https://github.com/rolanddelmundo/MembraneForger.git
cd MembraneForger
pip install ".[figures]"          # or pip install -e ".[figures]" to keep running from the clone
membraneforger --version
```

After this, `membraneforger` and `python -m membraneforger` are the same command and work from any directory.

!!! tip "Python needs the real `gmx` executable"
    A shell alias or function is invisible to Python, so put the executable itself on `PATH` or pass
    `--gmx /path/to/gmx`. To point at PPM, set one variable, with `res.lib` in the same folder:
    `export MEMBRANEFORGER_PPM=/path/to/ppm3_code/immers`.

## Check the installation

The pre-flight check reads your inputs as the build will, checks GROMACS, mstool and PPM, and writes nothing. It takes
seconds and ends with `CHECK PASS` or `CHECK FAIL`.

```bash
membraneforger --check --aa membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_prot-lig.pdb --orientation none --out example_out
```

Fix every `FAIL` line. `WARN` lines don't stop a build.

## On a Slurm cluster

`bash membraneforger/slurm/setup.sh` builds a self-contained environment with GROMACS, and
`membraneforger/slurm/run_membraneforger.sbatch` runs one structure per array task. Fill in its three marked blocks
(scheduler, software, this run). See
[membraneforger/slurm/README.md](https://github.com/rolanddelmundo/MembraneForger/blob/main/membraneforger/slurm/README.md).

**Next:** [Your first build](first_build.md).
