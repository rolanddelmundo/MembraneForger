# MembraneForger

MembraneForger builds an all-atom membrane–protein system for GROMACS from a single structure file.

You give it an all-atom protein (or protein–ligand) PDB. It places the protein in a pre-equilibrated Martini 3
coarse-grained cell membrane, converts the membrane to all atoms with [mstool](https://github.com/ksy141/mstool), and
returns a solvated, energy-minimized CHARMM36 system (`em.gro`) that has passed a set of automated checks.

## How it works

1. **Orient** the protein in the membrane (from an OPM entry, PPM 3.0, or your own coordinates).
2. **Check** the coarse-grained membrane frame (bundled, or your own Martini 3 frame) before anything is placed in it.
3. **Embed** the protein in it.
4. **Cut** the membrane to the protein plus a margin, and check that the cut kept the lipid packing.
5. **Backmap** the lipids to all atoms around the fixed protein.
6. **Build** the CHARMM36 topology, add water and 0.15 M NaCl, and minimize.
7. **Audit** the result. `em.gro` is written only when every check passes.

The output is minimized, not equilibrated. Run restrained equilibration before production.

## Requirements

- Python 3.10 or later with numpy, scipy, networkx, pandas < 3, openmm and mstool 0.3.9 or 0.3.10.
  matplotlib is optional and draws the validation figures.
- GROMACS (tested with 2023.3 and 2025.3). If it is not on your `PATH` as `gmx`, pass it with `--gmx`.
- PPM 3.0, only if your structure is not already oriented and has no matching OPM entry
  (see [Advanced options](membraneforger/docs/advanced.md#ppm-30)).

## Installation

Clone the repository and install the Python packages. MembraneForger runs from the clone; there is nothing else to
install.

```bash
git clone https://github.com/rolanddelmundo/MembraneForger.git
cd MembraneForger
pip install numpy scipy networkx "pandas<3" openmm "mstool>=0.3.9,<=0.3.10" matplotlib
```

On a Slurm cluster, `bash membraneforger/slurm/setup.sh` builds a self-contained environment with GROMACS instead, and
`membraneforger/slurm/run_membraneforger.sbatch` runs one structure per array task; its three marked blocks (scheduler, software,
this run) are all that needs filling in for any cluster. See [membraneforger/slurm/README.md](membraneforger/slurm/README.md).

## Quick start

This builds the bundled example, a glucagon receptor–tirzepatide–Gs complex, in its own coarse-grained membrane.
The structure is already oriented, so PPM is not needed.

```bash
python -m membraneforger \
    --aa membraneforger/examples/preeq_cg_cellmem/6WHC_MTZP_prot-lig.pdb \
    --cg membraneforger/examples/preeq_cg_cellmem/6WHC_MTZP_cg_cellmem.gro \
    --orientation none \
    --out example_out
```

Backmapping takes most of the run time, about 45 minutes on 8 cores. When the build finishes, the last line of
`example_out/membranebuilder.log` says `PASS` and `example_out/em.gro` exists.

For your own protein in a bundled membrane, name the membrane-spanning chain and the side of its N terminus:

```bash
python -m membraneforger --aa complex.pdb --orient-chain R --nterm-side out --out my_build
```

`--nterm-side` refers to residue 1 as it appears in your file. A receptor that still carries its signal peptide
usually needs `in`. The [tutorial](membraneforger/docs/membraneforger_tutorial.md) walks through a full build.

To build the same protein in five different membranes, add `--runs 5` with `--cg 2`: one system per bundled kappa opioid receptor frame, in
`my_build/KOR1/`, `my_build/KOR2/`, `my_build/KOR3/`, `my_build/KOR5/` and `my_build/KOR6/` (see [Advanced options](membraneforger/docs/advanced.md#several-membranes-for-one-complex)).

## Inputs and outputs

You provide an all-atom structure (`--aa complex.pdb`: every chain you want simulated, each with a unique chain ID;
no water, ions, crystallization additives or insertion codes; ligands need a topology in `membraneforger/forcefield/toppar/`) and,
optionally, a membrane: `--cg 1` (default, a bundled GPR139 membrane), `--cg 2` (kappa opioid receptor), `--cg NAME`
such as `KOR5` (`--list-membranes` lists them), or `--cg FILE` for a Martini 3 frame of your own complex. The
orientation options (`--orient-chain`, `--nterm-side`, `--orientation`) say which chain sits in the membrane. Every
option is in the [command-line reference](membraneforger/docs/cli_reference.md).

The main files in the `--out` directory:

| File | Content |
|---|---|
| `em.gro` | minimized all-atom system, ready for equilibration |
| `topol.top`, `toppar/` | GROMACS topology and force-field files |
| `index_ini.ndx` | index groups `Protein_LIG`, `MEMB` and `SOL_ION` |
| `membranebuilder.log` | build log; the last line says `PASS` or `FAIL` |
| `membrane_validation.md` | membrane quality report, with figures |
| `audit.json`, `run_manifest.json` | audit result, settings, versions and checksums |
| `int/` | every intermediate file, kept for inspection and restarts |

Each stage checks its own result, and a failed check stops the build and names the stage. The checks before and after
backmapping and minimization, and the thresholds, are described in the
[tutorial](membraneforger/docs/membraneforger_tutorial.md#5-validation). After you equilibrate,
`python -m membraneforger.equilibration` grades convergence and membrane properties on the trajectory.

## Documentation

- [Tutorial](membraneforger/docs/membraneforger_tutorial.md): build stages, validation and thresholds, limitations.
- [Command-line reference](membraneforger/docs/cli_reference.md): every option and its default.
- [Advanced options](membraneforger/docs/advanced.md): orientation modes, single-pass proteins, your own membrane, box size,
  lipid composition, restarts.
- [Troubleshooting](membraneforger/docs/troubleshooting.md): common failures and what to check.
- [Running on Slurm](membraneforger/slurm/README.md) and [the bundled membranes](membraneforger/examples/preeq_cg_cellmem/README.md).
- [Module map](membraneforger/docs/modules.md) and [Python example](membraneforger/examples/build_from_python.py), for developers.
- [Changelog](membraneforger/CHANGELOG.md).

## Tests

```bash
python -m unittest discover -s membraneforger/tests -t .
```

The command-line tests need GROMACS on the `PATH`.

## Citation and license

If you use MembraneForger, cite this repository ([CITATION.cff](CITATION.cff)). Also cite the tools and force fields
it relies on: Martini 3, mstool, CHARMM36 and CHARMM-GUI, CGenFF, GROMACS, and PPM 3.0 or OPM when you use them for
orientation.

MembraneForger is released under the MIT license ([LICENSE.txt](LICENSE.txt)). The force-field files in `membraneforger/forcefield/`
and the mapping data in `membraneforger/backmap_data/` come from third-party projects and remain under their own terms.
