# MembraneForger

A tested workflow for converting a Martini 3 coarse-grained protein–membrane simulation into an all-atom
CHARMM36 / GROMACS system.

MembraneForger takes one frame of your coarse-grained simulation and the all-atom structure of the same protein
complex, and returns a solvated, neutralized, energy-minimized all-atom system that has passed an independent audit.

```bash
python -m membraneforger --all-atom prot-lig.pdb --coarse-grain system.gro --out output_directory
```

If any check fails, the build stops with an error and no final structure is written.

## Contents

1. [What you need](#what-you-need)
2. [Before MembraneForger: the coarse-grained simulation](#before-membraneforger-the-coarse-grained-simulation)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Quick start](#quick-start)
6. [What MembraneForger does](#what-membraneforger-does)
7. [Outputs](#outputs)
8. [Running on a cluster](#running-on-a-cluster)
9. [Tests](#tests)
10. [Stage renders](#stage-renders)
11. [Supported systems and known limits](#supported-systems-and-known-limits)
12. [Repository layout](#repository-layout)
13. [Version history](#version-history)
14. [Citation and licensing](#citation-and-licensing)

## What you need

| Input | Option | Description |
|---|---|---|
| Coarse-grained frame | `--coarse-grain` | One frame (`.gro` or `.pdb`) of a Martini 3 simulation of the protein complex in its membrane. It supplies the periodic cell, the membrane composition and geometry, and where the protein sits. |
| All-atom complex | `--all-atom` | A PDB file of the same protein complex with its ligands. It is the only source of protein and ligand chemistry and is moved as one rigid body; none of its atoms are rebuilt. |

The two files do not need matching chain IDs or residue numbers: chains are matched by sequence. Neither input file
is ever modified.

## Before MembraneForger: the coarse-grained simulation

MembraneForger starts from a finished coarse-grained simulation. It does not build or run one. The steps below are
how I prepared mine; they are my own usage, not part of MembraneForger, and must be done before running it.

1. **Coarse-grain the protein** with Martinize2 (Vermouth), using the Martini 3 force field.
2. **Build the membrane** around the coarse-grained protein with INSANE, choosing the lipid composition.
3. **Minimize, equilibrate and run** the coarse-grained system in GROMACS.
4. **Take one frame** of the trajectory as a `.gro` file. That file is the `--coarse-grain` input.

The scripts I used for steps 1 to 3 are the first two stages of the earlier MembraneForger release, kept under
`archive/` with their own documentation (`archive/README.md`, `archive/docs/stage1.md`, `archive/docs/stage2.md`).
They are provided as a record and are no longer maintained.

Any Martini 3 simulation prepared another way works too, as long as it uses the lipids listed under
[Supported systems and known limits](#supported-systems-and-known-limits).

## Requirements

- Python 3.10 or newer with numpy, scipy and networkx:

```bash
pip install numpy scipy networkx
```

- GROMACS (tested with 2025.3). Pass the command with `--gmx` if it is not called `gmx`.
- [mstool](https://github.com/ksy141/mstool) 0.3.9 or 0.3.10. Other versions are refused. If mstool lives in a
  different Python environment, pass that interpreter with `--mstool-python`.
- Optional, only for the stage renders: a Python that can import `vmd` (vmd-python), Pillow, and the `tachyon` ray
  tracer.

## Installation

There is no installation step. Clone the repository and run the package from the repository folder:

```bash
git clone https://github.com/rolanddelmundo/MembraneForger.git
cd MembraneForger
python -m membraneforger --help
```

## Quick start

The repository includes one complete input pair, a glucagon receptor–tirzepatide–Gs complex in a ten-lipid membrane:

```bash
python -m membraneforger \
    --all-atom examples/6WHC_MTZP_run1/prot-lig.pdb \
    --coarse-grain examples/6WHC_MTZP_run1/system.gro \
    --out example_out --ntomp 8
```

This takes about an hour on 8 cores; membrane backmapping is most of it. When the build finishes, the last line of
`example_out/membranebuilder.log` starts with `PASS` and `example_out/em.gro` exists.

Common options:

| Option | Meaning |
|---|---|
| `--out DIR` | output directory |
| `--gmx CMD` | GROMACS command |
| `--mstool-python PATH` | Python interpreter that has mstool |
| `--ntomp N` | threads for backmapping and minimization |
| `--toppar DIR` | your own force-field folder (default `forcefield/`; or set `MEMBRANEFORGER_TOPPAR`) |
| `--data DIR` | your own backmapping data folder (default `backmap_data/`; or set `MEMBRANEFORGER_DATA`) |
| `--membrane FILE` | finish a system that is already all-atom instead of giving the two inputs |

Run `python -m membraneforger --help` for the full list.

## What MembraneForger does

1. **Reads and checks both inputs**, and records their checksums.
2. **Classifies every coarse-grained residue** as protein, lipid species, water or ion. Martini 2 files and unknown
   residues are refused, never silently dropped.
3. **Places the all-atom complex.** Each all-atom chain is matched to its coarse-grained chain by sequence, and the
   whole complex is fitted onto the coarse-grained backbone as one rigid body.
4. **Backmaps the membrane** around the placed complex with mstool.
5. **Reviews the backmapped lipids** for bonds threaded through rings and for wrong stereochemistry. If it finds
   one, it backmaps again with a new seed (up to five attempts).
6. **Builds the CHARMM36 topology**: `pdb2gmx` for the protein, molecule topologies for lipids and ligands.
7. **Rebuilds the box, adds water and 0.15 M NaCl**, and writes the index groups.
8. **Energy-minimizes** with GROMACS.
9. **Audits the result** from the output files alone (`grompp -maxwarn 0`, `gmx check`, energies, geometry, ring
   threading). Only after the audit passes is the final structure published as `em.gro`.

MembraneForger stops at the energy-minimized system. Equilibration and production molecular dynamics are up to you.

## Outputs

| File | Meaning |
|---|---|
| `em.gro` | the validated, energy-minimized system; present only after a passing audit |
| `topol.top`, `toppar/` | topology and the molecule files it includes |
| `index_ini.ndx` | index groups `System`, `Protein_LIG`, `MEMB`, `SOL_ION` |
| `membrane.pdb` | the placed complex with the backmapped membrane, before topology building |
| `aa_cg_mapping.tsv` | which all-atom chain was matched to which coarse-grained chain, and how well |
| `run_manifest.json` | inputs and their checksums, settings, per-stage results and timings, output checksums |
| `audit.json` | the independent audit, check by check |
| `membranebuilder.log` | the full transcript; its last line starts with `PASS` or `FAIL` |

## Running on a cluster

`membraneforger/cpu.submit` is a Slurm array template that runs one build per task from a tab-separated case table
(name, all-atom PDB, coarse-grained file):

```bash
sbatch --array=0-7 \
    --export=ALL,REPO=/path/to/MembraneForger,CASES=cases.tsv,OUTROOT=/path/to/outputs \
    membraneforger/cpu.submit
```

Edit the partition, memory and time limit at the top of the file for your cluster. To drive a build from Python with
your own acceptance settings, see `membraneforger/example.py`.

## Tests

```bash
python -m unittest discover -s membraneforger/tests -t .
```

Set `MEMBRANEFORGER_GMX` if GROMACS is not on the path as `gmx`.

What has been tested:

- 111 unit, negative and invariance tests (the command above).
- Complete builds of nine receptor–peptide–Gs systems (glucagon, GLP-1 and GIP receptors), each passing the audit.
- The bundled example, built from a fresh clone of this layout.
- Two deliberate failures, a truncated minimization and a membrane with a threaded ring, both of which are refused.

Other receptors, lipids and force-field residues have not been exercised.

## Stage renders

`membraneforger/visualization.py` renders every stage of a finished build with VMD and Tachyon: fixed camera, ambient
occlusion, and one colour per lipid species and per protein chain.

```bash
python membraneforger/visualization.py example_out --out renders --mode publication \
    --all-atom examples/6WHC_MTZP_run1/prot-lig.pdb \
    --coarse-grain examples/6WHC_MTZP_run1/system.gro
```

## Supported systems and known limits

- **Martini 3 only.** Supported lipids: POPC, DOPC, POPE, DOPE, POPS, DOPS, PSM, cholesterol, PIP2 (`SAP6`) and GM3.
- **Rigid placement.** The all-atom complex is placed as one rigid body. A chain that moved relative to the others
  during the coarse-grained run keeps its all-atom pose, and the build is refused if it no longer fits.
- **Box and membrane.** Orthorhombic cells with one planar bilayer. Triclinic cells and more than one bilayer are
  refused.
- **Protein.** Inter-chain disulfides and residue insertion codes are refused.
- **GM3 stereochemistry.** Six GM3 stereocentres have malformed definitions in the mstool mapping and are not
  controlled. Every run with GM3 reports this as a warning.
- **Cholesterol and PIP2.** mstool's mappings for these lipids do not use the Martini 3 beads `R6` and `C4`. Those
  beads are dropped before backmapping and reported.
- **Reproducibility.** Backmapping is not bit-reproducible: two runs with the same seed give slightly different
  lipid coordinates, and therefore different water and ion counts. Both pass the same checks.

## Repository layout

| Path | Contents |
|---|---|
| `membraneforger/` | the Python package, one module per pipeline stage (listed in `membraneforger/README.md`) |
| `membraneforger/tests/` | unit, negative and invariance tests |
| `membraneforger/example.py` | running a build from Python |
| `membraneforger/cpu.submit` | Slurm array template |
| `forcefield/` | CHARMM36 force field for GROMACS, lipid and ligand topologies, and the parameter packages of the modified residues used in the example |
| `backmap_data/` | mstool mapping additions (`map.dat`) and extra force-field XML files |
| `examples/6WHC_MTZP_run1/` | the example input pair: `prot-lig.pdb` and `system.gro` |
| `examples/preequilibrated_gpcr_cellmembrane/` | eighteen Martini 3 frames (30 µs) of the kappa opioid receptor and GPR139 in a cell-membrane model; coarse-grained inputs only |
| `archive/` | the earlier MembraneForger (v0.2.0) exactly as released, including the coarse-grained setup stages |

## Version history

- **1.0.0** — the two-input pipeline described here.
- **0.2.0, 0.1.0** — a four-stage workflow that started from a PDB alone and also built the coarse-grained system.
  It is kept unchanged under `archive/` and is no longer maintained.

## Citation and licensing

If you use MembraneForger, cite this repository and the tools and force fields it relies on: Martini 3, Martinize2
and Vermouth, mstool, CHARMM36 and CGenFF, and GROMACS.

MembraneForger's own source files are released under the MIT license (`LICENSE.txt`). The force-field files under
`forcefield/` and the mapping data under `backmap_data/` are third-party data and remain under the terms of their
own projects. `archive/` carries its own license files.
