# MembraneForger

A tested workflow for converting a Martini 3 coarse-grained protein–membrane simulation into an all-atom
CHARMM36 / GROMACS system.

MembraneForger takes one frame of your coarse-grained simulation and the all-atom structure of the same protein
complex, and returns a solvated, neutralized, energy-minimized all-atom system that has passed an independent audit.

```bash
python -m membraneforger --all-atom prot-lig.pdb --coarse-grain system.gro --orient-chain R --out output_directory
```

`--orient-chain` names the chain that actually spans or associates with the membrane (the receptor, transporter or
channel subunit). MembraneForger determines the membrane orientation of that chain and moves the complete complex,
with its partners, peptides, ligands, ions and modified residues, as one rigid object. With a single protein chain the
option can be left out; with several chains it is required, and the program says so.

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
| Anchor chain | `--orient-chain` | The chain that spans or associates with the membrane. Its membrane orientation positions the whole complex. Optional only when the input has one protein chain. |

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
- For membrane orientation: PPM 3.0, the standalone program behind the OPM/PPM web server. Its Fortran source
  (`ppm3_code/`, with `res.lib`) is distributed by the OPM team; compile it with `make` (needs `gfortran`), which
  produces the executable `immers`. Point MembraneForger at it with `--ppm-exe /path/to/immers` or
  `MEMBRANEFORGER_PPM`, or put it on `PATH`; `res.lib` must sit next to the executable. PPM is not bundled here.
  When the structure is a PDB entry that OPM holds, its orientation is taken from OPM instead (one download, cached)
  and PPM is not needed for that build.
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
    --orient-chain R --nterm-side out --ppm-exe /path/to/immers \
    --out example_out --ntomp 8
```

The receptor is chain R and its N terminus is extracellular, which PPM needs because the file carries no PDB ID. With
`--pdb-id 6WHC` the orientation is taken from the OPM entry instead, and neither PPM nor `--nterm-side` is needed.

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
| `--orient-chain C`, `--orient-chains A,B` | anchor chain(s) that define the membrane orientation |
| `--nterm-side in\|out` | side of the membrane of the first anchor chain's N terminus (needed by PPM) |
| `--orientation auto\|ppm\|opm\|none` | orientation source; `none` uses the coordinates as given |
| `--pdb-id ID` | exact PDB ID, for the OPM reference (default: `HEADER` record, never the file name) |
| `--box auto\|X,Y,Z` | `auto` (default: slice the membrane around the complex), or an opt-in box in nm |
| `--xy-buffer NM` | membrane kept around the complex in x and y with `--box auto` (default 1.0) |

Run `python -m membraneforger --help` for the full list.

## Membrane orientation

```bash
python -m membraneforger ... --orient-chain R                                          # auto
python -m membraneforger ... --orientation ppm --orient-chain R --nterm-side out       # force PPM 3.0
python -m membraneforger ... --orientation opm --pdb-id 7F6G --orient-chain R          # force the OPM entry
python -m membraneforger ... --orientation none                                         # coordinates as given
```

- `auto` (the default) uses the exact OPM/OPRLM entry of the structure when a PDB ID is known (`--pdb-id`, else the
  `HEADER` record of the input) and the anchor chain matches it by sequence and by its membrane-embedded backbone;
  otherwise it runs PPM 3.0 on the anchor chain's own coordinates. Custom models, predicted and mutated structures
  therefore work without a PDB ID.
- `--nterm-side in|out` is the side of the membrane the N terminus of the first anchor chain lies on (`in` =
  cytoplasmic). PPM needs it. It is read from the OPM entry when there is one; otherwise it must be given, and the
  build stops before anything expensive runs and says so. It is never guessed from a protein family or name.
- `--orient-chains A,B` orients from several chains at once. Only the anchor chain(s) are submitted to PPM; every
  other chain, ligand and ion is moved with the derived transform.
- `--orientation none` keeps the previous behaviour for inputs that are already oriented and placed: the all-atom
  complex is fitted onto the coarse-grained protein with a free rigid fit.

Orientation moves the complex into the OPM/PPM frame (normal +z, midplane z = 0, IN negative z, OUT positive z). The
coarse-grained frame is then registered to it laterally: the anchor's membrane-embedded residues set a rotation about
z and an xy translation, and z = 0 is put on the coarse-grained bilayer midplane (halfway between the phosphate
planes). The tilt and depth by which the equilibrated coarse-grained pose differs are measured and reported
(`orientation_report.json`, `run_manifest.json`) and never adopted; a difference above 20 degrees or 0.5 nm is
refused, because the lipid cavity of the coarse-grained membrane would no longer fit the oriented complex.

**Orientation does not make a coarse-grained membrane generic.** PPM makes the all-atom side work for any protein.
The `--coarse-grain` frame, including the pre-equilibrated KOR and GPR139 frames, is an equilibrated system that
contains its own receptor and a lipid cavity shaped around it; it is usable only with an all-atom structure of that
same protein. A different protein needs a coarse-grained system of its own.

Advanced options: `--ppm-membrane CODE` selects a PPM 3.0 membrane model (default: PPM's undefined flat bilayer, which
assumes nothing about the biological membrane); `--ppm-heteroatoms` submits the anchor chains' heteroatoms to PPM;
`--opm-file` uses an already downloaded OPM/OPRLM coordinate file (offline); `--opm-cache` sets the download cache.
PPM is always run in its planar single-membrane mode.

## Box

`BOX = auto` (the default) cuts the coarse-grained membrane before it is backmapped, so a large frame costs no more
than the complex needs. This is the slicing of the earlier workflow (`archive/scripts/aa_stage4.py`, `slice_membrane`,
buffer 1.0 nm), applied to the coarse-grained lipids instead of the finished all-atom system:

- x and y: the window is the extent of the placed complex plus `--xy-buffer` (default 1.0 nm) on each side. An axis
  whose window is as wide as the cell is not cut and keeps its periodicity.
- A lipid is kept only if **all** of its beads lie inside the window, taken in its periodic image nearest the complex;
  glycolipids such as GM3 are kept or dropped as one molecule.
- The window becomes the new periodic cell. Bead pairs of different lipids closer than 0.30 nm that only the new
  periodicity brings together are removed (the lipid with most clashes first) and reported; contacts that were already
  close in the source frame are data, are left alone and are counted separately.
- z is sized around the bilayer midplane (halfway between the phosphate planes) with 1.5 nm of water above and below
  the complex and membrane, never less than 2 nm over the full extent.
- Lipid counts, per-leaflet composition, area kept and the seam report are logged and written to `run_manifest.json`
  (`slice`, `box`). Slicing changes the composition: leaflet ratios are those of the kept lipids.

The 6WHC example keeps 313 of 392 lipids (its cell is cut in x only, because the complex nearly spans y). For the
18 nm GPCR frames the default keeps about 75 to 90 of about 1,300 lipids, which leaves the protein only 2 nm of lipid
from its own periodic image; use `--xy-buffer 1.5` or `2.0` for production systems:

| `--xy-buffer` | KOR run1 cell (nm) | KOR lipids | GPR139 run3 cell (nm) | GPR139 lipids |
|---|---|---|---|---|
| 1.0 | 6.2 x 5.9 | 75 | 6.1 x 7.0 | 88 |
| 1.5 | 7.2 x 6.9 | 120 | 7.1 x 8.0 | 138 |
| 2.0 | 8.2 x 7.9 | 171 | 8.1 x 9.0 | 194 |
| 3.0 | 10.2 x 9.9 | 307 | 10.1 x 11.0 | 339 |

(Measured by cutting the frames with the coarse-grained protein standing in for the complex; no frame has been
backmapped.) A box of your own is opt-in: `--box x,y,z` in nm. x and y may not exceed the coarse-grained cell and must
leave at least 0.5 nm of membrane on each side of the complex (the membrane is cut to that size around the complex;
the full cell width means no cut), and z must be at least the automatic minimum.

## What MembraneForger does

1. **Reads and checks both inputs**, and records their checksums.
2. **Orients the all-atom complex** in the membrane frame from the anchor chain and verifies that every atom received
   the same rigid transform.
3. **Classifies every coarse-grained residue** as protein, lipid species, water or ion. Martini 2 files and unknown
   residues are refused, never silently dropped.
4. **Registers the all-atom complex** to the coarse-grained frame. Each all-atom chain is matched to its
   coarse-grained chain by sequence; the coarse-grained frame supplies only the lateral position (a rotation about the
   membrane normal and a translation), so the orientation is never tilted.
5. **Slices the coarse-grained membrane** to the complex plus a buffer in x and y (`BOX = auto`), keeping whole lipids
   and removing seam overlaps, so that only those lipids are backmapped.
6. **Backmaps the membrane** around the placed complex with mstool.
7. **Reviews the backmapped lipids** for bonds threaded through rings and for wrong stereochemistry. If it finds
   one, it backmaps again with a new seed (up to five attempts).
8. **Builds the CHARMM36 topology**: `pdb2gmx` for the protein, molecule topologies for lipids and ligands.
9. **Sizes the box in z**, adds water and 0.15 M NaCl, and writes the index groups.
10. **Energy-minimizes** with GROMACS.
11. **Audits the result** from the output files alone (`grompp -maxwarn 0`, `gmx check`, energies, geometry, ring
   threading). Only after the audit passes is the final structure published as `em.gro`.

MembraneForger stops at the energy-minimized system. Equilibration and production molecular dynamics are up to you.

## Outputs

| File | Meaning |
|---|---|
| `em.gro` | the validated, energy-minimized system; present only after a passing audit |
| `oriented.pdb`, `orientation_report.json` | the complete complex in the membrane frame, and the orientation record: provider, anchor, N-terminal side, PPM/OPM parameters and hashes, fit RMSD, rotation, translation, validation and CG registration |
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

- 187 unit, negative, invariance, orientation and slicing tests (the command above; three need real PPM, the network or an mstool-free interpreter and are skipped otherwise).
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
- **Orientation.** Needs an OPM entry of the structure or a PPM 3.0 installation; PPM's N-terminal side must be
  supplied when no OPM entry gives it. The coarse-grained frame must contain the same protein complex, because its
  membrane is an equilibrated, protein-specific patch rather than a generic template.
- **Rigid placement.** The all-atom complex is placed as one rigid body. A chain that moved relative to the others
  during the coarse-grained run keeps its all-atom pose, and the build is refused if it no longer fits.
- **Box and membrane.** Orthorhombic cells with one planar bilayer. Triclinic cells and more than one bilayer are
  refused. With `BOX = auto` the membrane is cut around the complex, so the lipid composition of the build is that of
  the kept lipids, and the cut edges are a seam that was never equilibrated; minimization and the audit are the checks
  on it.
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
