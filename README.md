## MembraneForger: all-atom membrane systems from Martini 3 coarse-grained frames

MembraneForger takes two files, an all-atom protein/ligand complex and one frame of a Martini 3 coarse-grained
simulation of the same complex in a membrane, and builds a validated, energy-minimized all-atom CHARMM36 / GROMACS
system from them.

```bash
python -m membraneforger --all-atom prot-lig.pdb --coarse-grain system.gro --out output_directory
```

What it does, in order:

1. checks both inputs and records their checksums (the inputs are never modified);
2. classifies every coarse-grained residue (protein, lipid species, water, ions) and refuses Martini 2 files;
3. matches each all-atom chain to its coarse-grained chain by sequence and places the all-atom complex on the
   coarse-grained one with a rigid fit;
4. backmaps the membrane around the placed complex with [mstool](https://github.com/ksy141/mstool);
5. checks the backmapped lipids for bonds threaded through rings and for wrong stereochemistry, and backmaps again
   with a new seed when it finds one;
6. builds the CHARMM36 topology (`pdb2gmx` for the protein, molecule topologies for lipids and ligands);
7. rebuilds the box, adds water and 0.15 M NaCl, and writes the index groups;
8. energy-minimizes with GROMACS;
9. audits the result from the output files alone (`grompp -maxwarn 0`, `gmx check`, energies, geometry) and only then
   publishes `em.gro`.

A build that fails any check stops with an error and does not write `em.gro`.

The coarse-grained input frames come from Martini 3 molecular dynamics; the coarse-grained systems for those
simulations were prepared with Martinize2 and Vermouth. MembraneForger starts from a frame of such a simulation and
does not run the coarse-grained simulation itself.

This is version 1.0.0. The earlier four-stage workflow (versions 0.1.0 and 0.2.0) is kept unchanged under `archive/`
and is no longer maintained.

**Requirements**

- Python 3.10 or newer with numpy, scipy and networkx:

```bash
pip install numpy scipy networkx
```

- GROMACS (tested with 2025.3). Pass the command with `--gmx` if it is not `gmx`.
- mstool 0.3.9 or 0.3.10 (other versions are refused). If mstool lives in a different Python environment, pass that
  interpreter with `--mstool-python`.
- Optional, for the stage renders only: a Python that can import `vmd` (vmd-python), Pillow, and the `tachyon` ray
  tracer.

**Installation**

No installation step is needed. Clone the repository and run the package from the repository folder:

```bash
git clone <repository URL>
cd membraneforger
python -m membraneforger --help
```

**Repository layout**

| Path | Contents |
|---|---|
| `membraneforger/` | the Python package (one module per pipeline stage, listed in `membraneforger/README.md`) |
| `membraneforger/tests/` | unit, negative and invariance tests |
| `membraneforger/example.py` | how to run a build from Python with your own acceptance settings |
| `membraneforger/cpu.submit` | Slurm array template, one build per task |
| `forcefield/` | CHARMM36 force field for GROMACS, lipid and ligand topologies, and the parameter packages of the modified residues used in the example |
| `backmap_data/` | mstool mapping additions (`map.dat`) and extra force-field XML files |
| `examples/6WHC_MTZP_run1/` | one complete input pair: `prot-lig.pdb` (all-atom complex) and `system.gro` (Martini 3 frame) |
| `archive/` | the previous MembraneForger (v0.2.0) exactly as released, with its own licence files; not used by v1.0.0 |

Use `--toppar` and `--data` (or `MEMBRANEFORGER_TOPPAR` and `MEMBRANEFORGER_DATA`) to point at your own force-field
and backmapping data folders.

**Example**

```bash
python -m membraneforger --all-atom examples/6WHC_MTZP_run1/prot-lig.pdb \
    --coarse-grain examples/6WHC_MTZP_run1/system.gro --out example_out --ntomp 8
```

The backmapping step dominates the run time (about 45 minutes on 8 cores for this example). The last line of
`example_out/membranebuilder.log` starts with `PASS` or `FAIL`.

**Outputs**

| File | Meaning |
|---|---|
| `em.gro` | the validated, energy-minimized system (present only after a passing audit) |
| `topol.top`, `toppar/`, `index_ini.ndx` | topology, included molecule files, index groups (`System`, `Protein_LIG`, `MEMB`, `SOL_ION`) |
| `membrane.pdb` | the placed complex with the backmapped membrane, before topology building |
| `run_manifest.json` | inputs and their checksums, settings, per-stage results and timings, output checksums |
| `audit.json` | the independent audit, check by check |
| `aa_cg_mapping.tsv` | which all-atom chain was matched to which coarse-grained chain, and how well |
| `membranebuilder.log` | the full transcript |

**Tests**

```bash
python -m unittest discover -s membraneforger/tests -t .
```

Set `MEMBRANEFORGER_GMX` if GROMACS is not on the path as `gmx`.

**Stage renders**

`membraneforger/visualization.py` renders every stage of a finished build with VMD and Tachyon (fixed camera,
ambient occlusion, one colour per lipid species and per protein chain):

```bash
python membraneforger/visualization.py example_out --out renders --mode publication \
    --all-atom examples/6WHC_MTZP_run1/prot-lig.pdb --coarse-grain examples/6WHC_MTZP_run1/system.gro
```

**Supported systems and known limits**

- Martini 3 only. Lipids: POPC, DOPC, POPE, DOPE, POPS, DOPS, PSM, cholesterol, PIP2 (`SAP6`) and GM3.
- The all-atom complex is placed as one rigid body.
- Orthorhombic boxes with one bilayer. Triclinic cells, more than one bilayer, inter-chain disulfides and residue
  insertion codes are refused.
- Six GM3 stereocentres have malformed definitions in the mstool mapping and are not controlled; every run reports
  them as a warning.
- mstool's cholesterol and PIP2 mappings do not use the Martini 3 beads `R6` and `C4`; these beads are dropped
  before backmapping and reported.
- Only the systems listed in the tutorial were exercised.

**Third-party components**

The CHARMM36 force-field files under `forcefield/` and the mapping data under `backmap_data/` come from their
respective projects (CHARMM36 / CHARMM-GUI, CGenFF, mstool). Cite those projects, Martini 3 and GROMACS when you
publish results obtained with this package.
