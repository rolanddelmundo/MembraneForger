## MembraneForger: all-atom membrane systems from a single PDB

MembraneForger takes an all-atom protein (or protein/ligand) PDB, embeds it in a pre-equilibrated 40 × 40 nm
Martini 3 coarse-grained membrane, backmaps the membrane around it, and returns a validated, energy-minimized
all-atom CHARMM36 / GROMACS system (`em.gro`).

**You only need one input file: the all-atom structure.** The membrane ships with the package.

```bash
python -m membraneforger --all-atom protein.pdb --out output_directory
```

The box is sized automatically from the protein and membrane. To set it yourself, give the three edge lengths
in Å:

```bash
python -m membraneforger --all-atom protein.pdb --box 80 80 100 --out output_directory
```

If you already have a Martini 3 frame of your complex in its own membrane, pass it with `--coarse-grain` and that
membrane is used instead of the bundled one. In both cases the pipeline:

1. checks the inputs and records their checksums (inputs are never modified);
2. places the all-atom complex in the coarse-grained membrane (rigid fit to the matching coarse-grained chain
   when a frame is given, centred in the bundled membrane otherwise);
3. backmaps the membrane with [mstool](https://github.com/ksy141/mstool), re-running with a new seed when it
   finds a bond threaded through a ring or a wrong stereocentre;
4. builds the CHARMM36 topology (`pdb2gmx` for the protein, molecule topologies for lipids and ligands);
5. sets the box, adds water and 0.15 M NaCl, writes the index groups, and energy-minimizes with GROMACS;
6. audits the result from the output files alone (`grompp -maxwarn 0`, `gmx check`, energies, geometry) and
   only then publishes `em.gro`.

**Bundled membrane**

The pre-equilibrated membrane is a plasma-membrane mimic built with INSANE and run with Martini 3. Final leaflet
compositions, in mole percent of each leaflet (`examples/leaflet_composition.py`, average over the equilibrated
frames, rounded to whole numbers; cholesterol flip-flops, so its two numbers differ from the 25:25 it was built with):

| Lipid | Outer leaflet | Inner leaflet |
|---|---|---|
| Cholesterol (CHOL) | 27 | 23 |
| POPC | 19 | 5 |
| DOPC | 19 | 5 |
| POPE | 5 | 21 |
| DOPE | 5 | 21 |
| Sphingomyelin (PSM) | 15 | 0 |
| GM3 | 10 | 0 |
| POPS | 0 | 8 |
| DOPS | 0 | 7 |
| PIP2 (SAP6) | 0 | 10 |

A build that fails any check stops with an error and does not write `em.gro`. The last line of
`membranebuilder.log` starts with `PASS` or `FAIL`. Backmapping dominates the run time (about 45 minutes on
8 cores for the bundled example).

**Requirements**

- Python 3.10 or newer with numpy, scipy and networkx (`pip install numpy scipy networkx`).
- GROMACS (tested with 2025.3); pass the command with `--gmx` if it is not `gmx`.
- mstool 0.3.9 or 0.3.10; if it lives in another Python environment, pass that interpreter with `--mstool-python`.

No installation step is needed: clone the repository and run `python -m membraneforger --help` from its folder.

**Outputs**

| File | Meaning |
|---|---|
| `em.gro` | the validated, energy-minimized system (present only after a passing audit) |
| `topol.top`, `toppar/`, `index_ini.ndx` | topology, included molecule files, index groups (`System`, `Protein_LIG`, `MEMB`, `SOL_ION`) |
| `membrane.pdb` | the placed complex with the backmapped membrane, before topology building |
| `run_manifest.json`, `audit.json` | inputs, settings, per-stage results and checksums; the independent audit, check by check |
| `membranebuilder.log` | the full transcript |

**Repository layout**

`membraneforger/` holds the package (one module per stage, see `membraneforger/README.md`), its tests, a Python
usage example and a Slurm array template. `forcefield/` is the CHARMM36 force field with lipid and ligand
topologies, `backmap_data/` the mstool mapping additions, and `examples/preeq_cg_cellmem/` the pre-equilibrated
coarse-grained membranes (`<system>_cg_cellmem.gro`) with a matching all-atom example (`6WHC_MTZP_prot-lig.pdb`). Use `--toppar` and `--data` (or `MEMBRANEFORGER_TOPPAR` and `MEMBRANEFORGER_DATA`) to point at your own
folders. The earlier workflow (v0.1.0 and v0.2.0) is kept unchanged under `archive/` and is no longer maintained.

Tests: `python -m unittest discover -s membraneforger/tests -t .` (set `MEMBRANEFORGER_GMX` if GROMACS is not
on the path as `gmx`). Stage renders with VMD and Tachyon: `python membraneforger/visualization.py --help`.

**Supported systems and known limits**

- Martini 3 only. Lipids: POPC, DOPC, POPE, DOPE, POPS, DOPS, PSM, cholesterol, PIP2 (`SAP6`) and GM3.
- The all-atom complex is placed as one rigid body; orthorhombic boxes with one bilayer only. Triclinic cells,
  inter-chain disulfides and residue insertion codes are refused.
- Six GM3 stereocentres are malformed in the mstool mapping and are not controlled; every run reports them.
  The cholesterol and PIP2 beads `R6` and `C4` are not used by mstool and are dropped before backmapping.

**Third-party components**

The force-field files under `forcefield/` and the mapping data under `backmap_data/` come from CHARMM36 /
CHARMM-GUI, CGenFF and mstool. Cite those projects, Martini 3 and GROMACS when you publish results obtained with
this package.
