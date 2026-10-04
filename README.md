## MembraneForger: all-atom membrane systems from Martini 3 coarse-grained frames

MembraneForger takes two files, an all-atom protein/ligand complex and one frame of a Martini 3 coarse-grained
simulation of the same complex in a membrane, and builds a validated, energy-minimized all-atom CHARMM36 / GROMACS
system from them.

```bash
python -m membraneforger --all-atom prot-lig.pdb --coarse-grain system.gro --orient-chain R --out output_directory
```

`--orient-chain` names the chain that actually spans or associates with the membrane (the receptor, the
transporter, the channel subunit). MembraneForger determines the membrane orientation of that chain and moves the
complete complex, with its partners, peptides, ligands, ions and modified residues, as one rigid object. If the input
has a single protein chain the option can be left out; with several chains it is required, and the program says so.

What it does, in order:

1. checks both inputs and records their checksums (the inputs are never modified);
2. orients the all-atom complex in the membrane frame from the anchor chain (membrane normal along z, bilayer
   midplane at z = 0, extracellular side at positive z) and proves that every atom received the same rigid transform;
3. classifies every coarse-grained residue (protein, lipid species, water, ions) and refuses Martini 2 files;
4. matches each all-atom chain to its coarse-grained chain by sequence; the coarse-grained frame then supplies only
   the lateral position (a rotation about the membrane normal and a translation): the orientation is never tilted;
5. backmaps the membrane around the placed complex with [mstool](https://github.com/ksy141/mstool);
6. checks the backmapped lipids for bonds threaded through rings and for wrong stereochemistry, and backmaps again
   with a new seed when it finds one;
7. builds the CHARMM36 topology (`pdb2gmx` for the protein, molecule topologies for lipids and ligands);
8. sizes the box (`BOX = auto`: the membrane cell in x and y, z centred on the bilayer midplane with 1.5 nm of water
   beyond the complex on both sides), adds water and 0.15 M NaCl, and writes the index groups;
9. energy-minimizes with GROMACS;
10. audits the result from the output files alone (`grompp -maxwarn 0`, `gmx check`, energies, geometry) and only
    then publishes `em.gro`.

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
- For membrane orientation: PPM 3.0, the standalone program behind the OPM/PPM web server. Its Fortran source
  (`ppm3_code/`, with `res.lib`) is distributed by the OPM team; compile it with `make` (needs `gfortran`), which
  produces the executable `immers`. Point MembraneForger at it with `--ppm-exe /path/to/immers` or
  `MEMBRANEFORGER_PPM`, or put it on `PATH`; `res.lib` must sit next to the executable. PPM is not bundled here.
  When the structure is a PDB entry that OPM holds, its orientation can be taken from OPM instead (one download,
  cached), and PPM is not needed for that build.
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
    --coarse-grain examples/6WHC_MTZP_run1/system.gro --orient-chain R --nterm-side out \
    --ppm-exe /path/to/immers --out example_out --ntomp 8
```

The receptor is chain R; its N terminus is extracellular (`--nterm-side out`), which PPM needs because the file
carries no PDB ID. With `--pdb-id 6WHC` the orientation is taken from the OPM entry instead and neither PPM nor
`--nterm-side` is needed. The backmapping step dominates the run time (about 45 minutes on 8 cores for this
example). The last line of `example_out/membranebuilder.log` starts with `PASS` or `FAIL`.

**Membrane orientation**

```bash
python -m membraneforger --all-atom complex.pdb --coarse-grain system.gro --orient-chain R      # auto
python -m membraneforger ... --orientation ppm --orient-chain R --nterm-side out                # force PPM 3.0
python -m membraneforger ... --orientation opm --pdb-id 7F6G --orient-chain R                   # force the OPM entry
python -m membraneforger ... --orientation none                                                  # coordinates as given
```

- `auto` (the default) uses the exact OPM/OPRLM entry of the structure when a PDB ID is known (`--pdb-id`, else the
  `HEADER` record of the input; never the file name) and the anchor chain matches it by sequence and by its
  membrane-embedded backbone; otherwise it runs PPM 3.0 on the anchor chain's own coordinates. Custom models,
  predicted and mutated structures therefore work without a PDB ID.
- `--nterm-side in|out` is the side of the membrane the N terminus of the first anchor chain lies on (`in` =
  cytoplasmic). PPM needs it. It is read from the OPM entry when there is one; otherwise it must be given, and the
  build stops before anything expensive runs and says so. It is never guessed from a protein family or name.
- `--orient-chains A,B` orients from several chains at once (a dimer, a multi-subunit channel). Only the anchor
  chain(s) are submitted to PPM; every other chain, ligand and ion is moved with the derived transform.
- `--orientation none` keeps the previous behaviour for inputs that are already oriented and placed: the all-atom
  complex is then fitted onto the coarse-grained protein with a free rigid fit.

Orientation moves the complex into the OPM/PPM frame (normal +z, midplane z = 0, IN negative z, OUT positive z).
The coarse-grained frame is then registered to it laterally: the anchor's membrane-embedded residues set a rotation
about z and an xy translation, and z = 0 is put on the coarse-grained bilayer midplane (halfway between the two
phosphate planes). The tilt and depth by which the equilibrated coarse-grained pose differs are measured and
reported (`orientation_report.json`, `run_manifest.json`) and never adopted; a difference above 20 degrees or 0.5
nm is refused, because the lipid cavity of the coarse-grained membrane would no longer fit the oriented complex.
PPM makes the orientation step general; whether a given coarse-grained membrane fits the oriented complex is a
separate question that the registration check answers.

Advanced options: `--ppm-membrane CODE` selects a PPM 3.0 membrane model (default: PPM's undefined flat bilayer,
which assumes nothing about the biological membrane); `--ppm-heteroatoms` submits the anchor chains' heteroatoms to
PPM; `--opm-file` uses an already downloaded OPM/OPRLM coordinate file (offline), `--opm-cache` sets the download
cache. PPM is always run in its planar single-membrane mode; curved membranes are never requested.

**Box**

`BOX = auto` (the default) keeps the membrane's x and y cell, as it must for a periodic membrane, and sizes z around
the bilayer midplane so that the complex and the membrane are covered by 1.5 nm of water above and below (never less
than 2 nm over the full extent). A box of your own is opt-in: `--box x,y,z` in nm, where x and y must equal the
membrane cell and z must be at least the automatic minimum; anything else is refused with the value that would work.

**Outputs**

| File | Meaning |
|---|---|
| `em.gro` | the validated, energy-minimized system (present only after a passing audit) |
| `oriented.pdb`, `orientation_report.json` | the complete complex in the membrane frame, and the orientation record (provider, anchor, N-terminal side, PPM/OPM parameters and hashes, fit RMSD, rotation, translation, validation, CG registration) |
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
- The coarse-grained input must contain the same complex (every coarse-grained protein segment needs an all-atom
  counterpart). Its membrane is an equilibrated, protein-specific patch, not a generic template: orientation makes
  the all-atom side general, but a different protein needs a coarse-grained system of its own.
- Orientation needs either an OPM entry of the structure or a PPM 3.0 installation; PPM's N-terminal side must be
  supplied when no OPM entry provides it.
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
