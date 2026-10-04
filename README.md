## MembraneForger: all-atom membrane systems from a single all-atom PDB

MembraneForger: software package that embeds an all-atom protein (or protein–ligand) structure into a
pre-equilibrated Martini 3 cell-membrane model, backmaps the membrane with
[mstool](https://github.com/ksy141/mstool), and returns a validated, energy-minimized all-atom
CHARMM36 / GROMACS system (`em.gro`). You only need one input file: your all-atom structure, oriented with the
membrane normal along z (as from OPM or CHARMM-GUI).

**Installing Requirements**

Following Python packages are required: numpy, scipy, networkx, pandas, openmm, mstool (0.3.9 or 0.3.10).
We recommend using pip to install them on your local machine:

```
pip install numpy
pip install scipy
pip install networkx
pip install pandas
pip install openmm
pip install mstool==0.3.9
```

GROMACS is also required (tested with 2023.3 and 2025.3). If it is not on your path as `gmx`, pass it with `--gmx`.

**Installation**

No particular installation procedure is necessary. Clone the repository and run the package from its folder:

```bash
python -m membraneforger --aa protein.pdb --out output_directory
```

| Option | Meaning |
|---|---|
| `--aa PDB` | your all-atom structure (required) |
| `--cg 1` / `--cg 2` / `--cg FILE` | the membrane: `1` = bundled kappa opioid receptor frame (default), `2` = bundled GPR139 frame, or a Martini 3 frame of your own complex (add `--embed` to use only its membrane) |
| `--box X Y Z` | box edges in Å; default: the membrane's x and y, z from the protein height. Smaller x and y trim the membrane around the protein |
| `--dellipid LIPID` | remove a lipid species (repeatable): CHOL, POPC, DOPC, POPE, DOPE, POPS, DOPS, PSM, DPG3, SAP6 |
| `--addlipid LIPID` | turn the removed lipids into this one instead (POPC, DOPC, POPE, DOPE, POPS, DOPS, PSM) |
| `--bilayer-z Z` | z of the bilayer centre in your structure, in Å (default: found from the hydrophobic belt) |
| `--ntomp N`, `--gmx CMD`, `--mstool-python PATH` | threads, GROMACS command, Python with mstool |

Run `python -m membraneforger --help` for the rest. The protein is moved as one rigid body onto the spot the
frame's own receptor occupied, lipids overlapping it are removed, and the membrane is backmapped around it; the
build then adds water and 0.15 M NaCl, minimizes, and audits the result (`grompp -maxwarn 0`, `gmx check`,
energies, geometry). `em.gro` is written only when every check passes, with `topol.top`, `toppar/`,
`index_ini.ndx` and `run_manifest.json`; the last line of `membranebuilder.log` says `PASS` or `FAIL`.
Backmapping dominates the run time (about 45 minutes on 8 cores). Equilibrate the system before production.

**Bundled membranes**

Eighteen frames (30 µs, Martini 3) of a GPCR in an asymmetric ten-species cell-membrane model, in
`examples/preeq_cg_cellmem/`: `KOR1`–`KOR9` (kappa opioid receptor) and `GPR1`–`GPR9` (GPR139), about
18.3 × 18.3 × 19.2 nm each. Leaflet composition in mole percent, averaged over the 18 frames
(`examples/leaflet_composition.py`):

| Leaflet | CHOL | POPC | DOPC | POPE | DOPE | PSM | DPG3 | POPS | DOPS | SAP6 |
|---|---|---|---|---|---|---|---|---|---|---|
| Outer | 28 | 19 | 19 | 5 | 5 | 14 | 10 | 0 | 0 | 0 |
| Inner | 22 | 5 | 5 | 21 | 21 | 0 | 0 | 8 | 7 | 10 |

Tests: `python -m unittest discover -s membraneforger/tests -t .`

The force field under `forcefield/` and the mapping data under `backmap_data/` come from CHARMM36 / CHARMM-GUI,
CGenFF and mstool; cite those projects, Martini 3 and GROMACS when you publish results obtained with this package.
