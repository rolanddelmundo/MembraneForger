## MembraneForger: all-atom membrane systems from a single all-atom PDB

MembraneForger: software package that embeds an all-atom protein (or protein–ligand) structure into a
pre-equilibrated 40 × 40 nm Martini 3 plasma-membrane mimic, backmaps the membrane with
[mstool](https://github.com/ksy141/mstool), and returns a validated, energy-minimized all-atom
CHARMM36 / GROMACS system (`em.gro`). You only need one input file: your all-atom structure.

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
python -m membraneforger --all-atom protein.pdb --out output_directory
```

The box is sized automatically. To set it yourself, give the three edge lengths in Å: `--box 80 80 100`.
To use your own Martini 3 membrane instead of the bundled one, add `--coarse-grain your_frame.gro`.

**Bundled membrane**

Leaflet composition in mole percent (`examples/leaflet_composition.py`, averaged over the equilibrated frames in
`examples/preeq_cg_cellmem/`, rounded to whole numbers):

| Lipid | CHOL | POPC | DOPC | POPE | DOPE | PSM | GM3 | POPS | DOPS | PIP2 (SAP6) |
|---|---|---|---|---|---|---|---|---|---|---|
| Outer leaflet | 27 | 19 | 19 | 5 | 5 | 15 | 10 | 0 | 0 | 0 |
| Inner leaflet | 23 | 5 | 5 | 21 | 21 | 0 | 0 | 8 | 7 | 10 |

**Output**

`em.gro` is written only after every check passes (`grompp -maxwarn 0`, `gmx check`, energies, geometry);
`topol.top`, `toppar/` and `index_ini.ndx` go with it, and `membranebuilder.log` ends with `PASS` or `FAIL`.
Backmapping dominates the run time (about 45 minutes on 8 cores).

Tests: `python -m unittest discover -s membraneforger/tests -t .`

The force field under `forcefield/` and the mapping data under `backmap_data/` come from CHARMM36 / CHARMM-GUI,
CGenFF and mstool; cite those projects, Martini 3 and GROMACS when you publish results obtained with this package.
