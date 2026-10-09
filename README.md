## MembraneForger: all-atom preequilibrated membrane

MembraneForger: software package that embeds an all-atom protein (or protein–ligand) structure into a
pre-equilibrated Martini 3 coarse-grained cell-membrane, backmaps the membrane with
[mstool](https://github.com/ksy141/mstool), and returns a validated, energy-minimized all-atom
CHARMM36 / GROMACS system (`em.gro`). You only need one input file: your all-atom structure. Name the chain that
spans or associates with the membrane (`--orient-chain R`) and MembraneForger orients the complete complex in the
membrane as one rigid object (see Membrane orientation); with `--orientation none` your structure is used as given and
must already have its membrane normal along z (as from OPM or CHARMM-GUI).

**Installing Requirements**

Following Python packages are required: numpy, scipy, networkx, pandas (below 3, for mstool), openmm, mstool (0.3.9
or 0.3.10); matplotlib is optional (it draws the validation figures).
We recommend using pip to install them on your local machine:

```
pip install numpy
pip install scipy
pip install networkx
pip install "pandas<3"
pip install openmm
pip install mstool==0.3.9
pip install matplotlib
```

GROMACS is also required (tested with 2023.3 and 2025.3). If it is not on your path as `gmx`, pass it with `--gmx`.

**Installation**

No particular installation procedure is necessary. Clone the repository and run the package from its folder:

```bash
python -m membraneforger --aa complex.pdb --orient-chain R --out output_directory
```

| Option | Meaning |
|---|---|
| `--aa PDB` | your all-atom structure (required) |
| `--orient-chain C` | the chain that spans or associates with the membrane (optional when the input has one protein chain); `--orient-chains A,B` for several |
| `--orient-residues FIRST-LAST` | orient on these residues of the anchor chain only (e.g. a transmembrane helix, `343-363`); PPM sees them alone and the whole complex follows. With `--orientation none` (input already along z) it names the membrane-embedded segment instead: nothing is rotated, and the bilayer centre goes at the midpoint of that segment's CA z range |
| `--nterm-side in\|out` | side of the membrane the N terminus of that chain lies on; required by PPM unless an OPM entry gives it |
| `--orientation auto\|ppm\|opm\|none` | orientation source (default `auto`); `none` uses your coordinates as given |
| `--pdb-id ID` | exact PDB ID, to use its OPM orientation (default: the `HEADER` record of the input, never the file name) |
| `--ppm-exe PATH` | the compiled PPM 3.0 program (`immers`), or set `MEMBRANEFORGER_PPM` |
| `--cg 1` / `--cg 2` / `--cg FILE` | the membrane: `1` = one of the bundled GPR139 frames at random (default), `2` = one of the bundled kappa opioid receptor frames at random, or `custom=FILE`, a Martini 3 frame of your own complex (add `--embed` to use only its membrane) |
| `--embed-site hole\|free` | where the protein goes in a bundled (or `--embed`) frame: `hole` (default) = where the frame's receptor was; `free` = unbroken membrane, placed so that the membrane cut around the protein stays at least 1 nm from that receptor, for a protein much smaller than it, such as a single transmembrane helix (a lone helix wants a larger `--xy-buffer`, e.g. 3) |
| `--box X Y Z` | opt-in box edges in Å. Default (`BOX = auto`): the membrane is cut to the complex plus `--xy-buffer` (1.0 nm) in x and y, z from the protein height. A smaller x and y cuts the membrane around the complex; x and y may not exceed the membrane patch |
| `--xy-buffer NM` | membrane kept around the complex on each side in x and y when the box is auto (default 1.0) |
| `--apl-validate yes\|no`, `--apl-tolerance PERCENT` | the slice gate: stop before backmapping when a leaflet's area per lipid changes by more than the tolerance (default 5 %, warning above 3 %) against the membrane it was cut from |
| `--rdf-validate yes\|no`, `--no-slice-offset` | compare the lateral headgroup RDF of the slice with its source (default yes); keep the slice window centred instead of choosing the offset that best preserves lipid density |
| `--dellipid LIPID` | remove a lipid species (repeatable): CHOL, POPC, DOPC, POPE, DOPE, POPS, DOPS, PSM, DPG3, SAP6 |
| `--addlipid LIPID` | turn the removed lipids into this one instead (POPC, DOPC, POPE, DOPE, POPS, DOPS, PSM) |
| `--bilayer-z Z` | z of the bilayer centre in your structure, in Å, with `--orientation none` (default: found from the hydrophobic belt; an oriented complex has it at 0) |
| `--ntomp N`, `--gmx CMD`, `--mstool-python PATH` | threads, GROMACS command, Python with mstool |

Run `python -m membraneforger --help` for the rest. The protein is oriented, moved as one rigid body onto the spot
the frame's own receptor occupied, the lipids it touches are pushed aside at the coarse-grained level (only those still
overlapping it afterwards are removed, and an empty pocket in a leaflet stops the build), the membrane is cut around it and backmapped;
the build then adds water and 0.15 M NaCl, minimizes, and audits the result (`grompp -maxwarn 0`, `gmx check`,
energies, geometry, lipid stereochemistry including every GM3 sugar centre). `em.gro` is written only when every check passes, with `topol.top`, `toppar/`,
`index_ini.ndx` and `run_manifest.json`; the last line of `membranebuilder.log` says `PASS` or `FAIL`. Disulfides are
taken from the input geometry (two cysteine SG atoms within 3.0 A, within one chain; SSBOND records are not read),
listed in the manifest and re-measured after minimization (`docs/membraneforger_tutorial.md` 4.1).
Nothing the build makes on the way is deleted: at the end of every run, passed or failed, the intermediates
(`oriented.pdb`, `membrane.pdb`, `prot-memb.pdb`, `boxed.gro`, `solv.gro`, `solv_ions.gro`, the `.mdp`/`.tpr`/`.log`/`.edr`/`.trr`
files of both minimizations, and `em.unverified.gro` when the build fails before publishing) are moved into `<out>/int/`,
next to the working files (`int/work`: CG cut, mstool steps, earlier backmapping attempts) and the audit's scratch
runs (`int/audit`). To redo the steps after backmapping, start from the kept membrane:
`python -m membraneforger --membrane <out>/int/membrane.pdb --out <out>`.
Backmapping dominates the run time (about 45 minutes on 8 cores). Equilibrate the system before production.

**Membrane orientation**

Name the chain that spans or associates with the membrane. MembraneForger determines its membrane orientation and moves
the complete complex, with its partners, peptides, ligands, ions and modified residues, as ONE rigid body (membrane
normal along +z, bilayer midplane at z = 0, cytoplasmic side negative z):

- `auto` (default): when a PDB ID is known and the anchor chain matches that exact OPM/OPRLM entry by sequence and by
  its membrane-embedded backbone, the entry's orientation is used (one download, cached; nothing of the entry is
  copied into your structure). Otherwise PPM 3.0 is run on the anchor chain's own coordinates, so models, predicted and
  mutated structures work without a PDB ID.
- PPM needs the side of the N terminus of the first anchor chain. It is read from the OPM entry when there is one;
  otherwise give `--nterm-side in` or `out`. It is never guessed from a protein family or name, and the build stops
  before anything expensive runs if it is missing.
- With several protein chains `--orient-chain` is required; the error lists the chains.
- `--orient-residues 343-363` submits only those residues of the anchor chain(s) to PPM (numbers as in the input),
  so a single-pass protein is oriented on its transmembrane helix rather than on its ectodomain. The rest of the
  complex follows the helix rigidly. `--nterm-side` is then required and gives the side of the first selected residue
  (`out` for a type I protein such as RAGE); no OPM entry is used. The report lists every protein residue outside the
  selection whose CA lies in the hydrophobic slab (`validation.frame.outside_segment_residues_inside_slab`, a `WARN`
  in the log), which for a single-pass protein should be at most the residues flanking the helix:
  `python -m membraneforger --aa rage.pdb --orientation ppm --orient-residues 343-363 --nterm-side out --out rage_out`.
- With `--orientation none` the input is used as given, and the bilayer centre is found by a hydrophobic-belt search
  over the whole protein, which a domain with a hydrophobic surface can pull off the helix. `--orient-residues` then
  names the embedded segment instead (no PPM, no `--nterm-side`): the bilayer centre is the midpoint of that
  segment's CA z range, and protein CA outside it within 15 A of that centre are listed (`WARN`). `--nterm-side`,
  when given, is checked on the anchor chain's own N terminus relative to that centre; the wrong side stops the build. For a single-pass
  protein in a bundled frame, add `--embed-site free`:
  `python -m membraneforger --aa sp.pdb --orientation none --orient-residues 152-177 --nterm-side out --embed-site free --out sp_out`.
- PPM 3.0 is not bundled. Compile its Fortran source (distributed by the OPM team, `ppm3_code/`) with `make` and pass the
  `immers` executable; `res.lib` must sit next to it. Advanced: `--ppm-membrane CODE`, `--ppm-heteroatoms`,
  `--opm-file FILE`, `--opm-cache DIR`.
- The orientation is never undone afterwards: placement into the membrane is a rotation about z and a translation only.
  The orientation, placement and box are recorded in `orientation_report.json` and `run_manifest.json`.

A bundled membrane is a patch equilibrated around a GPCR. Orientation makes the all-atom side general, but a membrane
frame of your own complex (`--cg FILE`) is a patch around that protein: it fits that protein only.

**Box**

`BOX = auto` cuts the membrane before backmapping, so a large frame costs no more than the complex needs: the window
is the complex's extent plus 1.0 nm in x and y; an axis as wide as the cell is not cut. The cut is made lipid by
lipid under periodic boundaries: each lipid is made whole, one headgroup anchor represents it (PO4, cholesterol ROH,
the GM3 sugar centroid; `membraneforger/lipids.py`), and the complete lipid is kept when the periodic image of its
anchor falls in the half-open window, whether or not a tail bead crosses the edge. The window may shift by up to
0.5 nm to the position that best preserves the source's lipid density and composition, and the new periodic seam is
relaxed in place at the coarse-grained level instead of deleting lipids. (The earlier rule kept a lipid only when all
its beads were inside, which deleted a band of edge lipids and left the cut membrane 12 to 30 % short of lipids per
area; see `docs/membraneforger_tutorial.md` 5.3.)

Before anything is backmapped the slice is validated against the membrane it was cut from (the embedded membrane):
per leaflet, the area per lipid from a periodic Voronoi tessellation of the headgroup anchors with the protein
footprint discarded, the composition, the lateral headgroup RDF, the bilayer thickness and the protein orientation,
plus molecular integrity. A leaflet whose area per lipid changes by more than 5 % stops the build (`--apl-validate`,
`--apl-tolerance`). The backmapped membrane and the minimized system are measured the same way. Everything is
written to `membrane_validation.md` / `.json`, `membrane_validation_lipids.tsv` (one row per lipid with its Voronoi
area) and four figures; `python -m membraneforger.equilibration` adds the convergence and equilibrium-property gates
once a trajectory exists. The tutorial describes every measurement and threshold.

**Bundled membranes**

Eighteen frames (30 µs, Martini 3) of a GPCR in an asymmetric ten-species cell-membrane model, in
`examples/preeq_cg_cellmem/`: `KOR1`–`KOR9` (kappa opioid receptor) and `GPR1`–`GPR9` (GPR139), about
18.3 × 18.3 × 19.2 nm each. Leaflet composition in mole percent, averaged over the 18 frames
(`examples/leaflet_composition.py`):

| Leaflet | CHOL | POPC | DOPC | POPE | DOPE | PSM | DPG3 | POPS | DOPS | SAP6 |
|---|---|---|---|---|---|---|---|---|---|---|
| Outer | 28 | 19 | 19 | 5 | 5 | 14 | 10 | 0 | 0 | 0 |
| Inner | 22 | 5 | 5 | 21 | 21 | 0 | 0 | 8 | 7 | 10 |

The bundled frames had been written with their coordinates rotated about z relative to the box (a rotational fit
saved without its box), each by its own angle; they have been rotated back in place (`examples/align_frames.py`), and
MembraneForger checks every frame it reads, rotating a frame with the same artefact back onto its periodic cell and
logging the angle (`docs/membraneforger_tutorial.md` 1.1).

**On a Slurm cluster**: `bash slurm/setup.sh` installs everything above (no Anaconda or modules needed) and
`slurm/run_membraneforger.sbatch` builds one structure per array task; see `slurm/README.md`, which also covers PPM 3.0.

Tests: `python -m unittest discover -s membraneforger/tests -t .` (the command-line tests need GROMACS on the path).

Documentation: `docs/membraneforger_tutorial.md` (build stages, validation, thresholds, limitations).

The force field under `forcefield/` and the mapping data under `backmap_data/` come from CHARMM36 / CHARMM-GUI,
CGenFF and mstool; cite those projects, Martini 3 and GROMACS when you publish results obtained with this package.
