# Advanced options

This page collects the details that used to live in the top-level README. For the full option list see the
[command-line reference](cli_reference.md). For the science behind each stage and every threshold see the
[tutorial](membraneforger_tutorial.md).

## Membrane orientation

Name the chain that spans or associates with the membrane. MembraneForger determines its membrane orientation and
moves the complete complex, with its partners, peptides, ligands, ions and modified residues, as one rigid body. After
orientation the membrane normal is +z, the bilayer midplane is at z = 0 and the cytoplasmic side is at negative z.

- **`auto` (default).** When a PDB ID is known and the anchor chain matches that exact OPM/OPRLM entry by sequence and
  by its membrane-embedded backbone, the entry's orientation is used (one download, cached; nothing of the entry is
  copied into your structure). Otherwise PPM 3.0 is run on the anchor chain's own coordinates, so models, predicted
  and mutated structures work without a PDB ID.
- **N-terminus side.** PPM needs the side of the N terminus of the first anchor chain. It is read from the OPM entry
  when there is one; otherwise give `--nterm-side in` or `out`. It is never guessed from a protein family or name,
  and the build stops before anything expensive runs if it is missing. The answer refers to residue 1 *as it appears
  in your file*: a receptor that still carries its signal peptide usually needs `in`, not the textbook `out`.
- **Several chains.** With several protein chains `--orient-chain` is required; the error lists the chains.
- **Orienting on one segment.** `--orient-residues 343-363` submits only those residues of the anchor chain(s) to PPM
  (numbers as in the input), so a single-pass protein is oriented on its transmembrane helix rather than on its
  ectodomain. The rest of the complex follows the helix rigidly. `--nterm-side` is then required and gives the side of
  the first selected residue (`out` for a type I protein such as RAGE); no OPM entry is used. Protein residues outside
  the selection whose CA lies in the hydrophobic slab are listed as a `WARN`
  (`validation.frame.outside_segment_residues_inside_slab`):

  ```bash
  python -m membraneforger --aa rage.pdb --orientation ppm --orient-residues 343-363 --nterm-side out --out rage_out
  ```

- **`--orientation none`.** The input is used as given and must already have its membrane normal along z (as from
  OPM or CHARMM-GUI). The bilayer centre is found by a hydrophobic-belt search over the whole protein, which a domain
  with a hydrophobic surface can pull off the helix; `--bilayer-z` sets it explicitly. `--orient-residues` then names
  the embedded segment instead (no PPM, no `--nterm-side`): the bilayer centre is the midpoint of that segment's CA z
  range, and protein CA outside it within 15 Å of that centre are listed (`WARN`). `--nterm-side`, when given, is
  checked on the anchor chain's own N terminus relative to that centre; the wrong side stops the build.
- **The orientation is never undone.** Placement into the membrane is a rotation about z and a translation only. The
  orientation, placement and box are recorded in `orientation_report.json` and `run_manifest.json`.

### PPM 3.0

PPM 3.0 is not bundled. Compile its Fortran source (distributed by the OPM team, `ppm3_code/`) with `make` and pass the
`immers` executable with `--ppm-exe` or `MEMBRANEFORGER_PPM`; `res.lib` must sit next to it. On a cluster,
`slurm/setup.sh` compiles it for you ([slurm/README.md](../slurm/README.md#ppm-30)). PPM is not needed with
`--orientation none` or when an exact OPM entry matches (`--pdb-id`). Advanced: `--ppm-membrane CODE`,
`--ppm-heteroatoms`, `--opm-file FILE`, `--opm-cache DIR`.

## Choosing the membrane

| You have | Use |
|---|---|
| Only an all-atom structure | `--cg 1` (GPR139 frames, default) or `--cg 2` (kappa opioid receptor frames). The complex is embedded where the frame's receptor was. |
| A protein much smaller than a GPCR (one transmembrane helix) | add `--embed-site free` and a larger `--xy-buffer`, e.g. 3 |
| A Martini 3 frame of your own complex | `--cg FILE`: the all-atom structure is fitted onto the frame's coarse-grained protein and no lipid is removed |
| A Martini 3 membrane you want to use with a different protein | `--cg FILE --embed` |
| A specific bundled frame, for a repeatable build | `--cg GPR1` (any of `GPR1` … `GPR9`, `KOR1` … `KOR9`; `--list-membranes` lists them) |
| The same complex in five distinct membranes | `--runs 5`: one build per frame, GPR1 to GPR5 (or KOR1 to KOR5 with `--cg 2`), see [several membranes for one complex](#several-membranes-for-one-complex) |

A bundled membrane is a patch equilibrated around a GPCR. Orientation makes the all-atom side general, but a membrane
frame of your own complex is a patch around that protein and fits only that protein.

For a single-pass protein in a bundled frame:

```bash
python -m membraneforger --aa sp.pdb --orientation none --orient-residues 152-177 --nterm-side out \
    --embed-site free --xy-buffer 3 --out sp_out
```

When the complex is embedded, the lipids it touches are pushed aside at the coarse-grained level. Only those still
overlapping it afterwards are removed, and an empty pocket in a leaflet larger than 1.0 nm³ stops the build. The
tutorial ([section 2](membraneforger_tutorial.md#2-placement-and-embedding)) explains the push, the pocket gate and
how `--embed-site free` chooses its spot.

### Bundled membranes

Eighteen frames (30 µs, Martini 3) of a GPCR in an asymmetric ten-species cell-membrane model, in
`examples/preeq_cg_cellmem/`: `KOR1`–`KOR9` (kappa opioid receptor) and `GPR1`–`GPR9` (GPR139), about
18.3 × 18.3 × 19.2 nm each. Leaflet composition in mole percent, averaged over the 18 frames
(`examples/leaflet_composition.py`):

| Leaflet | CHOL | POPC | DOPC | POPE | DOPE | PSM | DPG3 | POPS | DOPS | SAP6 |
|---|---|---|---|---|---|---|---|---|---|---|
| Outer | 28 | 19 | 19 | 5 | 5 | 14 | 10 | 0 | 0 | 0 |
| Inner | 22 | 5 | 5 | 21 | 21 | 0 | 0 | 8 | 7 | 10 |

Martini and CHARMM36 names differ for two species. GM3 is `DPG3` in the coarse-grained frames (CER16, BGLC, BGAL and
ANE5A merged into one molecule) and PIP2 is `SAP6` there and `SAPI25` in the all-atom files.
[examples/preeq_cg_cellmem/README.md](../examples/preeq_cg_cellmem/README.md) lists the frames in detail.

### Several membranes for one complex

The bundled frames are independent 30 µs replicates of one cell-membrane model, so building the same complex in
several of them gives replicate systems that differ in their lipid arrangement, not in composition (each frame is
within 1–2 mole percent of the averages above). `--runs N` does that in one command: the complex is built N times,
each time in the next frame of the chosen series (GPR1, GPR2, … for `--cg 1`; KOR1, KOR2, … for `--cg 2`), one
complete, validated system per run in `<out>/<frame>/`, and the command's last line counts the runs that passed.

```bash
python -m membraneforger --aa complex.pdb --orient-chain R --nterm-side out --runs 5 --out complex_runs
# complex_runs/GPR1/em.gro ... complex_runs/GPR5/em.gro
```

The runs are sequential, so five runs take about five times as long as one (`RUNS=5` in
`slurm/run_membraneforger.sbatch`; raise its time limit). `--cg NAME` builds in one particular frame
(`--list-membranes` lists them), which is the way to spread the runs over separate jobs.

The bundled frames had been written with their coordinates rotated about z relative to the box, each by its own
angle. They have been rotated back in place (`examples/align_frames.py`), and MembraneForger checks every frame it
reads, rotating a frame with the same artefact back onto its periodic cell and logging the angle
([tutorial 1.1](membraneforger_tutorial.md#11-a-frame-must-be-periodic-in-its-box)).

### Changing the lipid composition

`--dellipid LIPID` removes every molecule of a species (repeatable). `--addlipid LIPID` turns the removed molecules
into another species with the same twelve-bead Martini 3 layout (`POPC`, `DOPC`, `POPE`, `DOPE`, `POPS`, `DOPS`,
`PSM`) instead of deleting them.

## Box size

By default (`BOX = auto`) the membrane is cut before backmapping, so a large frame costs no more than the complex
needs. The window is the complex's extent plus `--xy-buffer` (1.0 nm) in x and y; an axis as wide as the cell is not
cut. z is sized from the protein height around the bilayer midplane.

The cut is made lipid by lipid under periodic boundaries. Each lipid is made whole, one headgroup anchor represents it
(PO4, cholesterol ROH, the GM3 sugar centroid; `membraneforger/lipids.py`), and the complete lipid is kept when the
periodic image of its anchor falls in the window. The window may shift by up to 0.5 nm to the position that best
preserves the source's lipid density and composition (`--no-slice-offset` turns this off), and the new periodic seam
is relaxed in place at the coarse-grained level instead of deleting lipids.

`--box X Y Z` (Å) sets the box instead. x and y may only be smaller than the membrane patch; z must leave the default
water padding. A z too short for the complex is refused before backmapping.

The default 1.0 nm buffer is small for large complexes. If the slice gate stops the build, rerun with a larger
`--xy-buffer` (2.0 is a common choice for production systems).

## Outputs in detail

The top level of the output directory holds the final system and the reports:

| File | Content |
|---|---|
| `em.gro` | minimized coordinates; written only when every check passes |
| `topol.top`, `toppar/` | GROMACS topology and force-field includes |
| `index_ini.ndx` | groups `System`, `Protein_LIG`, `MEMB`, `SOL_ION` |
| `membranebuilder.log` | build log; its last line says `PASS` or `FAIL` |
| `audit.json` | independent audit of the minimized system |
| `membrane_validation.md`, `.json` | stage-by-stage membrane validation |
| `membrane_validation_lipids.tsv` | one row per lipid with its Voronoi area |
| `membrane_validation_*.png` | validation figures (needs matplotlib) |
| `orientation_report.json` | orientation source, transform and checks |
| `run_manifest.json` | inputs, settings, software versions, checksums, stage times, disulfides |

Nothing the build makes on the way is deleted. At the end of every run, passed or failed, the intermediates
(`oriented.pdb`, `membrane.pdb`, `prot-memb.pdb`, `boxed.gro`, `solv.gro`, `solv_ions.gro`, the `.mdp`/`.tpr`/`.log`/
`.edr`/`.trr` files of both minimizations, and `em.unverified.gro` when the build fails before publishing) are moved
into `<out>/int/`. Next to them are the working files (`int/work`: coarse-grained cut, mstool steps, earlier
backmapping attempts) and the audit's scratch runs (`int/audit`).

### Restarting after backmapping

Backmapping dominates the run time. To redo the steps after it (topology, water, ions, minimization, audit), start
from the kept membrane:

```bash
python -m membraneforger --membrane <out>/int/membrane.pdb --out <out>
```

Rerunning into the same `<out>` replaces that build's other outputs and intermediates and keeps only
`int/membrane.pdb`. Give a new `--out` to keep the earlier build as it is.

## Disulfides

Disulfides are taken from the input geometry: two cysteine SG atoms within 3.0 Å, within one chain. `SSBOND` records
are not read. The pairs are listed in `run_manifest.json` (`disulfides`) and re-measured after minimization. A
disulfide between two chains is refused. See [tutorial 4.1](membraneforger_tutorial.md#41-disulfides).

## Using MembraneForger from Python

`membraneforger/example.py` shows how to call the build from Python, and
[membraneforger/README.md](../membraneforger/README.md) maps each module to its role.
