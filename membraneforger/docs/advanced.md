# Advanced options

This page collects the details behind the top-level README. For the full option list see the
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
`membraneforger/slurm/setup.sh` compiles it for you ([membraneforger/slurm/README.md](../slurm/README.md#ppm-30)). PPM is not needed with
`--orientation none` or when an exact OPM entry matches (`--pdb-id`). Advanced: `--ppm-membrane CODE`,
`--ppm-heteroatoms`, `--opm-file FILE`, `--opm-cache DIR`.

## Choosing the membrane

| You have | Use |
|---|---|
| Only an all-atom structure | `--cg 1` (GPR139 frame, default) or `--cg 2` (kappa opioid receptor frames). The complex is embedded where the frame's receptor was. |
| A protein much smaller than a GPCR (one transmembrane helix) | add `--embed-site free` and a larger `--xy-buffer`, e.g. 3 |
| A Martini 3 frame of your own complex | `--cg FILE`: the all-atom structure is fitted onto the frame's coarse-grained protein and no lipid is removed |
| A Martini 3 membrane you want to use with a different protein | `--cg FILE --embed` |
| A specific bundled frame, for a repeatable build | `--cg GPR139_Membrane5` (or `KOR_Membrane1`, …; `--list-membranes` lists the names that exist) |
| The same complex in five distinct membranes | `--runs 5` with `--cg 2`: one build per frame, KOR_Membrane1, KOR_Membrane2, KOR_Membrane3, KOR_Membrane5, KOR_Membrane6 (the GPR139 series holds one frame), see [several membranes for one complex](#several-membranes-for-one-complex) |

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

Eight frames (50 µs, Martini 3) of a GPCR in an asymmetric ten-species cell-membrane model, in
`membraneforger/examples/preeq_cg_cellmem/`: `KOR_Membrane1`, `KOR_Membrane2`, `KOR_Membrane3`, `KOR_Membrane5`, `KOR_Membrane6`, `KOR_Membrane7`, `KOR_Membrane9` (kappa opioid receptor) and `GPR139_Membrane5`
(GPR139), about 18.3 × 18.3 × 19.2 nm each (the number is the simulation run the frame comes from). Leaflet
composition in mole percent, averaged over the 8 frames
(`membraneforger/examples/leaflet_composition.py`):

| Leaflet | CHOL | POPC | DOPC | POPE | DOPE | PSM | DPG3 | POPS | DOPS | SAP6 |
|---|---|---|---|---|---|---|---|---|---|---|
| Upper | 28 | 19 | 19 | 5 | 5 | 14 | 10 | 0 | 0 | 0 |
| Lower | 22 | 5 | 5 | 21 | 21 | 0 | 0 | 8 | 7 | 10 |

Martini and CHARMM36 names differ for two species. GM3 is `DPG3` in the coarse-grained frames (CER16, BGLC, BGAL and
ANE5A merged into one molecule) and PIP2 is `SAP6` there and `SAPI25` in the all-atom files.
[membraneforger/examples/preeq_cg_cellmem/README.md](../examples/preeq_cg_cellmem/README.md) lists the frames in detail.

Lipid contacts around the receptor over the 50 µs of the simulations behind `KOR_Membrane5` and `GPR139_Membrane5` (lipids within 5 and
5-10 Å of the receptor, per leaflet): [KOR_Membrane5](../examples/preeq_cg_cellmem/figures/KOR_Membrane5_lipid_contacts_0-50us.svg),
[GPR139_Membrane5](../examples/preeq_cg_cellmem/figures/GPR139_Membrane5_lipid_contacts_0-50us.svg).

### Several membranes for one complex

The bundled frames are independent 50 µs replicates of one cell-membrane model, so building the same complex in
several of them gives replicate systems that differ in their lipid arrangement, not in composition (each frame is
within 1–2 mole percent of the averages above). `--runs N` does that in one command: the complex is built N times,
each time in the next frame of the chosen series (`--cg 2`: KOR_Membrane1, KOR_Membrane2, KOR_Membrane3, KOR_Membrane5, … ; the GPR139 series holds one frame, so `--runs` above 1 needs `--cg 2`), one
complete, validated system per run in `<out>/<frame>/`, and the command's last line counts the runs that passed.

```bash
python -m membraneforger --aa complex.pdb --orient-chain R --nterm-side out --cg 2 --runs 5 --out complex_runs
# complex_runs/KOR_Membrane1/em.gro ... complex_runs/KOR_Membrane6/em.gro (KOR_Membrane1, KOR_Membrane2, KOR_Membrane3, KOR_Membrane5, KOR_Membrane6)
```

The runs are sequential, so five runs take about five times as long as one (`RUNS=5` in
`membraneforger/slurm/run_membraneforger.sbatch`; raise its time limit). `--cg NAME` builds in one particular frame
(`--list-membranes` lists them), which is the way to spread the runs over separate jobs.

MembraneForger checks every frame it reads, rotating a frame whose coordinates were rotated about z relative to its box
back onto its periodic cell and logging the angle
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
preserves the source's lipid density and composition when `--slice-offset-search` is given (the default window is centred on the complex), and the new periodic seam
is relaxed in place at the coarse-grained level instead of deleting lipids.

`--box-nm X Y Z` (nm) sets the box instead. x and y may only be smaller than the membrane patch; z must leave the default
water padding (2.25 nm above the highest and below the lowest atom). A z too short for the complex is refused before backmapping.

The automatic z is the full extent of the complex and membrane plus the padding at both ends, so a receptor that reaches
further on one side of the bilayer does not add water on the other side. The bilayer is then not at the middle of the box;
`run_manifest.json` records the layout, the clearances and `midplane_fraction_z`. The 2.25 nm default is the CHARMM-GUI water
thickness (22.5 A) and can give a larger box than the earlier 1.5 nm: `--z-pad-nm 1.5` keeps the smaller one, and
`--symmetric-box-z --z-pad-nm 1.5` restores the earlier box centred on the midplane.

After slicing, the log and the slice report give the lipid shells around the transmembrane region and the distance of
the complex to its own image in x and y (`xy_check`): `estimated_shell_equivalents_per_side` is a geometric proxy (half the closest approach
of the transmembrane region to its image, in 0.8 nm units), and `voronoi_shells` counts shells by molecular adjacency on
the lipid coordinates of each leaflet (a periodic Voronoi tessellation of the lipid anchors and the receptor atoms near
the leaflet; first-shell lipids share a cell edge with the receptor, each next shell is adjacent to the previous one;
`voronoi_shells_per_side` is half the lipid rows between the receptor and its image, the lower leaflet). A build with fewer
than two shells by either measure (at least 2-3 are recommended, Feng et al. 2023, section 10.1) or a complex closer
than the cutoff to its image gets a WARN; nothing is resized, so widen the membrane with `--xy-buffer` if you see it.

The default 1.0 nm buffer is small for large complexes. If the slice gate stops the build, or the xy check warns, rerun
with a larger `--xy-buffer` (1.6 is two 0.8 nm lipid shells, 2.4 three; 2.0 is common for production systems).

## Outputs in detail

The top level of the output directory holds the final system and the reports:

| File | Content |
|---|---|
| `em.gro` | minimized coordinates; written only when every check passes |
| `topol.top`, `toppar/` | GROMACS topology and force-field includes |
| `index_ini.ndx` | groups `System`, `Protein_LIG`, `MEMB`, `SOL_ION`, and `LIG` when the system has a ligand |
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

Rerunning into an `<out>` that holds an earlier build is refused unless `--overwrite` is given; checkpoints and trajectories
there also need `--overwrite-restart` and are never deleted otherwise. With `--overwrite` the other outputs and intermediates are
replaced and only `int/membrane.pdb` is kept. Give a new `--out` to keep the earlier build as it is.

## Disulfides

Disulfides are taken from the input geometry: two cysteine SG atoms within 3.0 Å, within one chain. `SSBOND` records
are not read. The pairs are listed in `run_manifest.json` (`disulfides`) and re-measured after minimization. A
disulfide between two chains is refused. See [tutorial 4.1](membraneforger_tutorial.md#41-disulfides).

A pair close enough to be bonded that should nevertheless stay reduced — a contact a predicted model produces where
the chemistry does not support a disulfide, such as two cysteines of a cytosolic protein — is named with
`--exclude-disulfide CHAIN:RESID-CHAIN:RESID`, in the input's own numbering, repeatable:

    --exclude-disulfide A:860-A:1252

That pair is then not bonded, not repaired to 2.03 Å and not expected in the topology: `pdb2gmx` is asked before
linking (`-ss`) and told to leave it, so both cysteines are built as reduced CYS carrying their HG1 thiol hydrogen.
Every other pair of the same chain is still bonded.
Coordinates are never changed, so the two sulfurs start where the model put them and the minimization relaxes the
contact. The exclusions appear in the build log and in `run_manifest.json` (`excluded_disulfides`).

## Using MembraneForger from Python

`membraneforger/examples/build_from_python.py` shows how to call the build from Python, and
[modules.md](modules.md) maps each module to its role.

## Ligands during minimization

A ligand here is any non-protein residue of the `--aa` input (GDP, Mg2+, a small molecule), as listed in
`run_manifest.json`. Both minimizations (`emres` and `em`) freeze the ligand heavy atoms in x, y and z
(`freezegrps = LIG`, `freezedim = Y Y Y` in `int/emres.mdp` and `int/em.mdp`; the `LIG` group is the last group of
`index_ini.ndx`, a subset of `Protein_LIG`). A ligand therefore keeps its placed coordinates, and the protein, lipids and
solvent relax around it. Ligand hydrogens are not frozen, so a hydrogen that was missing from the input and was placed by
the build can still be corrected. Peptide residues and lipidated residues are protein, not ligands, and are not frozen.
Systems without a ligand have no `LIG` group and no freeze lines.

Nothing moves a frozen ligand out of a clash. A deep overlap with a ligand shows up as an Fmax or contact failure naming
the ligand pair, and the fix is in the input pose or in the neighbouring residue, not in the minimization. After EM the
build checks that no frozen atom moved (beyond 0.002 nm, periodic image allowed for), and the audit checks that `LIG`
lies inside `Protein_LIG`. The build does not repack or rebuild protein side chains; side chains are only minimized.

## Lipidated residues

The force field carries three lipidated protein residues that the input may use by name, with every heavy atom of the
lipid present in the input (the build never places a lipid chain itself): `CYSG` (S-geranylgeranyl cysteine, the Ggamma
C terminus), `CYSP` (S-palmitoyl cysteine) and `GLYP` (N-palmitoyl glycine, a Galpha N terminus after Met1 removal; its
backbone N carries the chain, so a chain that starts with `GLYP` gets the `None` N terminus). Hydrogens are rebuilt by
pdb2gmx from `merged.hdb`. `GLYP` is CHARMM's `GLYM` (N-myristoyl glycine, `toppar_all36_lipid_prot.str`) with the chain
extended from C14 to C16 using the palmitate atom types and charges of `CYSP`; no new parameters are involved. The
lipidated lysines of the peptide agonists (`KTZ`, `KRT`, `KSM`) are different: their CGenFF side chains are taken from the
ligand packages in `forcefield/toppar/` and the input must carry all their hydrogens.
