# MembraneForger tutorial: building and validating an all-atom membrane system

MembraneForger takes an all-atom protein (or protein-ligand) structure and a pre-equilibrated Martini 3 membrane
and returns a validated, energy-minimized CHARMM36/GROMACS system. This document explains what each stage does,
how the build is validated, and how to read the validation report (`membrane_validation.md`, `.json`, the figures)
that every build writes next to `em.gro`.

```
equilibrated Martini membrane (frame)
    -> orientation of the complex (OPM / PPM) and placement into the frame            (orient, mapping / embed)
    -> embedded membrane: the lipids overlapping the complex are removed              (the REFERENCE for slicing)
    -> PBC-aware, whole-lipid, anchor-based crop to the complex + buffer               (slice)
    -> gate 1: APL, composition, lateral RDF, thickness, protein orientation, integrity (slice_check)
    -> mstool backmapping around the rigid complex                                    (backmap, assemble)
    -> gate 2: the all-atom membrane (membrane.pdb), then the minimized system (em.gro)
    -> topology, box, water, ions, restrained + free minimization, audit               (publish em.gro)
    -> gate 3 and 4 on the equilibration / production trajectory                       (python -m membraneforger.equilibration)
```

Sections 1-4 are the build; section 5 is the validation, with 5.3 (slicing and global APL), 5.5 (species-specific
APL) and 5.12 (lateral RDF) the parts most readers come for; section 6 collects the thresholds; section 7 the
limitations; section 8 the settings.

## 1. Inputs

- `--aa complex.pdb`: the all-atom complex, the only source of protein and ligand chemistry.
- `--cg 1 | 2 | FILE`: a bundled GPR139 or kappa opioid receptor frame (18.3 x 18.3 nm, about 1,300 lipids, ten
  species, 30 us of Martini 3), or your own Martini 3 frame.
- `--orient-chain R` names the chain that spans the membrane; `--nterm-side in|out` is needed by PPM when no OPM
  entry exists. The complete complex moves as one rigid body.

## 2. Placement and embedding

With a bundled membrane (or `--embed`) the oriented complex is placed where the frame's own receptor was, its
bilayer centre on the midplane found between the two PO4 planes, and every lipid with a bead within 0.40 nm of a
heavy atom of the complex is removed (`embed_complex`). With your own frame of the same complex the all-atom
structure is fitted onto the frame's coarse-grained protein instead and nothing is removed.

The membrane that exists after this step is the **embedded membrane**. It is the reference every later stage is
compared with: it is the membrane the slice is actually cut from. The uncut frame is also measured ("CG frame"
row) so that the effect of embedding itself is visible: it raises the leaflet APL by about 1-2 A^2 (the cavity the
removed lipids leave is assigned to their neighbours), which is reported and not graded.

## 3. Slicing (BOX = auto)

The window is the complex's x, y extent plus `--xy-buffer` (1.0 nm) on each side; an axis whose window would reach
the cell is not cut. The cut is made lipid by lipid:

1. every lipid is made whole under the periodic boundaries of the source cell;
2. one **anchor** represents the lipid in the plane (`membraneforger/lipids.py`, `LIPID_SLICE_ANCHORS`):

   | Lipid | Martini 3 anchor | atomistic anchor (CHARMM36) |
   |---|---|---|
   | POPC, DOPC, POPE, DOPE, POPS, DOPS, PSM | `PO4` | `P` |
   | PIP2 (`SAP6` / `SAPI25`) | `PO4` (the diester phosphate, not P4/P5) | `P` |
   | cholesterol (`CHOL` / `CHL1`) | `ROH` | `O3` |
   | GM3 (`DPG3` / `GLPA`) | centroid of the three sugar residues' beads (GLC, GAL, NMC) | `NF` (ceramide amide N) |
   | anything else | centre of geometry of the whole residue | centre of geometry of the heavy atoms |

3. the periodic image of the anchor that falls into the **half-open** window `lower <= a < lower + width` is found.
   The window is never wider than the cell, so exactly one image qualifies: no lipid is counted twice and a lipid on
   a boundary is counted on one side only. This is the decision one gets from tiling the membrane 3 x 3 and cutting
   the window out of the tiling, without building the tiling;
4. if the anchor is in the window the **complete residue** is kept and moved rigidly by the same shift; a tail bead
   that crosses the new cell edge stays where it is (mstool runs with `pbc=True`, GROMACS wraps);
5. the window may be shifted on a cropped axis by up to 0.5 nm (never past the 0.5 nm minimum margin) in steps of
   0.1 nm; the offset whose lipid counts best match the reference leaflet density and composition is used
   (`--no-slice-offset` keeps the centred window). On the bundled frames the gain is small: the within-leaflet
   scatter against the whole cell changed from 5.9 / 4.2 % (lower / upper, standard deviation over 18 frames) to
   5.5 / 3.8 %;
6. the new periodic seam is relaxed in place at the coarse-grained level: a short steepest descent pushes apart
   beads of different lipids that only the cut brought within 0.30 nm of each other, with the complex as a fixed
   repulsive wall, every bead held near its position and every lipid held in shape by intramolecular distance
   restraints. Pairs that were already that close in the source frame are data and are left alone. Only a lipid
   still in a hard-core overlap (< 0.15 nm) afterwards is removed; on the 18 bundled frames none was. The
   relaxation moves about 40 % of the lipids by 0.02 nm on average and 0.3 nm at most.

### Why the previous rule was wrong

The earlier rule kept a lipid only when **every** bead lay inside the window. Lipid tails reach a nanometre
sideways, so a band of lipids along every cut edge lost one tail bead to the outside and was deleted, while the new
cell kept the full window area. The lipid number fell by 15-30 % per unit area and the two leaflets, with different
compositions and tail disorder, lost different amounts. Energy minimization cannot create lipids, and NPT
relaxation cannot repair it either: both leaflets share one x, y cell, so no single lateral contraction can restore
two different deficits. (Earlier report drafts quoted 59 and 70 A^2 per lipid for the two leaflets from an
area-over-count estimate and predicted a 24 % NPT shrink; those numbers were approximate and the interpretation was
wrong. The Voronoi measurement below replaces them.)

## 4. Backmapping and the rest of the build

mstool backmaps the sliced membrane around the complex held as a rigid obstacle, the structure is checked for ring
threading, stereochemistry and clashes (another seed is tried when needed), the topology is built, the box is set in
z, water and 0.15 M NaCl are added, the system is minimized with and then without the lipid dihedral restraints,
and an independent audit must pass before `em.gro` appears. Nothing in this part changed; the measurements of
gate 2 were added around it.

## 5. Validation

### 5.1 What is measured, and where

Every stage is reduced to one record per lipid (anchor position, leaflet, species) and the heavy atoms of the
complex, so the same code measures the coarse-grained and the atomistic stages. The report's stage matrix:

| Validation | CG frame | Embed | Slice | Backmap | Minimized | Equilibrated | Final |
|---|---|---|---|---|---|---|---|
| Ring penetration / clashes (build) | | | ✓ | ✓ | ✓ | | |
| Voronoi APL (global and per species) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Lateral headgroup RDF | ✓ | ✓ | ✓ | | | | |
| Bilayer thickness (global, local map) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Composition / leaflet asymmetry | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Z-density / core hydration | | | | | ✓ | ✓ | ✓ |
| Tail interdigitation | | | | ✓ | ✓ | ✓ | ✓ |
| Protein tilt / insertion depth | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Tail order S_CD | | | | | | ✓ | ✓ |
| Area / volume convergence | | | | | | ✓ | ✓ |
| Lateral diffusion, K_A, leaflet tension | | | | | | | ✓ |

The build runs the first five columns. The last two need a trajectory and are run afterwards with
`python -m membraneforger.equilibration` (section 5.9). The lateral RDF is confined to the coarse-grained stages
because the all-atom anchors (P, O3, NF) are not the Martini beads: a cross-resolution curve is reported in the
JSON but not graded. Everything ends in one record format (`membraneforger/qc.py`): metric, stage, leaflet, value,
units, reference, uncertainty, N, deviation, status (PASS / WARNING / FAIL / INSUFFICIENT SAMPLING / REFERENCE
NEEDED / NOT RUN) and the reason. Poor sampling is never turned into a PASS or a FAIL.

### 5.2 Two different areas per lipid

Keep these apart; the report does.

- **Global leaflet APL**: `APL_global = (A_xy - A_protein) / N_lipids`. It tests overall packing.
- **Local, per-lipid area**: the cell of each lipid in a 2D periodic Voronoi tessellation of its leaflet, generated
  from the headgroup anchors with the heavy atoms of the complex that lie between the midplane and the headgroup
  plane as additional generators whose cells are discarded. The protein footprint is therefore never assigned to a
  lipid, `sum(A_i) = A_xy - A_protein` exactly, and the global APL is the mean of the per-lipid areas.
- **Species APL**: the mean (and median, SD, IQR) of the per-lipid areas of that species, per leaflet.
  "Available area divided by the number of POPC" is not the APL of POPC in a mixed bilayer and is not computed.
- **Head-plane area**: a second tessellation of the phosphate/amide anchors alone, sterols left out because their
  ROH sits about 0.5 nm below the phosphate plane. This is the area per lipid head in the headgroup plane
  (about 65 A^2 in these membranes against 48-52 A^2 for the all-lipid value, the difference being the cholesterol).

### 5.3 Slicing and the global APL

The slice is compared with the embedded membrane in two ways, and only the first is a pass/fail criterion:

- **Construction check**: the leaflet APL of the slice against the Voronoi areas of the **same lipids in the
  uncut embedded membrane**. The cut keeps those lipids with their neighbours; only the new seam can change their
  areas, so this isolates what slicing did. PASS within 3 %, WARNING to 5 %, FAIL beyond (and the build stops
  before backmapping unless `--apl-validate no`; `--apl-tolerance` moves the FAIL ceiling).
- **Representativeness**: the slice against the whole embedded cell, and the lipid count against the count the
  reference APL predicts for the accessible area of the new cell. A finite crop around a protein samples the
  membrane next to that protein, not the cell average, so this is graded against the empirical distribution of
  equal-size windows cut everywhere in the embedded membrane (central 95 % PASS, 95-99 % WARNING) and never stops
  the build.

Measured on the bundled kappa opioid receptor frame `KOR1` with the oriented glucagon receptor complex (6WHC)
embedded, cut to 12.7 x 10.6 nm (the values come from `membrane_validation.md`; the old-rule row is the lipid
count the all-beads-inside rule keeps in the same window, with the APL that count implies):

| Stage / algorithm | Upper lipids | Lower lipids | Upper APL (A^2) | Lower APL (A^2) | dAPL upper vs embed | dAPL lower vs embed | Verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| CG frame (frame's own receptor) | 684 | 618 | 47.0 | 52.0 | | | |
| Embedded CG membrane (reference) | 656 | 560 | 48.0 | 52.6 | reference | reference | |
| Old all-beads-inside slice | 211 | 158 | 54.9 | 60.5 | +14.4 % | +15.0 % | FAIL |
| New anchor slice, centred window | 249 | 199 | 46.0 | 47.7 | -4.2 % | -9.3 % | (whole cell) |
| New anchor slice, offset + relaxed seam | 249 | 193 | 46.5 | 49.6 | -3.0 % | -5.8 % | PASS |
| same lipids in the uncut membrane | 249 | 193 | 47.5 | 50.8 | construction: -2.0 % | -2.4 % | PASS |

The old rule deleted 73 lipids whose anchor was inside the window (35 lower, 38 upper) because one tail bead
crossed the edge. The new slice keeps them; the seam relaxation removed none. The same build on the 11.2 nm `6WHC`
frame (fit, not embedded; one axis cut): old rule +11.5 / +13.6 %, new slice +1.8 / -0.7 % against the whole cell
and -0.1 / +2.5 % against the same lipids uncut.

Over all 18 bundled frames (36 leaflets) the construction deviation of the new slice has a median of 1.6 %, a 95th
percentile of 3.8 % and a largest value of 5.7 %, which is where the 3 / 5 % ceilings come from; the old rule gave
+12 to +32 % on the same frames. The whole-cell deviation of a correct crop scatters more (standard deviation 4-6 %,
lower leaflet systematically 3-4 % denser near the protein than the cell mean), which is why it is reported as
representativeness and not used to stop the build. A 50 -> 56 and 51 -> 59 A^2 change (the GIPR `gipsi_v3` build
that exposed the defect) is +12 and +16 %: FAIL on both leaflets under either comparison.

The repaired slice therefore preserves the embedded membrane's packing to within 2-3 % per leaflet on the
representative cases and within 5 % on every bundled frame tested; it does not preserve it "perfectly", and the
report says by how much.

### 5.4 Composition and leaflet asymmetry

For each leaflet and species the report gives N and the mole fraction in the reference and in the slice, the
difference in mole percent, the sampling standard deviation a random crop of that many lipids would have, and the
composition distance `D = 0.5 * sum |x_slice - x_reference|`. D is graded against the equal-size-window distribution
(and against three times the multinomial expectation when fewer than 20 windows exist). Rare species (PIP2, GM3,
10-20 molecules in a cut) are shown with their N and flagged "(few)" below 10; their mole fractions move by one or
two points between reference and slice and their species APL has a standard deviation of 20-25 A^2 across molecules,
so a difference of a few percent in their mean carries no weight.

### 5.5 Species-specific APL

The species table (and `membrane_validation_species.png`) lists, per stage, leaflet and species: N, mean, median,
SD and head-plane area of the per-lipid Voronoi cells, the mole fraction, and the change against the reference.
In the `6WHC` build the well-sampled species (N >= 30) change by -2.4 to +3.8 % between the embedded membrane and
the slice, the same order as the global check; species with fewer than 10 molecules (DOPC and POPC in the lower
leaflet, DOPE and POPE in the upper leaflet of that frame) move by up to 6.5 % and are marked. Cholesterol's cell
(about 45 A^2 in the upper leaflet) is smaller than a phospholipid's (51-64 A^2) because its anchor shares the plane
with the phosphates it sits between; the head-plane column, which excludes it, is the area a phospholipid head
occupies among heads (66-84 A^2). Per-lipid areas are written to `membrane_validation_lipids.tsv`.

### 5.6 Bilayer thickness and protein orientation

Thickness is the distance between the median anchor z of the two leaflets (non-sterol lipids), plus a 1 nm local
map whose spread is reported. Protein orientation is the tilt of the principal axis of the complex's heavy atoms
within 1.5 nm of the midplane, its insertion depth (mean z of those atoms relative to the midplane) and an inversion
flag. Across slicing (gate 1) thickness may change by 3 / 5 % (PASS / WARNING), tilt by 1 / 3 degrees and depth by
0.5 / 1.5 A; in the examples all three change by 0.0-1.4 % and 0.0 degrees, as expected from a rigid translation.

### 5.7 Integrity

Every kept lipid is complete, no lipid index appears twice, every anchor lies inside the new cell (within the seam
relaxation's displacement), no hard-core seam overlap remains, and no lipid bead ends closer to the complex than the
closest one was before the cut. Ring threading and the all-atom lipid-solute scan belong to the build itself and
run after backmapping and minimization as before.

### 5.8 Gate 2: backmapping integrity

The backmapped membrane (`membrane.pdb`) and the minimized system (`em.gro`) are measured with the atomistic
anchors: leaflet lipid counts and Voronoi APL (reported against the embedded reference, not graded across
resolutions), thickness (PASS within 10 %, WARNING to 15 %), tilt (2 / 5 degrees) and depth (1 / 2.5 A), tail
interdigitation (the overlap of the two leaflets' tail z-density profiles; the minimized value is graded against the
backmapped one, 10 / 20 %), and on `em.gro` the Z-density profiles and the core hydration: water oxygens within
0.8 nm of the midplane that are not within 0.5 nm of a protein heavy atom, as a percentage of the bulk water density
(PASS below 1 %, WARNING to 5 %), with a connected-component test for a water column spanning the core. A single
structure cannot tell a transient from a persistent path; the record says so. The GIPR data that motivated this work
(backmap 57 / 59, minimized 57 / 59 against a slice at 56 / 59 A^2) are consistent with what the report measures:
backmapping and minimization keep the lipid count and the cell, hence the leaflet mean; what they change is local.

### 5.9 Gates 3 and 4: equilibration and equilibrium properties

```
python -m membraneforger.equilibration --out BUILD_DIR --energy ener.xvg --frames frame*.gro \
    [--frame-dt PS] [--production-fraction 0.5] [--temperature 310] [--stress profile.dat] \
    [--reference-ka MN_PER_M] [--reference-d UM2_PER_S]
```

Convergence (gate 3) is judged on the final analysis window of each series: box area and volume from the energy
file; APL per leaflet, thickness, core hydration, protein tilt and depth from the frames. A series passes when its
fitted drift over the window is below 1 % of the mean and the first-half / second-half means differ by less than
2 %; WARNING to 3 / 5 %; FAIL beyond or when the block means are monotonic with more than 1 % drift; block standard
errors are reported. S_CD profiles per species and chain are computed from the last frame and reported as
REFERENCE NEEDED unless a matched profile is supplied (then RMSE 0.03 / 0.05).

Gate 4 uses the production fraction only: lateral diffusion from the drift-corrected x, y MSD of the anchors
(reported only when `MSD ~ t^alpha` with alpha in 0.9-1.1 and R^2 >= 0.98; WARNING for alpha in 0.8-0.9 or 1.1-1.2;
otherwise "not in a diffusive regime"), K_A from the equilibrium area fluctuations with a block-bootstrap
uncertainty (relative uncertainty above 40 % is INSUFFICIENT SAMPLING), and the leaflet tensions from a lateral
pressure profile (`gamma = 0.1 * sum (P_N - P_L) dz` in mN/m per leaflet; the difference passes within 5 mN/m,
warns to 10, fails beyond when the uncertainty supports it; a single leaflet above 10 mN/m is flagged). Without a
matched reference D and K_A are REFERENCE NEEDED: a universal target is not enforced. All these thresholds live in
`membraneforger/trajectory.py` (`TENSION_THRESHOLDS_MN_PER_M`) and `membraneforger/membrane_report.py` so they can
be revised as reference systems accumulate.

### 5.10 Overall status

PASS when every graded record passes; PASS WITH WARNINGS when nothing fails but something warns, is poorly sampled
by design (rare species) or lacks a reference; FAIL when any record fails (a leaflet beyond 5 %, a broken or
duplicated lipid, a water column through the core, a clearly non-converged series, ...); INSUFFICIENT SAMPLING when
nothing is wrong but a required equilibrium quantity could not be classified. A build stops at the slice gate when
its verdict is FAIL and `--apl-validate yes` (the default).

### 5.11 Figures

`membrane_validation_apl.png` (global APL by leaflet and stage, with the old rule's estimate as a dashed line),
`membrane_validation_species.png` (median and interquartile range of the per-lipid areas per species, leaflet and
stage, N on every bar), `membrane_validation_rdf.png` (embedded against sliced headgroup g(r), all-atom curves
dashed), `membrane_validation_composition.png` (counts and mole percent per species and leaflet). Areas are in A^2,
distances in A.

### 5.12 Lateral RDF

The in-plane pair distribution of the headgroup anchors of one leaflet, g(r), is computed with the periodic minimum
image of the stage's own cell in 1 A bins out to half the cell (at most 25 A), normalized by the leaflet's
lipid-accessible area so that g -> 1 whatever fraction of the cell the protein covers; for every species with at
least 20 molecules in the leaflet its same-species curve is added. Each curve is summarized by its first-shell peak
(the maximum of the 5-bin smoothed curve between 3 and 12 A), the following minimum, the first-shell centroid and
the coordination number up to the minimum.

The decisive comparison is **sliced CG against the same lipids in the uncut embedded membrane**: same
representation, same force field, same molecules, so the curves must agree within counting noise. The first-shell
peak may move by 0.5 A (PASS) or 1.0 A (WARNING); the RMS difference of the curves between 3 and 20 A may be at most
three times their combined Poisson noise. The slice's deviation from the whole-cell curve is in addition graded
against the empirical distribution of the same deviation over equal-size windows cut everywhere in the embedded
membrane (central 95 % PASS, 95-99 % WARNING, beyond FAIL): the slice must look like one of the parent's own windows. In the `KOR1` build the peak (7.5 A in both leaflets, the CHOL ROH-PO4
contact dominating the first shell; the POPC-POPC curve peaks near 8-9 A) does not move and the RMS difference is
1.5 noise units in both leaflets; in the `6WHC` build 0.7 / 0.8. The slice therefore keeps the lateral organization
of the membrane it was cut from. The whole-cell curve is reported alongside; it differs more (up to 3 noise units
on `KOR1`) because a cell with a larger protein fraction has a different accessible geometry, not because of the cut.

Across resolutions the anchors change from Martini beads to atoms, and the atomistic force field has its own
excluded-volume distances. Backmapping preserves the lateral topology inherited from the equilibrated Martini
membrane; subsequent atomistic minimization relaxes short-range contacts toward the preferred packing distribution
of the atomistic force field, which can move the first peak (the 8.5 A to about 6.5 A shift seen earlier) without
changing the leaflet mean APL. Such a shift is expected physics and is not graded; the test of slicing is the
CG-to-CG comparison above.

## 6. Thresholds in one place

The twelve validation tests, the stage each is best read at, and the criteria the code applies (`membraneforger/qc.py`
turns every one into a PASS / WARNING / FAIL / INSUFFICIENT SAMPLING / REFERENCE NEEDED record):

| # | Validation test | Best stage | What it validates | Primary output | PASS | WARNING | FAIL / investigate |
|---|---|---|---|---|---|---|---|
| 1 | Lateral headgroup RDF | CG frame -> embed -> slice | lateral lipid organization survives slicing | g(r) per leaflet, first-shell peak and minimum, coordination number, curve deviation | vs the same lipids uncut: peak shift <= 0.5 A and RMS <= 3 x counting noise; vs the whole cell: within the central 95 % of equal-size parent windows | 0.5-1.0 A shift, or 95-99 % of the window distribution | > 1.0 A shift, RMS > 3 x noise, or outside 99 % |
| 2 | Bilayer thickness | all stages | expansion / compression, hydrophobic thickness | D_HH from the anchor planes, 1 nm local map | slice <= 3 % from embed; backmap / minimized <= 10 %; equilibrium <= 3 % from a matched reference | 3-5 % (CG), 10-15 % (AA), 3-5 % (reference) | beyond; REFERENCE NEEDED without a reference |
| 3 | Lipid-tail order parameters | equilibrated AA | acyl-chain ordering, phase | S_CD per species, chain (sn-1, sn-2) and leaflet | profile RMSE <= 0.03 vs a matched profile | 0.03-0.05 | > 0.05; REFERENCE NEEDED without a profile |
| 4 | Z-density profiles | backmap -> minimized -> equilibrated | bilayer architecture along the normal | rho(z) of heads, tails, cholesterol, protein, ions, water | heads outside tails and water outside heads on both sides, cholesterol between; first-to-last change reported | cholesterol not between tails and heads | head density inside the tail region |
| 5 | Composition and leaflet asymmetry | frame -> embed -> slice -> final | the intended composition is kept | N and mole fraction per species and leaflet, distance D | within the central 95 % of equal-size parent windows (or 3 sigma of a random crop with < 20 windows) | 95-99 % | outside 99 % for species with N >= 10; rarer species are flagged "(few)", never failed alone |
| 6 | Leaflet differential tension | final equilibrated AA | mechanical balance of the asymmetric leaflets | gamma_upper, gamma_lower, delta gamma from a lateral pressure profile | abs(delta gamma) <= 5 mN/m with compatible uncertainty | 5-10 mN/m, or one leaflet above 10 mN/m | > 10 mN/m with the uncertainty supporting it; INSUFFICIENT SAMPLING otherwise |
| 7 | Lipid lateral diffusion | final equilibrated AA | fluidity, lipid mobility | drift-corrected x, y MSD, D, exponent alpha | alpha 0.9-1.1 with R^2 >= 0.98; D within a factor 2 of a matched reference | alpha 0.8-0.9 or 1.1-1.2; factor 2-5 | not diffusive (not graded) or > 5-fold |
| 8 | Area compressibility K_A | final equilibrated AA | elasticity | K_A from equilibrium area fluctuations, block-bootstrap uncertainty | <= 20 % from a matched reference, uncertainty <= 20 % | 20-40 % | > 40 %; uncertainty > 40 % = INSUFFICIENT SAMPLING |
| 9 | Hydrophobic-core water | minimized -> equilibrated | pores, cavities, packing defects | core water / bulk water (protein-associated water excluded), spanning-column test, persistence over frames | < 1 %; no column in any frame | 1-5 %; a column in fewer than half the frames (transient) | > 5 %; a column in half the frames or more (persistent) |
| 10 | Tail interdigitation | backmap -> minimized -> equilibrated | overlap of the opposing leaflets | overlap integral of the two tail-density profiles | <= 10 % from the reference (backmapped stage, then a matched reference) | 10-20 % | > 20 % |
| 11 | Protein insertion and orientation | embed -> slice -> backmap -> equilibration | placement and hydrophobic matching | tilt, insertion depth, inversion, embedded fraction | slicing: tilt <= 1 deg, depth <= 0.5 A; AA stages 2 deg / 1 A; equilibrium drift < 3 deg and < 1 A | 1-3 deg / 0.5-1.5 A (slice); 2-5 deg / 1-2.5 A (AA); 3-7 deg / 1-2 A drift | beyond, inversion, or a monotonic drift |
| 12 | Box-area / volume convergence | AA equilibration | stationarity of NPT equilibration | L_x L_y, V, APL and thickness vs time, block SEM | final-window drift < 1 % and half-window means within 2 % | 1-3 % or 2-5 % | > 3 %, > 5 %, or a persistent monotonic trend |

Plus the hard construction checks that are PASS or FAIL only: complete residues, no duplicate lipid, every anchor
inside the cell, no hard-core seam overlap, no new protein-lipid clash (gate 1); ring threading, stereochemistry
and lipid-solute clashes (gate 2, the build's own checks). Where the thresholds live: `Settings` in `config.py`
(APL and RDF ceilings, window samples), `membrane_report.CG_GATE` / `AA_GATE` (thickness, tilt, depth),
`structure_metrics` (core water, interdigitation, layering), `trajectory` (convergence, S_CD, diffusion, K_A,
`TENSION_THRESHOLDS_MN_PER_M`) and `equilibration.ORIENTATION_DRIFT`.

## 7. Limitations and assumptions

- The seam of a cut membrane was never equilibrated: the relaxation removes overlaps, it does not equilibrate the
  interface. Equilibrate the built system before production, as before.
- The construction check isolates slicing; the whole-cell comparison carries the heterogeneity of the parent
  membrane around the protein (the lower leaflet of the bundled frames is 3-4 % denser there than on average).
  The offset search reduces that only slightly within the margin it is allowed.
- Protein generators in the tessellation take area from the lipids next to them by the perpendicular-bisector rule;
  lipids whose heads sit above a protein loop at the headgroup level get small cells. The slab is bounded by the
  anchor plane (plus 0.2 nm) so that atoms in the water beyond the heads do not count.
- GM3's slicing anchor is the sugar centroid (as requested) while its atomistic anchor is the amide nitrogen; the two
  differ by the tilt of the sugars, which affects only the cross-resolution rows.
- Species RDFs and species APLs for 10-20 molecules are reported with N and should be read as such.
- Gates 3 and 4 have been tested on synthetic data (random walks, Gaussian area series, antisymmetric stress
  profiles) and on files, not yet on a production trajectory of a MembraneForger build; `MDAnalysis` is optional
  and only needed to read xtc/trr directly (write frames with `gmx trjconv -sep` otherwise).
- No GIPR build or trajectory is part of this repository; the before/after numbers above are from the bundled
  frames and the 6WHC example, which show the same defect (+12 to +32 % with the old rule) and the same repair.

## 8. Settings

`membraneforger/config.py` carries the defaults, each with the reasoning: `box_xy_buffer_nm` (1.0), the slice
offset search (`slice_optimize_offset`, `slice_offset_search_nm` 0.5, `slice_offset_step_nm` 0.1), the seam
(`seam_min_bead_nm` 0.30, `seam_hard_core_nm` 0.15, `seam_relax_steps` 400), the gate (`APL_VALIDATE`,
`APL_SLICE_WARNING_PERCENT` 3, `APL_SLICE_TOLERANCE_PERCENT` 5, `RDF_VALIDATE`, the RDF ceilings and
`slice_window_samples` 12 per axis). On the command line: `--apl-validate yes|no`, `--apl-tolerance PERCENT`,
`--rdf-validate yes|no`, `--no-slice-offset`, `--xy-buffer NM`.
