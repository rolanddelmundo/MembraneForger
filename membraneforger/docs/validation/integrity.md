# Integrity checks and gates

The construction checks, the backmapping gate, the equilibration gates, the overall status, the figures, the lateral RDF and the definitions of every graded quantity.

## 5.7 Integrity

Every kept lipid is complete, no lipid index appears twice, every anchor lies inside the new cell (within the seam
relaxation's displacement), no hard-core seam overlap remains, and no lipid bead ends closer to the complex than the
closest one was before the cut. Ring threading and the all-atom lipid-solute scan belong to the build itself and
run after backmapping and minimization.

## 5.8 Gate 2: backmapping integrity

The backmapped membrane (`membrane.pdb`) and the minimized system (`em.gro`) are measured with the atomistic
anchors: leaflet lipid counts and Voronoi APL (reported against the embedded reference, not graded across
resolutions), thickness (PASS within 10 %, WARNING to 15 %), tilt (2 / 5 degrees) and depth (1 / 2.5 A), tail
interdigitation (the overlap of the two leaflets' tail z-density profiles; the minimized value is graded against the
backmapped one, 10 / 20 %), and on `em.gro` the Z-density profiles and the core hydration: water oxygens within
0.8 nm of the midplane that are not within 0.5 nm of a protein heavy atom, as a percentage of the bulk water density
(PASS below 1 %, WARNING to 5 %), with a connected-component test for a water column spanning the core. A single
structure cannot tell a transient from a persistent path; the record says so.

The build's own audit of `em.gro` was extended as well. Its dihedral-restraint check reads the DIHRES rows of the
lipid topologies, and for GM3 (GLPA) those 24 rows cover only the ceramide C2/C3 centres, the C4=C5 bond, the ceramide
amide C2S-NF-C1F-C2F (held trans during `emres`; its CHARMM36 torsion has no force constant), the three
ring chairs and sialic acid C7/C8; a sugar carbon with an inverted configuration can sit in a correct chair and pass
them. The audit therefore also tests, geometrically, all 16 sugar stereocentres (glucose C1-C5, galactose C1-C5,
sialic acid C2 and C4-C8) and the two ceramide trans bonds of every GM3 with the chirality definitions of
`membraneforger/backmap_data/map.dat` (`membraneforger/stereo.py`, the same test as `membraneforger/examples/check_gm3_stereo.py`), and fails the
build on any inverted centre or cis bond (a ceramide trans bond under 90 degrees): a configuration never corrects itself in MD.
The report also counts each trans bond as trans (150 degrees or more), twisted (30 to 150) or cis (30 or less). The six inositol stereocentres of every SAPI25 are tested the same way with the SAPI25 chirality definitions of `map.dat`,
because its six chair rows stay satisfied when a hydrogen and a hydroxyl swap places at one carbon. The same test runs on
`membrane.pdb` right after backmapping, where a wrong GM3 counts as a wrong configuration in the backmap verdict and
another seed is tried, so a build only reaches the audit with a GM3 that was wrong on every attempt.

Two lipids threaded through each other, or an atom placed on top of another molecule's atom, are the other defects
minimization cannot repair (a frozen ligand cannot move away either, so an atom overlapping a ligand fails the same way and the finding names the ligand pair): steepest descent ends at a finite Fmax on the atom the bonds hold in place (one build ended
at Fmax 20878 kJ/mol/nm on a DOPC carbon 2.4 A from a POPC carbon, unchanged through 50000 steps). The lipid scan and
the start-contact check only look at lipids against the solute, so the build also tests the structure the way
GROMACS sees it (`membraneforger/contacts.py`): every non-bonded pair (different molecules, or the same molecule more
than `nrexcl` bonds apart) against `contact_rmin_fraction` (0.6) of its own Lennard-Jones Rmin, taken from the NBFIX
table or the combination rule, hydrogens included, under periodic boundaries; and every pair of bonds of two molecules
for a crossing (closest approach below `bond_crossing_nm`, 0.12, inside both bonds). At 0.6 Rmin the repulsion of a
heavy-atom pair is thousands of kJ/mol/nm, while hydrogen bonds and packed tails stay above about 0.7 Rmin. Pairs of
two solvent molecules are not judged: `gmx genion` puts an ion on a water's site, whose neighbours then sit at
water-water distances (0.6-0.7 of the chloride-oxygen Rmin), and a free water or ion relaxes in the first EM steps; a
solvent molecule caught in the protein or a lipid is what the gate is for. The gates: `boxed.gro` right after
backmapping (a membrane finding is logged as a warning and no seed is retried for it: a fit-mode complex such as the quick
start has 20 to 600 membrane pairs below 0.6 Rmin and up to 30 crossings at this stage whatever the seed, and plain EM removes
them; `membrane_overlap_before_em = "stop"` refuses
the seed and finally the build instead), the water and ion placement (a
water is dropped, an ion placement retried, when any of its atoms is inside the limit of a protein, ligand or lipid
atom), `solv_ions.gro` before EM (a solvent finding stops the build; a membrane finding or crossing only when `membrane_overlap_before_em` is `stop`), and `em.unverified.gro`
after EM (nothing may remain). Every gate's findings are written to `contacts.json`,
by stage, and an Fmax failure names the culprit's element, its force-field overlaps before and after EM and the stage a
lipid first clashed in (`em_clash_trace.json`).

Five builds on the code of this release, default settings, show the behaviour: no seed is retried for the membrane pairs, and
every build recovers under plain EM.

| Build | Attempts | Membrane pairs and crossings warned before EM | Fmax after EM (steps) | Audit | Closest membrane pair after EM (Rmin) |
|---|---|---|---|---|---|
| KOR (GCGR complex) | 1 | 3 | 430 (6) | PASS | 0.686 |
| GPR139 (GCGR complex) | 1 | 3 | 319 (9) | PASS | 0.685 |
| prot-cgmemb (the README example, own CG membrane) | 2 (seed 1 had a ring threading) | 144 | 287 (13) | PASS | 0.694 |
| V6G (GLP1R with orforglipron) | 1 | 163 | 424 (6) | PASS | 0.680 |
| TIR (GLP1R with tirzepatide) | 1 | 122 | 417 (7) | PASS | 0.687 |

After EM no build has a contact or crossing of any class. The closest pair of any class is the GDP O3B to Mg2+ coordination
contact at 0.601 to 0.610 of Rmin, just above the 0.6 limit.

Known limitation of the gate: atoms within `nrexcl` bonds of each other are all treated as excluded, but force-field
topologies reintroduce some of them as explicit 1-4 interactions (`[ pairs ]`, `[ pairtypes ]`). A badly overlapping
intramolecular 1-4 pair is therefore not examined, and the gate relies on backmapping and EM for it. No validated build
has shown such a case. A follow-up would audit real 1-4 distances and energies on validated builds, and only if that shows a
gap add a topology-aware check, benchmarked on valid and deliberately distorted lipids rather than a fixed cutoff.

The salt is counted in the water, not in the box. A membrane system is about 30 % protein and lipid by volume, so 0.15 M of
the whole cell would be about 0.21 M in the water. `gmx solvate` fills the box from `spc216.gro` (33.46 waters per nm^3), so
the number of waters that survive solvation gives the aqueous volume, and the target is N = round half up (0.15 M x 0.6022 nm^-3 per M x V_aq), with V_aq = N_water / 33.46 nm^-3
(`Settings.salt_volume` = "water"; "box" restores the cell-volume rule). Neutralization follows
`Settings.counterion_mode` (`--counterion-mode`). "additive" (default) adds N NaCl pairs and then |Q| counterions of the
balancing species on top (Q = -217, N = 176: 393 Na+ and 176 Cl-), so the requested salt is reached at any box height.
"inclusive" (the CHARMM-GUI 2017 convention) makes the counterion species equal N and the other ion N - |Q|, so a Q = -217
system at N = 176 would need a negative Cl-: it is set to 0, only the 217 Na+ are added, and a WARN line states the shortfall,
how much taller the box must be (at fixed xy) and that `--counterion-mode additive` reaches the target at any height. The ion
categories are reported apart, in the `ions:` log line and in `run_manifest.json` (`ions`): the neutralizing ions
(species and |Q|), the NaCl pairs, the Na+ and Cl- totals, and their concentrations in the water, c = n / (N_A V_aq) with
V_aq = N_water / 33.46 nm^-3 and N_A = 0.602214 nm^-3 per M. The ionic strength is I = 1/2 sum c_i z_i^2 = (c_Na + c_Cl) / 2
for the monovalent ions, which is the salt concentration plus half the neutralizing concentration: it is larger than the added
salt whenever the system carries charge (Q = -217: 0.150 M salt, 0.242 M ionic strength). The ions are placed once, by `gmx genion -np
-nn` with these counts (`-conc` would measure the whole cell); they equilibrate with the membrane afterwards, which
pulls counterions to the surface and leaves the bulk below the nominal value.

Backmapping and minimization keep the lipid count and the cell, hence the leaflet mean; what they change is local.

## 5.9 Gates 3 and 4: equilibration and equilibrium properties

```
python -m membraneforger.equilibration --out BUILD_DIR --energy ener.xvg --frames frame*.gro \
    [--frame-dt PS] [--production-fraction 0.5] [--temperature 310] [--stress profile.dat] \
    [--reference-ka MN_PER_M] [--reference-d UM2_PER_S]
```

Convergence (gate 3): box area and volume come from the energy file; APL per leaflet, thickness, core hydration,
protein tilt and depth from the frames. A series passes when its drift over the final window (the last half of the run)
is below 1 % of the window mean and the first-half / second-half means of the whole series differ by less than 2 % of
its mean (definitions in 5.13); WARNING to 3 / 5 %; FAIL beyond or when the block means are monotonic with more than 1 % drift; block standard
errors are reported. Core hydration is the exception: its counts are too small for a drift fit, so it is graded by occupancy (the mean
percent of bulk water with the single-frame limits, 1 / 5 %) and a trend is tested only with at least 20 frames and at least 10 free waters
per frame; with less, a passing level is INSUFFICIENT SAMPLING. S_CD profiles per species and chain are computed from the last frame and reported as
REFERENCE NEEDED unless a matched profile is supplied (then RMSE 0.03 / 0.05).

Gate 4 uses the production fraction only: lateral diffusion from the drift-corrected x, y MSD of the anchors
(reported only when alpha, the log-log slope of MSD against lag, is in 0.9-1.1 and the linear fit of MSD against lag has R^2 >= 0.98, both over lags between 20 % and 50 % of the longest lag; WARNING for alpha in 0.8-0.9 or 1.1-1.2;
otherwise "not in a diffusive regime"), K_A from the equilibrium area fluctuations with a block-bootstrap
uncertainty (relative uncertainty above 40 % is INSUFFICIENT SAMPLING), and the leaflet tensions from a lateral
pressure profile (`gamma = 0.1 * sum (P_N - P_L) dz` in mN/m for each side of the midplane; the difference passes within 5 mN/m,
warns to 10, fails beyond when the uncertainty supports it; a single leaflet above 10 mN/m is flagged). Without a
matched reference D and K_A are REFERENCE NEEDED: a universal target is not enforced. All these thresholds live in
`membraneforger/trajectory.py` (`TENSION_THRESHOLDS_MN_PER_M`) and `membraneforger/membrane_report.py` so they can
be revised as reference systems accumulate.

## 5.10 Overall status

PASS when every graded record passes; PASS WITH WARNINGS when nothing fails but something warns, is poorly sampled
by design (rare species) or lacks a reference; FAIL when any record fails (a leaflet beyond 5 %, a broken or
duplicated lipid, a water column through the core, a clearly non-converged series, ...); INSUFFICIENT SAMPLING when
nothing is wrong but a required equilibrium quantity could not be classified. A build stops at the slice gate when
its verdict is FAIL and `--apl-validate yes` (the default).

## 5.11 Figures

`membrane_validation_apl.png` (global APL by leaflet and stage, with the old rule's estimate as a dashed line),
`membrane_validation_species.png` (median and interquartile range of the per-lipid areas per species, leaflet and
stage, N on every bar), `membrane_validation_rdf.png` (embedded against sliced headgroup g(r), all-atom curves
dashed), `membrane_validation_composition.png` (counts and mole percent per species and leaflet). Areas are in A^2,
distances in A.

## 5.12 Lateral RDF

The in-plane pair distribution of the headgroup anchors of one leaflet, g(r), is computed with the periodic minimum
image of the stage's own cell in 1 A bins out to half the cell (at most 25 A), normalized by the leaflet's
lipid-accessible area so that g -> 1 whatever fraction of the cell the protein covers; for every species with at
least 20 molecules in the leaflet its same-species curve is added. Each curve is summarized by its first-shell peak
(the maximum of the 5-bin smoothed curve between 3 and 12 A), the following minimum, the first-shell centroid and
the coordination number up to the minimum.

The decisive comparison is **sliced CG against the same lipids in the uncut embedded membrane**: same
representation, same force field, same molecules, so the curves must agree within counting noise. The first-shell
peak may move by 0.5 A (PASS) or 1.0 A (WARNING); the RMS difference of the curves between 3 and 20 A may be at most
three times their combined Poisson noise. In these ten-species membranes the all-anchor first shell is weak (its
height is within 5 % of the uniform value, because cholesterol ROH, PO4 and the GM3 sugar centroid sit at different
distances from each other), so its position is not graded and the curve difference decides; the species curves
carry the structure (POPC-POPC first peak 7.5 A, DOPC-DOPC 9.5 A, CHOL-CHOL 11.5 A in the KOR frame). The slice's
deviation from the whole-cell curve is in addition graded
against the empirical distribution of the same deviation over equal-size windows cut everywhere in the embedded
membrane (below its 95th percentile PASS, below the 99th WARNING, above FAIL; a slice that deviates less than the
parent's windows is never graded down): the slice must look like one of the parent's own windows. In the KOR build (50 us) the RMS difference against the same lipids uncut is 1.6
(upper) and 1.7 (lower) noise units, and the POPC-POPC and DOPC-DOPC first peaks of
the slice lay within one 1 A bin of the reference's. The slice therefore keeps the lateral organization of the
membrane it was cut from. Against the whole embedded cell the KOR slice sits at the 0th (upper) and 18th (lower) percentile of
the 144 equal-size windows (PASS); a window around the protein is not a typical window, which is representativeness,
not a construction error.

Across resolutions the anchors change from Martini beads to atoms, and the atomistic force field has its own
excluded-volume distances. Backmapping preserves the lateral topology inherited from the equilibrated Martini
membrane; subsequent atomistic minimization relaxes short-range contacts toward the preferred packing distribution
of the atomistic force field, which can move the first peak without
changing the leaflet mean APL. Such a shift is expected physics and is not graded; the test of slicing is the
CG-to-CG comparison above.

## 5.13 Definitions of the graded quantities

Every formula below is what the code computes (module in brackets). Symbols: A_xy cell area, A_acc lipid-accessible
area (the cell minus the protein footprint), N lipids in the leaflet, x_s mole fraction of species s.

```
Voronoi area [packing]     A_i = area of the periodic Voronoi cell of lipid i (shoelace formula);  sum A_i = A_acc = A_xy - A_protein
APL                        APL = A_acc / N
composition distance       D_comp = 1/2 sum_s |x_s(slice) - x_s(reference)|;  random crop of n lipids:
                           E[D_comp] ~ 1/2 sum_s sqrt(2/pi) sqrt(x_s (1 - x_s) / n)
g(r) [rdf]                 g(r_k) = A_acc n_k / [ N(N-1)/2 * pi (r_{k+1}^2 - r_k^2) ],  n_k = pairs in shell k
coordination number        n = rho * sum_k g(r_k) pi (r_{k+1}^2 - r_k^2) up to the first minimum,  rho = (N-1) / A_acc
curve difference           rms of g_sample - g_reference over 3-20 A, in units of the counting noise
                           sigma^2 = sum over both curves of g_ref / E_k  (Poisson counts, E_k = counts expected at g = 1)
lipid density near the     rho~(R) = [N(<R) / A(<R)] / [N / A_open], R = in-plane distance of a lipid anchor to the nearest
protein [rdf]              protein heavy atom within 0.6 nm of the leaflet's median anchor z (the protein's cross-section
                           at headgroup depth); A(<R) and A_open are areas on a 0.05 nm grid, A_open the part farther
                           than 0.3 nm from that cross-section (not the Voronoi A_acc). First-shell value: R = 1.0 nm
thickness                  D_HH = median z(upper anchors) - median z(lower anchors), non-sterol lipids
protein tilt               theta = arccos |e_z|, e = leading eigenvector of the covariance of the protein atoms within
                           1.5 nm of the midplane
core water                 C = 100 * [N_free / (2 h A_xy)] / rho_bulk,  h = half-height of the core, N_free = core waters
                           not associated with the protein, rho_bulk = bulk water number density
core water [trajectory]    level = mean of C over the frames, graded like one frame (PASS <= 1 %, WARNING <= 5 %, FAIL beyond);
                           trend (>= 20 frames and >= 10 free waters per frame only) = mean N_free of the second half minus
                           the first half, WARNING if above 2 sigma (Poisson error); otherwise INSUFFICIENT SAMPLING
tail interdigitation       Omega = sum_z min(rho_u, rho_l) / sum_z 1/2 (rho_u + rho_l)  (0 = disjoint, 1 = identical tail
                           profiles). Not the local overlap integral of 4 rho_u rho_l / (rho_u + rho_l)^2 used elsewhere
ion count [solvation]      N = round half up (c * 0.6022 nm^-3 per M * V_aq),  V_aq = N_water / 33.46 nm^-3;
                           ionic strength I = 1/2 sum_i c_i z_i^2
drift [trajectory]         delta_drift = 100 * slope * L_w / <s>_w  (linear fit over the final window of length L_w)
half difference            100 * (<s>_second half - <s>_first half) / <s>_whole series
block SEM                  std(block means, ddof = 1) / sqrt(B), B = 5 blocks
S_CD                       S_CD = < (3 cos^2 theta - 1) / 2 >, theta = angle of the C-H vector to z; negative for ordered
                           chains. NMR reports |S_CD| (GROMACS g_order prints -S_CD); compare with -S or |S|
lateral diffusion          MSD(tau) = 4 D tau + c in two dimensions; D from the linear fit, alpha from log MSD vs log tau
K_A                        K_A = k_B T <A> / var(A)  (sample variance, A the box area), uncertainty from a block bootstrap
                           (B = 5 blocks, 200 resamples)
leaflet tension            gamma_side = 0.1 * sum over the side (P_N - P_L) dz  [mN/m, P in bar, z in nm],
                           P_N = P_zz, P_L = (P_xx + P_yy) / 2;  delta = gamma_upper - gamma_lower
```

D from the MSD is the value in the periodic cell: it is not corrected for the finite box, and in a membrane it depends
on the box width logarithmically (Vogele, Koefinger and Hummer, PRL 2018), so compare it with a reference of a similar
size. K_A from the box-area fluctuations is an apparent modulus; it depends on the system size and on undulations, and
it needs a barostat that samples the NPT fluctuations (not Berendsen).
