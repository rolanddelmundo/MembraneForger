# Thresholds { #6-thresholds-in-one-place }

Every validation test in one table, with the stage it is best read at and its PASS, WARNING and FAIL limits.


Gate 0 (the frame, below the table) and the twelve validation tests, the stage each is best read at, and the criteria the code applies (`membraneforger/qc.py`
turns every one into a PASS / WARNING / FAIL / INSUFFICIENT SAMPLING / REFERENCE NEEDED record):

| # | Validation test | Best stage | What it validates | Primary output | PASS | WARNING | FAIL / investigate |
|---|---|---|---|---|---|---|---|
| 1 | Lateral headgroup RDF | CG frame -> embed -> slice | lateral lipid organization survives slicing | g(r) per leaflet, first-shell peak and minimum, coordination number, curve deviation | vs the same lipids uncut: peak shift <= 0.5 A and RMS <= 3 x counting noise; vs the whole cell: below the 95th percentile of equal-size parent windows | 0.5-1.0 A shift, or 95th-99th percentile of the window distribution | > 1.0 A shift, RMS > 3 x noise, or above the 99th percentile |
| 2 | Bilayer thickness | all stages | expansion / compression, hydrophobic thickness | D_HH from the anchor planes, 1 nm local map | slice <= 3 % from embed; backmap / minimized <= 10 %; equilibrium <= 3 % from a matched reference | 3-5 % (CG), 10-15 % (AA), 3-5 % (reference) | beyond; REFERENCE NEEDED without a reference |
| 3 | Lipid-tail order parameters | equilibrated AA | acyl-chain ordering, phase | S_CD per species, chain (sn-1, sn-2) and leaflet | profile RMSE <= 0.03 vs a matched profile | 0.03-0.05 | > 0.05; REFERENCE NEEDED without a profile |
| 4 | Z-density profiles | backmap -> minimized -> equilibrated | bilayer architecture along the normal | rho(z) of heads, tails, cholesterol, protein, ions, water | heads outside tails and water outside heads on both sides, cholesterol between; first-to-last change reported | cholesterol not between tails and heads | head density inside the tail region |
| 5 | Composition and leaflet asymmetry | frame -> embed -> slice -> final | the intended composition is kept | N and mole fraction per species and leaflet, distance D_comp | below the 95th percentile of shell-stratified draws around the complex (or 3 sigma of a random crop with < 20 draws) | 95th-99th percentile | above the 99th percentile for species with N >= 10; rarer species are flagged "(few)", never failed alone |
| 6 | Leaflet differential tension | final equilibrated AA | mechanical balance of the asymmetric leaflets | gamma_upper, gamma_lower, delta gamma from a lateral pressure profile | abs(delta gamma) <= 5 mN/m with compatible uncertainty | 5-10 mN/m, or one leaflet above 10 mN/m | > 10 mN/m with the uncertainty supporting it; INSUFFICIENT SAMPLING otherwise |
| 7 | Lipid lateral diffusion | final equilibrated AA | fluidity, lipid mobility | drift-corrected x, y MSD, D (periodic cell), exponent alpha | alpha 0.9-1.1 with R^2 >= 0.98; D within a factor 2 of a matched reference | alpha 0.8-0.9 or 1.1-1.2; factor 2-5 | not diffusive (not graded) or > 5-fold |
| 8 | Area compressibility K_A | final equilibrated AA | elasticity | K_A from equilibrium area fluctuations, block-bootstrap uncertainty | <= 20 % from a matched reference, uncertainty <= 20 % | 20-40 % | > 40 %; uncertainty > 40 % = INSUFFICIENT SAMPLING |
| 9 | Hydrophobic-core water | minimized -> equilibrated | pores, cavities, packing defects | core water / bulk water (protein-associated water excluded), spanning-column test, persistence over frames | < 1 %; no column in any frame | 1-5 %; a column in fewer than half the frames (transient) | > 5 %; a column in half the frames or more (persistent) |
| 10 | Tail interdigitation | backmap -> minimized -> equilibrated | overlap of the opposing leaflets | overlap integral of the two tail-density profiles | <= 10 % from the reference (backmapped stage, then a matched reference) | 10-20 % | > 20 % |
| 11 | Protein insertion and orientation | embed -> slice -> backmap -> equilibration | placement and hydrophobic matching | tilt, insertion depth, inversion, embedded fraction | slicing: not graded (the complex does not move; reported as a diagnostic); AA stages 2 deg / 1 A; equilibrium drift < 3 deg and < 1 A | 2-5 deg / 1-2.5 A (AA); 3-7 deg / 1-2 A drift | beyond, inversion, or a monotonic drift |
| 12 | Box-area / volume convergence | AA equilibration | stationarity of NPT equilibration | L_x L_y, V, APL and thickness vs time, block SEM | final-window drift < 1 % and half-window means within 2 % | 1-3 % or 2-5 % | > 3 %, > 5 %, or a persistent monotonic trend |

The coarse-grained frame the build starts from is graded on its own before anything is embedded (gate 0,
`MembraneValidation.validate_cg_frame`; there is no reference yet, so the limits are absolute and sit in
`membrane_report.CG_FRAME_GATE`, set just outside what the original eight bundled 50 us frames show). PASS /
WARNING / FAIL: both leaflets populated (empty = FAIL); leaflet APL 42-60 A^2 (WARNING 30-100); headgroup-to-headgroup
thickness 38-47 A (WARNING 28-56); PO4-normal tilt from z <= 5 deg (WARNING <= 10); non-sterol lipids whose headgroup and
beads lie on opposite sides of the midplane <= 0.5 %; longest bead-to-bead link inside a lipid <= 0.8 nm (WARNING <= 1.0;
Martini bonds are <= 0.6 nm); lipid-protein bead pairs closer than 0.30 nm at most max(5, beads/2000) (WARNING up to ten
times that). `--cg-validate no` reports a FAIL without stopping the build. A membrane model unlike the bundled ones (the GCGR example
has APL 54-59 A^2 and D_HH 37 A) may sit outside the PASS windows; that is a WARNING, not a defect.

Plus the hard construction checks that are PASS or FAIL only: complete residues, no duplicate lipid, every anchor
inside the cell, no hard-core seam overlap, no new protein-lipid clash (gate 1); ring threading, the CHARMM-GUI
dihedral restraints, the 16 GM3 sugar stereocentres and 2 ceramide trans bonds, and lipid-solute clashes (gate 2,
the build's own checks). Where the thresholds live: `Settings` in `config.py`
(APL and RDF ceilings, window samples), `membrane_report.CG_GATE` / `AA_GATE` (thickness, tilt, depth),
`structure_metrics` (core water, interdigitation, layering), `trajectory` (convergence, S_CD, diffusion, K_A,
`TENSION_THRESHOLDS_MN_PER_M`) and `equilibration.ORIENTATION_DRIFT`.
