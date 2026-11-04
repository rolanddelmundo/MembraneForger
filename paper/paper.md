---
title: 'MembraneForger: validated all-atom membrane–protein systems from equilibrated Martini 3 membranes'
tags:
  - Python
  - molecular dynamics
  - membrane proteins
  - GPCR
  - Martini 3
  - CHARMM36
  - GROMACS
  - backmapping
authors:
  - name: Roland Del Mundo
    # orcid: add the author's ORCID here before submitting (a placeholder fails the PDF build)
    affiliation: 1
affiliations:
  - name: City of Hope, Duarte, California, United States
    index: 1
date: 4 November 2026
bibliography: paper.bib
---

# Summary

MembraneForger builds an energy-minimized, all-atom CHARMM36 [@Klauda2010; @Huang2017] membrane–protein system for
GROMACS [@Abraham2015] from a single protein or protein–ligand structure file. Instead of packing lipids around the
protein from a library, it places the protein into a Martini 3 [@Souza2021] cell-membrane model that has already been
simulated for 50 µs around a G protein-coupled receptor (GPCR), cuts the membrane to the protein plus a margin, converts
the lipids to all atoms with mstool [@Kim2023], and then builds the topology, water, 0.15 M NaCl and minimization. Every
stage is followed by a quantitative check, and the final structure (`em.gro`) is written only when all of them pass.
The package ships 18 such membranes (nine around the kappa opioid receptor and nine around GPR139; ten lipid species,
asymmetric leaflets including GM3 and PIP2, about 1,300 lipids each), the force-field files, restrained equilibration
templates and a tool that grades the equilibration trajectory. One build of a GPCR–G-protein complex takes 15 to 30
minutes on eight cores.

# Statement of need

Lipids are not a passive solvent for membrane proteins. Cholesterol, anionic lipids, PIP2 and gangliosides gather at
specific sites on GPCRs and other membrane proteins, and that local composition takes microseconds of lateral diffusion
to form. All-atom simulations usually start from a membrane packed at random around the protein, so whatever lipid
environment they sample has to emerge within a trajectory that is too short for it. Coarse-grained Martini 3
simulations reach that timescale, but converting such a membrane back to all atoms around an all-atom protein is a long
manual procedure: orient and place the protein, remove the coarse-grained protein, cut a box, backmap, repair, solvate,
ionize, minimize, and check that nothing broke along the way. Each step can silently produce a structure that GROMACS
accepts but that is physically wrong.

MembraneForger is for computational biophysicists and structural biologists who want all-atom membrane–protein systems
that start from an equilibrated, protein-adapted lipid environment, with a record of how each system was made and
checked. It turns that procedure into one command (`membraneforger --aa complex.pdb`), runs on a laptop or a Slurm
cluster, and writes a validation report, an audit file and a manifest of settings, software versions and input
checksums with every build, so that several replicate systems in different equilibrated membranes can be built and
compared in the same way.

# State of the field

CHARMM-GUI Membrane Builder [@Jo2009; @Wu2014] is the most widely used route to all-atom membrane–protein systems. It
packs lipids around the protein from a library of conformations and provides equilibration inputs for several MD
engines, but the starting lipid arrangement has no memory of protein–lipid interactions, and it is used mainly through a
web interface. The Martini Maker of CHARMM-GUI [@Qi2015] and insane [@Wassenaar2015] build coarse-grained membranes, and
backward [@Wassenaar2014] converts coarse-grained structures to all atoms with per-molecule mapping files.

mstool [@Kim2023] is the closest tool and the one MembraneForger uses for the conversion itself. It backmaps Martini
lipids to all atoms around a fixed protein and uses mapping and isomeric information to keep chirality and cis/trans
bonds correct. mstool is a backmapping library: it does not orient or embed a protein in an existing membrane, cut a
periodic patch from a larger one, build the GROMACS topology, solvent and ions, or decide whether a finished system is
fit to simulate. Building and testing systems with it showed three problems that a user of mstool alone would have to
find and fix:

- **Configurations can change after backmapping.** mstool holds trans double bonds and amides with restraints only in
  its first refinement stage and its structure check only reports. The ceramide amide of GM3 and PSM has a CHARMM36
  torsion with a zero force constant, so a later restrained minimization could twist it (about −71°) after mstool had
  placed it trans. MembraneForger adds dihedral restraints for these bonds during minimization and tests every
  stereocentre of GM3 and PIP2 geometrically after backmapping and after minimization.
- **Some defects are invisible until `mdrun`.** Two lipids threaded through each other, or atoms placed on top of
  another molecule's atoms, pass backmapping, solvation and ionization and are first noticed when steepest descent stalls
  (one build stopped at Fmax 20,878 kJ mol⁻¹ nm⁻¹ after 50,000 steps on a DOPC carbon 2.4 Å from a POPC carbon).
  MembraneForger checks every non-bonded pair against its own Lennard-Jones Rmin from the force field and every pair of
  bonds of two molecules for a crossing, under periodic boundaries, before minimization; on failure it retries
  backmapping with another seed.
- **Input frames are not always periodic in their box.** A frame written after a rotational fit (`gmx trjconv -fit
  rotxy+transxy`) folds its corners onto the box edges. MembraneForger detects the rotation from overlapping beads,
  rotates the frame back and refuses frames that no rotation fixes.

Beyond these, cutting a patch out of a larger membrane changes the lipid packing at the new edge. MembraneForger cuts
whole lipids by their anchors under periodic boundaries, relaxes the seam, and grades the cut against the uncut membrane
(area per lipid, leaflet composition, lateral radial distribution functions, thickness and protein orientation) before
any atom is placed.

# Software design

MembraneForger is a Python package (about 20,000 lines in about 40 modules) that drives GROMACS and mstool as external
programs. A session object runs the stages in order (orient, frame check, embed, slice, backmap, topology, solvate,
minimize, audit), keeps every intermediate file in `<out>/int/` for inspection and restarts, and can stop after any
stage: `membraneforger slice`, `backmap` and `solvate` return the cut coarse-grained membrane, the backmapped membrane,
or a solvated and minimized system built from any assembled protein–membrane PDB. mstool is imported only by a worker
process, so it can live in its own Python environment. Orientation uses an OPM entry [@Lomize2012], PPM 3.0
[@Lomize2022] or the user's coordinates. Ligands take CGenFF [@Vanommeslaeghe2010] or user topologies, are matched to
them by bond graph, and are frozen during minimization so that docked or experimental poses are kept.

Validation is organized as gates: gate 0 grades the coarse-grained frame, gate 1 the cut membrane, gate 2 the
backmapped membrane and the minimized system (thickness, tilt, interdigitation, core hydration, stereochemistry, ring
threading, contacts and bond crossings, backbone geometry), and gates 3 and 4 the equilibration and production
trajectories. Thresholds are collected in one table of the documentation. `membraneforger --check` runs the pre-flight
part of a build in seconds, reporting what `pdb2gmx` will do with each chain, terminus, histidine, disulfide and ligand.

The test suite has about 590 unit, negative and invariance tests, run in continuous integration together with linting
and a pre-flight run on the bundled example. Full builds, equilibration templates and speed changes were tested on an
HPC cluster with GROMACS 2023.3, 2025.3 and 2026.3. Documentation, including a first-build walkthrough, an input
preparation guide, a validation reference, a command-line reference and troubleshooting, is published at
<https://membraneforger.github.io>.

# Research impact statement

MembraneForger was developed for simulations of class A and class B GPCRs, including the glucagon receptor with
tirzepatide and Gs (the bundled example), GPR139 with a small-molecule agonist, and lipidated receptor constructs
(S-geranylgeranyl and S-palmitoyl cysteine, N-palmitoyl glycine). The 18 bundled membranes come from 18 independent
50 µs Martini 3 simulations, so a user can build the same complex in several independently equilibrated membranes and
treat them as replicates. Since its first release (v0.1.0, July 2026) the project has gone through more than 60 pull
requests and three releases, adding the validation gates described above, pip installation, stage commands,
restrained equilibration templates and the documentation site.

# AI usage disclosure

Claude Code (Anthropic) was used to assist with code optimization and packaging. The author designed the method,
reviewed all changes, ran the validation builds and takes full responsibility for the software and this paper.

# Acknowledgements

The author thanks the developers of Martini, mstool, CHARMM-GUI, GROMACS, OpenMM [@Eastman2017] and OPM/PPM, on whose work
MembraneForger builds. Simulations were run on the City of Hope high-performance computing resources.

# References
