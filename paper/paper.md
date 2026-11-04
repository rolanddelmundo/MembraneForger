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
    orcid: 0009-0007-1702-2846
    affiliation: 1
affiliations:
  - name: City of Hope, Duarte, California, United States
    index: 1
date: 4 November 2026
bibliography: paper.bib
---

# Summary

MembraneForger builds all-atom membrane–protein systems for GROMACS [@Abraham2015]. The input is one protein or
protein–ligand structure file. The output is an energy-minimized CHARMM36 [@Klauda2010; @Huang2017] system.

The protein is placed into a Martini 3 [@Souza2021] membrane that was already simulated for 50 µs around a G
protein-coupled receptor (GPCR). The membrane is cut to the protein plus a margin. The lipids are converted to all
atoms with mstool [@Kim2023]. MembraneForger then builds the topology, adds water and 0.15 M NaCl, and minimizes.

Each stage ends with a quantitative check. The final structure is written only when every check passes.

The package includes 18 membranes: nine around the kappa opioid receptor and nine around GPR139. Each has about 1,300
lipids of ten species, with asymmetric leaflets that include GM3 and PIP2. The package also includes the force-field
files, restrained equilibration templates and a tool that grades equilibration trajectories. A GPCR–G-protein build
takes 15 to 30 minutes on eight cores.

# Statement of need

Lipids are not a passive solvent. Cholesterol, anionic lipids, PIP2 and gangliosides bind specific sites on GPCRs.
This local composition takes microseconds to form.

All-atom simulations usually start from lipids packed at random around the protein. The protein-adapted lipid
environment must then form during a trajectory that is too short.

Martini 3 simulations reach the needed timescale. Converting such a membrane to all atoms is a long manual task. The
user must orient and place the protein, cut a box, backmap, repair, solvate, ionize and minimize. Any step can give a
structure that GROMACS accepts but that is physically wrong.

MembraneForger does this in one command (`membraneforger --aa complex.pdb`). It runs on a laptop or a Slurm cluster.
Every build writes a validation report, an audit file and a manifest of settings, software versions and input
checksums. Replicate systems in different equilibrated membranes are built and checked the same way. The target users
are computational biophysicists and structural biologists.

# State of the field

CHARMM-GUI Membrane Builder [@Jo2009; @Wu2014] is the most used tool for all-atom membrane–protein systems. It packs
lipids from a conformer library and writes equilibration inputs for several MD engines. Its starting lipid arrangement
has no memory of protein–lipid interactions. It is used mainly through a web interface.

CHARMM-GUI Martini Maker [@Qi2015] and insane [@Wassenaar2015] build coarse-grained membranes. backward
[@Wassenaar2014] converts coarse-grained structures to all atoms with per-molecule mapping files.

mstool [@Kim2023] is the closest tool. MembraneForger uses it for the conversion step. mstool backmaps Martini lipids
around a fixed protein and uses isomeric information to keep chirality and cis/trans bonds correct. It is a
backmapping library. It does not embed a protein in an existing membrane, cut a periodic patch, build the GROMACS
topology, solvent and ions, or judge whether the final system is fit to simulate.

# Contributions

## Backmapping of GM3, PIP2 and PSM

Cell-membrane models need GM3, PIP2 and sphingomyelin. Their sugars, inositol ring and amide are hard to backmap correctly.
To build these lipids reliably, we corrected and extended their Martini 3 to CHARMM36 backmapping.

- **GM3.** In Martini 3, GM3 is one molecule (`DPG3`). In CHARMM36 it is four residues: ceramide, glucose, galactose
  and sialic acid. MembraneForger supplies the mapping and the OpenMM template that join them. The starting mapping
  had glucose and galactose beads swapped, six wrong chirality definitions, no sialic acid C7 definition and no
  ceramide trans definitions. All were corrected. All 16 sugar stereocentres and both ceramide trans bonds of every
  GM3 are now tested after backmapping and after minimization.
- **PIP2 (`SAP6` in Martini 3, `SAPI25` in CHARMM36).** The ring-chair restraints of the inositol stay satisfied when a
  hydrogen and a hydroxyl swap places at one carbon. MembraneForger tests all six inositol stereocentres directly. This
  found a real inverted C1 in a build whose chair restraints were all satisfied.
- **PSM and the GM3 ceramide.** The CHARMM36 amide torsion of both lipids has a zero force constant, so nothing holds
  the amide trans. In one build, restrained minimization twisted a GM3 amide that mstool had placed trans to about
  −71°. In 349 archived builds, about 3 % of PSM amides ended below 150°. MembraneForger adds a trans restraint for
  this amide during restrained minimization only. With the restraint, all 27 GM3 amides of a test build stayed trans (172.6–179.2°), compared with four below 150° without it.

These checks matter because a wrong configuration never corrects itself in MD. A build with one inverted centre or cis amide is refused, and
backmapping is retried with a new seed.

## Improvements over mstool alone

- **Configurations are held after backmapping.** mstool restrains trans bonds only in its first refinement stage, and
  its structure check only reports. MembraneForger restrains them through minimization and tests them at the end.
- **Defects that appear only in `mdrun` are caught first.** Threaded lipids and overlapping atoms pass backmapping,
  solvation and ionization. One build stalled at Fmax 20,878 kJ mol⁻¹ nm⁻¹ after 50,000 steps, on a DOPC carbon 2.4 Å
  from a POPC carbon. MembraneForger tests every non-bonded pair against its Lennard-Jones Rmin and every pair of bonds
  for a crossing, under periodic boundaries. A failed test triggers a new backmapping seed.
- **Frames are checked before use.** A frame saved after a rotational fit folds its corners onto the box edges.
  MembraneForger detects the rotation and reverses it. It refuses a frame that no rotation fixes.
- **The cut membrane is graded.** Cutting a patch changes lipid packing at the new edge. MembraneForger cuts whole lipids
  under periodic boundaries and relaxes the seam. It compares the cut with the uncut membrane before any atom is placed:
  area per lipid, leaflet composition, lateral radial distribution functions, thickness and protein orientation.
- **The full system is built.** mstool returns a membrane. MembraneForger returns a solvated, ionized, minimized and
  audited GROMACS system with its topology.

## Simple and feasible to run

- **One input and one command.** `membraneforger --aa complex.pdb` is a full build. No coarse-grained simulation is
  needed: the 18 bundled 50 µs membranes are ready to use.
- **Modest hardware.** A GPCR–G-protein build takes 15 to 30 minutes on eight cores and about 1.5 GB of memory. It runs
  on a laptop or a Slurm cluster.
- **Standard installation.** `pip install .` installs the package and its Python requirements. GROMACS is the only
  other requirement. PPM 3.0 is needed only for structures without an OPM entry.
- **Problems found in seconds.** `membraneforger --check` tests the software, the input and the options before a build
  starts.
- **Clear results.** Each build ends in PASS or FAIL. A failure names the stage and the check.

# Software design

MembraneForger is a Python package of about 20,000 lines in 40 modules. It runs GROMACS and mstool as external
programs. mstool runs in a separate worker process, so it can use its own Python environment.

A session object runs the stages in order: orient, frame check, embed, slice, backmap, topology, solvate, minimize and
audit. Every intermediate file is kept in `<out>/int/` for inspection and restarts. A run can stop after any stage.
`membraneforger slice`, `backmap` and `solvate` return the cut membrane, the backmapped membrane, or a minimized
system from any assembled protein–membrane PDB.

Orientation uses an OPM entry [@Lomize2012], PPM 3.0 [@Lomize2022] or the user's coordinates. Ligands use CGenFF
[@Vanommeslaeghe2010] or user topologies, matched by bond graph. Ligands are frozen during minimization, so docked or
experimental poses are kept.

Validation is organized in gates. Gate 0 grades the coarse-grained frame. Gate 1 grades the cut membrane. Gate 2 grades
the backmapped membrane and the minimized system: thickness, tilt, interdigitation, core hydration, stereochemistry,
ring threading, contacts, bond crossings and backbone geometry. Gates 3 and 4 grade equilibration and production
trajectories. All thresholds are in one documentation table. `membraneforger --check` runs the pre-flight checks in
seconds.

The test suite has about 590 unit, negative and invariance tests. Continuous integration runs them with linting and a
pre-flight run on the bundled example. Full builds were tested on an HPC cluster with GROMACS 2023.3, 2025.3 and
2026.3. Documentation is at <https://membraneforger.github.io>. It includes a first-build walkthrough, an input guide,
a validation reference, a command-line reference and troubleshooting.

# Research impact statement

MembraneForger has been used in the author's laboratory at City of Hope since November 2025. It was released publicly
in July 2026.

It was developed for simulations of class A and class B GPCRs. These include the glucagon receptor with tirzepatide
and Gs (the bundled example), GPR139 with a small-molecule agonist, and lipidated receptors (S-geranylgeranyl and
S-palmitoyl cysteine, N-palmitoyl glycine).

The 18 bundled membranes come from 18 independent 50 µs Martini 3 simulations. A user can build one complex in several
independently equilibrated membranes and treat them as replicates.

Since v0.1.0 the project has had more than 60 pull requests and three releases. These added the validation gates, pip
installation, stage commands, equilibration templates and the documentation site.

# AI usage disclosure

Claude Code (Anthropic) was used for code optimization and packaging. The author designed the method, reviewed all
changes, ran the validation builds, and takes full responsibility for the software and this paper.

# Acknowledgements

The author thanks the developers of Martini, mstool, CHARMM-GUI, GROMACS, OpenMM [@Eastman2017] and OPM/PPM.
Simulations were run on City of Hope high-performance computing resources. Funding: [to be added by the author; state
whether the sponsor had any role].

# Conflicts of interest

The author declares no conflicts of interest.

# References
