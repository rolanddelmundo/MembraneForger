# MembraneForger

MembraneForger builds an all-atom membrane–protein system for GROMACS from a single structure file. It places your
all-atom protein (or protein–ligand) PDB in a pre-equilibrated Martini 3 coarse-grained cell membrane, converts the
membrane to all atoms with [mstool](https://github.com/ksy141/mstool), and returns a solvated, energy-minimized CHARMM36
system (`em.gro`) that has passed a set of automated checks.

[Your first build](first_build.md){ .md-button .md-button--primary }
[Source on GitHub](https://github.com/rolanddelmundo/MembraneForger){ .md-button }

## How it works

1. **Orient** the protein in the membrane (from an OPM entry, PPM 3.0, or your own coordinates).
2. **Check** the coarse-grained membrane frame (bundled, or your own Martini 3 CG membrane) before anything is placed in it.
3. **Embed** the protein in it.
4. **Cut** the membrane to the protein plus a margin, and check that the cut kept the lipid packing.
5. **Backmap** the lipids to all atoms around the fixed protein.
6. **Build** the CHARMM36 topology, add water and 0.15 M NaCl, and minimize.
7. **Audit** the result. `em.gro` is written only when every check passes.

The output is minimized, not equilibrated; [run restrained equilibration](equilibration.md) before production.

## Install

```bash
git clone https://github.com/rolanddelmundo/MembraneForger.git
cd MembraneForger
pip install ".[figures]"
membraneforger --check --aa membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_prot-lig.pdb --orientation none --out example_out
```

GROMACS (2023.3 or 2025.3 tested) must be installed separately. A GPCR complex builds in about 15 to 30 minutes on 8
cores. The [README](https://github.com/rolanddelmundo/MembraneForger#requirements) lists the requirements in full.

## Where to go next

| I want to | Read |
|---|---|
| Run the bundled example and read the result | [Your first build](first_build.md) |
| Prepare my own PDB | [Preparing the input](preparing_inputs.md), [Ligands](ligands.md) |
| Equilibrate the output | [Equilibration](equilibration.md) |
| Find an option | [Command-line reference](cli_reference.md), [Advanced options](advanced.md) |
| Understand what the checks measure | [Tutorial and validation reference](membraneforger_tutorial.md) |
| Fix a failed build | [Troubleshooting](troubleshooting.md) |
| Work on the code | [Module map](modules.md) |

## Citation and license

If you use MembraneForger, cite the repository ([CITATION.cff](https://github.com/rolanddelmundo/MembraneForger/blob/main/CITATION.cff))
and the tools and force fields it relies on: Martini 3, mstool, CHARMM36 and CHARMM-GUI, CGenFF, GROMACS, and PPM 3.0 or
OPM when used for orientation. MembraneForger is released under the MIT license.
