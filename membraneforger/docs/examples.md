# Examples

Each example is one command. Run `membraneforger --check` with the same options first. It takes seconds and tells you
what is missing. Every option is in the [command-line reference](cli_reference.md).

## The bundled example

A glucagon receptor and G protein complex (GCGR-Gs) with its peptide agonist, GDP and Mg2+, already oriented, with a
matching coarse-grained membrane. Run it from the repository folder.

```bash
membraneforger \
    --aa membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_prot-lig.pdb \
    --cg membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_cg_cellmem.gro \
    --orientation none --out example_out
```

[Your first build](first_build.md) walks through the output.

## Your protein in a bundled membrane

Name the membrane-spanning chain and the side its N terminus is on. MembraneForger orients the complex (OPM entry or
PPM 3.0) and places it in one of the bundled 18 x 18 nm plasma membranes.

```bash
membraneforger --aa complex.pdb --orient-chain R --nterm-side out --out build
```

`--cg 1` (the default) and `--cg 2` pick a random frame from the GPR139 and kappa opioid receptor series. To make a
build repeatable, name the frame:

```bash
membraneforger --list-membranes
membraneforger --aa complex.pdb --cg KOR_Membrane3 --orient-chain R --nterm-side out --out build
```

## Several independent starting membranes

The bundled frames come from independent 50 µs coarse-grained simulations of the same membrane model. `--runs N`
builds the protein in the first N frames of a series, so the systems have the same composition but different lipid
arrangements. One system is written per frame, in `<out>/<frame>/`.

```bash
membraneforger --aa complex.pdb --cg 2 --runs 5 --orient-chain R --nterm-side out --out runs
```

These are independent starting configurations. Replicates of one system in the usual MD sense come from running the
same `em.gro` several times with new starting velocities.

## A small or single-pass membrane protein

A single transmembrane helix goes on intact bilayer, away from the hole the bundled receptor left, with a wider margin
of lipids around it.

```bash
membraneforger --aa sp.pdb --orientation none --orient-residues 152-177 --nterm-side out \
    --embed-site free --xy-buffer 3 --out build
```

See [Choosing the membrane](advanced.md#choosing-the-membrane).

## Your own coarse-grained membrane

Any Martini 3 membrane frame works: a different composition, a larger patch, or a frame of your own complex.

```bash
membraneforger --aa complex.pdb --cg my_frame.gro --embed --out build
```

## A different lipid composition

Remove a lipid type from the bundled membrane, or convert it to another twelve-bead Martini 3 lipid (POPC, DOPC, POPE,
DOPE, POPS, DOPS or PSM). The change is reported per leaflet.

```bash
membraneforger --aa complex.pdb --dellipid DOPC --addlipid POPC --orient-chain R --nterm-side out --out build
```

## One stage on its own

The same pipeline can stop after a stage or start from your own assembled system. The options are the same as for a
build.

| Command | Stops after | Main output |
|---|---|---|
| `membraneforger slice ...` | cutting the membrane to the protein plus a buffer | `slice_cg.gro`, `placed_complex.pdb` |
| `membraneforger backmap ...` | converting the membrane to all atoms | `membrane.pdb` |
| `membraneforger solvate --membrane system.pdb ...` | starts from your protein + membrane PDB, then solvates and minimizes | `em.gro`, `topol.top` |

```bash
membraneforger slice --aa complex.pdb --xy-buffer 2.0 --out cut
membraneforger backmap --aa complex.pdb --cg my_frame.gro --out aa
membraneforger solvate --membrane system.pdb --out build
```

`solvate` needs a `CRYST1` record and CHARMM36 lipid names. See
[Restarting after backmapping](advanced.md#restarting-after-backmapping).

## Equilibrate and grade the result

```bash
MF="$(python -c 'import membraneforger, os; print(os.path.dirname(membraneforger.__file__))')"
bash "$MF/mdp/run_equilibration.sh" build
```

This runs six restrained steps (1.875 ns) and writes the production input. After a production run,
`membraneforger equilibration` grades box convergence, area per lipid, thickness and hydration. See
[Equilibration](equilibration.md).
