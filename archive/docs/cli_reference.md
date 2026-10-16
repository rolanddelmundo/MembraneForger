# Command-line reference

Run MembraneForger from the repository folder:

```bash
python -m membraneforger --aa complex.pdb --orient-chain R --nterm-side out --out my_build
```

`python -m membraneforger --help` prints the same options with their current defaults. Lengths on the command line
are in Å unless the option says nm.

## Inputs and output

| Option | Meaning | Default |
|---|---|---|
| `--aa PDB` (`--all-atom`) | All-atom protein or protein–ligand structure. The only source of protein and ligand chemistry. | required, unless `--membrane` |
| `--cg 1\|2\|NAME\|FILE` (`--coarse-grain`) | The membrane. `1` picks a bundled GPR139 frame at random, `2` a bundled kappa opioid receptor frame. A name (`GPR1` … `GPR9`, `KOR1` … `KOR9`, any case) picks that bundled frame, so the build is repeatable. A file (or `custom=FILE`) is a Martini 3 frame of your own complex. | `1` |
| `--runs N` | Build the complex N times, each time in a different bundled membrane of the `--cg 1` or `--cg 2` series (its first N frames, so `--runs 5` with `--cg 1` is GPR1 to GPR5), one system per run in `<out>/<frame>/`. Needs `--cg 1` or `2`; at most 9. | `1` |
| `--list-membranes` | Print the bundled frames by name and exit. | — |
| `--embed` | With `--cg FILE`, embed the protein into that frame's membrane instead of fitting it onto the frame's own protein. Always on for the bundled frames (`1`, `2` or a name). | off |
| `--embed-site hole\|free` | Where an embedded complex goes. `hole` is where the frame's receptor was. `free` is unbroken membrane at least 1 nm from that receptor, for a much smaller protein such as one transmembrane helix. | `hole` |
| `-o`, `--out DIR` | Output directory (with `--runs`, the parent of one directory per run). | `./<membrane file stem>_membraneforger`, with `--runs` `./<structure stem>_membraneforger` |
| `--name NAME` | System name used in logs and in `[ system ]`. | output directory name |
| `--membrane PDB` | Skip orientation and backmapping and start from an assembled protein + membrane PDB with a `CRYST1` record, for example `<out>/int/membrane.pdb` of an earlier build. Replaces `--aa` and `--cg`. | — |

## Membrane orientation

The complete complex (partners, peptides, ligands, ions) is moved as one rigid body. See
[Advanced options](advanced.md#membrane-orientation) for how each mode decides.

| Option | Meaning | Default |
|---|---|---|
| `--orientation auto\|ppm\|opm\|none` | `auto` uses the exact OPM entry when a PDB ID is known and matches, otherwise local PPM 3.0. `none` uses your coordinates as given; the membrane normal must already be along z. | `auto` |
| `--orient-chain C` | Chain that spans or associates with the membrane. Required when the input has more than one protein chain. | the only protein chain |
| `--orient-chains A,B` | Several anchor chains; the first one sets the N-terminus side. | — |
| `--orient-residues FIRST-LAST` | Orient on these residues of the anchor chain only, for example a transmembrane helix `343-363` (several: `343-363,370-380`). With `--orientation none` it names the membrane-embedded segment instead, and the bilayer centre goes at the midpoint of that segment's CA z range. | whole anchor chain |
| `--nterm-side in\|out` | Side of the membrane on which residue 1 of the first anchor chain, as it appears in your file, lies. Needed by PPM when no OPM entry gives it. | read from OPM, else required |
| `--pdb-id ID` | Exact PDB ID for the OPM reference. | `HEADER` record of the input (never the file name) |
| `--ppm-exe PATH` | PPM 3.0 executable (`immers`) with `res.lib` next to it. | `$MEMBRANEFORGER_PPM`, then `immers`/`ppm3` on `PATH` |
| `--bilayer-z Z` | With `--orientation none`, z of the bilayer centre in your structure, in Å. | found from the hydrophobic belt |
| `--ppm-membrane CODE` | Advanced: PPM 3.0 membrane code such as `PMm` or `GnI`. | undefined flat bilayer |
| `--ppm-heteroatoms` | Advanced: also submit the anchor chains' heteroatoms to PPM. | off |
| `--opm-file PDB` | An already downloaded OPM/OPRLM coordinate file, for offline use. | — |
| `--opm-cache DIR` | Cache for downloaded OPM files. | `~/.cache/membraneforger/opm` |

## Box and membrane composition

| Option | Meaning | Default |
|---|---|---|
| `--xy-buffer NM` | Membrane kept around the complex on each side in x and y when the box is automatic. Large complexes and single helices often need 2–3. | `1.0` nm |
| `--box X Y Z` | Opt-in box edges in Å. x and y may only be smaller than the membrane patch, which is then cut around the complex. z must leave the default water padding. | automatic |
| `--dellipid LIPID` | Remove every molecule of this lipid (repeatable): `CHOL`, `POPC`, `DOPC`, `POPE`, `DOPE`, `POPS`, `DOPS`, `PSM`, `DPG3`, `SAP6`. The aliases `GM3` (for `DPG3`) and `PIP2`/`SAPI25` (for `SAP6`) are accepted. | — |
| `--addlipid LIPID` | Turn the lipids removed by `--dellipid` into this lipid instead: `POPC`, `DOPC`, `POPE`, `DOPE`, `POPS`, `DOPS`, `PSM`. | — |

## Slice validation

The slice gate compares the cut membrane with the membrane it was cut from, before the expensive backmapping step.

| Option | Meaning | Default |
|---|---|---|
| `--apl-validate yes\|no` | Stop before backmapping when a leaflet's area per lipid changes by more than `--apl-tolerance`. | `yes` |
| `--apl-tolerance PERCENT` | Largest accepted change of the leaflet area per lipid across slicing (a warning is logged above 3 %). | `5.0` |
| `--rdf-validate yes\|no` | Compare the lateral headgroup RDF of the slice with its source. Warns only. | `yes` |
| `--no-slice-offset` | Keep the slice window centred on the complex instead of choosing the offset (up to 0.5 nm) that best preserves lipid density. | off |

## Programs, resources and run control

| Option | Meaning | Default |
|---|---|---|
| `--gmx CMD` | GROMACS command, for example `gmx_mpi` or a container wrapper. | `gmx`, `gmx_mpi` or `gmx_d` on `PATH` |
| `--mstool-python PATH` | Python interpreter that can import mstool. | the interpreter running MembraneForger |
| `--ntomp N` | Threads. | all cores |
| `--nsteps N` | Maximum energy-minimization steps. | `50000` |
| `--toppar DIR` | Force-field directory with `charmm36.ff/`, `toppar/` and `protein_parameters/`. | `$MEMBRANEFORGER_TOPPAR`, else `forcefield/` |
| `--data DIR` | Backmapping data directory with `map.dat` and extra force-field XML such as `GM3.xml`. | `$MEMBRANEFORGER_DATA`, else `backmap_data/` |
| `--fit-max-core-rmsd A` | Largest accepted core RMSD of the rigid fit onto your own frame's protein. | `5.0` Å |
| `--fit-min-core-fraction F` | Smallest accepted fraction of backbone pairs kept by that fit. | `0.5` |

## Environment variables

| Variable | Used for |
|---|---|
| `MEMBRANEFORGER_PPM` | default for `--ppm-exe` |
| `MEMBRANEFORGER_TOPPAR` | default for `--toppar` |
| `MEMBRANEFORGER_DATA` | default for `--data` |

## Equilibration analysis

After you have equilibrated the system, `python -m membraneforger.equilibration --help` lists the options of the
trajectory gates (convergence and equilibrium properties). The [tutorial](membraneforger_tutorial.md#59-gates-3-and-4-equilibration-and-equilibrium-properties)
describes what they measure.
