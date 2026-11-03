# Ligands, ions and other non-protein residues

A ligand is every residue of the `--aa` input whose name is not an amino acid (`config.AMINO`): a small molecule, a
nucleotide such as GDP, an ion such as Mg2+. The protein chemistry comes from CHARMM36 through `pdb2gmx`; a ligand's
chemistry comes only from a GROMACS topology file that you place in the force-field directory. Without it the build
stops at the `inputs` stage:

```
ERROR: inputs: all-atom input complex.pdb: no topology in .../membraneforger/forcefield for residue(s) ['LIG']
```

`python -m membraneforger --check --aa complex.pdb ...` reports the same in seconds, with the ligands it found.

## What is bundled

`membraneforger/forcefield/toppar/` holds the lipids, water and ions of the membrane model, and these ligands:

| Residue | Molecule | Source |
|---|---|---|
| `GDP` | guanosine diphosphate | CHARMM36 (CHARMM-GUI `psf2itp`) |
| `MG` | Mg2+ | CHARMM36 ions |
| `7V7`, `JNJ` | small-molecule ligands of the lab's receptors | CGenFF through CHARMM-GUI FF-Converter |
| `KTZ`, `KRT`, `KSM` | lipidated lysines of peptide agonists | treated as protein residues, see [Advanced options](advanced.md#lipidated-residues) |

Water, Na+, Cl- and K+ must not be in the input; the build adds its own (`SOL`, `SOD`, `CLA`).

## Adding your own ligand

1. **Get a CHARMM-compatible topology.** Run the ligand through [CGenFF](https://cgenff.com) (or CHARMM-GUI's Ligand
   Reader & Modeler, which calls it) with all hydrogens present, then convert the `.str`/`.rtf` to GROMACS with
   CHARMM-GUI's FF-Converter or `cgenff_charmm2gmx.py`. The result is an `.itp` with `[ moleculetype ]`, `[ atoms ]`,
   `[ bonds ]`, `[ pairs ]`, `[ angles ]`, `[ dihedrals ]` and CGenFF atom types (`CG331`, `NG2S1`, ...), as in
   `toppar/JNJ.itp`.
2. **Name it after the residue.** The file is `toppar/<RESNAME>.itp` and the `[ moleculetype ]` name is `<RESNAME>`,
   the residue name used in your PDB (three characters, upper case). `GDP.itp` declares `GDP`; a ligand named `LIG`
   in the PDB needs `LIG.itp` declaring `LIG`.
3. **Put it in the force-field directory.** Either copy it into `membraneforger/forcefield/toppar/` of your clone, or
   copy the whole `forcefield/` directory somewhere, add the file there and pass `--toppar /that/directory` (or set
   `MEMBRANEFORGER_TOPPAR`). `--toppar` must still contain `charmm36.ff/`, `toppar/` and `protein_parameters/`.
4. **Extra parameters.** If CGenFF produced parameters that are not in CHARMM36 (the `.prm` part of its output, with
   penalty scores), they must be visible to `grompp`: append them, converted to GROMACS units, to
   `charmm36.ff/ffbonded.itp` in your `--toppar` copy (`[ angletypes ]`, `[ dihedraltypes ]`, ...). `toppar/JNJ.prm.itp`
   shows that format; a `.prm.itp` next to the ligand is **not** read automatically.
5. **Check the pose.** Ligand heavy atoms are frozen in both minimizations
   ([Advanced options](advanced.md#ligands-during-minimization)), so a ligand that overlaps a lipid or a residue in the
   input stays there and the build fails the Fmax or contact check naming the pair. Fix the pose or the residue in
   the input.

### What the build does with the ligand

- **Atom names do not have to match the `.itp`.** The ligand is matched to its topology by the bond graph of its heavy
  atoms (`topology.match_by_graph`): element and connectivity, computed from the coordinates. A residue whose heavy-atom
  graph is not isomorphic to the `.itp` (a different tautomer or protonation, a missing atom) stops the build with
  `heavy-atom graph is not isomorphic to its ITP`.
- **Hydrogens may be missing.** A hydrogen the `.itp` has and the input lacks is placed on its parent atom at a fixed
  direction and corrected by minimization (ligand hydrogens are not frozen). A hydrogen the input has and the `.itp`
  lacks stops the build (`source hydrogens have no topology counterpart`): the protonation state of the input and of
  the topology must agree.
- **Charge.** The ligand charge is part of the system charge, which the ions neutralize (0.15 M NaCl on top,
  `--counterion-mode`).
- **Restraints for equilibration.** If the `.itp` has no `#ifdef POSRES` block, the build writes `toppar/posre_<RESNAME>.itp`
  restraining every heavy atom with `POSRES_FC_BB`, so the ligand follows the backbone schedule of
  [the equilibration](equilibration.md). An `.itp` that brings its own block (as `JNJ.itp` does, also with
  `POSRES_FC_BB`) keeps it.
- **Index groups.** Ligands are in `Protein_LIG`; their heavy atoms form the `LIG` group that the minimizations
  freeze. `run_manifest.json` lists the ligands of the build.

## Ions and cofactors bound to the protein

A structural ion (Mg2+, Zn2+, Ca2+) is a ligand like any other: keep it in the input with its CHARMM residue name
(`MG`, `ZN`, `CAL`) and make sure `toppar/<NAME>.itp` exists; `MG.itp` is bundled, others can be written in the same
one-atom form. Do not keep bulk ions or crystallization additives: the input check refuses `NA`, `CL`, `K`, water and
their CHARMM names, and anything else without a topology.

## Covalently attached groups

A group bonded to the protein (a prenyl, palmitoyl or myristoyl chain, a lipidated lysine) is not a ligand but a
modified residue, which `pdb2gmx` must know. The bundled ones are `CYSG`, `CYSP`, `GLYP`, `KTZ`, `KRT` and `KSM`
([Advanced options](advanced.md#lipidated-residues)). Another modified residue needs an RTP entry in `charmm36.ff/merged.rtp`
and hydrogens in `merged.hdb`; see `topology.add_lipidated_residue` for how the bundled ones are added.
