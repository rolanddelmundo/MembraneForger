# Preparing the all-atom input

Everything in the `--aa` PDB is simulated. Nothing missing is added except hydrogens, and wrongly named residues
are refused. Check a file before building:

```bash
membraneforger --check --aa complex.pdb --orient-chain R --nterm-side out --out my_build
```

It lists the chains, termini, histidines, disulfides and ligands it found, checks GROMACS, mstool and PPM, and ends
with `CHECK PASS` or `CHECK FAIL`.

## Checklist

| Keep | Remove |
|---|---|
| every protein chain you want simulated, each with its own chain ID | water (`HOH`, `SOL`, `TIP3`, `WAT`), bulk ions (`NA`, `CL`, `K`, `SOD`, `CLA`, `POT`) |
| peptides, nucleotides, structural ions and small-molecule ligands with a topology ([Ligands](ligands.md)) | crystallization additives, detergents, lipids of the crystal (the membrane is built around the protein) |
| hydrogens if you have them (they are not required; see below) | alternate locations other than A (altloc B is dropped with a note), extra `MODEL` records |
| the `HEADER` record if the file is a deposited structure (its PDB ID selects the OPM reference) | residue insertion codes (`52A`): renumber the chain |

The file must have one model, no duplicate atoms (same chain, residue and atom name twice) and finite coordinates.

### Chains and termini

- A chain is everything with one chain ID (or segment ID if the chain column is empty). A missing loop is joined
  across the gap by a peptide bond, so model the loop or split the fragments into separate chains.
- The C terminus of every chain must be recognizable: a terminal carboxylate with `OXT` (or `OT1`/`OT2`), or a
  CHARMM cap (`CT3`, the N-methylamide `NT`, `CAT`, ...; `CT2`, the amide `NT` with `HT1`/`HT2`). A chain whose last
  residue has no `OXT` stops the build with `unrecognized C terminus`; add the atom (any PDB tool that completes
  terminal oxygens will do) or cap the chain. The N terminus is `NH3+` unless the first residue carries an `ACE` cap
  (`CAY`, `CY`, `OY`) or is an N-acylated residue (`GLYP`).
- Peptide ligands (an agonist bound to a receptor) are protein chains like any other and get their own termini.

### Residue names and protonation

- Histidines named `HIS` are built as `HSD` (neutral, proton on ND1). To choose, name them `HSD`/`HID`, `HSE`/`HIE`
  (proton on NE2) or `HSP`/`HIP` (doubly protonated). `--check` counts the `HIS` residues it will rename.
- Aspartate, glutamate, lysine, arginine and the termini are built in their standard charge states; CHARMM36 names for
  other states are not supported by the renaming table.
- Cysteines whose SG atoms are within 3.0 A are built as a disulfide (the bond is set to 2.03 A if it was stretched);
  `--exclude-disulfide A:860-A:1252` keeps a pair reduced. A cysteine within 3.0 A of two partners stops the build.
  Inter-chain disulfides are not supported (each chain has its own topology).
- Lipidated residues (`CYSG`, `CYSP`, `GLYP`, `KTZ`, `KRT`, `KSM`) are recognized by name and need every heavy atom of
  the lipid present ([Advanced options](advanced.md#lipidated-residues)). The S-geranylgeranyl group must be all-E;
  `--prenyl-configuration fix` repairs a Z bond.

### Hydrogens and missing atoms

- Protein hydrogens are not needed: `pdb2gmx` is run with `-ignh` and rebuilds them for the chosen protonation states.
  Ligand hydrogens may be missing too; they are placed and relaxed in minimization. A ligand hydrogen the topology does
  not have is an error ([Ligands](ligands.md)).
- Missing heavy atoms are **not** rebuilt, and a truncated side chain (common in cryo-EM models) is built
  incomplete without an error. Complete side chains first (PDBFixer, Modeller, ChimeraX or CHARMM-GUI PDB Reader).
- Missing residues are not added: each chain ends where the file ends.

### Orientation inputs

The whole complex moves as one rigid body, so ligands and partner proteins keep their positions.

- `--orient-chain R` names the chain that spans the membrane; required when the input has several protein chains.
  `--orient-chains A,B` for a dimer.
- `--nterm-side in|out` is the side of the membrane on which residue 1 of that chain, *as it appears in your file*,
  lies. A GPCR whose signal peptide was removed is `out`; a receptor model that still carries its signal peptide or
  starts with an N-terminal tag is usually `in`. Needed by PPM 3.0 unless an exact OPM entry gives it.
- `--orientation none` takes the coordinates as given: the membrane normal must be along z and the extracellular side
  at +z, as in files downloaded from OPM or exported by CHARMM-GUI. A PPM web-server result works after its `DUM` atoms
  are removed.
- `--orientation auto` (default) downloads the OPM entry when the structure's PDB ID is known (`--pdb-id` or the
  `HEADER` record) and matches the anchor chain; otherwise it runs PPM 3.0 locally, which must be installed
  ([Advanced options](advanced.md#ppm-30)). Offline, pass `--opm-file` or use `--orientation none`.

### Size and the bundled membranes

The bundled membranes are 18.3 x 18.3 nm and have a hole where a GPCR was.

- GPCR-sized protein: goes into the hole (`--embed-site hole`, the default).
- Single transmembrane helix: use `--embed-site free` and a larger `--xy-buffer`.
- Much larger protein (transporter dimer, channel): use your own Martini 3 frame (`--cg FILE --embed`,
  [Advanced options](advanced.md#choosing-the-membrane)).

Soluble partners such as a G protein or antibody fragment are fine; the box height is sized to the complex.

## A worked example

The bundled GCGR complex (`membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_prot-lig.pdb`) is a prepared input:
five chains (receptor R with ACE/CT3 caps, peptide P with a CT2 amide, Gs alpha A, beta B and gamma C), the ligands
GDP and MG, four disulfides in R, a `CYSG` on the gamma C terminus, 28 histidines left as `HIS`, and no water or bulk
ions. `membraneforger --check --aa membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_prot-lig.pdb --orientation none --out x`
shows how each is handled.
