# Stage 4

Stage 4 consumes the canonical Stage 3 structure `outputs/<run_id>/stage3/final_all_atom.pdb`, Stage 3 diagnostic DMS, nonprotein resources, the normalized all-atom protein input, and configured ligand parameters.

Stage 4 owns `gmx pdb2gmx` protein-chain ITP generation and final AA topology
assembly. Rosetta and PyRosetta plugins are optional and disabled in the
generic workflow.

For CHARMM36 terminal capping, Stage 4 uses CHARMM terminal patches rather than
standalone cap residues. The N-terminal acetyl patch is `ACE`; the C-terminal
N-methylamide patch is `CT3`. If `final_all_atom.pdb`, the protein input, or
`pdb2gmx` output contains explicit `ACE`/`NME`/`NMA`/`CT3` cap residues as
separate residues, Stage 4 fails with an actionable error instead of silently
reinterpreting that chemistry.

In the current implementation, terminal patch application is performed by
`gmx pdb2gmx -ter`; MembraneForger answers the terminal prompts with `ACE` and
`CT3`. For that reason, the selected GROMACS CHARMM36 force-field directory must
expose those terminal options in its `pdb2gmx` terminal databases. The selected
CHARMM toppar tree must also define `PRES ACE` and `PRES CT3` in its protein
topology. `forcefield.itp` proves a force-field file is present, but does not by
itself prove terminal patch availability or compatibility between resource
trees.

Before real Stage 4 construction, MembraneForger runs a CHARMM resource
preflight that checks required GROMACS force-field files, required toppar files,
terminal patch definitions, useful release/provenance metadata, and whether the
resource pair exposes a shared release/date marker. If compatibility cannot be
established, Stage 4 fails before using the resources.

The public Stage 4 workflow produces an energy-minimized all-atom system:

```text
minimized_all_atom.gro
minimized_all_atom.top
index_ini.ndx
```

It does not claim restrained all-atom equilibration or production readiness.
Protein atom matching includes chain, residue identity, insertion code, and atom
name. If pdb2gmx output lacks chain identifiers and duplicate residue numbers
make a multichain mapping ambiguous, Stage 4 fails instead of assigning atoms to
the first chain.
