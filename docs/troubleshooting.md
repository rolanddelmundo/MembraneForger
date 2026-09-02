# Troubleshooting

## Missing GROMACS

`doctor` reports a GROMACS failure when no usable `gmx` executable is found.
Install GROMACS or set `GMX_BIN=/path/to/gmx`.

## Missing CHARMM36

Stage 4 needs a local CHARMM36 GROMACS force-field directory containing
`forcefield.itp`, `ffbonded.itp`, `ffnonbonded.itp`, residue and terminal
databases, water, and ion parameters. This resource is not bundled. Set:

```bash
export MEMBRANEFORGER_CHARMM36_ROOT=/path/to/charmm36.ff
```

or configure `stage4.resources.charmm36_forcefield`.

## Missing Membrane Parameters

Stage 4 also needs membrane topology/parameter files for the all-atom membrane
molecules, including compatible CHARMM protein topology/parameter files with
`PRES ACE` and `PRES CT3`. Set:

```bash
export MEMBRANEFORGER_MEMBRANE_TOPPAR=/path/to/charmm36_membrane
```

or configure `stage4.resources.membrane_toppar`.

## Unsupported Molecular Component

Selected unsupported `HETATM` or nonprotein components fail instead of being
removed. Remove the component from the selected system, provide explicit
validated parameters, or choose a supported chemistry path.

## Unknown Config Key

Unknown keys fail because misspelled scientific options are dangerous. Compare
the key with `config/workflow.yaml` or regenerate a run config with
`./run_pipeline.sh init`.

## Conflicting Config Values

Public config sections are authoritative. If a public value conflicts with an
explicit legacy `stage1` through `stage4` duplicate, validation fails. Keep one
value and remove the conflicting duplicate.

## Ambiguous Stage 4 Atom Matching

Stage 4 matches protein atoms by chain, residue number, insertion code, and atom
name. If `pdb2gmx` output lacks chain identifiers and residue numbers are
duplicated across chains, Stage 4 fails rather than guessing.

## `run_cg_production: true`

CG production is not implemented in this public workflow. Set:

```yaml
simulation:
  run_cg_production: false
```

## Cleanup Refusal

`clean --confirm` only removes marked MembraneForger run directories beneath
`runs/`. It refuses the repository root, parent directories, arbitrary paths,
unmarked directories, and symlinks.

## CHARMM Terminal Capping

MembraneForger uses CHARMM's native terminal-patch system. An acetylated
N-terminus uses the `ACE` patch, while an N-methylamide C-terminus uses the
`CT3` patch. These modify the terminal amino acid in the CHARMM topology rather
than being treated as ordinary extra protein residues. Do not manually add an
`NME` residue for the CHARMM36 workflow.

In the current implementation, `gmx pdb2gmx -ter` applies these terminal patch
choices. If Stage 4 reports that `n_terminal_patch=ACE` or
`c_terminal_patch=CT3` is not available in the selected `charmm36.ff` terminal
database, use a CHARMM36 GROMACS force-field directory that actually provides
those `pdb2gmx` terminal options. A directory containing `forcefield.itp` alone
is not enough to prove terminal patch compatibility. A standalone `[ NME ]`
residue definition inside a force-field file is not itself an error.

## CHARMM Resource Pair Not Compatible

If preflight reports that resource compatibility is not established, the
configured GROMACS CHARMM36 tree and toppar tree do not expose a shared
release/date marker. Provide a paired CHARMM36 GROMACS force-field directory and
matching toppar resources from an established distribution, or perform and
document a local scientific compatibility audit before using them.
