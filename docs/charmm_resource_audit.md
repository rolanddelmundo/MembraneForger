# CHARMM Resource Audit

Audit date: 2026-09-02

This audit inspected the exact local resource trees:

- `resources/external/charmm36.ff`
- `resources/forcefields/charmm36/toppar`

The parameter files were not replaced, regenerated, upgraded, normalized, or
modified.

## Public Release Decision

Neither tree is vendored in the public repository. Both remain ignored and
must be provided by users or maintainers who have compatible local resources and
appropriate redistribution/use rights.

The local `charmm36.ff` tree is closest to a November 2018 CHARMM36 GROMACS
force-field tree, but it has substantial local differences and additions. The
local `toppar` tree is a mixed CHARMM/CHARMM-GUI-derived parameter tree, with
CHARMM36/CHARMM36m protein markers and CHARMM-GUI FF-Converter output.

Neither exact tree matched the authoritative upstream distributions checked.
Because redistribution terms for these exact local trees were not established,
they are blocked from public vendoring.

## Compatibility Status

The two local trees cannot currently be treated as a canonical paired force
field. Their compatibility is not proven by path existence or by the presence of
individual parameter files. Stage 4 preflight therefore requires a paired,
compatible CHARMM36 GROMACS force-field directory and matching CHARMM
topology/parameter resource tree before real all-atom construction proceeds.

Local/private scientific validation may still use these resources, but any such
run must record the exact paths, resource metadata, tree hashes, GROMACS
version, and MembraneForger commit. Use:

```bash
./run_pipeline.sh scientific-integration \
  --config runs/<run>/config/config.yaml \
  --run-dir runs/<run> \
  --preflight-only
```

Remove `--preflight-only` only when you intend to run the real stage-gated
scientific workflow.

## Terminal Patches

The inspected `toppar` tree defines the CHARMM terminal patches required by
MembraneForger:

- N-terminus acetyl patch: `PRES ACE`
- C-terminus N-methylamide patch: `PRES CT3`

These are CHARMM topology patches. They modify the terminal amino-acid residue
in the topology; they are not ordinary extra protein residues added as
standalone `ACE` or `NME` residues.

In the current Stage 4 implementation, `gmx pdb2gmx -ter` applies the protein
terminal choices. MembraneForger feeds `ACE` for the N-terminus and `CT3` for
the C-terminus into that prompt sequence. Therefore the GROMACS `charmm36.ff`
terminal databases must expose those options for this implementation. The
presence of a standalone `[ NME ]` residue definition inside a force-field file
is not itself an error; the prohibited behavior is MembraneForger synthesizing
standalone cap residues in the protein coordinate handoff.

## Public User Requirement

Public users must provide their own mutually compatible CHARMM36 resources:

```bash
export MEMBRANEFORGER_CHARMM36_ROOT=/path/to/charmm36.ff
export MEMBRANEFORGER_MEMBRANE_TOPPAR=/path/to/toppar
./run_pipeline.sh doctor --config runs/<run>/config/config.yaml
```

Only a configuration-specific doctor or `scientific-integration --preflight-only`
run can report whether those local resources are ready for Stage 4.
