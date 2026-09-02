# Software Dependencies

Executable resolution order is:

1. documented environment variable;
2. active `PATH` when appropriate;
3. actionable mode-specific failure.

Supported variables:

- `MEMBRANEFORGER_MSTOOL_ROOT`
- `MEMBRANEFORGER_CACHE_DIR`
- `MEMBRANEFORGER_CHARMM36_ROOT`
- `MEMBRANEFORGER_MEMBRANE_TOPPAR`
- `MEMBRANEFORGER_CGENFF_ROOT`
- `MEMBRANEFORGER_LIGAND_PARAMS_ROOT`
- `DSSP_BIN`
- `PYROSETTA_PYTHON`
- `ROSETTA_BIN`
- `ROSETTA_DATABASE`
- `MOLFILE_TO_PARAMS`
- `GMX_BIN`

`mstool` is resolved from the vendored source tree at
`resources/vendor/mstool` by default. `MEMBRANEFORGER_MSTOOL_ROOT` is an
explicit maintainer override. Compiled extensions are built into
`${MEMBRANEFORGER_CACHE_DIR}` or the platform cache, not into the repository.

Stage 2 uses the bundled exact INSANE-derived script at
`scripts/insane_M3_lipids_new.py`; MembraneForger uses that repository script
for this path.

DSSP is required only for DSSP-dependent modes. PyRosetta and Rosetta are
required only for modes that explicitly enable those integrations. Minimal
public dry-runs do not require Rosetta, PyRosetta, CHARMM/CGenFF, GLPA
parameters, or private regression data.

Stage 4 resolves CHARMM36 resources with documented environment variables first,
then normalized config paths. `MEMBRANEFORGER_CHARMM36_ROOT` must point to a
`charmm36.ff` directory containing the required GROMACS force-field files and
ACE/CT3 terminal database entries. `MEMBRANEFORGER_MEMBRANE_TOPPAR` must point
to a compatible CHARMM topology/parameter tree whose protein topology defines
`PRES ACE` and `PRES CT3`.

Use `./run_pipeline.sh doctor --config CONFIG` or
`./run_pipeline.sh scientific-integration --config CONFIG --run-dir RUN_DIR
--preflight-only` to check configuration-specific Stage 4 readiness.

Run:

```bash
python scripts/dependency_resolver.py
bash setup.sh --check
```
