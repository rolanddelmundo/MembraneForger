External/user-supplied resources are not redistributed in this directory.

Use this directory only for local, ignored copies if your local policy permits, or
configure absolute paths/environment variables in a run config:

- `MEMBRANEFORGER_CHARMM36_ROOT`
- `MEMBRANEFORGER_MEMBRANE_TOPPAR`
- `MEMBRANEFORGER_CGENFF_ROOT`
- `MEMBRANEFORGER_LIGAND_PARAMS_ROOT`

For Stage 4, provide a mutually compatible CHARMM36 GROMACS force-field tree
and CHARMM topology/parameter tree, then run:

```bash
./run_pipeline.sh doctor --config runs/<run>/config/config.yaml
./run_pipeline.sh scientific-integration --config runs/<run>/config/config.yaml --run-dir runs/<run> --preflight-only
```
