Runtime force-field policy:

- Redistributable Martini resources are retained under `resources/forcefields/martini/`.
- CHARMM36 and membrane toppar payloads are not redistributed here.
- Configure external all-atom resources with:
  - `MEMBRANEFORGER_CHARMM36_ROOT=/path/to/charmm36.ff`
  - `MEMBRANEFORGER_MEMBRANE_TOPPAR=/path/to/membrane/toppar`

Run `./run_pipeline.sh doctor --config runs/<run>/config/config.yaml` to verify them.
