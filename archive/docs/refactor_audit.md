# Portability Refactor Audit

This audit records the repository state addressed by the public workflow
refactor.

## Original Execution Architecture

- `run_pipeline.sh` delegated directly to `scripts/membraneforger_workflow.py`.
- Public staging/dry-run behavior existed, but there were no first-class
  `doctor`, `init`, `run`, `resume`, or `status` commands.
- Production stage execution was implemented in
  `scripts/test_membrane_workflow.py`, despite the test-oriented filename.
- Stage 2 real execution was scaffold-replacement oriented; de novo INSANE
  membrane construction was not routed from the public workflow.
- Stage 4 outputs under `outputs/<run>/stage4` were not standalone because
  topology includes referenced force-field files outside the published output.

## Dependency And Resource Findings

- Active Python packages on the audit machine did not match the pinned
  `environments/environment.yml` versions.
- Repository-local mstool is vendored under `resources/vendor/mstool`; the
  pinned lock and tree checksum are recorded in
  `resources/vendor/mstool.lock.yaml`.
- GROMACS, martinize2/Vermouth, OpenMM, and Martini resources were detected.
- An external INSANE executable failed because `pkg_resources` was missing, so
  the bundled INSANE-compatible script is the default checked path.
- CHARMM36, membrane toppar, GLPA parameters, and private regression outputs
  were present locally but are not redistributed.

## Cleanup Actions

The following were moved to `.local/quarantine/20260901T194240Z/`:

- Large private run fixtures under `inputs/`, `outputs/`, `work/`, and `logs/`.
- Restricted/user-supplied CHARMM36, membrane toppar, and GLPA parameter files.
- Local template payloads whose redistribution status was not required for the
  public workflow.
- Editor/cache/generated artifacts such as `.DS_Store`, `__pycache__`,
  `.pytest_cache`, backup files, and stale generated MDP artifacts.

No raw scientific input was edited in place.

## New Public Workflow

- `./run_pipeline.sh doctor`
- `./run_pipeline.sh init --pdb receptor.pdb --output-dir runs/receptor`
- `./run_pipeline.sh validate --config runs/receptor/config/config.yaml`
- `./run_pipeline.sh dry-run --config runs/receptor/config/config.yaml`
- `./run_pipeline.sh run --config runs/receptor/config/config.yaml`
- `./run_pipeline.sh resume --run-dir runs/receptor`
- `./run_pipeline.sh status --run-dir runs/receptor`

Stage states are recorded in `status.json` and final aliases are created under
`final/cg` and `final/all_atom` when validated stage outputs exist.

## Remaining Scientific Blockers

- Full end-to-end execution still requires representative scientific inputs and
  per-stage review.
- Stage 4 real execution is blocked until compatible external CHARMM36 and
  membrane toppar resources are configured.
- A scientifically meaningful tiny end-to-end fixture is still needed; the
  existing minimal PDB is only suitable for public CLI smoke tests.
