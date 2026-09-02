# MembraneForger

MembraneForger prepares protein-in-membrane structures through four checked
stages. The input is a PDB file containing atomic coordinates for the protein or
protein complex you want to place in a membrane.

Stage 1 inspects and prepares the molecular input. It records what chains,
residues, water molecules, alternate locations, and unsupported nonprotein
components are present before any scientific assumptions are made.

Stage 2 builds a simplified membrane representation. This uses Martini
coarse graining: several atoms are grouped into larger particles so large
molecular systems are easier to construct and manipulate.

Stage 3 converts the simplified model back to individual atoms. This is called
backmapping. The canonical Stage 3 output is:

```text
outputs/<run_id>/stage3/final_all_atom.pdb
```

Stage 4 prepares the atom-by-atom system using CHARMM36-compatible molecular
parameters. CHARMM36 is a molecular force field: a set of parameters that tells
simulation software how atoms interact. MembraneForger does not bundle the
development CHARMM resource set; Stage 4 needs compatible external/user-supplied
CHARMM36 resources whose exact redistribution rights and compatibility must be
established. Stage 4 prepares and checks an energy-minimized all-atom system; it
does not claim production molecular dynamics.

Green software tests mean the workflow code and safeguards passed. They do not
prove that every biological system, protonation state, ligand, membrane
composition, or force-field choice is scientifically valid.

## Beginner Workflow

```bash
git clone <repository-url> MembraneForger
cd MembraneForger

./setup.sh --verify
./run_pipeline.sh doctor
./run_pipeline.sh dry-run --config config/workflow.yaml
```

That checks the local software environment, confirms that the public workflow
can be planned, and does not require external CHARMM36 resources.

To start a real run, stage a PDB and review the generated config:

```bash
./run_pipeline.sh init \
  --pdb receptor.pdb \
  --output-dir runs/receptor

./run_pipeline.sh doctor --config runs/receptor/config/config.yaml
./run_pipeline.sh validate --config runs/receptor/config/config.yaml
./run_pipeline.sh dry-run --config runs/receptor/config/config.yaml
```

Then fill in the membrane composition, box/orientation assumptions, and any
external Stage 4 resource paths. When the configuration-specific doctor passes,
launch the supported stage-gated workflow:

```bash
./run_pipeline.sh run --config runs/receptor/config/config.yaml --through stage4
```

If a stage requires scientific review, the run stops with `REVIEW_REQUIRED` and
prints the exact `resume` command.

## Public Commands

```bash
./setup.sh --sync
./setup.sh --verify

./run_pipeline.sh doctor [--config CONFIG]
./run_pipeline.sh init --pdb receptor.pdb --output-dir runs/receptor
./run_pipeline.sh validate --config runs/receptor/config/config.yaml
./run_pipeline.sh dry-run --config runs/receptor/config/config.yaml
./run_pipeline.sh run --config runs/receptor/config/config.yaml --through stage4
./run_pipeline.sh resume --run-dir runs/receptor --accept-review
./run_pipeline.sh stage stage3 --run-dir runs/receptor
./run_pipeline.sh status --run-dir runs/receptor
./run_pipeline.sh clean --run-dir runs/receptor --confirm
```

Legacy stage wrappers under `stages/` remain available for advanced users, but
new runs should use `run_pipeline.sh`.

## Inputs

Required input is a readable PDB with `ATOM` or `HETATM` records. `init` copies
the original file to `runs/<name>/input/original.pdb` and records its SHA256.
The original PDB is never modified.

Before computation, MembraneForger reports chain IDs, residue counts, alternate
locations, waters, heteroatoms, and missing backbone atoms. Scientifically
ambiguous choices such as membrane orientation, termini, nonstandard residues,
ligands, glycans, metals, and missing heavy atoms must be resolved by the user.

## Configuration

`runs/<name>/config/config.yaml` is organized around user concepts:

```text
input
membrane
coarse_grained
simulation
backmapping
all_atom
compute
review
advanced
```

The generated file also contains legacy `stage1` through `stage4` sections for
compatibility. Public sections are normalized once before execution and are
authoritative; conflicting explicit legacy duplicates fail validation.

## Stages

- `stage1`: prepare the protein for Martini CG with Martinize2.
- `stage2`: build and equilibrate the CG membrane system. De novo INSANE-driven
  setup is the normal path; scaffold replacement remains an advanced mode.
- `stage3`: backmap the accepted CG structure with mstool and validate mappings.
- `stage4`: prepare CHARMM36 all-atom topology, solvent, ions, index, and
  minimization-ready files.

For CHARMM36 terminal capping, MembraneForger uses CHARMM's native terminal
patch system. An acetylated N-terminus uses the `ACE` patch, while an
N-methylamide C-terminus uses the `CT3` patch. These modify the terminal amino
acid in the CHARMM topology rather than being treated as ordinary extra protein
residues. Do not add standalone `NME` residues for the CHARMM36 workflow.

Stage states are written to `status.json`:

```text
NOT_STARTED READY RUNNING COMPLETE FAILED BLOCKED REVIEW_REQUIRED
```

## Outputs

Run output is predictable:

```text
runs/my_system/
  input/
  config/
  logs/
  provenance/
  stage1_cg_setup/
  stage2_cg_simulation/
  stage3_backmapping/
  stage4_all_atom/
  final/
    cg/
    all_atom/
  status.json
```

Stable final aliases are created only after the corresponding stage validates:

```text
final/cg/backmap_input.gro
final/cg/structure.gro
final/cg/topology.top
outputs/<run_id>/stage3/final_all_atom.pdb
final/all_atom/backmapped.pdb
final/all_atom/backmapped.dms
final/all_atom/minimized_all_atom.gro
final/all_atom/topology.top
```

## Dependencies

`./run_pipeline.sh doctor --config CONFIG` distinguishes globally required,
configuration-required, optional, disabled, licensed, and missing-but-not-needed
dependencies.

Core open dependencies:

```text
Python, PyYAML, NumPy, OpenMM, Vermouth/martinize2, GROMACS, bundled INSANE
script, Martini resources, vendored mstool source
```

Optional or external/user-supplied dependencies:

```text
DSSP/mkdssp, CHARMM36, membrane toppar resources, CGenFF, ligand parameters,
PyRosetta, Rosetta, Rosetta database, molfile_to_params.py
```

Current generic Stage 4 terminates at an energy-minimized all-atom system.
Restrained all-atom equilibration and production are not claimed as completed
outputs by the public workflow.

External/user-supplied resources are not bundled. Configure them with
environment variables or absolute paths in the run config:

```bash
export MEMBRANEFORGER_CHARMM36_ROOT=/path/to/charmm36.ff
export MEMBRANEFORGER_MEMBRANE_TOPPAR=/path/to/membrane/toppar
export MEMBRANEFORGER_CGENFF_ROOT=/path/to/cgenff
export MEMBRANEFORGER_LIGAND_PARAMS_ROOT=/path/to/ligand_params
```

## Containers

Containers provide the open-source runtime only. External/user-supplied
resources must be mounted and configured at run time.

```bash
docker build -f containers/Dockerfile -t membraneforger:portable .
docker run --rm -it \
  -v "$PWD:/workspace/MembraneForger" \
  -v "$PWD/runs:/workspace/MembraneForger/runs" \
  membraneforger:portable ./run_pipeline.sh doctor

apptainer build MembraneForger.sif containers/Apptainer.def
apptainer exec --bind "$PWD:/workspace/MembraneForger" \
  MembraneForger.sif ./run_pipeline.sh doctor
```

## Reproducibility

Each run records:

```text
MembraneForger commit, command line, input checksum, effective config,
dependency paths/versions, resource checksums, stage timestamps, warnings,
stage states, and final output checksums
```

Resource provenance is tracked in `resources/RESOURCE_MANIFEST.tsv`,
`config/resources.lock.yaml`, `docs/third_party_inventory.tsv`, and
`docs/runtime_provenance.tsv`.

Maintainers with local CHARMM resources can record Stage 4 resource preflight
metadata and tree hashes before a real scientific run:

```bash
./run_pipeline.sh scientific-integration \
  --config runs/receptor/config/config.yaml \
  --run-dir runs/receptor \
  --preflight-only
```

This command does not copy restricted resources into the repository or into
public artifacts.

## Troubleshooting

See `docs/troubleshooting.md` for common failures, including missing GROMACS,
missing CHARMM36 or membrane parameters, unsupported chemistry, unknown or
conflicting config keys, ambiguous Stage 4 atom matching, unsupported
`run_cg_production: true`, and cleanup refusals.

## Tests

```bash
bash -n run_pipeline.sh setup.sh stages/stage1_setup_cg/run.sh \
  stages/stage2_run_cg/run.sh stages/stage4_prepare_aa/run.sh
python -m compileall -q membraneforger scripts stages tests
python -m pytest -q
```

Real Stage 1-4 and full CG-to-AA validation must be reported as `BLOCKED`, not
`PASS`, when required external scientific resources are unavailable.

## Citation And Licensing

If you use MembraneForger, cite this repository and the upstream tools listed
in `THIRD_PARTY_NOTICES.md`.

MembraneForger-authored source files remain MIT-licensed. The distributed
bundled pipeline also includes GPL-covered scientific components:
`resources/vendor/mstool` is GPL-3.0-only, and
`scripts/insane_M3_lipids_new.py` is GPL-2.0-or-later. Because the runnable
bundled pipeline imports and uses these components, the combined source
distribution is provided under GPL-3.0-only terms while preserving the MIT grant
for MembraneForger-authored files used separately. See `LICENSE`,
`LICENSE_DEPENDENCIES.md`, and `THIRD_PARTY_NOTICES.md`.
