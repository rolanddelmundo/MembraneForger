# What is in `archive/`

`archive/` is the unmaintained 0.1.0 / 0.2.0 four-stage workflow (setup, CG run, backmapping, AA preparation),
kept unchanged when 1.0.0 introduced the two-input pipeline. Nothing in `membraneforger/`, `slurm/`, `examples/`,
`.github/workflows/` or the tests imports or reads it (checked by searching those trees for `archive`). Its own two
workflows under `archive/.github/` never run, because GitHub only reads the root `.github/`. The current guides and
`CITATION.cff` used to be stored there and now live in `docs/` and the repository root.

| Path | Size | Content | Candidate for archival outside the repository |
|---|---|---|---|
| `resources/vendor/` | 86 MB | tracked copy of the mstool source used by the legacy Stage 3, with its policy files | yes, largest item; the current pipeline uses a pip-installed mstool |
| `resources/forcefields/` | 16 MB | Martini force-field files and their origin notes | yes, but see the licence note below |
| `resources/aa_mdp/`, `cg_mdp/`, `mappings/`, `templates/`, `ligand_params/` | about 0.3 MB | legacy GROMACS parameter files, an older `map.dat` and `GM3.xml` | `mappings/map.dat` differs from `backmap_data/map.dat` (swapped GM3 glucose/galactose assignment fixed in 1.x); keep it separate and never copy it over |
| `membraneforger/`, `stages/`, `scripts/`, `run_pipeline.sh`, `setup.sh` | about 0.6 MB | the legacy code | yes, together with the vendored copy |
| `config/`, `containers/`, `environments/`, `inputs/`, `examples/`, `work/`, `logs/`, `outputs/` | small | legacy configuration, Dockerfile, environment files, placeholders | yes |
| `tests/` | 0.1 MB | legacy tests, not run in CI | yes |
| `docs/` | 0.1 MB | legacy stage guides, audits and inventories | yes |
| `LICENSES/`, `LICENSE*`, `THIRD_PARTY_NOTICES.md`, `SECURITY.md`, `CONTRIBUTING.md`, `.editorconfig`, `.gitattributes` | small | licence texts and notices for the vendored and force-field material | no, not before the notices are reviewed |

Licence note: the licence texts and third-party notices for the vendored mstool and the Martini files exist only here.
Before anything is removed or relocated, decide where those notices will live.

Possible later step, only with the owner's approval: tag the current commit (for example `legacy-0.2.0`), then
delete or move `archive/` to a separate repository. The history keeps every file, so nothing is lost.

Legacy tests, inventories and `archive/README.md` still mention paths such as `CITATION.cff` and
`docs/third_party_inventory.tsv` relative to the archive root; they were not updated because the archive is
unmaintained.
