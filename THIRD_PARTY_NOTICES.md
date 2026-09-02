# Third-Party Notices

Detailed file-level records are in `docs/third_party_inventory.tsv` and the
concise runtime inventory is in `docs/runtime_provenance.tsv`.

## Effective Combined Pipeline License

MembraneForger-authored files remain MIT-licensed. The bundled runnable
pipeline imports GPL-3.0-only mstool and includes GPL-2.0-or-later
INSANE-derived code. The combined bundled MembraneForger pipeline is therefore
distributed under GPL-3.0-only terms. Full texts are retained in `LICENSES/`.

## Martini Force Fields

- Copyright: upstream Martini contributors as recorded in the official source
- License: Apache-2.0
- Version/commit: `784591ebdc91d762ed4df986c4650546c938f776`
- Source: `https://github.com/marrink-lab/martini-forcefields.git`
- Required notice: retain Apache-2.0 license text in `LICENSES/Apache-2.0.txt`
- Required citations: Martini 3 general citation and molecule-specific
  citations from upstream documentation
- Local modifications: none
- Status: vendored public resource

## mstool

- License: GPL-3.0-only
- Version/commit: `2d37f9d3e89279ddd9125cc74da1f5e01153586c`
- Source: `https://github.com/ksy141/mstool.git`
- Local modifications: no upstream source-code changes; MembraneForger adds
  `MODIFICATIONS.md` and copies the upstream license into the package
  directory. Local generated C extension outputs are ignored, and native
  extensions build from vendored `.pyx` source into an external cache.
- Status: vendored source under `resources/vendor/mstool/`

See `resources/vendor/mstool/MODIFICATIONS.md`.

## INSANE-Derived Membrane Builder

- License: GPL-2.0-or-later
- Script metadata: `previous = "20140603.11.TAW"`
- Source lineage: INSert membrAN by Tsjerk A. Wassenaar, with local
  lipid-template modifications marked in the script header
- Local path: `scripts/insane_M3_lipids_new.py`
- Modification notice: script header and `docs/insane_modifications.md`
- Status: bundled exact script used by Stage 2 de novo membrane construction

## Restricted External Dependencies

CHARMM36, CGenFF, Rosetta, PyRosetta, Rosetta database files,
`molfile_to_params.py`, and user ligand parameters are external user-supplied
dependencies and are not redistributed.
