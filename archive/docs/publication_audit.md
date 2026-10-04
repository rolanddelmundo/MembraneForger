# Publication Audit

Audit date: 2026-07-15.

The public snapshot contains project-authored workflow code, documentation,
minimal synthetic examples, verified Apache-2.0 Martini files, vendored mstool
source, and the bundled INSANE-derived membrane builder. External resources
with unestablished redistribution terms are excluded.

| path or pattern | disposition | public status |
|---|---|---|
| `resources/vendor/mstool/` | retained exact GPL-3.0-only source with checksum lock | VENDORED-MODIFIED |
| `resources/forcefields/martini/` | retained byte-identical from official Apache-2.0 upstream | KEEP_VERIFIED |
| `resources/forcefields/charmm36.ff/` | removed; user supplies external directory with established local rights | EXTERNAL - NOT REDISTRIBUTED |
| `resources/forcefields/toppar/` | removed; user supplies CHARMM/CGenFF/toppar files with appropriate local rights | EXTERNAL - NOT REDISTRIBUTED |
| `resources/ligand_params/GLPA/` parameter assets | removed; user supplies reproducible external parameters | EXTERNAL - NOT REDISTRIBUTED |
| `examples/legacy_kor/` | removed from public snapshot | EXTERNAL - NOT REDISTRIBUTED |
| `examples/test_membrane/` | removed from public snapshot | EXTERNAL - NOT REDISTRIBUTED |
| Rosetta, PyRosetta, Rosetta database, `molfile_to_params.py` | user-supplied installations only | EXTERNAL - NOT REDISTRIBUTED |
| DSSP/mkdssp | external executable, installable in public CI | EXTERNAL - NOT REDISTRIBUTED |
| `examples/minimal/` | project-authored synthetic public dry-run fixture | KEEP_VERIFIED |

All public-snapshot entries are resolved by removal, verified retention, or
user-supplied external configuration.
