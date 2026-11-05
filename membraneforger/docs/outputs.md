# Output files

A successful build writes these files to the `--out` folder.

| File | What it is |
|---|---|
| `em.gro` | The finished system: protein, membrane, water and ions, energy-minimized. It is written only if every check passed. Start equilibration from this file. |
| `topol.top`, `toppar/` | The force-field description of every molecule (CHARMM36). GROMACS needs these together with `em.gro`. |
| `index_ini.ndx` | Atom groups (protein and ligands, membrane, water and ions) used by the equilibration files. |
| `membrane_validation.md` | A readable report on the membrane: packing, composition and thickness, with figures and an overall verdict (PASS, PASS WITH WARNINGS or FAIL). |
| `audit.json` | Independent checks on the final system: topology, net charge, energy, ring threading and lipid geometry. |
| `run_manifest.json` | The exact settings, software versions and checksums. Attach it when you report a problem. |
| `membranebuilder.log` | A line-by-line log. Its last line is `PASS` or `FAIL`. Read every `WARN:` line. |
| `int/` | Intermediate files from every stage, for inspection or restarts. A failed build leaves `int/em.unverified.gro` here. |

What each check in `membrane_validation.md` and `audit.json` measures is explained under [Validation](validation/index.md).
The full list of files, including the validation figures, is in [Outputs in detail](advanced.md#outputs-in-detail).
