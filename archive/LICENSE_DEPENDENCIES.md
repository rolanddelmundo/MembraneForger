# License Structure And External Dependencies

MembraneForger distinguishes the license of project-authored files from the
effective license of the combined bundled pipeline.

## Effective Source Distribution License

MembraneForger-authored files remain MIT-licensed. The public source
distribution also bundles and directly uses GPL-covered scientific components:

- `resources/vendor/mstool/`: GPL-3.0-only
- `scripts/insane_M3_lipids_new.py`: GPL-2.0-or-later

Because MembraneForger imports vendored GPL-3.0-only `mstool` into the Python
process, the combined bundled MembraneForger pipeline is distributed under
GPL-3.0-only terms. The INSANE-derived GPL-2.0-or-later component is compatible
with that combined GPL-3.0-only distribution through its "or later" option.

The root `LICENSE` records this combined-distribution notice. The MIT license
text for MembraneForger-authored files is retained in `LICENSES/MIT.txt`.

## Open Resources

- Martini 3 parameter files retained under `resources/forcefields/martini/`
  are byte-identical to the official Apache-2.0 upstream repository at commit
  `784591ebdc91d762ed4df986c4650546c938f776`.
- `mstool` source is vendored under `resources/vendor/mstool/` at upstream
  commit `2d37f9d3e89279ddd9125cc74da1f5e01153586c`. It retains its
  GPL-3.0-only license. Local packaging/modification details are recorded in
  `resources/vendor/mstool/MODIFICATIONS.md`.
- `scripts/insane_M3_lipids_new.py` is a bundled INSANE-derived script based on
  the script metadata `previous = "20140603.11.TAW"` with local lipid-template
  modifications. It retains GPL-2.0-or-later terms; the full GPL-2.0 text is
  retained in `LICENSES/GPL-2.0-or-later.txt`. Modification details are recorded
  in `docs/insane_modifications.md`.

## EXTERNAL - NOT REDISTRIBUTED

- CHARMM36 GROMACS force-field resources: set
  `MEMBRANEFORGER_CHARMM36_ROOT`.
- Compatible CHARMM membrane topology/parameter resources: set
  `MEMBRANEFORGER_MEMBRANE_TOPPAR`.
- CGenFF/toppar resources: set `MEMBRANEFORGER_CGENFF_ROOT`.
- GLPA or other ligand parameters: set `MEMBRANEFORGER_LIGAND_PARAMS_ROOT`.
- DSSP/mkdssp: install a package that provides `mkdssp` or set `DSSP_BIN`.
- PyRosetta: set `PYROSETTA_PYTHON` to an authorized Python environment.
- Rosetta utilities: set `ROSETTA_BIN`.
- Rosetta database: set `ROSETTA_DATABASE`.
- Rosetta `molfile_to_params.py`: set `MOLFILE_TO_PARAMS`.

Rosetta and PyRosetta are never downloaded, installed, cached, vendored,
uploaded, or redistributed by public CI.

The local development CHARMM resource trees are not vendored because their
exact redistribution rights and paired compatibility have not been established.
See `docs/charmm_resource_audit.md`.
