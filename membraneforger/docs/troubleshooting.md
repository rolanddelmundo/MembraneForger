# Troubleshooting

When a build stops, `membranebuilder.log` ends with `FAIL`, and the first `ERROR:` line names the stage that failed.
Read the `WARN:` lines just before it too. Nothing is deleted on failure: the intermediates are in `<out>/int/` and
the stage is recorded in `run_manifest.json` (`error.stage`).

## Log-line prefixes

| Prefix | Meaning |
|---|---|
| `CHECK:` | a stage starts and its checks run |
| `PASS:` | a check succeeded |
| `WARN:` | something to review; the build continues |
| `ERROR:` | a check failed; the build stopped |
| `INFO:` | a measurement or a count |

Running `membraneforger --check` with the same options before a build catches most of the failures below in seconds.

## Common failures

| Symptom | Likely cause | What to do |
|---|---|---|
| `cannot import numpy` (or scipy, networkx, pandas, openmm) before anything runs | the Python requirements are not installed in this interpreter | `pip install .` from the repository, or the `pip install` line the message prints. |
| `ERROR` in stage `inputs` | water, bulk ions, an insertion code, duplicate atoms, or a residue with no topology | Remove the offending records ([Preparing the input](preparing_inputs.md)). |
| `no topology in ... for residue(s) ['XYZ']` (stage `inputs`) | a ligand, ion or modified residue without an `.itp` | Add `toppar/XYZ.itp` to the force-field directory, see [Ligands](ligands.md). |
| `heavy-atom graph is not isomorphic to its ITP` or `source hydrogens have no topology counterpart` (stage `topology`) | the ligand in the input and its `.itp` differ in connectivity or protonation | Make the input and the topology describe the same molecule ([Ligands](ligands.md)). |
| `unrecognized C terminus` (stage `topology`) | the last residue of a chain has no `OXT`/`OT1`/`OT2` and no `CT2`/`CT3` cap | Add the terminal oxygen or cap the chain ([Preparing the input](preparing_inputs.md#chains-and-termini)). |
| OPM download fails (stage `orient`) | no network on the node | Pass `--opm-file` with the downloaded entry, or `--orientation ppm` or `none`. |
| `Multiple protein chains were detected` (stage `orient`) | several protein chains and no anchor chain named | Add `--orient-chain` with the chain the error suggests. |
| The build asks for `--nterm-side` | PPM is orienting the protein and no OPM entry gives the side | Add `--nterm-side in` or `out` for residue 1 of that chain as it is in your file. |
| Protein upside down in `int/oriented.pdb` | `--nterm-side` wrong for this file, often because a signal peptide is still attached | Rerun with the other side. The extracellular domain must end up at +z. |
| PPM fails or is not found | wrong path, or `res.lib` missing next to `immers` | Check `--ppm-exe` or `MEMBRANEFORGER_PPM`. See [Advanced options](advanced.md#ppm-30). Already-oriented inputs can use `--orientation none`. |
| Empty pocket stops the build after embedding | the protein is much smaller than the frame's receptor hole | Use `--embed-site free` with a larger `--xy-buffer`, or your own membrane frame. |
| `ERROR` in stage `slice_check` | the cut changed the lipid packing beyond `--apl-tolerance` | Read `membrane_validation.md`. Rerun with a larger `--xy-buffer` (2.0 or more). |
| `ERROR` in stage `mstool` | the interpreter cannot import mstool | Install mstool 0.3.9 or 0.3.10, or point `--mstool-python` at a Python that has it. |
| Backmapping fails after all seeds | ring threading or a stereochemistry error that every seed reproduces | Read the `backmap` lines of the log and `int/work/mstool/`. |
| `ERROR` in stage `validate` after minimization | the largest force stayed above the target, or a clash remained | The log names the atom and its neighbours. `em_clash_trace.json` traces a lipid culprit back to the first stage it clashes. Ligand heavy atoms are frozen in EM, so a clash with a ligand is not relieved by moving the ligand: fix the ligand pose or the clashing residue in the input. |
| `ERROR` in stage `validate`: `the frozen ligand moved during EM` | a ligand atom in `LIG` shifted by more than 0.002 nm, so the freeze did not take effect | Compare `int/solv_ions.gro` with `int/em.unverified.gro` and check that `int/emres.mdp` and `int/em.mdp` contain `freezegrps = LIG`. |
| GROMACS error | `gmx` not found, or a topology problem | Pass `--gmx`. Read the GROMACS message in the log. |
| `--runs` or `--cg` rejects a frame name | the name is not a bundled frame, or `--runs` is larger than the series (9 frames each) | `--list-membranes` prints the valid names; `--runs N` needs `--cg 1` or `2`. |
| The build is `PASS` but the log has `WARN:` lines, or the report says PASS WITH WARNINGS | a measurement is near a limit or poorly sampled (rare lipid species, Fmax, a contact) | The build is usable; read each warning and the matching section of the [tutorial](validation/membrane.md#5-validation) before production. Equilibrate and check the trajectory with `python -m membraneforger.equilibration`. |
| `input ... would be overwritten by the build` | the input file lives in the output directory under a generated name | Choose another `--out`. |

## Where to look for each stage

| Stage | Look at |
|---|---|
| `inputs` | the input file named in the message |
| `orient` | `orientation_report.json`, the chains of your input, and the PPM transcript in the log |
| `slice`, `slice_check`, `membrane_check` | `membrane_validation.md`, `membrane_validation.json`, the slice report in `run_manifest.json` |
| `mstool`, `backmap` | the mstool transcript at the end of the log and `int/work/mstool/` |
| `topology` | the `pdb2gmx` transcript in the log |
| `solvate`, `ions` | `int/boxed.gro` or `int/solv.gro`, and the GROMACS transcript in the log |
| `minimize`, `validate` | `int/em.log`, `int/em.unverified.gro`, `em_clash_trace.json` |
| `audit` | `audit.json` |

## Checking a finished build

Before equilibration, confirm that `audit.json` reports `"result": "PASS"` and that GROMACS accepts the files
without warnings. From the output directory:

```bash
# with a ligand, int/em.mdp freezes the LIG group of index_ini.ndx
gmx grompp -f int/em.mdp -c em.gro -p topol.top -n index_ini.ndx -o test.tpr -maxwarn 0
```

The [tutorial](validation/index.md) explains every check and threshold, and its
[limitations section](validation/limitations.md#7-limitations-and-assumptions) lists what the build does not do.
The next step is the restrained equilibration in [Equilibration](equilibration.md).
