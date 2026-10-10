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

## Common failures

| Symptom | Likely cause | What to do |
|---|---|---|
| `ERROR` in stage `inputs` | water, bulk ions, an insertion code, duplicate atoms, or a residue with no topology | Remove the offending records. A ligand needs its `.itp` in `forcefield/toppar/`. |
| `Multiple protein chains were detected` (stage `orient`) | several protein chains and no anchor chain named | Add `--orient-chain` with the chain the error suggests. |
| The build asks for `--nterm-side` | PPM is orienting the protein and no OPM entry gives the side | Add `--nterm-side in` or `out` for residue 1 of that chain as it is in your file. |
| Protein upside down in `int/oriented.pdb` | `--nterm-side` wrong for this file, often because a signal peptide is still attached | Rerun with the other side. The extracellular domain must end up at +z. |
| PPM fails or is not found | wrong path, or `res.lib` missing next to `immers` | Check `--ppm-exe` or `MEMBRANEFORGER_PPM`. See [Advanced options](advanced.md#ppm-30). Already-oriented inputs can use `--orientation none`. |
| Empty pocket stops the build after embedding | the protein is much smaller than the frame's receptor hole | Use `--embed-site free` with a larger `--xy-buffer`, or your own membrane frame. |
| `ERROR` in stage `slice_check` | the cut changed the lipid packing beyond `--apl-tolerance` | Read `membrane_validation.md`. Rerun with a larger `--xy-buffer` (2.0 or more). |
| `ERROR` in stage `mstool` | the interpreter cannot import mstool | Install mstool 0.3.9 or 0.3.10, or point `--mstool-python` at a Python that has it. |
| Backmapping fails after all seeds | ring threading or a clash that every seed reproduces | Read the `backmap` lines of the log and `int/work/mstool/`. |
| `ERROR` in stage `validate` after minimization | the largest force stayed above the target, or a clash remained | The log names the atom and its neighbours. `em_clash_trace.json` traces a lipid culprit back to the first stage it clashes. |
| GROMACS error | `gmx` not found, or a topology problem | Pass `--gmx`. Read the GROMACS message in the log. |
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
gmx grompp -f int/em.mdp -c em.gro -p topol.top -n index_ini.ndx -o test.tpr -maxwarn 0
```

The [tutorial](membraneforger_tutorial.md) explains every check and threshold, and its
[limitations section](membraneforger_tutorial.md#7-limitations-and-assumptions) lists what the build does not do.
