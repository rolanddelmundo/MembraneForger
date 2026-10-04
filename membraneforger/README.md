# MembraneForger

All-atom protein/ligand complex + Martini 3 coarse-grained membrane system -> validated, energy-minimized
all-atom CHARMM36/GROMACS system.

```bash
python -m membraneforger --all-atom prot-lig.pdb --coarse-grain system.gro --orient-chain R --out output_directory
```

| File | Role |
|---|---|
| `cli.py` | command line, resource discovery |
| `pipeline.py` | stage order, timing, fail-closed publication of `em.gro` |
| `config.py` | CHARMM/GROMACS constants and the tunable `Settings` |
| `martini.py` | the Martini 3 classification layer and residue-name aliases |
| `validation.py` | input validation, input/stale-output protection |
| `orientation.py` | membrane orientation from the anchor chain(s): OPM exact reference or local PPM 3.0, one rigid transform for the whole complex, validation, `orientation_report.json` |
| `alignment.py` | sequence correspondence and rigid placement of the all-atom complex (a free fit, or, after orientation, a rotation about z plus translation that keeps the orientation) |
| `backmapping.py`, `mstool_worker.py` | membrane backmapping; the worker is the only file that imports mstool |
| `topology.py` | structure repair and CHARMM36 topology construction |
| `solvation.py` | box (`BOX = auto` around the bilayer midplane, or an opt-in user box), water, ions, index groups |
| `minimization.py` | EM and its validation |
| `audit.py` | independent audit (own parsers, `grompp -maxwarn 0`, `gmx check`, `gmx energy`, geometry) |
| `reporting.py` | `run_manifest.json` |
| `visualization.py` | VMD + Tachyon stage renders (runs in a python that can import `vmd`) |
| `runtools.py`, `structio.py` | logging/subprocess wrapper, file readers and writers |
| `example.py`, `cpu.submit` | use from Python, and a Slurm array template |
| `tests/` | `python -m unittest discover -s membraneforger/tests -t .` from the repository root |

Requirements: Python >= 3.10 with numpy, scipy, networkx; GROMACS; an interpreter with mstool 0.3.9 or 0.3.10
(`--mstool-python`). Installation, example and outputs are described in the top-level `README.md`.
