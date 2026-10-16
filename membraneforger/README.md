# MembraneForger

All-atom protein/ligand complex (+ a bundled or own Martini 3 coarse-grained membrane) -> validated,
energy-minimized all-atom CHARMM36/GROMACS system.

```bash
python -m membraneforger --aa complex.pdb --orient-chain R --out output_directory     # bundled membrane
python -m membraneforger --aa prot-lig.pdb --cg system.gro --orient-chain R --out output_directory   # own frame
```

| File | Role |
|---|---|
| `cli.py` | command line, resource discovery |
| `pipeline.py` | stage order, timing, fail-closed publication of `em.gro` |
| `config.py` | CHARMM/GROMACS constants and the tunable `Settings` |
| `martini.py` | the Martini 3 classification layer and residue-name aliases |
| `validation.py` | input validation, input/stale-output protection |
| `orientation.py` | membrane orientation from the anchor chain(s): OPM exact reference or local PPM 3.0, one rigid transform for the whole complex, validation, `orientation_report.json` |
| `alignment.py` | sequence correspondence and rigid placement of the all-atom complex on its own CG frame (a free fit, or, after orientation, a rotation about z plus a translation that keeps the orientation) |
| `embedding.py` | placing the complex into a bundled membrane, `--dellipid`/`--addlipid` |
| `lipids.py` | the one table of lipid anchors (Martini bead / atomistic atom per species) used by slicing and every analysis |
| `slicing.py` | `BOX = auto`: PBC-aware, anchor-based, whole-lipid crop of the coarse-grained membrane to the complex plus a buffer (offset search, seam relaxation) before backmapping; also applies a `--box` x and y |
| `packing.py` | periodic Voronoi area per lipid, global and species APL, composition, equal-window distributions |
| `rdf.py` | lateral headgroup RDFs per leaflet and their comparison |
| `structure_metrics.py` | bilayer thickness, protein tilt and depth, Z-density profiles, core hydration, tail interdigitation |
| `trajectory.py` | convergence, S_CD order parameters, MSD and diffusion, K_A, leaflet tension (gates 3 and 4) |
| `qc.py` | the metric record and the PASS / WARNING / FAIL / INSUFFICIENT SAMPLING classification |
| `membrane_report.py` | the staged validation of a build: slice gate, all-atom stages, `membrane_validation.md/.json`, figures |
| `equilibration.py` | `python -m membraneforger.equilibration`: gates 3 and 4 on a finished trajectory |
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

Requirements: Python >= 3.10 with numpy, scipy, networkx (matplotlib optional, for the validation figures); GROMACS;
an interpreter with mstool 0.3.9 or 0.3.10 (`--mstool-python`). Installation, example and outputs are described in the top-level `README.md`.
