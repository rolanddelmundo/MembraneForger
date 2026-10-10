# Pre-equilibrated coarse-grained cell membranes (Martini 3)

Eighteen frames of a single GPCR in an asymmetric, ten-species cell-membrane model, each taken at 30 µs of a
Martini 3 simulation. They are the bundled membranes of MembraneForger: `--cg 1` picks one of the GPR139 frames at random, `--cg 2`
one of the kappa opioid receptor frames, and `--cg NAME` (for example `--cg KOR5`) a particular frame.
`python -m membraneforger --list-membranes` lists them, and `--runs N` builds one complex in the first N frames of a
series, giving N replicate systems that differ in their lipid arrangement.

| Files | Receptor | Replicates |
|---|---|---|
| `KOR1_cg_cellmem.gro` … `KOR9_cg_cellmem.gro` | kappa opioid receptor (PDB 8F7W) | 9 |
| `GPR1_cg_cellmem.gro` … `GPR9_cg_cellmem.gro` | GPR139 (PDB 7F6G) | 9 |

Each frame holds about 60,000 beads in a cell of about 18.3 × 18.3 × 19.2 nm: the receptor (about 290 residues),
about 1,300 lipids, water and ions. The membranes were built with INSANE at CHOL:POPC:DOPC:POPE:DOPE:PSM:DPG3 =
25:20:20:5:5:15:10 (upper leaflet) and CHOL:POPC:DOPC:POPE:DOPE:POPS:DOPS:SAP6 = 25:5:5:20:20:8:7:10 (lower
leaflet). Composition at 30 µs, mole percent of each leaflet averaged over the 18 frames
(`python examples/leaflet_composition.py examples/preeq_cg_cellmem/KOR?_cg_cellmem.gro examples/preeq_cg_cellmem/GPR?_cg_cellmem.gro`):

| Leaflet | CHOL | POPC | DOPC | POPE | DOPE | PSM | DPG3 | POPS | DOPS | SAP6 |
|---|---|---|---|---|---|---|---|---|---|---|
| Upper | 28 | 19 | 19 | 5 | 5 | 14 | 10 | 0 | 0 | 0 |
| Lower | 22 | 5 | 5 | 21 | 21 | 0 | 0 | 8 | 7 | 10 |

`6WHC_MTZP_cg_cellmem.gro` and `6WHC_MTZP_prot-lig.pdb` are a complete input pair (glucagon receptor–tirzepatide–Gs
in the same membrane model, 11.2 × 11.2 nm) used by the tests:

```bash
python -m membraneforger --aa examples/preeq_cg_cellmem/6WHC_MTZP_prot-lig.pdb \
    --cg examples/preeq_cg_cellmem/6WHC_MTZP_cg_cellmem.gro --orientation none --out example_out
```

History of the 18 KOR/GPR frames: as first committed, their coordinates were rotated about z relative to the box line
(by 3 to 83 degrees, a different angle per frame, from a rotational fit saved without its box), so they were not
periodic in their box and rendered as a "diamond" inside it. They were rotated back about the box centre in place,
coordinates and velocities, with `python examples/align_frames.py examples/preeq_cg_cellmem/*_cg_cellmem.gro`, which
reports every file as periodic (no bead pairs of different lipids within 0.30 nm once wrapped). MembraneForger runs
the same check on every frame it reads (`align_frame_to_box`) and corrects a rotated frame with a warning that names
the angle. `6WHC_MTZP_cg_cellmem.gro` was never affected.
