# Pre-equilibrated coarse-grained cell membranes (Martini 3)

Eighteen frames of a single GPCR in an asymmetric, ten-species cell-membrane model, each taken at 30 µs of a
Martini 3 simulation. They are the bundled membranes of MembraneForger and can also be given explicitly with
`--coarse-grain`.

| Files | Receptor | Replicates |
|---|---|---|
| `KOR1_cg_cellmem.gro` … `KOR9_cg_cellmem.gro` | kappa opioid receptor (PDB 8F7W) | 9 |
| `GPR1_cg_cellmem.gro` … `GPR9_cg_cellmem.gro` | GPR139 (PDB 7F6G) | 9 |

Each frame holds about 60,000 beads in a cell of about 18.3 × 18.3 × 19.2 nm: the receptor (about 290 residues),
about 1,300 lipids, water and ions. The membranes were built with INSANE at CHOL:POPC:DOPC:POPE:DOPE:PSM:GM3 =
25:20:20:5:5:15:10 (upper leaflet) and CHOL:POPC:DOPC:POPE:DOPE:POPS:DOPS:SAP6 = 25:5:5:20:20:8:7:10 (lower
leaflet). Composition at 30 µs, mole percent of each leaflet averaged over the 18 frames
(`python examples/leaflet_composition.py examples/preeq_cg_cellmem/KOR?_cg_cellmem.gro examples/preeq_cg_cellmem/GPR?_cg_cellmem.gro`):

| Leaflet | CHOL | POPC | DOPC | POPE | DOPE | PSM | GM3 | POPS | DOPS | SAP6 |
|---|---|---|---|---|---|---|---|---|---|---|
| Upper | 28 | 19 | 19 | 5 | 5 | 14 | 10 | 0 | 0 | 0 |
| Lower | 22 | 5 | 5 | 21 | 21 | 0 | 0 | 8 | 7 | 10 |

`6WHC_MTZP_cg_cellmem.gro` and `6WHC_MTZP_prot-lig.pdb` are a complete input pair (glucagon receptor–tirzepatide–Gs
in the same membrane model, 11.2 × 11.2 nm) used by the tests:

```bash
python -m membraneforger --all-atom examples/preeq_cg_cellmem/6WHC_MTZP_prot-lig.pdb \
    --coarse-grain examples/preeq_cg_cellmem/6WHC_MTZP_cg_cellmem.gro --out example_out
```
