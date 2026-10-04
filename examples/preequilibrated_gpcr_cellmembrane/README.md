# Pre-equilibrated GPCR cell-membrane frames (Martini 3)

Eighteen coarse-grained frames of a single GPCR in an asymmetric, ten-species cell-membrane model, each taken at
30 µs of a Martini 3 simulation. They are ready to use as the `--coarse-grain` input of MembraneForger.

| Files | Receptor | Replicates |
|---|---|---|
| `8F7W_KOR_run1_30us.gro` … `run9` | kappa opioid receptor (PDB 8F7W) | 9 |
| `GPR139_7F6G_run1_30us.gro` … `run9` | GPR139 (PDB 7F6G) | 9 |

Each frame holds about 60,000 beads in a cell of about 18.2 × 18.2 × 19.2 nm: the receptor (about 290 residues),
about 1,300 lipids, water and ions. The membranes were built with INSANE with these lipid ratios:

| Leaflet | Composition |
|---|---|
| Upper | CHOL 70 : POPC 56 : DOPC 56 : POPE 14 : DOPE 14 : PSM 42 : GM3 28 |
| Lower | CHOL 75 : POPC 15 : DOPC 15 : POPE 60 : DOPE 60 : POPS 24 : DOPS 21 : SAP6 30 |

Only the coarse-grained frames are provided. To build an all-atom system you also need an all-atom structure of the
same receptor for `--all-atom`:

```bash
python -m membraneforger --all-atom your_receptor.pdb \
    --coarse-grain examples/preequilibrated_gpcr_cellmembrane/8F7W_KOR_run1_30us.gro --out kor_run1
```

MembraneForger's reader and Martini 3 classification accept all eighteen files. No complete build from these frames
is included in the test results of this repository.
