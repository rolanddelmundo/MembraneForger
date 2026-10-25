# Pre-equilibrated coarse-grained cell membranes (Martini 3)

Thirteen frames of a single GPCR in an asymmetric, ten-species cell-membrane model, each taken at 50 µs of a
Martini 3 simulation (nine kappa opioid receptor frames and four GPR139 frames). They are the bundled membranes of MembraneForger: `--cg 1` picks one of the GPR139 frames, `--cg 2`
one of the kappa opioid receptor frames, and `--cg NAME` (for example `--cg KOR_Membrane5`) a particular frame.
`python -m membraneforger --list-membranes` lists them, and `--runs N` builds one complex in the first N frames of a
series (at most as many as the series holds: nine for KOR, four for GPR139), giving N replicate systems that differ in
their lipid arrangement.

| Files | Receptor | Replicates |
|---|---|---|
| `KOR_Membrane1` to `KOR_Membrane9` (`KOR_Membranen_cg_cellmem.gro`) | kappa opioid receptor (PDB 8F7W) | 9 |
| `GPR139_Membrane2` to `GPR139_Membrane5` (`GPR139_Membranen_cg_cellmem.gro`) | GPR139 (PDB 7F6G) | 4 |

The number in a name is the replicate (run) of the simulation the frame was taken from, so it is kept even though the
series has gaps: the frames are the 50 µs frame of runs 1 to 9 of 8F7W and of runs 2 to 5 of 7F6G, the
runs that reached 50 µs. Each is the unmodified `run50.gro` of its run. `--list-membranes` shows the names present.

Each frame holds about 60,000 beads in a cell of about 18.3 × 18.3 × 19.2 nm: the receptor (about 290 residues),
about 1,300 lipids, water and ions. The membranes were built with INSANE at CHOL:POPC:DOPC:POPE:DOPE:PSM:DPG3 =
25:20:20:5:5:15:10 (upper leaflet) and CHOL:POPC:DOPC:POPE:DOPE:POPS:DOPS:SAP6 = 25:5:5:20:20:8:7:10 (lower
leaflet). Composition at 50 µs, mole percent of each leaflet averaged over the 13 frames
(`python membraneforger/examples/leaflet_composition.py membraneforger/examples/preeq_cg_cellmem/KOR_Membrane?_cg_cellmem.gro membraneforger/examples/preeq_cg_cellmem/GPR139_Membrane?_cg_cellmem.gro`):

| Leaflet | CHOL | POPC | DOPC | POPE | DOPE | PSM | DPG3 | POPS | DOPS | SAP6 |
|---|---|---|---|---|---|---|---|---|---|---|
| Upper | 28 | 19 | 19 | 5 | 5 | 14 | 10 | 0 | 0 | 0 |
| Lower | 22 | 5 | 5 | 21 | 21 | 0 | 0 | 8 | 7 | 10 |

## Lipid contacts around the receptor

Number of lipids in contact with the receptor over the 50 µs of the simulation the frame comes from, per lipid species
and leaflet (the leaflet is assigned from the lipid's position relative to the bilayer midplane). A lipid counts when
any of its beads is within 5 Å (shell 0-5 Å) or 5-10 Å of any receptor bead, averaged over 501 frames (every 10th) of
each 5 µs segment; the first column, 0 µs, is the energy-minimized start. DPG3 counts as one molecule. A white cell
means no lipid of that species contacted the receptor in that leaflet.

| Simulation | Figure |
|---|---|
| 8F7W (kappa opioid receptor), run 5, the run of `KOR_Membrane5` | [KOR_Membrane5_lipid_contacts_0-50us.svg](figures/KOR_Membrane5_lipid_contacts_0-50us.svg) |
| 7F6G (GPR139), run 5, the run of `GPR139_Membrane5` | [GPR139_Membrane5_lipid_contacts_0-50us.svg](figures/GPR139_Membrane5_lipid_contacts_0-50us.svg) |

Each is a single replicate, so differences between species are not averaged over runs.

The figures come from the simulation trajectories, which are not part of the repository, with a MDAnalysis script (a
residue counts once per frame when any of its beads is within the cutoff of any receptor bead; segments `run5` to `run50`
at 5 µs each, every 10th frame). The text in the SVGs is editable Arial, 12 pt and larger.

`GCGR_MTZP_cg_cellmem.gro` and `GCGR_MTZP_prot-lig.pdb` are a complete input pair (glucagon receptor–tirzepatide–Gs
in the same membrane model, 11.2 × 11.2 nm) used by the tests:

```bash
python -m membraneforger --aa membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_prot-lig.pdb \
    --cg membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_cg_cellmem.gro --orientation none --out example_out
```

Every frame is periodic in its box: `python membraneforger/examples/align_frames.py membraneforger/examples/preeq_cg_cellmem/*_cg_cellmem.gro` finds
no bead pairs of different lipids within 0.30 nm once wrapped. MembraneForger runs the same check on every frame it
reads (`align_frame_to_box`) and corrects a rotated frame with a warning that names the angle.
