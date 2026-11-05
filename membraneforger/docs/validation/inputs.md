# Inputs { #1-inputs }

What the build reads, and the check the membrane frame must pass before anything is embedded.


- `--aa complex.pdb`: the all-atom complex, the only source of protein and ligand chemistry.
- `--cg 1 | 2 | NAME | FILE`: a bundled GPR139 or kappa opioid receptor frame, random (`1`, `2`) or by name (`--list-membranes` lists them) (18.3 x 18.3 nm, about 1,300 lipids, ten
  species, 50 us of Martini 3), or your own Martini 3 frame.
- `--orient-chain R` names the chain that spans the membrane; `--nterm-side in|out` is needed by PPM when no OPM
  entry exists. The complete complex moves as one rigid body.

## 1.1 A frame must be periodic in its box

Everything after this point (slicing images, the seam, Voronoi cells, RDFs, equal-size windows) uses the periodic
cell on the frame's box line. A frame whose coordinates were rotated about z after the simulation (the usual cause is
a rotational fit, `gmx trjconv -fit rotxy+transxy` or similar, written without transforming the box) is no longer
periodic in that box: wrapping folds the corners of the rotated square onto its edges, beads of different molecules
land on top of each other there, the opposite corners stay empty, and the "diamond" one sees when rendering the file
is real data, not drift. The bundled 50 us frames are periodic in their boxes (`membraneforger/examples/align_frames.py`
finds no rotation in any of them); the same script rotates a user frame back in place (coordinates and velocities).

MembraneForger detects this when it reads a frame (`validation.align_frame_to_box`): it counts bead pairs of
different lipids closer than 0.30 nm after wrapping (an equilibrated Martini membrane has none; a rotated frame has
hundreds to thousands), scans the 90 degrees a square box allows for the rotation that removes them, refines it to
0.01 degree, and rotates the whole frame, water and ions included, back about the box centre. The angle and the
overlap counts before and after are logged (`WARN`) and recorded in `run_manifest.json` (`coarse_grain.frame_alignment`).
A frame that no rotation makes periodic (a non-square box, or coordinates that simply do not belong to the box line)
stops the build with a message saying so. The check stays on for every frame, so a user frame with the same
artefact is corrected the same way.
