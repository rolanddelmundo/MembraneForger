# Limitations and settings

What the checks assume, what they cannot catch, and where the defaults live.

## 7. Limitations and assumptions

- The seam of a cut membrane was never equilibrated: the relaxation removes overlaps, it does not equilibrate the
  interface. Equilibrate the built system before production.
- The construction check isolates slicing; the whole-cell comparison carries the heterogeneity of the parent
  membrane around the protein (on the bundled frames the lower leaflet there has a 5-13 % smaller APL than the cell average, the upper leaflet -4 to +0.4 %).
  Neither the centred window (default) nor the offset search (`--slice-offset-search`) reduces that on these frames.
- Protein generators in the tessellation take area from the lipids next to them by the perpendicular-bisector rule;
  lipids whose heads sit above a protein loop at the headgroup level get small cells. The slab is bounded by the
  anchor plane (plus 0.2 nm) so that atoms in the water beyond the heads do not count.
- GM3's slicing anchor is the sugar centroid (as requested) while its atomistic anchor is the amide nitrogen; the two
  differ by the tilt of the sugars, which affects only the cross-resolution rows.
- Species RDFs and species APLs for 10-20 molecules are reported with N and should be read as such.
- Gates 3 and 4 have been tested on synthetic data (random walks, Gaussian area series, antisymmetric stress
  profiles) and on files, not yet on a production trajectory of a MembraneForger build; `MDAnalysis` is optional
  and only needed to read xtc/trr directly (write frames with `gmx trjconv -sep` otherwise).
- No GIPR build or trajectory is part of this repository; the before/after numbers above are from the bundled
  frames and the GCGR example, which show the same defect with the old rule and the same repair.
- The reader checks every frame, corrects a rotated user frame (`membraneforger/examples/align_frames.py` does the same in
  place) and refuses one whose coordinates do not fit its box at all rather than guessing.

## 8. Settings

`membraneforger/config.py` carries the defaults, each with the reasoning: `box_xy_buffer_nm` (1.2), the slice
the opt-in offset search (`slice_optimize_offset` False, `slice_offset_search_nm` 0.5, `slice_offset_step_nm` 0.1), the ion mode (`counterion_mode` additive), the seam
(`seam_min_bead_nm` 0.30, `seam_hard_core_nm` 0.15, `seam_relax_steps` 400), the force-field contact gate
(`contact_rmin_fraction` 0.6, `bond_crossing_nm` 0.12), the gate (`APL_VALIDATE`,
`APL_SLICE_WARNING_PERCENT` 3, `APL_SLICE_TOLERANCE_PERCENT` 5, `RDF_VALIDATE`, the RDF ceilings and
`slice_window_samples` 12 per axis). On the command line: `--apl-validate yes|no`, `--apl-tolerance PERCENT`,
`--rdf-validate yes|no`, `--cg-validate yes|no` (gate 0, `cg_validate`), `--slice-offset-search`, `--counterion-mode additive|inclusive`, `--xy-buffer NM`.
