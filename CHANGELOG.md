# Changelog

## 2.0.0 (2026-10-10)

MembraneForger now needs only an all-atom structure. It orients the complex in the membrane, embeds it in a bundled
pre-equilibrated Martini 3 cell membrane, cuts the membrane to size, and validates every stage before and after
backmapping.

This is a major version because two defaults changed. Commands written for 1.0.0 may need updating:

- **Orientation is on by default** (`--orientation auto`). An input with several protein chains now needs
  `--orient-chain`, and PPM needs `--nterm-side` when no OPM entry matches. Add `--orientation none` to use your
  coordinates as given, as 1.0.0 did.
- **Intermediates move to `<out>/int/`.** `oriented.pdb`, `membrane.pdb`, `boxed.gro`, the minimization files and
  the other steps on the way to `em.gro` are kept there instead of in the top-level output directory.

### New

- **Single-input builds.** `--cg 1` (default) and `--cg 2` embed the complex in a bundled GPR139 or kappa opioid
  receptor membrane: 18 frames of an asymmetric ten-species cell membrane, 30 µs of Martini 3 each, in
  `examples/preeq_cg_cellmem/`. `--cg FILE --embed` embeds into a frame of your own.
- **Membrane orientation** from an exact OPM/OPRLM entry or a local PPM 3.0 run, applied to the whole complex as one
  rigid body (`--orientation`, `--orient-chain(s)`, `--nterm-side`, `--pdb-id`, `--ppm-exe`).
- **`--orient-residues`** orients on a transmembrane segment only, or with `--orientation none` names the embedded
  segment that sets the bilayer centre.
- **`--embed-site free`** places a small protein, such as a single transmembrane helix, on unbroken membrane away
  from the frame's receptor hole.
- **Automatic box (`BOX = auto`).** The coarse-grained membrane is cut to the complex plus `--xy-buffer` before
  backmapping, lipid by lipid under periodic boundaries, with a density-preserving window offset and seam
  relaxation. `--box X Y Z` sets the box instead.
- **Slice gate** before backmapping: per-leaflet Voronoi area per lipid, composition, lateral headgroup RDF,
  thickness and protein orientation against the source membrane (`--apl-validate`, `--apl-tolerance`,
  `--rdf-validate`, `--no-slice-offset`).
- **Membrane validation report** for every build: `membrane_validation.md` / `.json`, a per-lipid table and figures,
  including per-species density around the protein.
- **Equilibration gates.** `python -m membraneforger.equilibration` grades convergence, order parameters, diffusion,
  area compressibility and leaflet tension on a finished trajectory.
- **Lipid composition edits** with `--dellipid` and `--addlipid`.
- **Slurm support.** `slurm/setup.sh` builds a self-contained environment (Python packages, GROMACS, PPM 3.0
  compilation) and `slurm/run_membraneforger.sbatch` runs one structure per array task.

### Changed

- Lipids that overlap the embedded complex are pushed aside at the coarse-grained level, and only those still
  overlapping are removed. An empty pocket larger than 1.0 nm³ in a leaflet stops the build.
- Minimization runs first with the lipid dihedral restraints, then without them, and both are audited.
- When minimization ends above the force target, the log names the atom and its neighbours, and
  `em_clash_trace.json` traces a lipid culprit back to the first stage at which it clashes.
- Disulfides are detected from the input geometry, listed in `run_manifest.json`, and their SG–SG distances are
  audited after minimization.
- Every coarse-grained frame is checked for coordinates rotated relative to its box and is rotated back when needed.
  The bundled frames were corrected on disk.

### Fixed

- GM3 backmapping: swapped glucose and galactose bead assignments, six wrong chirality definitions, the missing sialic
  acid C7 chirality definition and the missing ceramide trans definitions. Every GM3 sugar stereocentre is now
  checked in the built membrane.
- GM3 molecules split into CHARMM sugar residues, and GM3 and PIP2 in `.gro` stages, are read correctly.
- A `--box` z too short for the complex is refused before backmapping.

### Documentation

- The README is rewritten as a short getting-started guide. The full option table moved to
  `docs/cli_reference.md`, and orientation, membrane, box and output details moved to `docs/advanced.md`.
- New `docs/troubleshooting.md`.
- The bundled-example command in `examples/preeq_cg_cellmem/README.md` now includes `--orientation none`, which the
  multi-chain example needs under the new orientation default.

## 1.0.0 (2026-10-03)

Two-input pipeline: an all-atom protein–ligand complex plus one Martini 3 coarse-grained frame of the same complex
become a validated, energy-minimized CHARMM36/GROMACS system. The earlier four-stage workflow (0.1.0, 0.2.0) is kept
unchanged under `archive/` and is no longer maintained. Tagged as `v1.0.0` without a GitHub release.

## 0.2.0 (2026-09-02) and 0.1.0 (2026-07-15)

The earlier four-stage workflow. See the GitHub releases and `archive/`.
