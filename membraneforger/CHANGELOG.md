# Changelog

## 1.1.0 (2026-10-18)

MembraneForger now needs only an all-atom structure. It orients the complex in the membrane, embeds it in a bundled
pre-equilibrated Martini 3 cell membrane, cuts the membrane to size, and validates every stage before and after
backmapping.

Two defaults changed since 1.0.0, and the review decisions below change more. Commands written for 1.0.0 may need updating:

- **Orientation is on by default** (`--orientation auto`). An input with several protein chains now needs
  `--orient-chain`, and PPM needs `--nterm-side` when no OPM entry matches. Add `--orientation none` to use your
  coordinates as given, as 1.0.0 did.
- **Intermediates move to `<out>/int/`.** `oriented.pdb`, `membrane.pdb`, `boxed.gro`, the minimization files and
  the other steps on the way to `em.gro` are kept there instead of in the top-level output directory.

### Review decisions, 50 µs frames and repository layout

- **Review decisions (defaults changed).**
  - The slice window is centred by default; the offset search is opt-in (`--slice-offset-search`). Construction deviation median/95th/max on the eight bundled frames: 0.71/1.88/2.19 % centred, 1.55/3.39/4.87 % searched.
  - Counterions default to `additive` (neutralizing ions plus the requested salt pairs); `--counterion-mode inclusive` restores the old behaviour. The manifest and log report neutralizing ions, salt pairs, Na/Cl totals and ionic strength separately.
  - The slice-stage thickness record was dropped from the report by a name-collision in `all_records()`; it is now graded (3 % PASS / 5 % WARNING). Protein tilt and insertion depth at the slice are diagnostics (NOT RUN), since slicing does not move the complex.
  - A build into an `--out` that already holds an earlier build is refused unless `--overwrite` is given; checkpoints and trajectories also need `--overwrite-restart`. Slurm reruns into an existing folder must pass `--overwrite` (for example through `EXTRA`).
  - The default `--out` is `./<all-atom structure stem>_membraneforger` for a single build and for `--runs`; duplicate run folders are refused.
  - Geometry options are in nm: `--box-nm`, `--bilayer-z-nm`, `--fit-max-core-rmsd-nm`. The Å spellings still work with a conversion warning. `--no-slice-offset` is a hidden no-op.
  - New accessible-shell-area normalization of the lateral headgroup RDF (`rdf.lateral_rdf_shell`), reported separately and not used by any gate; `examples/validate_rdf_normalization.py` compares it with the existing curve and with random geometry-matched references.

- **Definitions of the graded quantities** (tutorial 5.13): the formulas the code computes for the Voronoi area, g(r), the
  protein-density profile, thickness, tilt, core water, interdigitation, ion count, drift, S_CD, MSD, K_A and leaflet
  tension, with the conventions that were implicit (S_CD sign, fit windows, periodic-cell D, anchor-based APL). The
  composition distance is written D_comp. The seam-repulsion docstring now states the force and the potential.
- **The coarse-grained frame is graded before embedding** (gate 0, `membrane_validation.md`): both leaflets populated, absolute
  leaflet APL and bilayer thickness, bilayer tilt, lipids on the wrong leaflet, broken molecules and lipid-protein bead
  overlaps, as PASS / WARNING / FAIL records. The eight bundled frames pass. A FAIL stops the build
  unless `--cg-validate no`.
- **Repository layout.** The repository root is now `README.md`, `LICENSE.txt`, `AUTHORS`, `CITATION.cff` and the
  `membraneforger/` folder, as in BaNDyT. `forcefield/`, `backmap_data/`, `examples/`, `docs/`, `slurm/`,
  `CHANGELOG.md` and `ruff.toml` moved into `membraneforger/` (`python -m membraneforger` still runs from the root;
  Slurm: `bash membraneforger/slurm/setup.sh`). The command-line, tutorial, advanced and troubleshooting guides are back in
  `docs/` and `CITATION.cff` is back in the root; both had been moved into `archive/`, which broke the README links.
  `archive/` (the unmaintained 0.1/0.2 workflow, with its licence texts and notices for the vendored mstool and Martini
  files) is removed from main; commit `7e26de9` still holds it. `membraneforger/cpu.submit` is now
  `slurm/examples/array_cases.sbatch`, `membraneforger/example.py` is `examples/build_from_python.py`, the module map
  is `docs/modules.md` and `test_contact_gates.sbatch` is in `slurm/examples/`. The README is shorter. No code logic,
  option, output or data file changed; only the folders the program looks in for its bundled files.
- **Bundled frames and example files have clearer names.** `KOR1` to `KOR9` are now `KOR_Membrane1` to `KOR_Membrane9`,
  `GPR5` is `GPR139_Membrane5`, and `6WHC_MTZP_*` is `GCGR_MTZP_*`; the tutorial and guides call the example complex
  GCGR. `--cg` still accepts the old names (`KOR5`, `GPR5`). Output folders of `--runs` take the new names.
- **Bundled membranes are now 50 µs frames.** The 18 frames at 30 µs are replaced by the 50 µs frame of every
  simulation that reached 50 µs: `KOR1`, `KOR2`, `KOR3`, `KOR5`, `KOR6`, `KOR7`, `KOR9` and `GPR5` (the number is the
  simulation run, so the numbering has gaps). `KOR4`, `KOR8` and `GPR1` to `GPR4`, `GPR6` to `GPR9` are gone: the
  GPR139 series now holds one frame, so `--cg 1 --runs N` accepts N = 1 only. The earlier frames stay in the git
  history. The leaflet composition is unchanged at the precision of the table.
- **Tutorial numbers re-measured on the 50 µs frames.** Lipids removed by embedding, pocket volumes, single-helix
  pockets, seam removals, the slice-offset effect and the construction deviation now come from the eight bundled
  frames (`membraneforger/examples/measure_tutorial_numbers.py`); on these frames the offset search no longer reduces the whole-cell
  spread. The `KOR1` stage table and its RDF statements are re-measured as well. Remaining 30 µs-frame, GIPR and unmeasured statements in the tutorial were removed or replaced by 50 µs values; gate 0 is now in the tutorial's flow diagram and settings, `advanced.md` uses Upper/Lower like the other guides, and the README mstool pin matches the supported range (0.3.9 or 0.3.10). New tutorial figures and tables: the KOR1 membrane before and after embedding, the per-frame APL table (section 5.3) and the lipid density around the protein (section 5.6a), drawn by `examples/make_tutorial_figures.py`. Removed the unused `pipeline.refuse_wrong_stereochemistry` (the backmapping verdict already stops on a wrong configuration) and two unused variables in `audit.py`; the README steps and the module map now list gate 0, `contacts.py` and `stereo.py`.
- **Series and `--runs` limit come from the frames present.** The nine-frame limit is gone: `--runs` is bounded by
  the number of frames the series holds, and the error message names them; `--list-membranes` and the unknown-name
  error list the names that exist.
- **`--runs N`** builds the complex in N different bundled membranes, the first N frames of the `--cg 1` or `--cg 2`
  series (at most the number of frames the series holds), one validated system per run in `<out>/<frame>/`. `RUNS` in `membraneforger/slurm/run_membraneforger.sbatch` passes it.
- **Bundled frames by name.** `--cg NAME` builds in one particular bundled membrane (`--list-membranes` lists the names) instead of a random one. `--cg 1`, `--cg 2` and `--cg FILE`
  are unchanged.
- **Portable Slurm script.** `slurm/run_membraneforger.sbatch` keeps every cluster-specific value in three marked
  blocks (scheduler, software, this run) with working defaults, so only those are edited; submitted as is it builds
  the bundled 6WHC example. `slurm/examples/gemini_6WHC.sbatch` is the filled-in Gemini version.

### New

- **Single-input builds.** `--cg 1` (default) and `--cg 2` embed the complex in a bundled GPR139 or kappa opioid
  receptor membrane: the 50 µs Martini 3 frames of an asymmetric ten-species cell membrane (eight frames, see below), in
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

- A membrane contact or bond crossing found before EM no longer rejects a backmapping seed or stops the build: it is logged as a
  warning and the check after EM (0.6 of Rmin, no crossings) decides. Another seed leaves the same overlaps (the quick start and 28
  GLP1R builds start with 23 to 585 such pairs) and plain EM removes them, so retrying only added minutes.
  `membrane_overlap_before_em = "stop"` restores the earlier behaviour.
- **Ligand force-field packages trimmed.** The four packages under `forcefield/` (MTZP, MRTR, MSEM, V6G) keep only the files the build reads (the ligand `.top`, `charmm36.ff/ffbonded.itp`, and for V6G `V6G_gmx.top`, `V6G_base_ff_params.itp` and `charmm36.ff/V6G_ffbonded.itp`) plus their `.pdb` and `posre` files; the unused copies of `charmm36.ff` inside them are gone (`forcefield/` 64 MB to 22 MB). The build uses the shared `forcefield/charmm36.ff/`, so no output changes.
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

- The GM3 stereochemistry audit makes every molecule whole across the periodic boundary before testing, using the
  `.gro` box or the PDB `CRYST1` record. A GPU `mdrun` writes minimization output with molecules wrapped into the
  cell atom by atom, which the audit read as inverted sugar stereocentres — a false positive that failed correct
  builds.
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
unchanged under `archive/` (removed in 1.1.0, see git history) and is no longer maintained. Tagged as `v1.0.0` without a GitHub release.

## 0.2.0 (2026-09-02) and 0.1.0 (2026-07-15)

The earlier four-stage workflow. See the GitHub releases.
