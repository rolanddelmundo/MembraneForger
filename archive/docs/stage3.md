# Stage 3

Stage 3 consumes `outputs/<run_id>/stage2/cg_backmap_input.gro` and the normalized all-atom reference input. `inputs/opm_reference.pdb` is required only when `stage3.placement.use_opm: true`.

Stage 3 imports `mstool` from `resources/vendor/mstool` by default and records the resolved module path in provenance. It must not rely on `work/runtime_vendor/mstool`.

Before backmapping starts, Stage 3 introspects `mstool.Ungroup` and validates the exact arguments MembraneForger will pass to the pinned vendored API. GM3/DPG3 mapping and XML additions are loaded only when the Stage 2 topology or enabled glycolipid finalizer requires them.

`fixed` means the protein coordinates already present in the backmapped CG-to-AA structure. `moving` means the AA protein/ligand reference coordinates. The transform is applied to the appropriate moving AA atoms and is not silently swapped.

The canonical Stage 3 structure output is `outputs/<run_id>/stage3/final_all_atom.pdb`. Optional chemistry adapters may transform the raw `step4_final.pdb` before that canonical output is written.
