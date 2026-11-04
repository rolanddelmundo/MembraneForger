# Your first build

This page builds the bundled glucagon receptor (GCGR) example from start to finish. Nothing needs to be downloaded or
oriented. The [tutorial](membraneforger_tutorial.md) explains what each check measures, but you don't need it here.

## 1. Install and check

```bash
git clone https://github.com/rolanddelmundo/MembraneForger.git
cd MembraneForger
pip install ".[figures]"
membraneforger --version
```

GROMACS must be on your `PATH` (or passed with `--gmx`). On a cluster, load its module first. Then run the pre-flight
check. It takes seconds and writes nothing:

```bash
membraneforger --check \
    --aa membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_prot-lig.pdb \
    --cg membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_cg_cellmem.gro \
    --orientation none --out example_out
```

Fix every `FAIL` line before building. `WARN` lines don't stop the build.

## 2. Build

Run the same command without `--check`:

```bash
membraneforger \
    --aa membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_prot-lig.pdb \
    --cg membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_cg_cellmem.gro \
    --orientation none --out example_out
```

This takes 15 to 30 minutes on 8 cores. Progress is written to `example_out/membranebuilder.log`, and its last line is
`PASS` or `FAIL`. To restart an interrupted build, add `--overwrite`.

## 3. Read the result

In `example_out/`:

| File | What it is |
|---|---|
| `em.gro`, `topol.top`, `toppar/`, `index_ini.ndx` | The minimized system. Written only if every check passed. |
| `membranebuilder.log` | The log. Read every `WARN:` line. |
| `audit.json` | Final checks on the system. |
| `membrane_validation.md` | Membrane quality report. PASS WITH WARNINGS is normal for rare lipids. |
| `run_manifest.json` | Settings, versions and checksums. Attach it to bug reports. |
| `int/` | Intermediate files. A failed build leaves `int/em.unverified.gro` here. |

Open `em.gro` in a viewer. The protein should span the bilayer with its extracellular side at +z, the ligands should
be where you put them, and there should be no water inside the membrane.

## 4. Equilibrate

`em.gro` is minimized, not equilibrated. To run the six restrained steps (1.875 ns):

```bash
MF="$(python -c 'import membraneforger, os; print(os.path.dirname(membraneforger.__file__))')"
bash "$MF/mdp/run_equilibration.sh" example_out
```

To grade the result, see [Equilibration](equilibration.md#checking-the-result).

## 5. Your own protein

Prepare the PDB as in [Preparing the input](preparing_inputs.md), add topologies for any ligands ([Ligands](ligands.md)),
then run the check:

```bash
membraneforger --check --aa complex.pdb --cg KOR_Membrane1 --orient-chain R --nterm-side out --out my_build
```

- `--cg 1` or `--cg 2` picks a random frame. To get a repeatable build, name a frame instead (`--list-membranes`
  shows all 18).
- The bundled membranes are 18.3 x 18.3 nm, sized for a GPCR. Proteins of other sizes need
  [these options](preparing_inputs.md#size-and-the-bundled-membranes).

If a build stops, see [Troubleshooting](troubleshooting.md).
