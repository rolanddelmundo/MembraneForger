# Your first build

This page takes you from a fresh clone to a minimized, audited system and says what to look at on the way. It uses the
bundled glucagon receptor complex, so nothing has to be downloaded or oriented. The [tutorial](membraneforger_tutorial.md)
is the reference for what each check measures; you do not need it to get through a first build.

## 1. Install and check

```bash
git clone https://github.com/rolanddelmundo/MembraneForger.git
cd MembraneForger
pip install ".[figures]"
membraneforger --version
```

GROMACS has to be installed separately and be a real executable on your `PATH` (or passed with `--gmx`); on a cluster,
load its module first. Then run the pre-flight check, which writes nothing and takes seconds:

```bash
membraneforger --check \
    --aa membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_prot-lig.pdb \
    --cg membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_cg_cellmem.gro \
    --orientation none --out example_out
```

Every line starts with `OK`, `WARN`, `FAIL` or `INFO`, and the last says `CHECK PASS` or `CHECK FAIL`. Fix every `FAIL`
before building (the message says what is missing). A `WARN` is worth reading but does not stop the build.

## 2. Build

Run the same command without `--check`:

```bash
membraneforger \
    --aa membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_prot-lig.pdb \
    --cg membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_cg_cellmem.gro \
    --orientation none --out example_out
```

It takes about 15 to 30 minutes on 8 cores. Progress goes to the screen and to `example_out/membranebuilder.log`, one
line per check, starting with `CHECK:`, `PASS:`, `WARN:`, `INFO:` or `ERROR:`. The stages run in the order of the
[How it works](../../README.md#how-it-works) list in the README: orientation, the membrane-frame check, embedding,
slicing and its gate, backmapping, topology, solvation and ions, minimization, audit. The last line of the log is
`PASS` or `FAIL`. If you interrupt a build, run it again with `--overwrite`; MembraneForger refuses to write into a
folder that already holds a build.

## 3. Read the result

In `example_out/`:

1. `membranebuilder.log`: the last line must say `PASS`. Search it for `WARN:` and read each one.
2. `audit.json`: `"result": "PASS"`, with the independent checks that were run on the final system.
3. `membrane_validation.md` (with the `membrane_validation_*.png` figures): the quality report of the cut membrane. Its
   verdict is PASS, PASS WITH WARNINGS, FAIL or INSUFFICIENT SAMPLING; the meaning of each is in the
   [tutorial](membraneforger_tutorial.md#510-overall-status). A PASS WITH WARNINGS is normal for rare lipid species that
   have only a few molecules in a small patch.
4. `em.gro`, `topol.top`, `toppar/` and `index_ini.ndx`: the system. Nothing is written as `em.gro` unless every check passed;
   a failed build leaves `int/em.unverified.gro` and the other intermediates for inspection, and
   [Troubleshooting](troubleshooting.md) says where to look for each stage.
5. `run_manifest.json`: the exact settings, software versions and checksums. Attach it when you report a problem.

Open `em.gro` in a viewer (VMD, PyMOL, ChimeraX) and look at it once: the protein should span the bilayer with its
extracellular side up (+z), the ligands should sit where they were in your input, and there should be no water inside
the membrane core.

## 4. Equilibrate

`em.gro` is minimized, not equilibrated:

```bash
bash membraneforger/mdp/run_equilibration.sh example_out
```

This runs the six restrained steps (1.875 ns) and writes `step7_production.tpr`; see [Equilibration](equilibration.md).
Afterwards `python -m membraneforger.equilibration` grades box convergence, area per lipid, thickness and the other
membrane properties; it needs `--out`, the box series and some frames from the production run, so use the full command in
[Equilibration](equilibration.md#checking-the-result).

## 5. Your own protein

Prepare the PDB following [Preparing the input](preparing_inputs.md) (one model, no water or bulk ions, complete heavy atoms,
termini, unique chain IDs; hydrogens are optional), add a topology for any ligand ([Ligands](ligands.md)), and run `--check` on it:

```bash
membraneforger --check --aa complex.pdb --cg KOR_Membrane1 --orient-chain R --nterm-side out --out my_build
```

Two choices affect reproducibility and fit:

- `--cg 1` and `--cg 2` draw a random frame of nine, so two identical commands give different membranes. Name a frame
  (`--cg GPR139_Membrane1`; `--list-membranes` shows all 18) for a repeatable build, or use `--runs N` to build the same
  protein in the first N frames and compare them.
- The bundled membranes are patches around a GPCR, 18.3 x 18.3 nm. A protein of GPCR size fits the hole the receptor left; a
  much smaller or larger one needs the options in [Preparing the input](preparing_inputs.md#size-and-the-bundled-membranes).

If a build stops, the log names the stage and the file to look at, and [Troubleshooting](troubleshooting.md) lists the
common causes. To report a problem, see [Contributing](../../CONTRIBUTING.md).
