# Equilibrating a build

`em.gro` is minimized, not equilibrated: the seam of the cut membrane, the water around the complex and the lipids
that were pushed aside have never moved at finite temperature. Run the restrained equilibration below before any
production simulation. The templates in `membraneforger/mdp/` follow the CHARMM-GUI membrane schedule (six steps,
1.875 ns in all) and set every restraint switch the topology expects, so the files of a build can be used as they are.

## In one command

```bash
bash membraneforger/mdp/run_equilibration.sh my_build            # gmx on PATH, 8 threads
bash membraneforger/mdp/run_equilibration.sh my_build gmx_mpi 16
```

The script runs `grompp` and `mdrun` for `step6.1_equilibration` to `step6.6_equilibration` in `my_build/`, each
restrained to the positions of `em.gro`, then writes `step7_production.tpr` and prints the `mdrun` command for it.
On a cluster, put the same lines in a job script; each step is one `grompp` plus one `mdrun`.

## What the steps do

| Step | Ensemble | dt (fs) | Length | `POSRES_FC_BB` | `POSRES_FC_SC` | `POSRES_FC_LIPID` | `DIHRES_FC` |
|---|---|---|---|---|---|---|---|
| 6.1 | NVT, velocities generated | 1 | 125 ps | 4000 | 2000 | 1000 | 1000 |
| 6.2 | NVT | 1 | 125 ps | 2000 | 1000 | 400 | 400 |
| 6.3 | NPT semi-isotropic | 1 | 125 ps | 1000 | 500 | 400 | 200 |
| 6.4 | NPT | 2 | 500 ps | 500 | 200 | 200 | 200 |
| 6.5 | NPT | 2 | 500 ps | 200 | 50 | 40 | 100 |
| 6.6 | NPT | 2 | 500 ps | 50 | 0 | 0 | 0 |
| 7 | NPT, no restraints | 2 | your choice | - | - | - | - |

Force constants are in kJ mol^-1 nm^-2 (positions) and kJ mol^-1 rad^-2 (dihedrals). Temperature 310 K with
v-rescale (Nose-Hoover in production), pressure 1 bar semi-isotropic with C-rescale (Parrinello-Rahman in production),
CHARMM36 non-bonded settings (1.2 nm force-switched LJ from 1.0 nm, PME), `h-bonds` constraints. The three coupling
groups are the index groups of `index_ini.ndx`: `Protein_LIG`, `MEMB` and `SOL_ION`. Change `ref_t` and `gen-temp` in
every file to simulate at another temperature.

## The restraint switches

The topology a build writes declares four position-restraint macros (`topol.top` defines them as 0 when nothing sets
them) and includes the restraint rows under `#ifdef POSRES`; the lipid topologies carry their CHARMM-GUI dihedral
restraints under `#ifdef DIHRES`. The mdp `define` line sets them:

| Macro | Restrains | Where it comes from |
|---|---|---|
| `POSRES_FC_BB` | protein backbone heavy atoms (N, CA, C, O and the cap atoms), and every heavy atom of a ligand | `toppar/posre_PRO*.itp`, `toppar/posre_<LIGAND>.itp` (a ligand whose `.itp` brought its own block uses the same macro) |
| `POSRES_FC_SC` | protein side-chain heavy atoms | `toppar/posre_PRO*.itp` |
| `POSRES_FC_LIPID` | lipid heavy atoms | `toppar/posre_<LIPID>.itp`, for lipids that came without a restraint block |
| `DIHRES` / `DIHRES_FC` | lipid dihedral restraints: ring chairs, double bonds, the GM3 and PSM amides, inositol | `[ dihedral_restraints ]` of the lipid `.itp` files |
| `BBDIHRES` | the EM-only backbone phi/psi restraints of `--backbone-restraints` | **never define it in MD**; the templates do not |
| `POSRES_FC_LIG` | nothing (kept at 0 for topologies of earlier releases) | - |

Ligand heavy atoms were frozen in the two minimizations (`freezegrps = LIG` in `int/em.mdp`). The templates do not
freeze anything: the ligand follows the backbone schedule through `POSRES_FC_BB`, so it relaxes gently with the
protein. If you write your own mdp files, keep `-DPOSRES` with the four force constants and `-DDIHRES -DDIHRES_FC=...`
in the equilibration steps, and leave every define out of production.

## Checking the result

A later rebuild into the same folder is refused once `step6.*.cpt` or `.xtc` files exist, unless
`--overwrite --overwrite-restart` is given; MembraneForger never deletes them.

`python -m membraneforger.equilibration` grades the equilibrated system (box convergence, area per lipid, thickness,
core hydration, protein tilt, lateral diffusion, area compressibility) from the energy file and evenly spaced frames:

```bash
cd my_build
echo "Box-X Box-Y Box-Z" | gmx energy -f step7_production.edr -o box.xvg
gmx trjconv -s step7_production.tpr -f step7_production.xtc -o frame.gro -sep -pbc mol -dt 1000   # one frame per ns
python -m membraneforger.equilibration --out . --energy box.xvg --frames frame*.gro --temperature 310
```

The report is `equilibration_validation.md`; the [tutorial](membraneforger_tutorial.md#59-gates-3-and-4-equilibration-and-equilibrium-properties)
explains what each record measures and its thresholds. The box z shrinks by 1.4 to 1.6 nm over the six steps
([Advanced options](advanced.md#choosing---z-pad-nm)); the default 2.25 nm water padding leaves about 2.7 nm between
periodic images of the solute afterwards, which the report checks.
