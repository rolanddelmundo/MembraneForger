# Validation

MembraneForger takes an all-atom protein (or protein-ligand) structure and a pre-equilibrated Martini 3 membrane and
returns a validated, energy-minimized CHARMM36/GROMACS system. These pages explain what each stage does, how the build
is checked, and how to read the validation report (`membrane_validation.md`, `.json` and the figures) that every build
writes next to `em.gro`. New to MembraneForger? Start with [Your first build](../first_build.md).

<div class="mf-flow" markdown>
[Inputs](inputs.md){ .mf-flow-step }
[Placement and slicing](placement.md){ .mf-flow-step }
[Backmapping](backmapping.md){ .mf-flow-step }
[Membrane checks](membrane.md){ .mf-flow-step }
[Integrity and gates](integrity.md){ .mf-flow-step }
</div>

## At a glance: one validated build

GPR139 with dynorphin A (1&ndash;8) and its Gq protein (GC8N), built in the bundled frame `GPR139_Membrane1` on 16
threads. Every number below is read from that build's `run_manifest.json`, `audit.json`, `membrane_validation.json` and
log.

<div class="mf-stats" markdown>
<div><b>391</b><span>lipids kept, 203 upper / 188 lower</span></div>
<div><b>20 / 20</b><span>audit checks passed</span></div>
<div><b>0 of 679</b><span>lipid rings pierced</span></div>
<div><b>442</b><span>largest force after minimization, kJ mol<sup>-1</sup> nm<sup>-1</sup> (target 500)</span></div>
<div><b>0</b><span>water molecules in the bilayer core</span></div>
<div><b>10 min</b><span>wall time, backmapping 4.7 min</span></div>
</div>

The final system has 197,829 atoms: 44,798 waters, 271 Na<sup>+</sup> and 122 Cl<sup>-</sup>, net charge 0. Leaflet
areas per lipid stay within 0.2 &Aring;<sup>2</sup> from the slice to the minimized system (45.4 / 41.0 &Aring;<sup>2</sup>
after slicing, 45.3 / 40.8 &Aring;<sup>2</sup> after minimization), and the receptor tilt moves from 65.29&deg; to 64.80&deg;.
The audit reports three warnings worth reading: a few sugar and inositol ring chairs outside their restraints (GM3 4 of
45, PIP2 1 of 24) and three backbone residues that moved more than 60&deg; (no backbone restraints were used).

## The stages and their checks

```
equilibrated Martini membrane (frame)
    -> orientation of the complex (OPM / PPM)                                          (orient)
    -> gate 0: the frame is graded on its own, absolute limits                         (validate_cg_frame)
    -> placement of the complex into the frame                                         (embed)
    -> embedded membrane: lipids pushed aside, only those it occupies removed, no pocket (the REFERENCE for slicing)
    -> PBC-aware, whole-lipid, anchor-based crop to the complex + buffer               (slice)
    -> gate 1: APL, composition, lateral RDF, thickness, protein orientation, integrity (slice_check)
    -> mstool backmapping around the rigid complex                                    (backmap, assemble)
    -> gate 2: the all-atom membrane (membrane.pdb), then the minimized system (em.gro)
    -> topology, box, water, ions, restrained + free minimization (ligands frozen), audit   (publish em.gro)
    -> gate 3 and 4 on the equilibration / production trajectory                       (python -m membraneforger.equilibration)
```

| Page | Covers |
|---|---|
| [Inputs](inputs.md) | the files the build reads; gate 0, the check on the membrane frame |
| [Placement and slicing](placement.md) | embedding the complex, cutting the membrane to the complex plus a buffer |
| [Backmapping and assembly](backmapping.md) | rebuilding the lipids with all atoms, disulfides, solvation and minimization |
| [Membrane quality checks](membrane.md) | areas per lipid, composition, leaflet asymmetry, thickness, lipids around the protein |
| [Integrity checks and gates](integrity.md) | construction checks, gate 2 (backmapping), gates 3 and 4 (equilibration), lateral RDF, definitions |
| [Thresholds](thresholds.md) | every test and its PASS, WARNING and FAIL limits in one table |
| [Limitations and settings](limitations.md) | what the checks assume, and where the defaults live |

Most readers come for [slicing and the global APL](membrane.md#53-slicing-and-the-global-apl),
[species-specific APL](membrane.md#55-species-specific-apl) and the [lateral RDF](integrity.md#512-lateral-rdf).
