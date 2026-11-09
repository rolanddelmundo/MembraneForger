# Questions

## Will my protein work?

If it is an all-atom structure with complete side chains, one model, a unique chain ID per chain and a known membrane
orientation, yes. Ligands need a CHARMM/CGenFF topology. Go through [Preparing the input](preparing_inputs.md), then run
`membraneforger --check`, which reports what the build will do with every chain.

## How long does a build take?

About 15 to 30 minutes for a GPCR complex on 8 cores, with about 1.5 GB of memory. Backmapping and minimization take
most of it.

## Is the output ready for production MD?

No. `em.gro` is energy-minimized, not equilibrated. Run the six restrained steps first (1.875 ns, templates included),
then production. See [Equilibration](equilibration.md).

## Does MembraneForger parameterize my ligand?

No. It uses the topology you provide. GDP and Mg<sup>2+</sup> are bundled. Check CGenFF penalty scores before trusting
a ligand. See [Ligands](ligands.md).

## My protein is not a GPCR. Can I still use it?

Yes. The bundled membranes are 18 x 18 nm and sized for a GPCR. Small and single-pass proteins can be placed on intact
bilayer with `--embed-site free`, and any Martini 3 membrane frame of your own works with `--cg my_frame.gro`. See
[Examples](examples.md).

## What do I get from several runs?

`--runs N` builds the protein in N different bundled membrane frames, so each system starts with a different lipid
arrangement. These are independent starting configurations, not MD replicates of one system.

## The build passed, but the log has warnings.

The build is usable. A measurement is near a limit or poorly sampled, often a rare lipid species. Read each warning and
the matching section of [Validation](validation/index.md). If the build stops, see [Troubleshooting](troubleshooting.md).
