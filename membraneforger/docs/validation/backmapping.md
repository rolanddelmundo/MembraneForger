# Backmapping and assembly { #4-backmapping-and-the-rest-of-the-build }

How the cut membrane is rebuilt with all atoms around the unchanged complex, then solvated and minimized.


mstool backmaps the sliced membrane around the complex held as a rigid obstacle, the structure is checked for ring
threading and stereochemistry (another seed is tried when needed; contact findings between the membrane and the complex are logged as a warning and left to the minimization), the topology is built, the box is set in
z, water and 0.15 M NaCl (in the water, see below) are added, the system is minimized with and then without the lipid dihedral restraints (the restrained run to Fmax < 1000 kJ/mol/nm, the free run that is published to Fmax < 500; the heavy atoms of every ligand, meaning each non-protein residue of the input, are frozen in both, so the receptor, lipids and solvent relax around it; ligand hydrogens and all protein atoms stay free), and
an independent audit must pass before `em.gro` appears. Every intermediate is kept: at the end of each build, passed
or failed, the steps on the way to `em.gro` (`config.INTERMEDIATES`) are moved into `<out>/int/`, beside the working
files (`int/work`) and the audit's scratch runs (`int/audit`); the top level holds the final outputs and the reports.
A retried backmapping seed keeps the earlier attempt (`int/work/mstool_before_seed<N>`,
`int/work/membrane_seed<N>.pdb`).

## 4.1 Disulfides

Disulfides are read from the all-atom input's geometry at the topology stage, which is the first stage with
anything to bond (orientation, placement, slicing and backmapping move the complex as one rigid body, so detecting
them earlier would give the same answer): two cysteine SG atoms closer than 3.0 A are bonded, a pair whose SG-SG
distance lies outside 1.8-2.2 A is rebuilt at 2.03 A before pdb2gmx, the topology is read back and must contain
exactly the detected pairs, and the pairs are listed in the log and in `run_manifest.json` (`disulfides`). The
audit measures every SG-SG bond of the topology on the minimized structure and fails the build if one is outside
1.8-2.2 A. Limitations: SSBOND records are not read and detection is by distance alone, so a reduced pair modelled
within 3.0 A becomes a bond and a pair stretched beyond 3.0 A in a model is missed (check the `disulfides` list in
the manifest against what you expect); a disulfide between two chains is refused, because each chain gets its own
topology; a cysteine within 3.0 A of two others is refused as ambiguous.

`--exclude-disulfide CHAIN:RESID-CHAIN:RESID` answers the first of those limitations: the named pair is left
reduced although its SG atoms are within 3.0 A, so a close contact a predicted model produces where the chemistry
does not support a bond does not become one. pdb2gmx is run with `-ss` for that chain and told not to link it, and
both cysteines are built as reduced CYS with their HG1 hydrogen; coordinates are untouched. See
[advanced options](../advanced.md#disulfides).
