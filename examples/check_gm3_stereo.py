#!/usr/bin/env python3
#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Check the sugar stereocentres and ceramide trans bonds of every GM3 (GLPA) in a built system.
#//=============================================================
"""Report GM3 ganglioside sugar stereocentres and ceramide trans bonds with the wrong configuration in a .gro or .pdb.

    python examples/check_gm3_stereo.py em.gro            # or any frame of a trajectory, e.g. from gmx trjconv
    python examples/check_gm3_stereo.py membrane.pdb

The same test the build's audit runs (membraneforger/stereo.py): all 16 sugar stereocentres (glucose C1-C5, galactose
C1-C5, sialic acid C2 and C4-C8) with the chirality definitions of backmap_data/map.dat, and the ceramide C4=C5 double
bond and amide with its [ trans ] definitions. Neither a stereocentre nor a C=C bond can invert in a classical MD
run, so the result for the first frame holds for the whole trajectory. Exit code 1 when any molecule is wrong.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from membraneforger.stereo import check_gm3_stereo  # noqa: E402


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    path = Path(sys.argv[1])
    report, fails = check_gm3_stereo(path, REPO / "backmap_data")
    if not report["molecules"]:
        print(f"{path}: no GLPA or GM3 residues")
        return 0
    n = report["molecules"]
    print(f"{path}: {n} GM3 molecules, {report['stereocentres_per_molecule']} sugar stereocentres and "
          f"{report['trans_bonds_per_molecule']} ceramide trans bonds each")
    for position, wrong in report["wrong_by_centre"].items():
        print(f"  {position:10s}: {wrong:3d} wrong ({100 * wrong / n:5.1f}%)")
    for bond, cis in report["cis_by_bond"].items():
        print(f"  {bond:16s}: {cis:3d} cis   ({100 * cis / n:5.1f}%)")
    affected = report["affected_molecules"]
    print(f"molecules with at least one wrong centre or bond: {len(affected)}/{n}"
          + ("" if not affected else " -> molecule(s) " + ", ".join(f"{m['molecule']} (residue {m['first_residue']})" for m in affected)))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
