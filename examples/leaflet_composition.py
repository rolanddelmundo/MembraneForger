"""Per-leaflet lipid composition of Martini 3 membrane frames, averaged over all the frames given.

    python examples/leaflet_composition.py examples/preeq_cg_cellmem/*_cg_cellmem.gro

Lipids are assigned to the upper or lower leaflet by the mean z of their beads relative to the bilayer midplane.
GM3 is split into consecutive CER/GLC/GAL/NMC residues in Martini 3 frames and is counted by its CER residue.
Percentages are mole percent of each leaflet, averaged over frames and rounded to whole numbers.
"""
import sys
from collections import Counter, defaultdict
from statistics import mean

LIPIDS = ("CHOL", "POPC", "DOPC", "POPE", "DOPE", "PSM", "GM3", "POPS", "DOPS", "SAP6")
ALIASES = {"CER": "GM3", "DPG3": "GM3", "DPSM": "PSM"}


def leaflets(path):
    """Return (upper, lower) Counters of lipid molecules for one .gro frame."""
    lines = open(path).read().splitlines()
    molecules = defaultdict(list)  # (residue number, name, chunk) -> bead z values
    previous, chunk = None, 0
    for line in lines[2:2 + int(lines[1])]:
        key = (int(line[0:5]), ALIASES.get(line[5:10].strip(), line[5:10].strip()))
        if key[1] not in LIPIDS:
            previous = None
            continue
        if key != previous:  # a new molecule, even when residue numbers repeat
            chunk += 1
            previous = key
        molecules[(key[0], key[1], chunk)].append(float(line[36:44]))
    midplane = mean(z for zs in molecules.values() for z in zs)
    upper, lower = Counter(), Counter()
    for (_, name, _), zs in molecules.items():
        (upper if mean(zs) > midplane else lower)[name] += 1
    return upper, lower


def main(paths):
    totals = {"upper": defaultdict(list), "lower": defaultdict(list)}
    for path in paths:
        for side, counts in zip(("upper", "lower"), leaflets(path)):
            n = sum(counts.values())
            for lipid in LIPIDS:
                totals[side][lipid].append(100.0 * counts[lipid] / n)
    print(f"{len(paths)} frame(s)")
    print(f"{'lipid':6} {'upper %':>8} {'lower %':>8}")
    for lipid in LIPIDS:
        up, lo = mean(totals["upper"][lipid]), mean(totals["lower"][lipid])
        print(f"{lipid:6} {round(up):8d} {round(lo):8d}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
