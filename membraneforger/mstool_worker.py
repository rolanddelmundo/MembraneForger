#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// mstool compatibility layer: the only file that imports mstool.
#// Runs as a child process:  python mstool_worker.py job.json
#//=============================================================
"""Report the installed mstool mapping, or backmap one membrane around a rigid all-atom complex."""
import functools
import json
import os
import sys
from pathlib import Path

__all__ = ['TESTED_MSTOOL', 'malformed_chirals', 'flipped_chirals', 'describe', 'forcefield_files', 'backmap', 'main']

TESTED_MSTOOL = ("0.3.9", "0.3.10")  # versions whose lipid mapping and Backmap(rock=...) path were exercised


def malformed_chirals(chirals: list, bonds: list) -> list:
    """List the chiral centres whose definition names an atom that is not bonded to the centre."""
    # A definition is [target, centre, C, D, E]; mstool restrains and reviews the sign of (target-centre).(CD x DE),
    # which only means "chirality" when all four atoms are substituents of the centre.
    bonded = {frozenset((str(a), str(b))) for a, b in bonds}
    return [c[1] for c in chirals if not all(frozenset((c[1], other)) in bonded for other in (c[0], c[2], c[3], c[4]))]


def flipped_chirals(atoms, chirals: list, bad_centres: list) -> dict:
    """Count, for one residue type, the centres whose handedness is opposite to its (well-formed) definition."""
    # atoms: mstool atom table of that residue type; same sign convention as mstool: (target-centre).(CD x DE) < 0.
    import numpy
    flips = {}
    for target, centre, c, d, e in (definition[:5] for definition in chirals):
        if centre in bad_centres:
            continue
        xyz = {name: atoms[atoms.name == name][["x", "y", "z"]].to_numpy() for name in (target, centre, c, d, e)}
        if len({len(v) for v in xyz.values()}) != 1:
            raise SystemExit(f"chirality review: atoms of centre {centre} are missing in some residues")
        sign = ((xyz[target] - xyz[centre]) * numpy.cross(xyz[d] - xyz[c], xyz[e] - xyz[d])).sum(axis=1)
        if (sign < 0).any():
            flips[centre] = int((sign < 0).sum())
    return flips


def describe(job: dict, mstool, mapping_add: list) -> dict:
    """Return the mstool version, location, residue mapping, and any chirality definitions that cannot be trusted."""
    maps = mstool.ReadMappings(mapping_add=mapping_add)
    ff, ff_add = forcefield_files(Path(job["data"]), Path(mstool.__file__).resolve().parent)
    xml = mstool.ReadXML(ff=ff, ff_add=ff_add)
    untrusted = {name: malformed_chirals(r["chiral"], xml.RESI[name]["bonds"])
                 for name, r in maps.RESI.items() if name in xml.RESI and r["chiral"]}
    return {"mstool_version": getattr(mstool, "__version__", "unknown"), "mstool_path": str(Path(mstool.__file__).parent),
            "mapping_files": mapping_add, "malformed_chirals": {name: bad for name, bad in untrusted.items() if bad},
            "residues": {name: {"beads": list(r["CGAtoms"]), "atoms": list(r["AAAtoms"])} for name, r in maps.RESI.items()}}


def forcefield_files(data: Path, package: Path) -> tuple[list, list]:
    """List the OpenMM force-field XML files: data-directory copies take precedence over the packaged ones."""
    local = lambda name, *alternatives: next((str(data / n) for n in alternatives + (name,) if (data / n).is_file()),
                                             str(package / "FF" / "charmm36" / name))
    ff = [local("charmm36.xml", "charmm36_local.xml"), local("pip.xml"), local("water.xml"), local("chyo.xml")]
    used = {Path(f).name for f in ff} | {"charmm36.xml"}
    return ff, [str(p) for p in sorted(data.glob("*.xml")) if p.name not in used]


def backmap(job: dict, mstool, mapping_add: list) -> None:
    """Backmap the membrane beads with the complex held as a rock and write the all-atom lipids as a table."""
    import random

    import numpy
    from openmm.app import ForceField
    random.seed(job["seed"])  # mstool places atoms around each bead at random; the seed makes an attempt repeatable
    numpy.random.seed(job["seed"])
    ff, ff_add = forcefield_files(Path(job["data"]), Path(mstool.__file__).resolve().parent)
    # Residues are matched to force-field templates by name, so look-alike templates (SAPI24/SAPI25) cannot be swapped.
    create = ForceField.createSystem

    @functools.wraps(create)
    def create_by_name(self, topology, *args, **kwargs):
        """Create the OpenMM system with every residue pinned to the template of its own name."""
        chosen = dict(kwargs.get("residueTemplates") or {})
        chosen.update({r: r.name for r in topology.residues() if r.name in self._templates and not r.name.startswith("ROCK")})
        return create(self, topology, *args, **dict(kwargs, residueTemplates=chosen))

    ForceField.createSystem = create_by_name
    # sanitizeMartini applies only our alias table: names longer than the 5-character .gro field travel as aliases.
    mstool.Backmap(structure=job["structure"], workdir=job["workdir"], rock=job["rock"], mapping_add=mapping_add,
                   ff=ff, ff_add=ff_add, use_AA_structure=True, AA_shrink_factor=0.8, pbc=True, turn_off_EMNVT=True,
                   nsteps=job["nsteps"], sanitizeMartini=True, changename=job["rename"])
    atoms = mstool.Universe(str(Path(job["workdir"]) / "step4_nonprotein.dms")).atoms
    with open(job["result"] + ".part", "w") as fh:
        fh.write("chain\tresid\tresname\tname\tx\ty\tz\n")
        for a in atoms.itertuples():
            fh.write(f"{a.chain}\t{a.resid}\t{a.resname}\t{a.name}\t{a.x:.3f}\t{a.y:.3f}\t{a.z:.3f}\n")
    maps = mstool.ReadMappings(mapping_add=mapping_add)
    xml = mstool.ReadXML(ff=ff, ff_add=ff_add)
    review = {}
    for name in sorted(set(atoms.resname)):
        chirals = maps.RESI.get(name, {}).get("chiral", [])
        if chirals and name in xml.RESI:
            flips = flipped_chirals(atoms[atoms.resname == name], chirals, malformed_chirals(chirals, xml.RESI[name]["bonds"]))
            if flips:
                review[name] = flips
    Path(job["result"] + ".chirality.json").write_text(json.dumps(review))
    os.replace(job["result"] + ".part", job["result"])


def main(argv: list) -> int:
    """Run one job file: mode 'describe' or 'backmap'."""
    job = json.loads(Path(argv[1]).read_text())
    import mstool
    data = Path(job["data"])
    mapping_add = [str(data / "map.dat")] if (data / "map.dat").is_file() else []
    if job["mode"] == "describe":
        Path(job["result"]).write_text(json.dumps(describe(job, mstool, mapping_add)))
    elif job["mode"] == "backmap":
        backmap(job, mstool, mapping_add)
    else:
        raise SystemExit(f"unknown worker mode {job['mode']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
