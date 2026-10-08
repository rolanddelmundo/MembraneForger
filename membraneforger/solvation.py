#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Box construction, solvation, ions and index groups.
#//=============================================================
"""Resize the box in z, add water and 0.15 M NaCl, and write the index groups."""
import math
import shutil
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .config import (
    EM_MDP,
    FMAX_TARGET,
    GENION_ATTEMPTS,
    INDEX_GROUPS,
    ION_RMIN_NM,
    MIN_Z_PAD_TOTAL_NM,
    SALT_M,
    SLAB_Z_PAD_NM,
    WATER_CLASH_NM,
    WATER_PROTECT_NM,
)
from .runtools import log, run_command
from .structio import element, read_gro, wrap, write_gro, xyz_nm

__all__ = ['rebox_system', 'solute_geometry', 'solvate_system', 'add_ions', 'make_index']

def rebox_system(system: dict, topology: dict, requested_z_nm: float | None = None, z_pad_nm: float = SLAB_Z_PAD_NM) -> list[float]:
    """Keep the membrane XY cell and size z around the bilayer midplane so water covers everything above and below."""
    # x and y are the membrane cell (the lipids are periodic in it; --box x and y were already applied by slicing). z is
    # centred on the bilayer midplane, found halfway between the two lipid phosphate planes (robust against a single
    # stray lipid or tall head groups, unlike the extent of all membrane atoms), and reaches z_pad_nm beyond the
    # farthest solute or membrane heavy atom on either side, but never less than MIN_Z_PAD_TOTAL_NM over the full
    # extent of all atoms. A requested z (from --box) is used when it is at least that tall, else refused.
    coords, cryst = topology["coords"], system["cryst1"].ljust(80)
    if not all(abs(float(cryst[i:i + 7]) - 90.0) < 0.02 for i in (33, 40, 47)):
        raise SystemExit("input CRYST1 is not orthorhombic")
    box_x, box_y = float(cryst[6:15]) / 10.0, float(cryst[15:24]) / 10.0
    heavy = [a for a in coords if element(a["atom"]) != "H"]
    memb_z = np.array([a["z"] for a in heavy if a["group"] == "MEMB"]) / 10.0
    solute_z = np.array([a["z"] for a in heavy if a["group"] == "Protein_LIG"]) / 10.0
    if not (len(memb_z) and len(solute_z)):
        raise SystemExit("rebox needs membrane lipids and a protein")
    phosphorus = np.array([a["z"] for a in heavy if a["group"] == "MEMB" and a["atom"] == "P"]) / 10.0
    middle = 0.5 * (memb_z.min() + memb_z.max())
    if len(phosphorus) >= 10 and (phosphorus < middle).any() and (phosphorus >= middle).any():
        planes = [float(np.median(phosphorus[phosphorus < middle])), float(np.median(phosphorus[phosphorus >= middle]))]
        center, centre_source = 0.5 * (planes[0] + planes[1]), "midpoint between the lipid phosphate planes"
    else:
        planes, center, centre_source = None, float(middle), "midpoint of the membrane heavy-atom extent (no phosphate planes)"
    top, bottom = float(solute_z.max()), float(solute_z.min())
    half = max(top - center, center - bottom, float(memb_z.max()) - center, center - float(memb_z.min())) + z_pad_nm
    all_z = np.array([a["z"] for a in coords]) / 10.0
    minimum = max(2.0 * half, float(all_z.max() - all_z.min()) + MIN_Z_PAD_TOTAL_NM)
    box_z, mode = math.ceil(minimum * 1000.0) / 1000.0, "auto"
    if requested_z_nm is not None:
        if requested_z_nm < box_z - 1e-6:
            raise SystemExit(f"requested box z {requested_z_nm * 10:.0f} A is below the {box_z * 10:.0f} A this system needs "
                             f"({z_pad_nm} nm of water above and below the complex and membrane)")
        box_z, mode = float(requested_z_nm), "user"
    shift = box_z / 2.0 - center
    if (all_z + shift).min() <= 0.0 or (all_z + shift).max() >= box_z:
        raise SystemExit("solute extends beyond the rebuilt z box")
    atoms = [dict(a, x=a["x"] / 10.0, y=a["y"] / 10.0, z=a["z"] / 10.0 + shift) for a in coords]
    write_gro(atoms, [box_x, box_y, box_z], system["out"] / "boxed.gro", f"{system['name']} boxed")
    system["box_report"] = {"mode": mode, "box_nm": [box_x, box_y, box_z], "xy_source": "membrane cell (CRYST1)",
                            "bilayer_midplane_nm": round(center, 4), "midplane_source": centre_source,
                            "phosphate_planes_nm": [round(p, 4) for p in planes] if planes else None,
                            "solute_above_midplane_nm": round(top - center, 4), "solute_below_midplane_nm": round(center - bottom, 4),
                            "membrane_above_midplane_nm": round(float(memb_z.max()) - center, 4),
                            "membrane_below_midplane_nm": round(center - float(memb_z.min()), 4),
                            "z_pad_nm": z_pad_nm, "minimum_z_nm": round(minimum, 4), "shift_nm": round(shift, 4),
                            "requested_z_nm": requested_z_nm}
    log(system["out"], f"rebox ({mode}): {box_x:.3f} x {box_y:.3f} x {box_z:.3f} nm; bilayer midplane z {center:.2f} nm "
                       f"({centre_source}), solute top {top - center:+.2f} nm, solute bottom {bottom - center:+.2f} nm, "
                       f"pad {z_pad_nm} nm, minimum z {minimum:.3f} nm")
    return [box_x, box_y, box_z]


def solute_geometry(topology: dict, solute: list[dict], box: list[float]) -> dict:
    """Build periodic lookup trees for the solute and find the two phosphate planes."""
    coords = topology["coords"]
    xyz = wrap(xyz_nm(solute), box)
    heavy = np.array([element(a["atom"]) != "H" for a in coords])
    protected = heavy & np.array([a["group"] == "Protein_LIG" for a in coords])
    z = xyz[np.array([a["group"] == "MEMB" and a["atom"] == "P" for a in coords])][:, 2]
    if not len(z):
        raise SystemExit("the membrane has no phospholipid (no lipid P atom) to define the bilayer planes")
    mid = float(np.median(z))
    if not (z < mid).any() or not (z >= mid).any():
        raise SystemExit("cannot split lipid phosphates into two leaflets")
    return {"heavy": cKDTree(xyz[heavy], boxsize=box), "protected": cKDTree(xyz[protected], boxsize=box),
            "planes": (float(np.median(z[z < mid])), float(np.median(z[z >= mid])))}


def solvate_system(system: dict, topology: dict, gmx: str, box: list[float]) -> None:
    """Solvate, then drop waters inside the bilayer core (except cavity waters) or touching the solute."""
    out = system["out"]
    run_command(out, gmx.split() + ["solvate", "-cp", "boxed.gro", "-cs", "spc216.gro", "-o", "solv_raw.gro"],
            produces=("solv_raw.gro",))
    atoms, solv_box = read_gro(out / "solv_raw.gro")
    if not np.allclose(solv_box[:3], box, atol=1e-3):
        raise SystemExit(f"gmx solvate changed the box: {solv_box}")
    nsolute = len(topology["names"])
    solute, water = atoms[:nsolute], atoms[nsolute:]
    if [a["atom"] for a in solute] != topology["names"] or [a["atom"] for a in water] != ["OW", "HW1", "HW2"] * (len(water) // 3):
        raise SystemExit("unexpected atom order after gmx solvate")
    geo = solute_geometry(topology, solute, box)
    oxygen = wrap(xyz_nm(water[0::3]), box)
    lower, upper = geo["planes"]
    in_core = (oxygen[:, 2] >= lower) & (oxygen[:, 2] <= upper)
    cavity = np.array([bool(hits) for hits in geo["protected"].query_ball_point(oxygen, WATER_PROTECT_NM)])
    clash = np.array([bool(hits) for hits in geo["heavy"].query_ball_point(oxygen, WATER_CLASH_NM)])
    keep = ~((in_core & ~cavity) | clash)
    kept = [dict(a, resname="TIP3", atom=name, resid=nsolute + i + 1)
            for i in np.flatnonzero(keep) for a, name in zip(water[3 * i:3 * i + 3], ("OH2", "H1", "H2"))]
    write_gro(solute + kept, box, out / "solv.gro", f"{system['name']} solvated")
    (out / "solv_raw.gro").unlink()
    nwater = len(kept) // 3
    topology["molecules"].append(("TIP3", nwater, topology["solvent_itps"]["TIP3"]["atoms"], "SOL_ION"))
    topology["names"] += ["OH2", "H1", "H2"] * nwater
    with (out / "topol.top").open("a") as fh:
        fh.write(f"{'TIP3':<12s} {nwater}\n")
    log(out, f"solvate: {len(oxygen)} waters added, {int((in_core & ~cavity).sum())} removed from the bilayer core "
             f"(phosphate planes {lower:.2f}/{upper:.2f} nm), {int(clash.sum())} within {WATER_CLASH_NM} nm of solute; "
             f"{nwater} kept ({int((in_core & cavity & keep).sum())} cavity waters near protein/ligand)")


def add_ions(system: dict, topology: dict, gmx: str, box: list[float]) -> None:
    """Add 0.15 M NaCl plus counterions, retrying genion until no ion lands in the bilayer or on the solute."""
    out, charge = system["out"], round(topology["charge"])
    if abs(topology["charge"] - charge) > 1e-3:
        raise SystemExit(f"non-integer system charge {topology['charge']:.4f}")
    pairs = int(math.floor(SALT_M * box[0] * box[1] * box[2] * 0.602214076 + 0.5))
    n_pos, n_neg = pairs + max(0, -charge), pairs + max(0, charge)
    (out / "ions.mdp").write_text(EM_MDP.format(emtol=FMAX_TARGET, nsteps=0, coulombtype="Cut-off", define=""))
    run_command(out, gmx.split() + ["grompp", "-f", "ions.mdp", "-c", "solv.gro", "-p", "topol.top",
                                "-o", "ions.tpr", "-po", "ions_mdout.mdp", "-maxwarn", "0"], produces=("ions.tpr",))
    nsolute = len(topology["names"]) - 3 * topology["molecules"][-1][1]
    water_ids = range(nsolute + 1, len(topology["names"]) + 1)
    (out / "genion.ndx").write_text("[ TIP3 ]\n" + "".join(" ".join(map(str, water_ids[i:i + 15])) + "\n"
                                                           for i in range(0, len(water_ids), 15)))
    nwater = topology["molecules"][-1][1] - n_pos - n_neg
    expected = topology["names"][:nsolute] + ["OH2", "H1", "H2"] * nwater + ["SOD"] * n_pos + ["CLA"] * n_neg
    shutil.copy(out / "topol.top", out / "topol.pre_genion.top")
    for attempt in range(1, GENION_ATTEMPTS + 1):
        shutil.copy(out / "topol.pre_genion.top", out / "topol.top")
        run_command(out, gmx.split() + ["genion", "-s", "ions.tpr", "-n", "genion.ndx", "-o", "solv_ions.gro", "-p", "topol.top",
                                    "-pname", "SOD", "-nname", "CLA", "-np", str(n_pos), "-nn", str(n_neg),
                                    "-rmin", str(ION_RMIN_NM), "-seed", str(1009 + 7919 * attempt)], stdin="TIP3\n",
                produces=("solv_ions.gro",))
        atoms, _ = read_gro(out / "solv_ions.gro")
        if [a["atom"] for a in atoms] != expected:
            raise SystemExit("unexpected atom order after gmx genion")
        geo = solute_geometry(topology, atoms[:nsolute], box)
        ions = wrap(xyz_nm(atoms[len(expected) - n_pos - n_neg:]), box)
        embedded = int(((ions[:, 2] >= geo["planes"][0]) & (ions[:, 2] <= geo["planes"][1])).sum())
        clashing = sum(bool(h) for h in geo["heavy"].query_ball_point(ions, WATER_CLASH_NM))
        if not embedded and not clashing:
            break
        log(out, f"genion attempt {attempt}: {embedded} ions inside the bilayer, {clashing} touching solute; retrying", "WARN")
    else:
        raise SystemExit(f"no valid ion placement after {GENION_ATTEMPTS} genion attempts")
    for stale in ("ions.tpr", "ions_mdout.mdp", "genion.ndx", "topol.pre_genion.top"):
        (out / stale).unlink(missing_ok=True)
    mol, _, water_atoms, group = topology["molecules"].pop()
    topology["molecules"] += [(mol, nwater, water_atoms, group),
                              ("SOD", n_pos, topology["solvent_itps"]["SOD"]["atoms"], "SOL_ION"),
                              ("CLA", n_neg, topology["solvent_itps"]["CLA"]["atoms"], "SOL_ION")]
    topology["names"] = expected
    topology["charge"] = sum(n * sum(a["charge"] for a in atoms) for _, n, atoms, _ in topology["molecules"])
    if abs(topology["charge"]) > 1e-3:
        raise SystemExit(f"system not neutral after genion: {topology['charge']:+.4f}")
    log(out, f"ions: {SALT_M} M NaCl in {box[0] * box[1] * box[2]:.1f} nm^3 -> {pairs} pairs; "
             f"{n_pos} SOD + {n_neg} CLA neutralize {charge:+d} e; {nwater} TIP3 remain")


def make_index(out: Path, topology: dict) -> dict:
    """Write index_ini.ndx with System, Protein_LIG, MEMB and SOL_ION."""
    groups, start = {name: [] for name in INDEX_GROUPS}, 1
    for _, count, atoms, group in topology["molecules"]:
        groups[group] += range(start, start + count * len(atoms))
        start += count * len(atoms)
    groups["System"] = list(range(1, start))
    with (out / "index_ini.ndx").open("w") as fh:
        for name in INDEX_GROUPS:
            ids = groups[name]
            fh.write(f"[ {name} ]\n" + "".join(" ".join(map(str, ids[i:i + 15])) + "\n" for i in range(0, len(ids), 15)) + "\n")
    return {name: len(ids) for name, ids in groups.items()}
