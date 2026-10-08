#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Build orchestration: stages, timing, fail-closed publication.
#//=============================================================
"""Run the two-input build stage by stage and publish em.gro only after validation and audit pass."""
import json
import os
import shutil
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from .alignment import map_all_atom_to_cg
from .audit import audit_run, check_ring_piercing, closest_contact_between
from .backmapping import assemble_membrane_pdb, backmap_membrane, read_mapping
from .config import AMINO, DEFAULT_LIGANDS, RENAME_MOLECULE, Settings
from .embedding import FRAME_PROTEIN_CLEARANCE_NM, check_box_z, edit_lipids, embed_complex, frame_protein_in_slice
from .martini import bilayer_midplane, classify_cg, make_membrane_whole
from .membrane_report import MembraneValidation, protein_of_cg_frame, stage_from_gro, stage_from_membrane_pdb
from .minimization import run_em, validate_em
from .orientation import OrientationRequest, check_orientation_preserved, orient_complex
from .reporting import write_run_manifest
from .runtools import LOG_NAME, gromacs_version, log, sha256
from .slicing import slice_membrane_cg
from .solvation import add_ions, make_index, rebox_system, solvate_system
from .stereo import check_gm3_stereo
from .structio import xyz_nm
from .topology import build_topology, lipid_clashes, prepare_structure, repair_structure, resolve_lipid_clashes
from .validation import align_frame_to_box, check_inputs_unchanged, clear_stale_outputs, protect_inputs, read_all_atom, read_cg

__all__ = ['Session', 'StageFailure', 'run_stage', 'prepare_inputs', 'backmap_and_assemble', 'heavy_atom_piercings', 'membrane_piercings',
           'backmap_verdict', 'refuse_frame_protein_in_slice',
           'build_topology_and_box', 'finish_system', 'build_from_two_inputs', 'build']

# What to look at when a stage fails; {out} and {work} are filled in at run time.
STAGE_HINTS = {
    "inputs": "the input file named in the message",
    "orient": "{out}/orientation_report.json, the chains of the all-atom input, and {work}/orientation/ppm/ with the PPM "
              "transcript in {out}/" + LOG_NAME,
    "slice": "the slice report in {out}/run_manifest.json and the log, the --box option and --xy-buffer",
    "slice_check": "{out}/membrane_validation.md and membrane_validation.json; --apl-tolerance and --apl-validate",
    "membrane_check": "{out}/membrane_validation.md and membrane_validation.json",
    "mstool": "the --mstool-python interpreter and the transcript at the end of {out}/" + LOG_NAME,
    "classify": "the residue and bead names of the coarse-grained input listed in the message",
    "mapping": "{out}/aa_cg_mapping.tsv if written, and the chain sequences of both inputs",
    "backmap": "the mstool transcript in {out}/" + LOG_NAME + " and {work}/mstool/",
    "assemble": "{work}/membrane_aa.tsv",
    "structure": "{out}/membrane.pdb",
    "topology": "the pdb2gmx transcript in {out}/" + LOG_NAME + " and {work}/topology_*/",
    "box": "{out}/prot-memb.pdb",
    "solvate": "{out}/boxed.gro and the gmx solvate transcript in {out}/" + LOG_NAME,
    "ions": "{out}/solv.gro and the gmx genion transcript in {out}/" + LOG_NAME,
    "index": "{out}/topol.top",
    "rings": "{out}/ring_piercing.json and the named atoms in {out}/boxed.gro",
    "minimize": "{out}/em.log and the grompp/mdrun transcript in {out}/" + LOG_NAME,
    "validate": "{out}/em.log and {out}/em.unverified.gro",
    "audit": "{out}/audit.json",
    "publish": "{out}/em.unverified.gro",
}


@dataclass
class Session:
    """Everything one build needs to know: where it runs, which tools it uses, and what it has recorded so far."""
    out: Path
    name: str
    gmx: str
    forcefield: Path
    data: Path | None
    python: str
    ntomp: int
    nsteps: int
    settings: Settings = field(default_factory=Settings)
    orientation: OrientationRequest = field(default_factory=lambda: OrientationRequest(mode="none"))
    embed: bool = False  # place the complex into the frame's membrane instead of fitting it onto the frame's protein
    embed_site: str = "hole"  # with embed: "hole" (where the frame's protein was) or "free" (farthest from it; embedding.free_site)
    box_a: tuple | None = None  # opt-in box edges (A): x and y cut the membrane, z sizes the water layers; None = BOX auto
    bilayer_z_a: float | None = None  # z of the bilayer centre in the all-atom input (A) when embedding without orientation
    delete_lipids: list = field(default_factory=list)
    add_lipid: str | None = None
    work: Path | None = None
    timings: dict = field(default_factory=dict)
    record: dict = field(default_factory=dict)
    analysis: object = None  # the MembraneValidation of this build (stage measurements; not written to the manifest)


class StageFailure(Exception):
    """A stage failed: carries the stage name, what failed and what to inspect next."""

    def __init__(self, stage: str, what: str, inspect: str):
        """Remember which stage failed, what failed and what to look at."""
        super().__init__(f"{stage}: {what}")
        self.stage, self.what, self.inspect = stage, what, inspect


def run_stage(session: Session, stage: str, function, *args, **kwargs):
    """Run one stage, time it, and turn any failure into a StageFailure naming the artifact to inspect."""
    start = time.time()
    if stage not in session.timings:
        log(session.out, f"stage {stage}", "CHECK")
    try:
        return function(*args, **kwargs)
    except StageFailure:
        raise
    except SystemExit as exc:
        raise StageFailure(stage, str(exc), STAGE_HINTS[stage].format(out=session.out, work=session.work)) from None
    except Exception as exc:
        raise StageFailure(stage, f"{type(exc).__name__}: {exc}",
                           STAGE_HINTS[stage].format(out=session.out, work=session.work)) from exc
    finally:
        session.timings[stage] = round(session.timings.get(stage, 0.0) + time.time() - start, 2)


def prepare_inputs(session: Session, all_atom: Path, coarse_grain: Path) -> dict:
    """Validate both inputs, classify the CG system and place the all-atom complex; everything backmapping needs."""
    out, work, record = session.out, session.work, session.record
    aa_atoms, notes = run_stage(session, "inputs", read_all_atom, all_atom, session.forcefield)
    for note in notes:
        log(out, note)
    # Orientation first: it needs only the all-atom input, and it must fail before anything expensive runs.
    (work / "orientation").mkdir(exist_ok=True)
    oriented = run_stage(session, "orient", orient_complex, aa_atoms, all_atom, session.orientation, session.settings,
                         out, work / "orientation")
    record["orientation"] = oriented["report"]
    write_orientation_report(out, oriented["report"])
    aa_atoms = oriented["oriented"]  # from here on the PPM/OPM frame is authoritative: nothing may tilt the complex
    cg_atoms, box = run_stage(session, "inputs", read_cg, coarse_grain)
    cg_atoms, frame_fix = run_stage(session, "inputs", align_frame_to_box, cg_atoms, box)
    if frame_fix["rotation_about_z_deg"]:
        log(out, f"{coarse_grain.name}: the frame's coordinates are rotated {frame_fix['rotation_about_z_deg']} degrees about z relative to its "
                 f"box ({frame_fix['overlapping_pairs_as_read']} overlapping bead pairs once wrapped, a rotational fit written without the box); "
                 f"rotated back about the box centre ({frame_fix['overlapping_pairs_after']} overlaps remain)", "WARN")
    described = run_stage(session, "mstool", read_mapping, out, work, session.python, session.data)
    mapping = described["residues"]
    record["mstool"] = {k: described[k] for k in ("mstool_version", "mstool_path", "mapping_files")}
    log(out, f"mstool {described['mstool_version']} at {described['mstool_path']}")
    cg = run_stage(session, "classify", classify_cg, cg_atoms, mapping)
    log(out, f"{coarse_grain.name}: {len(cg_atoms)} beads, box {box[0]:.3f} x {box[1]:.3f} x {box[2]:.3f} nm; Martini 3 "
             + ", ".join(f"{cls} {n}" for cls, n in sorted(cg["counts"].items())))
    log(out, "membrane composition: " + ", ".join(f"{name} {n}" for name, n in sorted(cg["composition"].items())))
    used = {m["aa"] for m in cg["membrane"]}
    record["mstool"]["malformed_chirals"] = {k: v for k, v in described.get("malformed_chirals", {}).items() if k in used}
    for name, centres in sorted(record["mstool"]["malformed_chirals"].items()):
        log(out, f"{name}: the mapping's chirality definitions for centre(s) {', '.join(centres)} name atoms that are not "
                 "bonded to the centre; stereochemistry there is neither restrained nor reviewed during backmapping", "WARN")
    for name, dropped in sorted(cg["dropped"].items()):
        log(out, f"{name}: bead(s) {', '.join(dropped)} carry no atoms in the all-atom mapping and are not used", "WARN")
    record["coarse_grain"] = {"beads": len(cg_atoms), "box_nm": box, "classes": dict(cg["counts"]), "frame_alignment": frame_fix,
                              "membrane_composition": dict(cg["composition"]), "unused_beads": cg["dropped"],
                              "discarded": "CG solvent and ions are not backmapped; atomistic water and 0.15 M NaCl "
                                           "are rebuilt after the box is resized"}
    membrane = cg["membrane"]
    if session.delete_lipids or session.add_lipid:
        membrane, edits = run_stage(session, "classify", edit_lipids, membrane, session.delete_lipids, session.add_lipid, mapping)
        record["lipid_edits"] = edits
        log(out, "lipid edits: " + "; ".join(f"{what} {', '.join(f'{n} {k}' for n, k in sorted(v.items())) or 'none'}"
                                             for what, v in edits.items()), "WARN")
    run_stage(session, "classify", require_membrane_topologies, membrane, session.forcefield)
    slab = run_stage(session, "classify", make_membrane_whole, membrane, box)
    midplane, midplane_how, leaflets = run_stage(session, "classify", bilayer_midplane, membrane, slab)
    record["coarse_grain"].update({"membrane_z_range_nm": [round(slab[0], 4), round(slab[1], 4)],
                                   "bilayer_midplane_nm": round(midplane, 4), "bilayer_midplane_source": midplane_how, **leaflets})
    log(out, f"CG bilayer midplane z {midplane:.3f} nm ({midplane_how})"
             + (f"; PO4 normal {leaflets['po4_normal_tilt_from_z_deg']} deg from z" if leaflets else ""))
    session.analysis = MembraneValidation(out, session.settings)
    run_stage(session, "membrane_check", session.analysis.add_cg, "cg_frame", membrane, [], box, midplane, protein_of_cg_frame(cg["protein"]))
    enabled, anchors = oriented["report"]["enabled"], tuple(oriented["report"]["anchor_chains"])
    if session.embed:
        # An oriented complex has its bilayer centre at z = 0 and goes straight in; the midplane is the CG bilayer's own.
        # bilayer centre of the input: 0 for an oriented complex; otherwise --bilayer-z, else the midpoint of the segment
        # --orient-residues names (orientation.embedded_segment), else the hydrophobic-belt search (None)
        centre_a = 0.0 if enabled else (session.bilayer_z_a if session.bilayer_z_a is not None
                                        else oriented["report"].get("segment_bilayer_centre_A"))
        settings = session.settings
        fit = run_stage(session, "mapping", embed_complex, aa_atoms, cg["protein"], membrane, box, slab,
                        centre_a, midplane if enabled else None, site=session.embed_site, window_buffer_nm=settings.box_xy_buffer_nm,
                        window_slack_nm=settings.slice_offset_search_nm if settings.slice_optimize_offset else 0.0,
                        requested_xy_nm=(session.box_a[0] / 10.0, session.box_a[1] / 10.0) if session.box_a else None)
        membrane = fit["membrane"]
        for note in fit["notes"]:
            log(out, f"embed: {note}", "WARN" if "removed" in note else "INFO")
        record["embedding"] = fit["metrics"]
        log(out, f"complex embedded: {fit['metrics']['membrane_molecules_kept']} membrane molecules kept around it", "PASS")
    else:
        fit = run_stage(session, "mapping", map_all_atom_to_cg, aa_atoms, cg["protein"], box, slab, session.settings,
                        midplane if enabled else None, anchors, oriented["report"].get("half_thickness_A"))
        for note in fit["notes"]:
            log(out, f"AA-CG: {note}", "WARN" if "keeps its all-atom pose" in note else "INFO")
        (out / "aa_cg_mapping.tsv").write_text("aa_chain\taa_residue\tcg_segment\tcg_residue\tdeviation_A\tfit\n"
                                               + "".join("\t".join(map(str, row)) + "\n" for row in fit["rows"]))
        record["aa_cg_mapping"] = fit["metrics"]
        log(out, f"AA-CG fit accepted: core RMSD {fit['metrics']['core_rmsd_A']} A <= {session.settings.fit_max_core_rmsd_a} A, "
                 f"core fraction {fit['metrics']['core_fraction']:.2f} >= {session.settings.fit_min_core_fraction}", "PASS")
    moved = xyz_nm(aa_atoms) @ fit["R"].T + fit["t"]
    placed = [dict(a, x=float(x), y=float(y), z=float(z)) for a, (x, y, z) in zip(aa_atoms, moved)]
    if enabled:  # prove the placement kept the orientation: normal still z, depth unchanged
        preserved = run_stage(session, "mapping", check_orientation_preserved, aa_atoms, placed, list(anchors), midplane * 10.0)
        registration = fit["metrics"].get("registration") or {"mode": "embed", "translation_A": fit["metrics"]["translation_A"]}
        oriented["report"]["registration"] = {**registration, "cg_midplane_source": midplane_how,
                                              **{f"cg_{k}": v for k, v in leaflets.items()}, "check": preserved,
                                              "final_membrane_normal": preserved["membrane_normal_after_registration"],
                                              "final_membrane_center_nm": [None, None, round(midplane, 4)]}
        write_orientation_report(out, oriented["report"])
        measured = registration.get("tilt_between_cg_anchor_pose_and_orientation_deg")
        log(out, f"orientation preserved through placement: normal tilt {preserved['normal_tilt_deg']:.1e} deg, anchor depth "
                 f"{preserved['anchor_ca_depth_A']:+.2f} A kept"
                 + (f"; CG pose differs by {measured} deg (measured, not adopted)" if measured is not None else ""), "PASS")
    embedded = run_stage(session, "membrane_check", session.analysis.add_cg, "embed", membrane, placed, box, midplane)
    for leaflet in ("upper", "lower"):
        e = embedded["leaflets"][leaflet]
        log(out, f"embedded membrane, {leaflet} leaflet: {e['lipids']} lipids, APL {e['apl_A2']} A^2 (Voronoi; protein footprint "
                 f"{e['protein_area_A2']:.0f} A^2); this is the reference the slice must preserve")
    requested = (session.box_a[0] / 10.0, session.box_a[1] / 10.0) if session.box_a else None
    cut = run_stage(session, "slice", slice_membrane_cg, membrane, placed, box, midplane, session.settings, requested,
                    session.analysis.slice_reference())
    record["slice"] = cut["report"]
    if session.embed and session.embed_site == "free":
        reached = run_stage(session, "slice", refuse_frame_protein_in_slice, cg["protein"], slab, box, placed, cut)
        record["slice"]["frame_protein_beads_in_window"] = reached
        log(out, f"slice window keeps {FRAME_PROTEIN_CLEARANCE_NM} nm from every bead of the frame's own protein: the hole it "
                 "leaves stays outside the cut membrane", "PASS")
    report, seam = cut["report"], cut["report"]["seam"]
    if requested:
        record["box_trim"] = {"requested_A": list(session.box_a), "see": "slice"}
    if report["cropped"]:
        log(out, f"slice ({report['mode']}): cell {box[0]:.2f} x {box[1]:.2f} -> {cut['box'][0]:.2f} x {cut['box'][1]:.2f} nm "
                 f"({report['area_fraction_kept']:.0%} of the area), {report['lipids_before']} -> {report['lipids_after']} lipids "
                 f"kept whole; leaflets {sum(report['leaflets_after']['lower'].values())}/"
                 f"{sum(report['leaflets_after']['upper'].values())}; seam: {seam['pairs_created_by_new_periodicity']} bead pairs "
                 f"closer than {seam['threshold_nm']} nm created by the new periodicity, {seam['lipids_removed']} lipids removed; "
                 "the cut edges were never equilibrated together, so equilibrate before production",
            "WARN" if seam["lipids_removed"] else "INFO")
    else:
        log(out, f"slice ({report['mode']}): the complex plus {session.settings.box_xy_buffer_nm} nm reaches the whole "
                 f"{box[0]:.2f} x {box[1]:.2f} nm cell; the membrane is not cut")
    selection = report["selection"]
    log(out, f"slice selection: {selection['kept_by_anchor_rule']} lipids by their headgroup anchor; the earlier all-beads-inside rule "
             f"would have kept {selection['kept_by_all_beads_inside_rule']} ({selection['rejected_by_old_rule_with_anchor_inside']} lipids "
             f"with the anchor inside lost to a tail bead across the edge: {selection['rejected_by_old_rule_with_anchor_inside_by_leaflet']})")
    if report["crop_offset"]["optimized"]:
        log(out, f"slice window offset {report['crop_offset']['offset_nm']} nm chosen from {report['crop_offset']['candidates']} candidates "
                 f"(score {report['crop_offset']['chosen']['score']} vs {report['crop_offset']['unshifted']['score']} centred)")
    run_stage(session, "membrane_check", session.analysis.add_cg, "slice", cut["membrane"], cut["placed"], cut["box"], midplane)
    try:
        record["slice_check"] = run_stage(session, "slice_check", session.analysis.validate_slice, cut, membrane, placed, box)
    finally:
        record["membrane_validation"] = session.analysis.record(cut)
        session.analysis.cut = cut
        session.analysis.write(cut)
    failed = ", ".join(k for k, v in record["slice_check"]["verdict"].items() if not v)
    log(out, f"slice gate: {'PASS' if record['slice_check']['pass'] else 'FAIL'} ({failed or 'all checks'}); membrane_validation.md",
        "PASS" if record["slice_check"]["pass"] else "WARN")
    placed, membrane, box = cut["placed"], cut["membrane"], cut["box"]
    if session.box_a:  # refuse a z that cannot hold the complex now, not after the slow backmapping
        run_stage(session, "mapping", check_box_z, placed, slab, session.box_a[2] / 10.0)
    return {"placed": placed, "membrane": membrane, "box": box, "mapping": mapping,
            "composition": Counter(RENAME_MOLECULE.get(m["aa"], m["aa"]) for m in membrane),
            "ligands": DEFAULT_LIGANDS | {a["resname"] for a in aa_atoms if a["resname"] not in AMINO}}


def refuse_frame_protein_in_slice(cg_protein: list, slab: tuple, box: list, placed: list, cut: dict) -> int:
    """With --embed-site free the cut membrane must not reach the frame's own protein, whose hole would come with it."""
    reached = frame_protein_in_slice(cg_protein, slab, box, placed, cut)
    if reached:
        raise SystemExit(f"the slice window ({cut['box'][0]:.2f} x {cut['box'][1]:.2f} nm) comes within "
                         f"{FRAME_PROTEIN_CLEARANCE_NM} nm of {reached} bead(s) of the frame's own protein, so the hole it leaves "
                         "would reach the cut membrane; this complex is too wide to sit clear of it "
                         "in this frame: use a smaller --xy-buffer or --box, or --embed-site hole")
    return reached


def backmap_and_assemble(session: Session, prepared: dict, seed: int) -> Path:
    """Backmap the membrane around the placed complex with one random seed and write membrane.pdb."""
    out, work = session.out, session.work
    lipids, isomers = run_stage(session, "backmap", backmap_membrane, prepared["membrane"], prepared["placed"],
                                prepared["box"], prepared["mapping"], out, work, session.python, session.data,
                                session.ntomp, session.settings.backmap_em_steps, seed)
    session.record["backmap_isomer_review"] = isomers
    counts = {k: v for k, v in isomers.items() if isinstance(v, int)}
    if any(counts.values()):
        log(out, "isomer review of the backmapped membrane: " + ", ".join(f"{k} {v}" for k, v in counts.items())
                 + f" (mstool list in {work}/mstool/log.txt; 'chiral' includes centres with malformed definitions)", "WARN")
    natoms = run_stage(session, "assemble", assemble_membrane_pdb, prepared["placed"], lipids, prepared["box"], out / "membrane.pdb")
    log(out, f"membrane.pdb (seed {seed}): {len(prepared['placed'])} protein/ligand atoms placed, {len(lipids)} membrane "
             f"molecules backmapped ({natoms} atoms); inventory matches the coarse-grained input", "PASS")
    # GM3's sugar stereocentres are not covered by its DIHRES rows, so they are tested geometrically here, before any
    # minimization, and a wrong one counts as a wrong configuration: the verdict then tries another seed.
    gm3, _ = run_stage(session, "backmap", check_gm3_stereo, out / "membrane.pdb", session.data)
    isomers["gm3_sugar_centres_or_cis_bonds"] = sum(gm3["wrong_by_centre"].values()) + sum(gm3["cis_by_bond"].values())
    isomers["gm3_stereochemistry"] = {k: gm3[k] for k in ("molecules", "wrong_by_centre", "cis_by_bond", "affected_molecules")}
    if isomers["gm3_sugar_centres_or_cis_bonds"]:
        log(out, f"GM3 stereochemistry of membrane.pdb: {len(gm3['affected_molecules'])} of {gm3['molecules']} molecules have an inverted "
                 f"sugar centre or a cis ceramide bond ({isomers['gm3_sugar_centres_or_cis_bonds']} in all); counted as wrong configurations",
            "WARN")
    elif gm3["molecules"]:
        log(out, f"GM3 stereochemistry of membrane.pdb: all 16 sugar stereocentres and both ceramide trans bonds correct in every one "
                 f"of the {gm3['molecules']} molecules", "PASS")
    if session.analysis is not None:
        summary = run_stage(session, "membrane_check", session.analysis.add, stage_from_membrane_pdb(out / "membrane.pdb"))
        log(out, "backmapped membrane: " + leaflet_apl_summary(summary) + " (cross-resolution anchors P/O3/NF; reported, not graded)")
        session.record["membrane_validation"] = session.analysis.record(session.analysis.cut)
        session.analysis.write(session.analysis.cut)
    return out / "membrane.pdb"


def leaflet_apl_summary(summary: dict) -> str:
    """One log phrase per leaflet: lipids, Voronoi APL and the change against the embedded reference."""
    return "; ".join(f"{l} leaflet {summary['leaflets'][l]['lipids']} lipids, APL {summary['leaflets'][l]['apl_A2']} A^2 "
                     f"({summary['vs_reference']['apl'][l]['delta_apl_percent']:+.1f} % vs embedded)" for l in ("upper", "lower"))


def write_orientation_report(out: Path, report: dict) -> None:
    """Write orientation_report.json (rewritten once the CG registration has been checked)."""
    (out / "orientation_report.json").write_text(json.dumps(report, indent=2, default=str) + "\n")


def require_membrane_topologies(membrane: list[dict], forcefield: Path) -> None:
    """Fail before backmapping if any membrane residue has no GROMACS topology to be built with."""
    missing = sorted({RENAME_MOLECULE.get(m["aa"], m["aa"]) for m in membrane}
                     - {p.stem for p in (forcefield / "toppar").glob("*.itp")})
    if missing:
        raise SystemExit(f"no GROMACS topology in {forcefield / 'toppar'} for membrane residue(s) {missing}")


def check_membrane_inventory(system: dict, expected: Counter | None) -> None:
    """Fail if the assembled structure does not hold exactly the membrane molecules the CG input had."""
    found = Counter(key[0] for key in system["others"] if key[0] not in system["ligands"])
    if expected is not None and found != expected:
        raise SystemExit(f"membrane.pdb lost membrane molecules: {dict(found)} != coarse-grained {dict(expected)}")


def heavy_atom_piercings(report: dict) -> list:
    """Pick the ring threadings made by a bond between two heavy atoms, which minimization cannot undo."""
    hydrogen = lambda label: label.split(":")[1].startswith("H")
    return [p for p in report["pierced"] if not any(hydrogen(label) for label in p["bond"])]


def membrane_piercings(pierced: list, solute_atoms: int) -> list:
    """Pick the ring threadings that involve a membrane atom; only those can change when the membrane is backmapped again."""
    return [p for p in pierced if max(p["ring_atoms"] + p["bond_atoms"]) > solute_atoms]


def review_ring_piercing(out: Path, structure: str) -> dict:
    """List every bond that threads a 5- or 6-membered ring in a structure and save the list."""
    report, _fails = check_ring_piercing(out, structure)
    (out / "ring_piercing.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def build_topology_and_box(session: Session, membrane: Path, ligands: set, expected: Counter | None, hashes: dict,
                           attempt: int = 0) -> dict:
    """Take an assembled protein + membrane PDB through repair, topology and box, and review ring threading."""
    out, gmx = session.out, session.gmx
    work = session.work / f"topology_{attempt}"  # pdb2gmx scratch; one per backmapping attempt
    work.mkdir()
    system = run_stage(session, "structure", prepare_structure, membrane, out, session.name, session.forcefield, ligands)
    system["inputs"] = {**hashes, membrane: system["sha"]}
    log(out, f"{membrane.name}: {system['natoms']} atoms, chains {','.join(system['chains'])}, sha256 {system['sha'][:12]}; "
             f"dropped altloc B of {', '.join(system['alternates']) or 'none'}")
    run_stage(session, "structure", check_membrane_inventory, system, expected)
    for repair in run_stage(session, "structure", repair_structure, system):
        log(out, f"repair: {repair}", "WARN")
    clashes = run_stage(session, "structure", lipid_clashes, system)
    notes = run_stage(session, "structure", resolve_lipid_clashes, system, clashes)
    session.record["removed_lipids"] = notes[1:]
    log(out, notes[0])
    for note in notes[1:]:  # the published system then holds one lipid fewer than the coarse-grained input
        log(out, note, "WARN")
    topology = run_stage(session, "topology", build_topology, system, work, gmx)
    for note in topology["notes"]:
        log(out, note)
    box = run_stage(session, "box", rebox_system, system, topology, session.box_a[2] / 10.0 if session.box_a else None,
                    session.settings.box_z_pad_nm)
    session.record["box"] = system["box_report"]
    rings = run_stage(session, "rings", review_ring_piercing, out, "boxed.gro")
    solute_atoms = sum(n * len(atoms) for _, n, atoms, group in topology["molecules"] if group == "Protein_LIG")
    contact = run_stage(session, "rings", closest_contact_between, out, "boxed.gro", solute_atoms)
    return {"system": system, "topology": topology, "box": box, "rings": rings, "solute_atoms": solute_atoms, "contact": contact}


def finish_system(session: Session, staged: dict) -> str:
    """Solvate, add ions, minimize, validate and audit a boxed system, then publish em.gro."""
    out, gmx, record = session.out, session.gmx, session.record
    system, topology, box, rings = staged["system"], staged["topology"], staged["box"], staged["rings"]
    record["ring_piercing_before_em"] = rings
    record["closest_membrane_solute_contact_before_em"] = staged["contact"]
    if rings["pierced"]:  # heavy-atom membrane threadings never get here; EM can undo the rest, and the audit re-checks
        inside = len(rings["pierced"]) - len(membrane_piercings(rings["pierced"], staged["solute_atoms"]))
        log(out, f"{len(rings['pierced'])} bond(s) pass through a ring before EM ({inside} within the all-atom input "
                 f"itself; ring_piercing.json), e.g. {'-'.join(rings['pierced'][0]['bond'])} through the ring of "
                 f"{rings['pierced'][0]['ring'][0]}; the audit refuses the build if any ring is still threaded after EM", "WARN")
    else:
        log(out, f"no bond threads any of {rings['rings_checked']} rings before EM", "PASS")
    record["disulfides"] = {label: pairs for label, pairs in system["disulfides"].items() if pairs}
    run_stage(session, "solvate", solvate_system, system, topology, gmx, box)
    run_stage(session, "ions", add_ions, system, topology, gmx, box)
    index = run_stage(session, "index", make_index, out, topology)
    log(out, "index_ini.ndx: " + ", ".join(f"{k} {v}" for k, v in index.items()))
    record["box_nm"] = box
    record["topology"] = {"molecules": [[mol, n] for mol, n, _, _ in topology["molecules"]],
                          "atoms": len(topology["names"]), "net_charge": round(topology["charge"], 4), "index_groups": index}
    run_stage(session, "minimize", run_em, out, topology, gmx, session.ntomp, session.nsteps)
    record["em"] = run_stage(session, "validate", validate_em, system, topology, index)
    log(out, f"EM validated: {record['em']['summary']}", "PASS")
    if session.analysis is not None:
        summary = run_stage(session, "membrane_check", session.analysis.add, stage_from_gro(out / "em.unverified.gro"))
        log(out, "minimized membrane: " + leaflet_apl_summary(summary))
        record["membrane_validation"] = session.analysis.record(session.analysis.cut)
        session.analysis.write(session.analysis.cut)
    record["audit"] = run_stage(session, "audit", audit_run, out, gmx, "em.unverified.gro", session.settings, session.data)
    log(out, f"independent audit: {len(record['audit']['checks'])} checks passed", "PASS")
    chairs = {mol: (e["rings_out_of_chair"], e["rings"]) for mol, e in record["audit"]["dihedral_restraints"].items() if "rings" in e}
    if chairs:
        log(out, "sugar/inositol rings outside the chair their topology restrains (MD can repair a ring, slowly): " + ", ".join(
            f"{mol} {bad}/{total}" for mol, (bad, total) in chairs.items()), "WARN" if any(b for b, _ in chairs.values()) else "PASS")
    run_stage(session, "publish", check_inputs_unchanged, system["inputs"])
    os.replace(out / "em.unverified.gro", out / "em.gro")
    return record["em"]["summary"]


def backmap_verdict(staged: dict, review: dict, settings: Settings, last: bool) -> dict:
    """Decide whether one backmapped membrane can go on to solvation, and say what is wrong if it cannot."""
    # Only defects that involve the membrane count: re-backmapping cannot change the all-atom input itself.
    pierced = membrane_piercings(staged["rings"]["pierced"], staged["solute_atoms"])
    heavy = heavy_atom_piercings({"pierced": pierced})
    wrong = review.get("chiral_well_formed", 0) + review.get("cistrans", 0) + review.get("gm3_sugar_centres_or_cis_bonds", 0)
    contact = staged["contact"]
    counts = {"heavy_atom_ring_threadings": len(heavy), "hydrogen_ring_threadings": len(pierced) - len(heavy),
              "solute_internal_ring_threadings": len(staged["rings"]["pierced"]) - len(pierced),
              "wrong_stereocentres_or_double_bonds": wrong, "closest_membrane_solute_contact_nm": contact["distance_nm"]}
    problem, stage = None, "rings"
    if heavy:
        problem = (f"{'-'.join(heavy[0]['bond'])} threads the ring of {heavy[0]['ring'][0]} "
                   f"({len(heavy)} heavy-atom threading(s) involving the membrane; minimization cannot undo this)")
    elif wrong:
        problem, stage = f"{wrong} lipid stereocentre(s) or double bond(s) have the wrong configuration", "backmap"
    elif contact["distance_nm"] < settings.min_start_contact_nm:
        problem = (f"{contact['atoms'][0]} starts {10 * contact['distance_nm']:.2f} A from {contact['atoms'][1]} "
                   f"(limit {10 * settings.min_start_contact_nm:.2f} A; minimization diverges from such an overlap)")
    elif pierced and not last:  # an X-H bond poking through a ring: worth another seed, tolerable on the last one
        problem = f"{'-'.join(pierced[0]['bond'])} pokes through the ring of {pierced[0]['ring'][0]}"
    return {"accept": problem is None, "problem": problem, "stage": stage, "counts": counts}


def refuse_backmap(problem: str, attempts: int) -> None:
    """Stop the build when the last backmapping attempt still has a defect minimization cannot repair."""
    raise SystemExit(f"after {attempts} backmapping attempt(s): {problem}")


def refuse_threaded_rings(heavy: list, attempts: int) -> None:
    """Stop the build when a heavy-atom bond still threads a ring and nothing more can be tried."""
    if heavy:
        raise SystemExit(f"{len(heavy)} ring(s) threaded by a heavy-atom bond after {attempts} backmapping attempt(s), e.g. "
                         f"{'-'.join(heavy[0]['bond'])} through the ring of {heavy[0]['ring'][0]}; minimization cannot undo this")


def refuse_wrong_stereochemistry(wrong: int, attempts: int) -> None:
    """Stop the build when backmapped lipids still have a wrong stereocentre or double-bond configuration."""
    if wrong:
        raise SystemExit(f"{wrong} lipid stereocentre(s) or double bond(s) have the wrong configuration after {attempts} "
                         "backmapping attempt(s); minimization cannot correct a configuration")


def build_from_two_inputs(session: Session, all_atom: Path, coarse_grain: Path, hashes: dict) -> dict:
    """Backmap and build up to the boxed system, backmapping again with a new seed while the membrane has a defect."""
    out, attempts = session.out, []
    if session.settings.backmap_attempts < 1:
        raise StageFailure("inputs", "Settings.backmap_attempts must be at least 1", "the settings passed to the build")
    prepared = prepare_inputs(session, all_atom, coarse_grain)
    for seed in range(1, session.settings.backmap_attempts + 1):
        membrane = backmap_and_assemble(session, prepared, seed)
        staged = build_topology_and_box(session, membrane, prepared["ligands"], prepared["composition"], hashes, seed)
        review = session.record["backmap_isomer_review"]
        last = seed == session.settings.backmap_attempts
        verdict = backmap_verdict(staged, review, session.settings, last)
        attempts.append({"seed": seed, **verdict["counts"]})
        session.record["backmap_attempts"] = attempts
        if verdict["accept"]:
            return staged
        if last:
            run_stage(session, verdict["stage"], refuse_backmap, verdict["problem"], len(attempts))
        log(out, f"backmap seed {seed}: {verdict['problem']}; backmapping again with another seed", "WARN")


def build(session: Session, all_atom: Path | None, coarse_grain: Path | None, membrane: Path | None, command: list) -> int:
    """Build one system end to end; em.gro appears only after validation and audit pass. Returns the exit code."""
    out = session.out
    inputs = [p for p in (all_atom, coarse_grain, membrane) if p]
    out.mkdir(parents=True, exist_ok=True)
    status, error, hashes = "FAIL", None, {}
    started = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    try:
        protect_inputs(inputs, out, keep=membrane, resources=(session.forcefield, session.data))
    except SystemExit as exc:  # nothing has been touched yet: report without clearing or logging into the directory
        print(f"[{out.name}] ERROR: inputs: {exc}. Inspect: the --out path")
        print(f"{session.name}\tFAIL\t{exc}")
        return 1
    removed = clear_stale_outputs(out, keep=membrane)
    for extra in ("run_manifest.json", "audit.json"):
        (out / extra).unlink(missing_ok=True)
    (out / LOG_NAME).write_text(f"INFO: membraneforger {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    if removed:
        log(out, f"removed stale outputs from an earlier run: {', '.join(removed)}", "WARN")
    session.work = out / "work"  # inside the output directory so failed intermediates survive on any node
    session.work.mkdir()
    try:
        hashes = {p: sha256(p) for p in inputs}
        for p, digest in hashes.items():
            log(out, f"input {p} sha256 {digest[:12]}")
        session.record["gromacs_version"] = run_stage(session, "inputs", gromacs_version, session.gmx)
        if membrane is None:
            staged = build_from_two_inputs(session, all_atom, coarse_grain, hashes)
        else:
            staged = build_topology_and_box(session, membrane, set(DEFAULT_LIGANDS), None, hashes)
            run_stage(session, "rings", refuse_threaded_rings, heavy_atom_piercings(staged["rings"]), 0)
            if staged["contact"]["distance_nm"] < session.settings.min_start_contact_nm:
                contact = staged["contact"]
                run_stage(session, "rings", refuse_backmap, f"{contact['atoms'][0]} starts {10 * contact['distance_nm']:.2f} A "
                                                            f"from {contact['atoms'][1]}", 0)
        summary = finish_system(session, staged)
        status = "PASS"
    except StageFailure as failure:
        error = {"stage": failure.stage, "what": failure.what, "inspect": failure.inspect}
        (out / "em.gro").unlink(missing_ok=True)
        log(out, f"work files kept in {session.work}")
        log(out, f"{failure.stage}: {failure.what}. Inspect: {failure.inspect}", "ERROR")
        summary = f"{failure.stage}: {failure.what}"
    except BaseException:
        (out / "em.gro").unlink(missing_ok=True)
        raise
    if status == "PASS":
        shutil.rmtree(session.work)
        session.work = None
        log(out, summary, "PASS")  # the last log line; the manifest written next hashes the finished log
    write_run_manifest(session, status, error, started, hashes, command)
    print(f"{session.name}\t{status}\t{summary}")
    return int(status != "PASS")
