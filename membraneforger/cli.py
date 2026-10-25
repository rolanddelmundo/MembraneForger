#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Command line: python membraneforger.py --all-atom AA.pdb --coarse-grain CG.gro --out DIR
#//=============================================================
"""Parse the command line, locate installation resources, and run one build."""
import argparse
import os
import random
import re
import shutil
import sys
from dataclasses import replace
from pathlib import Path

from .backbone_check import parse_segments
from .config import BOX, NSTEPS, NTERM_SIDE, ORIENT_CHAINS, ORIENT_RESIDUES, ORIENTATION, PDB_ID, PPM_MEMBRANE, Settings
from .embedding import CONVERTIBLE, EMBED_SITES, LIPID_NAMES, lipid_name
from .orientation import NTERM_SIDES, ORIENTATION_MODES, OrientationRequest, parse_residue_ranges
from .pipeline import Session, build
from .runtools import find_gromacs
from .validation import refuse_overwrite

__all__ = ['BUNDLED_MEMBRANES', 'bundled_frames', 'bundled_membrane', 'bundled_series', 'describe_bundled_membranes',
           'make_parser', 'locate_forcefield', 'locate_data', 'parse_chains', 'parse_box', 'parse_disulfide_pair',
           'orientation_request', 'resolve_lengths', 'run_sessions', 'angstrom_to_nm', 'nm_to_angstrom', 'default_out', 'main']

PACKAGE_DIR = Path(__file__).resolve().parent
BUNDLED_DIR = PACKAGE_DIR / "examples" / "preeq_cg_cellmem"
# --cg 1 / --cg 2: one frame drawn at random from the bundled pre-equilibrated membranes of that receptor
# (membraneforger/examples/preeq_cg_cellmem/README.md lists them); the chosen file is logged and recorded in run_manifest.json.
BUNDLED_MEMBRANES = {"1": "GPR*_cg_cellmem.gro", "2": "KOR*_cg_cellmem.gro"}
# --cg NAME: one particular bundled frame by name (--list-membranes lists them), for a repeatable build.
BUNDLED_FRAME_SUFFIX = "_cg_cellmem.gro"
RECEPTORS = {"GPR": "GPR139", "KOR": "kappa opioid receptor"}
# --runs N: the complex is built N times, each time in a different frame of the chosen series (the first N frames of
# the series, so the set is the same from one invocation to the next); the bound is the number of frames the series holds.


ANGSTROM_PER_NM = 10.0
# Every length on the command line is in nm and its option name says so (--box-nm, --bilayer-z-nm, --fit-max-core-rmsd-nm,
# --xy-buffer). The unsuffixed spellings of earlier releases (--box, --bilayer-z, --fit-max-core-rmsd) took angstroms; they still
# work, are converted here, and print which nm option they stand for. Giving both spellings of one option is an error.


def angstrom_to_nm(value: float) -> float:
    """Convert a length in angstrom to nm (1 nm = 10 A)."""
    return value / ANGSTROM_PER_NM


def nm_to_angstrom(value: float) -> float:
    """Convert a length in nm to angstrom (1 nm = 10 A)."""
    return value * ANGSTROM_PER_NM


def default_out(source: Path, cwd: Path | None = None) -> Path:
    """The default --out: ./<input structure stem>_membraneforger; the structure is the all-atom input (the assembled
    system with --membrane), the same for a single build and for --runs, whose frames each take a subfolder."""
    return Path(cwd or Path.cwd()) / f"{Path(source).stem}_membraneforger"


def bundled_frames() -> dict:
    """The bundled frames by name (for example KOR_Membrane1 or GPR139_Membrane5), each mapped to its file in
    membraneforger/examples/preeq_cg_cellmem/."""
    return {p.name[:-len(BUNDLED_FRAME_SUFFIX)]: p for p in sorted(BUNDLED_DIR.glob("*" + BUNDLED_FRAME_SUFFIX))
            if p.name[:-len(BUNDLED_FRAME_SUFFIX)][:3].upper() in RECEPTORS}


def bundled_series(code: str) -> list:
    """The bundled frames of --cg 1 (GPR139) or --cg 2 (kappa opioid receptor), in name order."""
    frames = sorted(BUNDLED_DIR.glob(BUNDLED_MEMBRANES[code]))
    if not frames:
        raise SystemExit(f"no bundled membrane matches membraneforger/examples/preeq_cg_cellmem/{BUNDLED_MEMBRANES[code]}")
    return frames


def bundled_membrane(code: str) -> Path | None:
    """Resolve a --cg value that names a bundled membrane; None when it is not one (a file path).

    "1" picks one of the GPR139 frames at random and "2" one of the kappa opioid receptor frames; a frame name such as
    "GPR139_Membrane5" or "KOR_Membrane5" (any case, with or without the _cg_cellmem.gro suffix) picks that frame; the
    names of earlier versions ("GPR5", "KOR5") still work.
    """
    if code in BUNDLED_MEMBRANES:
        return random.choice(bundled_series(code))
    if Path(code).is_file() or "/" in code or os.sep in code:
        return None  # a file of the user's own, or a path: --cg FILE
    name = code.strip().upper()
    for suffix in (BUNDLED_FRAME_SUFFIX.upper(), ".GRO"):
        name = name[:-len(suffix)] if name.endswith(suffix) else name
    old = re.fullmatch(r"(GPR|KOR)(\d+)", name)  # the names of earlier versions: KOR5 is KOR_Membrane5, GPR5 is GPR139_Membrane5
    if old:
        name = f"{'GPR139' if old.group(1) == 'GPR' else 'KOR'}_MEMBRANE{old.group(2)}"
    frames = bundled_frames()
    by_upper = {frame.upper(): path for frame, path in frames.items()}
    if name in by_upper:
        return by_upper[name]
    if re.fullmatch(r"(GPR|KOR)(139)?_MEMBRANE\d+", name):  # looks like a frame name, but no such frame is bundled
        raise SystemExit(f"no bundled membrane named {code}: the bundled frames are {', '.join(frames)} (--list-membranes)")
    return None  # a file path


def describe_bundled_membranes() -> str:
    """The text of --list-membranes: every bundled frame by name."""
    lines = ["bundled Martini 3 cell membranes (membraneforger/examples/preeq_cg_cellmem/, 50 us each; pass a name to --cg):"]
    for name, path in bundled_frames().items():
        lines.append(f"  {name:<18} {RECEPTORS[name[:3].upper()]:<22} {path.name}")
    lines.append("--cg 1 picks a GPR139 frame at random and --cg 2 a kappa opioid receptor frame; --runs N builds the complex in "
                 "the first N frames of that series, one system each")
    return "\n".join(lines)


def run_sessions(session: Session, all_atom: Path, frames: list, command: list) -> int:
    """--runs N: build the complex once per frame, each build in its own subdirectory of session.out, and sum up."""
    results = []
    tags = [frame.name[:-len(BUNDLED_FRAME_SUFFIX)] for frame in frames]
    if len(set(tags)) != len(tags):
        raise SystemExit(f"--runs: two frames would share one output folder ({', '.join(tags)})")
    for tag in tags:  # refuse before the first build rather than after an hour of building when a later run's folder is taken
        refuse_overwrite(session.out / tag, overwrite=session.overwrite, overwrite_restart=session.overwrite_restart)
    for frame in frames:
        tag = frame.name[:-len(BUNDLED_FRAME_SUFFIX)]
        run = replace(session, out=session.out / tag, name=f"{session.name}_{tag}", work=None, timings={}, record={}, analysis=None)
        print(f"[{session.name}] run {len(results) + 1} of {len(frames)}: membrane {frame.name} -> {run.out}")
        results.append((tag, build(run, all_atom, frame.resolve(), None, command)))
    failed = [tag for tag, code in results if code]
    print(f"[{session.name}] {len(frames) - len(failed)} of {len(frames)} runs passed"
          + (f"; failed: {', '.join(failed)}" if failed else "") + f" ({session.out})")
    return int(bool(failed))


def parse_chains(*values: str | None) -> tuple:
    """Turn --orient-chain / --orient-chains / ORIENT_CHAINS values ("R", "A,B", "A B") into a tuple of labels."""
    chains = []
    for value in values:
        for label in (value or "").replace(",", " ").split():
            if label not in chains:
                chains.append(label)
    return tuple(chains)


def parse_box(value: str | None) -> tuple | None:
    """Parse the BOX header value: 'auto' (None) or three positive edge lengths in nm as 'x,y,z'; returned in A (Session.box_a)."""
    text = (value or "auto").strip().lower()
    if text in ("", "auto"):
        return None
    try:
        parts = tuple(float(v) for v in text.replace("x", ",").split(","))
    except ValueError:
        parts = ()
    if len(parts) != 3 or not all(v > 0 for v in parts):
        raise SystemExit(f"BOX must be 'auto' or three positive edge lengths in nm as x,y,z, not {value!r}")
    return tuple(nm_to_angstrom(v) for v in parts)


def resolve_lengths(args, parser: argparse.ArgumentParser) -> None:
    """Fold the nm options and the angstrom spellings of earlier releases into args.box, args.bilayer_z and args.fit_max_core_rmsd (A).

    The session and Settings keep their internal lengths in A (the PDB unit); the command line is nm. An unsuffixed legacy option
    is read as A, converted, and reported once with the nm option it stands for; giving both spellings is an error.
    """
    for new, old, scale in (("box_nm", "box", 3), ("bilayer_z_nm", "bilayer_z", 1), ("fit_max_core_rmsd_nm", "fit_max_core_rmsd", 1)):
        given_new, given_old = getattr(args, new), getattr(args, old)
        flag_new, flag_old = "--" + new.replace("_", "-"), "--" + old.replace("_", "-")
        if given_new is not None and given_old is not None:
            parser.error(f"{flag_new} (nm) and {flag_old} (A, earlier releases) give the same length twice; use {flag_new} only")
        if given_old is not None:
            as_nm = [angstrom_to_nm(v) for v in (given_old if scale == 3 else [given_old])]
            print(f"WARNING: {flag_old} is in angstrom (earlier releases); read as {flag_new} {' '.join(f'{v:g}' for v in as_nm)} nm. "
                  f"Use {flag_new} (nm)", file=sys.stderr)
        elif given_new is not None:
            setattr(args, old, [nm_to_angstrom(v) for v in given_new] if scale == 3 else nm_to_angstrom(given_new))
    if args.fit_max_core_rmsd is None:
        args.fit_max_core_rmsd = Settings().fit_max_core_rmsd_a


def orientation_request(args) -> OrientationRequest:
    """Build the orientation request from the configuration header, overridden by the command line."""
    mode = (args.orientation or ORIENTATION).lower()
    if mode not in ORIENTATION_MODES:
        raise SystemExit(f"--orientation must be one of {', '.join(ORIENTATION_MODES)}, not {args.orientation!r}")
    nterm = (args.nterm_side or NTERM_SIDE).lower()
    if nterm not in NTERM_SIDES:
        raise SystemExit(f"--nterm-side must be one of {', '.join(NTERM_SIDES)}, not {args.nterm_side!r}")
    chains = parse_chains(args.orient_chain, args.orient_chains) or parse_chains(ORIENT_CHAINS)
    return OrientationRequest(mode=mode, chains=chains, nterm_side=nterm, pdb_id=(args.pdb_id or PDB_ID or None),
                              ppm_exe=args.ppm_exe, ppm_membrane=(args.ppm_membrane if args.ppm_membrane is not None else PPM_MEMBRANE),
                              ppm_heteroatoms=bool(args.ppm_heteroatoms),
                              residues=parse_residue_ranges(args.orient_residues if args.orient_residues is not None
                                                            else ORIENT_RESIDUES),
                              opm_file=args.opm_file.resolve() if args.opm_file else None,
                              opm_cache=args.opm_cache.resolve() if args.opm_cache else None)


def make_parser() -> argparse.ArgumentParser:
    """Build the argument parser; only --all-atom and --coarse-grain are per-system inputs."""
    defaults = Settings()
    parser = argparse.ArgumentParser(
        prog="membraneforger.py",
        description="Build a validated, energy-minimized all-atom CHARMM36 membrane system (em.gro) from an all-atom "
                    "protein/ligand PDB and a Martini 3 coarse-grained simulation frame.")
    parser.add_argument("--aa", "--all-atom", dest="all_atom", type=Path, metavar="PDB",
                        help="all-atom protein (or protein/ligand) structure with its membrane normal along z; "
                             "the only source of protein and ligand chemistry")
    parser.add_argument("--cg", "--coarse-grain", dest="coarse_grain", default="1", metavar="1|2|NAME|FILE",
                        help="membrane: 1 = a random bundled GPR139 frame, 2 = a random bundled kappa opioid receptor frame, "
                             "NAME = one particular bundled frame (--list-membranes lists the names; the protein is embedded "
                             "into any of these), or a Martini 3 frame of your own complex, given as FILE or custom=FILE (default: 1)")
    parser.add_argument("--runs", type=int, default=1, metavar="N",
                        help="build the complex N times, each in a different bundled membrane of the --cg 1 or 2 series (its "
                             "first N frames, at most the number the series holds), one system per run in <out>/<frame>/ (default: %(default)s)")
    parser.add_argument("--list-membranes", action="store_true", help="list the bundled membranes by name and exit")
    parser.add_argument("--embed", action="store_true",
                        help="embed the protein into the membrane of a --cg FILE instead of fitting it onto the frame's "
                             "protein (always the case for the bundled membranes)")
    parser.add_argument("--embed-site", choices=EMBED_SITES, default="hole",
                        help="where an embedded complex goes in the frame: hole = where the frame's own protein was; free = "
                             "the stretch of membrane farthest from that protein, an unbroken bilayer, for a complex much smaller "
                             "than the frame's receptor (e.g. one transmembrane helix) (default: %(default)s)")
    parser.add_argument("--box-nm", nargs=3, type=float, metavar=("X", "Y", "Z"),
                        help=f"opt-in box edges in nm (default: BOX = {BOX}: the membrane is sliced to the complex plus --xy-buffer "
                             "in x and y, and z is sized around the bilayer midplane); x and y may only be smaller than the "
                             "membrane patch, which is then cut around the complex, z must leave the default water padding")
    parser.add_argument("--box", nargs=3, type=float, metavar=("X", "Y", "Z"), help=argparse.SUPPRESS)  # earlier releases: edges in A
    parser.add_argument("--bilayer-z-nm", type=float, metavar="Z",
                        help="z of the bilayer centre in the all-atom input, in nm (default: found from the hydrophobic belt)")
    parser.add_argument("--bilayer-z", type=float, metavar="Z", help=argparse.SUPPRESS)  # earlier releases: A
    parser.add_argument("--dellipid", action="append", default=[], metavar="LIPID",
                        help=f"remove every molecule of this lipid (repeatable): {', '.join(LIPID_NAMES)}")
    parser.add_argument("--addlipid", metavar="LIPID",
                        help=f"turn the lipids removed by --dellipid into this lipid instead ({', '.join(CONVERTIBLE)})")
    parser.add_argument("--exclude-disulfide", action="append", default=[], metavar="CHAIN:RESID-CHAIN:RESID",
                        help="keep this cysteine pair reduced although its SG atoms are close enough to be bonded "
                             "(repeatable), in the input's own numbering, e.g. A:860-A:1252")
    parser.add_argument("-o", "--out", type=Path,
                        help="output directory (default: ./<all-atom input stem>_membraneforger, also with --runs, whose frames each get a "
                             "subfolder; with --membrane the stem of that file)")
    parser.add_argument("--overwrite", action="store_true",
                        help="build into an --out that already holds the outputs of an earlier build, replacing them (otherwise refused)")
    parser.add_argument("--overwrite-restart", action="store_true",
                        help="with --overwrite, also when --out holds checkpoints or trajectories of a simulation started from the build; "
                             "those files are never deleted, but em.gro, topol.top and the index files they depend on are replaced")
    parser.add_argument("--name", help="system name for logs and [ system ] (default: output directory name)")
    parser.add_argument("--toppar", type=Path, default=os.environ.get("MEMBRANEFORGER_TOPPAR"),
                        help="force-field directory with charmm36.ff/, toppar/, protein_parameters/ "
                             "(default: $MEMBRANEFORGER_TOPPAR, else membraneforger/forcefield/ in this repository)")
    parser.add_argument("--data", type=Path, default=os.environ.get("MEMBRANEFORGER_DATA"),
                        help="backmapping data directory with map.dat and extra force-field XML such as GM3.xml "
                             "(default: $MEMBRANEFORGER_DATA, else membraneforger/backmap_data/ in this repository)")
    parser.add_argument("--gmx", help="GROMACS command (default: gmx, gmx_mpi or gmx_d on PATH)")
    parser.add_argument("--mstool-python", help="python interpreter that can import mstool (default: the one running this)")
    parser.add_argument("--ntomp", type=int, default=os.cpu_count() or 1, help="threads (default: all cores)")
    parser.add_argument("--nsteps", type=int, default=NSTEPS, help="maximum EM steps (default: %(default)s)")
    parser.add_argument("--fit-max-core-rmsd-nm", type=float, default=None, metavar="NM",
                        help="largest accepted core RMSD of the AA-on-CG rigid fit, in nm "
                             f"(default: {angstrom_to_nm(defaults.fit_max_core_rmsd_a):g})")
    parser.add_argument("--fit-max-core-rmsd", type=float, default=None, metavar="A", help=argparse.SUPPRESS)  # earlier releases: A
    parser.add_argument("--fit-min-core-fraction", type=float, default=defaults.fit_min_core_fraction, metavar="F",
                        help="smallest accepted fraction of backbone pairs kept by the fit (default: %(default)s)")
    parser.add_argument("--membrane", type=Path, metavar="PDB",
                        help="skip backmapping and start from an already assembled protein + membrane PDB with CRYST1 "
                             "(replaces --all-atom/--coarse-grain)")
    orient = parser.add_argument_group("membrane orientation (the complete complex is moved as one rigid body)")
    orient.add_argument("--orientation", choices=ORIENTATION_MODES, default=None,
                        help=f"auto: exact OPM entry when a PDB ID is known and matches, else local PPM 3.0; ppm/opm: force "
                             f"one provider; none: use the input coordinates as given (default: {ORIENTATION})")
    orient.add_argument("--orient-chain", metavar="CHAIN", help="the chain that spans or associates with the membrane")
    orient.add_argument("--orient-chains", metavar="A,B", help="several anchor chains (the first one sets the N-terminus side)")
    orient.add_argument("--orient-residues", metavar="FIRST-LAST", default=None,
                        help="orient on these residues of the anchor chain(s) only, e.g. the transmembrane helix 343-363 "
                             "(several: 343-363,370-380; numbers as in the input). PPM is run on them alone and the whole "
                             "complex follows rigidly; needs --nterm-side for the first selected residue. With --orientation "
                             "none it names the membrane-embedded segment of an input already along z: nothing is rotated and "
                             "the bilayer centre goes at the midpoint of that segment's CA z range")
    orient.add_argument("--nterm-side", choices=NTERM_SIDES, default=None,
                        help=f"side of the membrane the N terminus of the first anchor chain lies on; needed by PPM when no "
                             f"exact OPM entry gives it (default: {NTERM_SIDE})")
    orient.add_argument("--pdb-id", metavar="ID", help="exact PDB ID of the structure, for the OPM reference (default: HEADER record)")
    orient.add_argument("--ppm-exe", metavar="PATH", help="PPM 3.0 executable (immers) with res.lib next to it "
                                                          "(default: $MEMBRANEFORGER_PPM, then immers/ppm3 on PATH)")
    orient.add_argument("--ppm-membrane", metavar="CODE", default=None,
                        help="advanced: PPM 3.0 membrane code such as PMm or GnI (default: undefined flat bilayer)")
    orient.add_argument("--ppm-heteroatoms", action="store_true", help="advanced: submit the anchor chains' heteroatoms to PPM too")
    orient.add_argument("--opm-file", type=Path, metavar="PDB", help="an already downloaded OPM/OPRLM coordinate file (offline use)")
    orient.add_argument("--opm-cache", type=Path, metavar="DIR", help="cache for downloaded OPM files (default: ~/.cache/membraneforger/opm)")
    parser.add_argument("--xy-buffer", "--xy-buffer-nm", dest="xy_buffer", type=float, default=defaults.box_xy_buffer_nm, metavar="NM",
                        help="membrane kept around the complex on each side in x and y when the box is auto (default: %(default)s nm)")
    parser.add_argument("--counterion-mode", choices=("additive", "inclusive"), default=defaults.counterion_mode,
                        help="additive: neutralizing counterions plus the requested 0.15 M NaCl pairs on top; inclusive: the "
                             "CHARMM-GUI 2017 rule, the counterion species equals the NaCl target and the other ion target - |charge| "
                             "(default: %(default)s)")
    check = parser.add_argument_group("slice validation (the sliced membrane against the embedded membrane it was cut from)")
    check.add_argument("--apl-validate", choices=("yes", "no"), default="yes" if defaults.apl_validate else "no",
                       help="stop before backmapping when a leaflet's Voronoi APL changes by more than --apl-tolerance (default: %(default)s)")
    check.add_argument("--apl-tolerance", type=float, default=defaults.apl_slice_tolerance_percent, metavar="PERCENT",
                       help="largest accepted change of the leaflet APL across slicing, in percent (default: %(default)s)")
    check.add_argument("--rdf-validate", choices=("yes", "no"), default="yes" if defaults.rdf_validate else "no",
                       help="compare the lateral headgroup RDF of the slice with the embedded membrane (warns only; default: %(default)s)")
    check.add_argument("--cg-validate", choices=("yes", "no"), default="yes" if defaults.cg_validate else "no",
                       help="stop before embedding when the coarse-grained frame has an empty leaflet, a broken molecule, a strongly "
                            "tilted bilayer or pervasive bead overlaps (default: %(default)s)")
    check.add_argument("--slice-offset-search", action="store_true",
                       help="opt in to shifting the slice window by up to 0.5 nm to the offset that best preserves lipid density "
                            "(default: the window is centred on the complex, which had the smaller APL deviation in the bundled frames)")
    check.add_argument("--no-slice-offset", action="store_true", help=argparse.SUPPRESS)  # kept from earlier releases: centring is now the default
    after = parser.add_argument_group("checks after minimization and handling of an earlier build")
    after.add_argument("--backbone-check", choices=("yes", "no", "report"), default=defaults.backbone_check,
                       help="compare the protein backbone after EM with the structure EM started from: yes stops the build on a peptide "
                            "bond EM twisted or flipped to another isomer, a C-alpha or Thr/Ile C-beta it inverted, or (with the EM backbone "
                            "restraints on) a phi/psi that moved more than 30 deg beyond the restraint width; report lists the same as "
                            "warnings; no skips it (default: %(default)s)")
    after.add_argument("--backbone-drift-margin", type=float, default=defaults.backbone_drift_margin_deg, metavar="DEG",
                       help="degrees added to the 30 deg restraint width before a phi/psi drift fails the backbone check (default: %(default)s)")
    after.add_argument("--strict-segments", metavar="CHAIN:FIRST-LAST[,...]", default="",
                       help="residues of the input that were modelled rather than observed (e.g. R:338-345,A:5-9): a Ramachandran outlier after EM "
                            "in one of them fails the backbone check, even if the residue was not an outlier before; chain letters and residue "
                            "numbers are those of the all-atom input")
    after.add_argument("--archive-previous", choices=("no", "reports", "all"), default="no",
                       help="with --overwrite, first copy the earlier build's reports (reports) or reports and structures (all: em.gro, "
                            "topol.top, index_ini.ndx) to <out>/previous_builds/<start time>_<status>/ so a rejected attempt survives its "
                            "rebuild (default: %(default)s)")
    return parser


def locate_forcefield(requested: Path | None) -> Path:
    """Return the force-field directory, checking it holds what the topology builder reads."""
    forcefield = requested or next((p for p in (PACKAGE_DIR / "forcefield", PACKAGE_DIR / "modified_pipeline") if p.is_dir()), None)
    if forcefield is None:
        raise SystemExit("cannot find the force-field directory; pass --toppar or set MEMBRANEFORGER_TOPPAR")
    forcefield = Path(forcefield).resolve()
    missing = [p for p in ("charmm36.ff/forcefield.itp", "charmm36.ff/merged.rtp", "toppar",
                           "protein_parameters/top_all36_prot.rtf") if not (forcefield / p).exists()]
    if missing:
        raise SystemExit(f"--toppar {forcefield} lacks {missing}")
    return forcefield


def locate_data(requested: Path | None) -> Path:
    """Return the backmapping data directory (map.dat and extra force-field XML files)."""
    candidates = (PACKAGE_DIR / "backmap_data", PACKAGE_DIR / "Incretin_backmap")
    return Path(requested or next((p for p in candidates if p.is_dir()), candidates[0])).resolve()


def parse_disulfide_pair(spec: str) -> tuple[str, int, int]:
    """Read one --exclude-disulfide CHAIN:RESID-CHAIN:RESID into (chain, lower resid, higher resid)."""
    match = re.fullmatch(r"\s*(\S):(\d+)\s*-\s*(\S):(\d+)\s*", spec)
    if not match:
        raise SystemExit(f"--exclude-disulfide {spec}: expected CHAIN:RESID-CHAIN:RESID, e.g. A:860-A:1252")
    left, a, right, b = match.group(1), int(match.group(2)), match.group(3), int(match.group(4))
    if left != right:
        raise SystemExit(f"--exclude-disulfide {spec}: both cysteines must be in the same chain")
    if a == b:
        raise SystemExit(f"--exclude-disulfide {spec}: names one residue twice")
    return (left, min(a, b), max(a, b))


def main(argv: list | None = None) -> int:
    """Run one build from command-line arguments and return its exit code."""
    parser = make_parser()
    args = parser.parse_args(argv)
    resolve_lengths(args, parser)
    if args.list_membranes:
        print(describe_bundled_membranes())
        return 0
    given_cg = "--cg" in sys.argv or "--coarse-grain" in sys.argv
    if args.membrane and (args.all_atom or given_cg):
        parser.error("--membrane replaces --aa and --cg; give one or the other")
    if not args.membrane and not args.all_atom:
        parser.error("--aa is required")
    embed, cg, frames = args.embed, str(args.coarse_grain), []
    if args.runs < 1:
        parser.error("--runs must be at least 1")
    if args.runs > 1 and (args.membrane or cg not in BUNDLED_MEMBRANES):
        parser.error("--runs builds in several bundled membranes of one series: use it with --cg 1 or --cg 2 (--aa, no --membrane)")
    if args.runs > 1:
        try:
            held = bundled_series(cg)
        except SystemExit as exc:
            parser.error(str(exc))
        if args.runs > len(held):
            parser.error(f"--runs must be between 1 and {len(held)}, the number of frames in the --cg {cg} series "
                         f"({', '.join(p.name[:-len(BUNDLED_FRAME_SUFFIX)] for p in held)})")
    try:
        bundled = bundled_series(cg)[:args.runs] if args.runs > 1 else [bundled_membrane(cg)]
    except SystemExit as exc:
        parser.error(str(exc))
    if bundled[0] is not None:
        args.coarse_grain, embed, frames = bundled[0], True, bundled
    else:
        args.coarse_grain = Path(cg[len("custom="):] if cg.lower().startswith("custom=") else cg)
    if args.membrane:
        args.coarse_grain = None
    for path in (args.all_atom, args.coarse_grain, args.membrane):
        if path and not path.is_file():
            parser.error(f"missing input file {path}")
    if args.box and (len(args.box) != 3 or min(args.box) <= 0):
        parser.error("--box-nm needs three positive edge lengths in nm")
    if args.membrane and (args.box or args.dellipid or args.addlipid or args.embed or args.embed_site != "hole"):
        parser.error("--box-nm, --dellipid, --addlipid, --embed and --embed-site need --aa and --cg, not --membrane")
    if args.embed_site != "hole" and not embed:
        parser.error("--embed-site applies when the protein is embedded: a bundled membrane (--cg 1, 2 or a frame name) or --embed")
    try:
        for name in args.dellipid + ([args.addlipid] if args.addlipid else []):
            lipid_name(name)
    except SystemExit as exc:
        parser.error(str(exc))
    try:
        exclude_disulfides = tuple(parse_disulfide_pair(spec) for spec in args.exclude_disulfide)
    except SystemExit as exc:
        parser.error(str(exc))
    try:
        forcefield = locate_forcefield(args.toppar)
        gmx = find_gromacs(args.gmx)
    except SystemExit as exc:
        parser.error(str(exc))
    data = locate_data(args.data)
    python = args.mstool_python or sys.executable
    if not args.membrane:
        if not (data / "map.dat").is_file():
            parser.error(f"--data {data} has no map.dat")
        if not shutil.which(python):
            parser.error(f"--mstool-python {python} not found")
    try:
        orientation = orientation_request(args)
        box = tuple(args.box) if args.box else (None if args.membrane else parse_box(BOX))
    except SystemExit as exc:
        parser.error(str(exc))
    if orientation.mode == "none" and orientation.residues and not args.membrane:
        if args.bilayer_z is not None:
            parser.error("--orient-residues with --orientation none sets the bilayer centre from that segment; give it or --bilayer-z-nm, not both")
        if not embed:
            parser.error("--orient-residues with --orientation none places an embedded protein (a bundled membrane or --embed); "
                         "a fitted frame takes the coarse-grained protein's own pose")
    if orientation.mode != "none" and args.bilayer_z is not None and not args.membrane:
        parser.error("--bilayer-z-nm applies with --orientation none; an oriented complex has its bilayer centre at z = 0")
    if args.membrane and (orientation.mode != "none" and (orientation.chains or orientation.residues or args.orientation or args.pdb_id)):
        parser.error("--membrane starts from an assembled system; membrane orientation options do not apply to it")
    if args.membrane:
        orientation = OrientationRequest(mode="none")
    if orientation.opm_file and not orientation.opm_file.is_file():
        parser.error(f"missing --opm-file {orientation.opm_file}")
    out = (args.out or default_out(args.all_atom or args.membrane)).resolve()
    if not args.xy_buffer > 0:
        parser.error("--xy-buffer must be positive")
    if not args.apl_tolerance > 0:
        parser.error("--apl-tolerance must be positive")
    if args.backbone_drift_margin < 0:
        parser.error("--backbone-drift-margin must not be negative")
    try:
        strict_segments = parse_segments(args.strict_segments)
    except ValueError as exc:
        parser.error(f"--strict-segments: {exc}")
    if args.archive_previous != "no" and not args.overwrite:
        parser.error("--archive-previous archives the build --overwrite replaces; give --overwrite too")
    settings = Settings(fit_max_core_rmsd_a=args.fit_max_core_rmsd, fit_min_core_fraction=args.fit_min_core_fraction,
                        box_xy_buffer_nm=args.xy_buffer, apl_validate=args.apl_validate == "yes", apl_slice_tolerance_percent=args.apl_tolerance,
                        rdf_validate=args.rdf_validate == "yes", cg_validate=args.cg_validate == "yes",
                        slice_optimize_offset=args.slice_offset_search and not args.no_slice_offset,
                        counterion_mode=args.counterion_mode, backbone_check=args.backbone_check,
                        backbone_drift_margin_deg=args.backbone_drift_margin)
    session = Session(out=out, name=args.name or out.name, gmx=gmx, forcefield=forcefield, data=data, python=python,
                      ntomp=args.ntomp, nsteps=args.nsteps, settings=settings, embed=embed, embed_site=args.embed_site,
                      box_a=box, bilayer_z_a=args.bilayer_z, delete_lipids=list(args.dellipid), add_lipid=args.addlipid,
                      exclude_disulfides=exclude_disulfides, orientation=orientation,
                      overwrite=args.overwrite, overwrite_restart=args.overwrite_restart,
                      archive_previous=args.archive_previous, strict_segments=strict_segments)
    resolved = [p.resolve() if p else None for p in (args.all_atom, args.coarse_grain, args.membrane)]
    if args.runs > 1:
        return run_sessions(session, resolved[0], frames, [sys.executable] + sys.argv)
    return build(session, resolved[0], resolved[1], resolved[2], [sys.executable] + sys.argv)


if __name__ == "__main__":
    sys.exit(main())
