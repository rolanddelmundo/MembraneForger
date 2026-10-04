#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Command line: python membraneforger.py --all-atom AA.pdb --coarse-grain CG.gro --out DIR
#//=============================================================
"""Parse the command line, locate installation resources, and run one build."""
import argparse
import os
import shutil
import sys
from pathlib import Path

from .config import BOX, NSTEPS, NTERM_SIDE, ORIENTATION, ORIENT_CHAINS, PDB_ID, PPM_MEMBRANE, Settings
from .orientation import NTERM_SIDES, ORIENTATION_MODES, OrientationRequest
from .pipeline import Session, build
from .runtools import find_gromacs

__all__ = ['make_parser', 'locate_forcefield', 'locate_data', 'parse_chains', 'parse_box', 'orientation_request', 'main']

REPOSITORY = Path(__file__).resolve().parents[1]


def parse_chains(*values: str | None) -> tuple:
    """Turn --orient-chain / --orient-chains / ORIENT_CHAINS values ("R", "A,B", "A B") into a tuple of labels."""
    chains = []
    for value in values:
        for label in (value or "").replace(",", " ").split():
            if label not in chains:
                chains.append(label)
    return tuple(chains)


def parse_box(value: str | None) -> tuple | None:
    """Parse BOX / --box: 'auto' (None) or three positive lengths in nm 'x,y,z'."""
    text = (value or "auto").strip().lower()
    if text in ("", "auto"):
        return None
    try:
        parts = tuple(float(v) for v in text.replace("x", ",").split(","))
    except ValueError:
        parts = ()
    if len(parts) != 3 or not all(v > 0 for v in parts):
        raise SystemExit(f"--box must be 'auto' or three positive lengths in nm as x,y,z, not {value!r}")
    return parts


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
                              opm_file=args.opm_file.resolve() if args.opm_file else None,
                              opm_cache=args.opm_cache.resolve() if args.opm_cache else None)


def make_parser() -> argparse.ArgumentParser:
    """Build the argument parser; only --all-atom and --coarse-grain are per-system inputs."""
    defaults = Settings()
    parser = argparse.ArgumentParser(
        prog="membraneforger.py",
        description="Build a validated, energy-minimized all-atom CHARMM36 membrane system (em.gro) from an all-atom "
                    "protein/ligand PDB and a Martini 3 coarse-grained simulation frame.")
    parser.add_argument("--all-atom", type=Path, metavar="PDB",
                        help="all-atom protein/ligand complex; the only source of protein and ligand chemistry")
    parser.add_argument("--coarse-grain", type=Path, metavar="PDB|GRO",
                        help="Martini 3 system with its box; supplies the membrane and where the protein sits")
    parser.add_argument("-o", "--out", type=Path, help="output directory (default: ./<coarse-grain stem>_membraneforger)")
    parser.add_argument("--name", help="system name for logs and [ system ] (default: output directory name)")
    parser.add_argument("--toppar", type=Path, default=os.environ.get("MEMBRANEFORGER_TOPPAR"),
                        help="force-field directory with charmm36.ff/, toppar/, protein_parameters/ "
                             "(default: $MEMBRANEFORGER_TOPPAR, else forcefield/ in this repository)")
    parser.add_argument("--data", type=Path, default=os.environ.get("MEMBRANEFORGER_DATA"),
                        help="backmapping data directory with map.dat and extra force-field XML such as GM3.xml "
                             "(default: $MEMBRANEFORGER_DATA, else backmap_data/ in this repository)")
    parser.add_argument("--gmx", help="GROMACS command (default: gmx, gmx_mpi or gmx_d on PATH)")
    parser.add_argument("--mstool-python", help="python interpreter that can import mstool (default: the one running this)")
    parser.add_argument("--ntomp", type=int, default=os.cpu_count() or 1, help="threads (default: all cores)")
    parser.add_argument("--nsteps", type=int, default=NSTEPS, help="maximum EM steps (default: %(default)s)")
    parser.add_argument("--fit-max-core-rmsd", type=float, default=defaults.fit_max_core_rmsd_a, metavar="A",
                        help="largest accepted core RMSD of the AA-on-CG rigid fit (default: %(default)s A)")
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
    parser.add_argument("--box", metavar="auto|X,Y,Z", default=None,
                        help=f"auto: x,y from the membrane cell and z sized around the bilayer midplane; or an opt-in box in nm "
                             f"whose x,y must equal the membrane cell (default: {BOX})")
    return parser


def locate_forcefield(requested: Path | None) -> Path:
    """Return the force-field directory, checking it holds what the topology builder reads."""
    forcefield = requested or next((p for p in (REPOSITORY / "forcefield", REPOSITORY / "modified_pipeline") if p.is_dir()), None)
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
    candidates = (REPOSITORY / "backmap_data", REPOSITORY / "Incretin_backmap")
    return Path(requested or next((p for p in candidates if p.is_dir()), candidates[0])).resolve()


def main(argv: list | None = None) -> int:
    """Run one build from command-line arguments and return its exit code."""
    parser = make_parser()
    args = parser.parse_args(argv)
    if args.membrane and (args.all_atom or args.coarse_grain):
        parser.error("--membrane replaces --all-atom and --coarse-grain; give one or the other")
    if not args.membrane and not (args.all_atom and args.coarse_grain):
        parser.error("--all-atom and --coarse-grain are both required")
    for path in (args.all_atom, args.coarse_grain, args.membrane):
        if path and not path.is_file():
            parser.error(f"missing input file {path}")
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
        box = parse_box(args.box if args.box is not None else BOX)
    except SystemExit as exc:
        parser.error(str(exc))
    if args.membrane and (orientation.mode != "none" and (orientation.chains or args.orientation or args.pdb_id)):
        parser.error("--membrane starts from an assembled system; membrane orientation options do not apply to it")
    if args.membrane:
        orientation = OrientationRequest(mode="none")
    if orientation.opm_file and not orientation.opm_file.is_file():
        parser.error(f"missing --opm-file {orientation.opm_file}")
    source = args.membrane or args.coarse_grain
    out = (args.out or Path.cwd() / f"{source.stem}_membraneforger").resolve()
    settings = Settings(fit_max_core_rmsd_a=args.fit_max_core_rmsd, fit_min_core_fraction=args.fit_min_core_fraction)
    session = Session(out=out, name=args.name or out.name, gmx=gmx, forcefield=forcefield, data=data, python=python,
                      ntomp=args.ntomp, nsteps=args.nsteps, settings=settings, orientation=orientation, box=box)
    resolved = [p.resolve() if p else None for p in (args.all_atom, args.coarse_grain, args.membrane)]
    return build(session, resolved[0], resolved[1], resolved[2], [sys.executable] + sys.argv)


if __name__ == "__main__":
    sys.exit(main())
