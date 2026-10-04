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

from .config import NSTEPS, Settings
from .pipeline import Session, build
from .runtools import find_gromacs

__all__ = ['make_parser', 'locate_forcefield', 'locate_data', 'main']

REPOSITORY = Path(__file__).resolve().parents[1]


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
    source = args.membrane or args.coarse_grain
    out = (args.out or Path.cwd() / f"{source.stem}_membraneforger").resolve()
    settings = Settings(fit_max_core_rmsd_a=args.fit_max_core_rmsd, fit_min_core_fraction=args.fit_min_core_fraction)
    session = Session(out=out, name=args.name or out.name, gmx=gmx, forcefield=forcefield, data=data, python=python,
                      ntomp=args.ntomp, nsteps=args.nsteps, settings=settings)
    resolved = [p.resolve() if p else None for p in (args.all_atom, args.coarse_grain, args.membrane)]
    return build(session, resolved[0], resolved[1], resolved[2], [sys.executable] + sys.argv)


if __name__ == "__main__":
    sys.exit(main())
