#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Example: build one system from Python instead of the command line.
#//=============================================================
"""Usage: python example.py ALL_ATOM.pdb COARSE_GRAIN.gro OUTPUT_DIR [GMX]"""
# Do not modify the next 5 lines
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import membraneforger as mf

# Modifications possible below
all_atom, coarse_grain, out = (Path(p).resolve() for p in sys.argv[1:4])
gmx = mf.find_gromacs(sys.argv[4] if len(sys.argv) > 4 else None)

# Acceptance criteria can be tightened or relaxed here; the defaults are documented in config.Settings.
settings = mf.Settings(fit_max_core_rmsd_a=5.0, fit_min_core_fraction=0.5)

session = mf.Session(out=out, name=out.name, gmx=gmx, forcefield=mf.locate_forcefield(None),
                     data=mf.locate_data(None), python=sys.executable,
                     ntomp=8, nsteps=mf.NSTEPS, settings=settings)
sys.exit(mf.build(session, all_atom, coarse_grain, None, [sys.executable] + sys.argv))
