"""Shared fixtures for the MembraneForger tests: repository example data and tiny synthetic systems."""
import json
import os
import shutil
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
import membraneforger as mf  # noqa: E402

first = lambda *paths: next((REPO / p for p in paths if (REPO / p).exists()), REPO / paths[0])
FORCEFIELD = first("membraneforger/forcefield", "modified_pipeline")
DATA = first("membraneforger/backmap_data", "Incretin_backmap")
AA_PDB = first("membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_prot-lig.pdb", "tests/glpa_sapi_6WHC_MTZP_run1/prot-lig.pdb")
CG_GRO = first("membraneforger/examples/preeq_cg_cellmem/GCGR_MTZP_cg_cellmem.gro", "final_gro_10us/6WHC_MTZP_run1.gro")
MAPPING = json.loads((Path(__file__).parent / "data/mapping_residues.json").read_text())
POPC = MAPPING["POPC"]["beads"]
CHOL3 = list(mf.MARTINI3_MEMBRANE["CHOL"]["beads"])


def bead(resid, resname, atom, x=1.0, y=1.0, z=1.0):
    """One coarse-grained bead record."""
    return {"resid": resid, "resname": resname, "atom": atom, "chain": "", "x": x, "y": y, "z": z}


def residue(resid, resname, names):
    """One coarse-grained residue with the given bead names."""
    return [bead(resid, resname, n, 1 + 0.01 * resid, 1 + 0.01 * i, 1.0) for i, n in enumerate(names)]


def small_system():
    """A minimal valid Martini 3 system: three protein residues, POPC, cholesterol, water and an ion."""
    return (residue(1, "ALA", ["BB", "SC1"]) + residue(2, "TRP", ["BB", "SC1", "SC2", "SC3", "SC4", "SC5"])
            + residue(3, "GLY", ["BB"]) + residue(4, "POPC", POPC) + residue(5, "CHOL", CHOL3)
            + residue(6, "W", ["W"]) + residue(7, "ION", ["NA"]))


def failure(function, *args, **kwargs):
    """Return the SystemExit message a call fails with, or None when it does not fail."""
    try:
        function(*args, **kwargs)
    except SystemExit as exc:
        return str(exc)
    return None


def ppm_skip_reason():
    """Why the real PPM 3.0 cannot be used here (None when it can): not found, or no res.lib beside the executable."""
    try:
        mf.find_ppm_executable(None)
    except mf.OrientationFailure as error:
        return f"real PPM 3.0 unusable ({error}); set MEMBRANEFORGER_PPM=/path/to/immers, with res.lib in the same folder"
    return None


def gmx_skip_reason():
    """Why GROMACS cannot be run from Python (None when it can). A shell alias is invisible here: use a real executable."""
    gmx = os.environ.get("MEMBRANEFORGER_GMX") or "gmx"
    if shutil.which(gmx.split()[0]) is None:  # like find_gromacs: the command may carry arguments (singularity run IMAGE gmx)
        return (f"no gmx executable found ({gmx!r}; shell aliases are invisible to Python); put a real gmx on PATH "
                "or set MEMBRANEFORGER_GMX=/path/to/gmx")
    return None


def skip_if(reason_function):
    """Decorator: skip a test or class, printing the reason, when reason_function() returns one."""
    return lambda target: unittest.skipIf(reason_function() is not None, reason_function() or "")(target)
