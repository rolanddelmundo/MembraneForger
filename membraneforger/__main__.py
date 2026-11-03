#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// python -m membraneforger  (or the `membraneforger` command after pip install)
#//=============================================================
"""Module entry point: check that the Python requirements are importable, then run the command line."""
import importlib
import sys

__all__ = ['REQUIRED', 'missing_dependencies', 'run']

# Import name -> pip name, for the message when one is missing (pyproject.toml holds the versions).
REQUIRED = {"numpy": "numpy", "scipy": "scipy", "networkx": "networkx", "pandas": "pandas<3", "openmm": "openmm"}


def missing_dependencies(required: dict = REQUIRED) -> list[str]:
    """The pip names of the required packages this interpreter cannot import."""
    missing = []
    for module, pip_name in required.items():
        try:
            importlib.import_module(module)
        except ImportError:
            missing.append(pip_name)
    return missing


def run(argv: list | None = None) -> int:
    """Run the command line; a missing requirement gives the install command instead of a traceback."""
    missing = missing_dependencies()
    if missing:
        quoted = " ".join(f'"{m}"' if "<" in m or ">" in m else m for m in missing)
        print(f"membraneforger: {sys.executable} cannot import {', '.join(missing)}.\n"
              f"Install the requirements with\n    {sys.executable} -m pip install {quoted}\n"
              "or install the package itself (pip install . from the repository), which pulls them in.", file=sys.stderr)
        return 2
    from .cli import main
    return main(argv)


if __name__ == "__main__":
    sys.exit(run())
