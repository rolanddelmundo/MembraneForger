# Vendor Resource Policy

`resources/vendor/mstool/` is tracked as exact source for the MembraneForger
Stage 3 backmapping path.

MembraneForger resolves `mstool` from this repository-specific path by default.
`MEMBRANEFORGER_MSTOOL_ROOT` is reserved for explicit maintainer override tests.
Compiled native extensions are built from this source into
`${MEMBRANEFORGER_CACHE_DIR}` or the platform cache, never into this directory.

```bash
python scripts/bootstrap_mstool.py --verify
```

The top-level MembraneForger MIT license does not relicense `mstool`; the
vendored mstool source retains its upstream GPL-3.0-only license.
