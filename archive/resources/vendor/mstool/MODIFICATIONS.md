# mstool Modifications

Vendored component: `mstool`

Upstream repository: `https://github.com/ksy141/mstool.git`

Base commit: `2d37f9d3e89279ddd9125cc74da1f5e01153586c`

License retained: GPL-3.0-only, with full text in `LICENSE`.

MembraneForger modification notice: modified by MembraneForger, 2026-09-02.

The vendored Python, Cython, XML, mapping, and data source files match the
upstream package tree at the base commit. The MembraneForger distribution adds
this modification/provenance notice and copies the upstream root `LICENSE` into
the vendored package directory so a clean source distribution carries the full
GPL-3.0-only text beside the imported package.

Local generated C extension outputs such as `lib/distancelib.c` and
`lib/qcprot.c`, native shared libraries, bytecode caches, and build directories
are ignored and must not be committed. MembraneForger builds platform-specific
mstool extensions from the vendored `.pyx` source into an external cache during
setup/doctor resolution.
