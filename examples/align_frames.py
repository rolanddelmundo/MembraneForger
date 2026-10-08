"""Rotate Martini frames back onto their periodic box, in place, keeping the .gro formatting, velocities and names.

    python examples/align_frames.py examples/preeq_cg_cellmem/*_cg_cellmem.gro

A frame whose coordinates were rotated about z after the simulation (a rotational fit saved without its box) is not
periodic in the box it declares. membraneforger.align_frame_to_box finds the angle that makes it periodic again; this
script applies that rotation about the box centre to every coordinate and velocity line of the file and rewrites it
with the same title, atom order, box line and column formats, so the corrected file replaces the original one for
one. Frames that are already periodic are reported and left untouched.
"""
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import membraneforger as mf  # noqa: E402


def rotate_file(path: Path) -> dict:
    """Rotate one .gro file in place by the angle align_frame_to_box finds; returns that function's report."""
    atoms, box = mf.read_cg(path)
    _, report = mf.align_frame_to_box(atoms, box)
    angle = report["rotation_about_z_deg"]
    if not angle:
        return report
    c, s = math.cos(math.radians(angle)), math.sin(math.radians(angle))
    cx, cy = 0.5 * box[0], 0.5 * box[1]
    lines = path.read_text().splitlines()
    n = int(lines[1])
    out = lines[:2]
    for line in lines[2:2 + n]:
        x, y = float(line[20:28]), float(line[28:36])
        rx, ry = cx + (x - cx) * c - (y - cy) * s, cy + (x - cx) * s + (y - cy) * c
        new = f"{line[:20]}{rx:8.3f}{ry:8.3f}{line[36:44]}"
        if len(line) >= 68 and line[44:68].strip():  # velocities: rotate their x, y components the same way
            vx, vy = float(line[44:52]), float(line[52:60])
            new += f"{vx * c - vy * s:8.4f}{vx * s + vy * c:8.4f}{line[60:68]}{line[68:]}"
        else:
            new += line[44:]
        out.append(new)
    out += lines[2 + n:]
    path.write_text("\n".join(out) + "\n")
    checked, verified = mf.align_frame_to_box(*mf.read_cg(path))
    if verified["rotation_about_z_deg"] or verified["overlapping_pairs_as_read"] > report["overlapping_pairs_after"]:
        raise SystemExit(f"{path.name}: rewritten file is still not periodic ({verified})")
    return report


def main(paths: list[str]) -> int:
    for name in paths:
        report = rotate_file(Path(name))
        print(f"{Path(name).name}: rotated {report['rotation_about_z_deg']} degrees about z; overlapping bead pairs "
              f"{report['overlapping_pairs_as_read']} -> {report['overlapping_pairs_after']}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1:]))
