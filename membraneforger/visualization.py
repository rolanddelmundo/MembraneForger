#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// Publication renders of every pipeline stage: VMD scenes ray-traced by Tachyon (headless).
#// Runs in a python that can import vmd (plus numpy and PIL):
#//   python visualization.py RUN_DIR --out DIR --mode preview|publication|hero [--all-atom PDB] [--coarse-grain GRO]
#//=============================================================
"""Render the stages of one build with one fixed camera, lipids shown atom by atom, and record how to reproduce each image."""
import argparse
import hashlib
import json
import math
import os
import shutil
import struct
import subprocess
import sys
import time
import zlib
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

__all__ = ['RenderMode', 'MODES', 'RENDER', 'LIGHTS', 'MATERIALS', 'PALETTE', 'CHAIN_COLOURS', 'LIPID_TYPES', 'ELEMENT_COLOURS',
           'COLOUR_IDS', 'VIEWS', 'Scene', 'find_tachyon', 'detect_renderer_capabilities', 'protein_chains', 'chain_colours',
           'assign_secondary_structure',
           'draw_complex_whole', 'draw_sticks', 'frame_scale', 'tune_scene_file', 'read_ppm48', 'write_png16', 'ray_trace',
           'verify_png', 'stage_scenes', 'render_run', 'main']



def run_path(run: Path, name: str) -> Path:
    """A file of a build directory: at the top level, else in int/ (runtools.run_path; this module runs standalone under VMD)."""
    top, kept = Path(run) / name, Path(run) / "int" / name
    return kept if not top.exists() and kept.exists() else top

@dataclass(frozen=True)
class RenderMode:
    """Pixel size, sampling and geometry detail of one render quality level."""
    name: str
    long_axis: int           # pixels on the long image axis
    aasamples: int           # Tachyon antialiasing samples per pixel
    aosamples: int           # Tachyon ambient-occlusion samples
    cartoon_resolution: int  # facets around a cartoon tube/ribbon (sticks and spheres are exact quadrics in Tachyon)
    key_lights: int          # directional lights the key light is split into (a soft penumbra instead of a hard edge)
    master: bool             # keep the 16-bit master next to the 8-bit PNG


# The one place render quality is defined. Sampling was chosen where 1:1 test crops stopped changing (see the tutorial).
MODES = {"preview": RenderMode("preview", 1500, 4, 8, 24, 1, False),
         "publication": RenderMode("publication", 5400, 12, 32, 40, 8, True),
         "hero": RenderMode("hero", 8100, 16, 32, 50, 12, True)}
# The one place the camera, sky light and output of every scene are defined.
RENDER = {"projection": "Orthographic", "background_rgb": (1.0, 1.0, 1.0), "axes": "Off", "depthcue": "off",
          "ambient_occlusion": "on", "ao_ambient": 0.90, "ao_direct": 0.30, "shadows": "on",
          "key_light_spread_deg": 7.0,  # half-angle of the cone the key light is spread over
          "tachyon_flags": ("-trans_vmd", "-shadow_filter_off"),  # translucent beads must not darken what is under them
          "dof": False, "dof_fnumber": 64.0,  # depth of field exists only in perspective; never used for publication
          "crop": None,  # (x, y, size): render only this square window of the full image, pixel for pixel (test crops)
          "fill": 0.90, "dpi": 600, "renderer": "Tachyon (external binary, VMD scene file)"}
# Studio lights as (position the light shines from, in screen axes x right / y up / z towards the viewer; intensity).
# VMD cannot set a light's intensity, so the lights are written into the Tachyon scene file (tune_scene_file).
LIGHTS = {"key": ((-0.35, 0.55, 1.0), 1.0), "fill": ((0.9, -0.15, 0.6), 0.35), "rim": ((0.0, 0.9, -0.45), 0.25)}
# Every material is created explicitly; "ambient" is light that reaches a surface however deeply it is buried, which
# keeps the colours of lipids inside the bilayer cross-section readable under ambient occlusion.
MATERIALS = {
    "MFProtein": {"ambient": 0.22, "diffuse": 0.80, "specular": 0.35, "shininess": 0.35, "opacity": 1.0},  # satin: broad highlight
    "MFLipid": {"ambient": 0.30, "diffuse": 0.78, "specular": 0.25, "shininess": 0.30, "opacity": 1.0},    # softer sheen
    "MFStick": {"ambient": 0.20, "diffuse": 0.78, "specular": 0.40, "shininess": 0.65, "opacity": 1.0},    # crisp
    "MFBead": {"ambient": 0.30, "diffuse": 0.78, "specular": 0.15, "shininess": 0.35, "opacity": 1.0},
    "MFGhostBead": {"ambient": 0.30, "diffuse": 0.70, "specular": 0.0, "shininess": 0.40, "opacity": 0.35},
}
# The molecular palette (user, 2026-10-03): ten lipid species and five protein chains, fifteen distinct colours.
PALETTE = {"Cyber Celadon": "#78D6B4", "Aqua Jelly": "#55E6D6", "Frosted Sky": "#80CFF4", "Digital Periwinkle": "#8FA9FF",
           "Dream Mauve": "#B88AD8", "Orchid Pop": "#DC75C8", "Electric Pool": "#22CFE5", "Lemon Flash": "#F4E76E",
           "Bubblegum Vinyl": "#F28FBC", "Peach Pearl": "#F3A58F", "Ice Chrome Blue": "#5DA7E8", "Popstar Pink": "#EE5FAA",
           "Cyber Mint": "#51C99D", "Digital Lavender": "#9E7DE2", "Candy Coral": "#F17782"}
# Colours of things that are not one of the fifteen molecules; none of them repeats a palette colour.
EXTRA_COLOURS = {"ligand carbon": "#8A8F98", "water": "#C9D6DF", "sodium": "#5E4B8B", "chloride": "#3F8F4F",
                 "CG protein bead": "#B4BCC6", "box": "#404040"}
# Heteroatoms keep conventional element colours; these VMD colour ids are the ones its by-name colouring already uses.
ELEMENT_COLOURS = {"O": (1, "#E6281E"), "N": (0, "#304FF0"), "S": (4, "#FFD530"), "P": (5, "#F58220"), "H": (8, "#FFFFFF")}
rgb_of = lambda code: tuple(round(int(code[k:k + 2], 16) / 255, 4) for k in (1, 3, 5))
rgb = lambda name: rgb_of({**PALETTE, **EXTRA_COLOURS}[name])
# Dedicated VMD colour ids: palette on 17-31 in palette order, the extra colours on free low ids.
COLOUR_IDS = {**{name: 17 + k for k, name in enumerate(PALETTE)},
              **dict(zip(EXTRA_COLOURS, (2, 3, 6, 7, 9, 11)))}
# Protein chains take these colours in file order; a sixth chain is an error, never a repeated colour.
CHAIN_COLOURS = ["Ice Chrome Blue", "Popstar Pink", "Cyber Mint", "Digital Lavender", "Candy Coral"]
# species: (all-atom residue names incl. 4/5-character truncations, CG residue names, palette colour of the carbons)
LIPID_TYPES = {
    "POPC": ("POPC", "POPC", "Cyber Celadon"),
    "DOPC": ("DOPC", "DOPC", "Aqua Jelly"),
    "POPE": ("POPE", "POPE", "Frosted Sky"),
    "DOPE": ("DOPE", "DOPE", "Digital Periwinkle"),
    "POPS": ("POPS", "POPS", "Dream Mauve"),
    "DOPS": ("DOPS", "DOPS", "Orchid Pop"),
    "PSM": ("PSM", "PSM", "Electric Pool"),
    "cholesterol": ("CHL1", "CHOL", "Lemon Flash"),
    "PIP2": ("SAPI SAPI2 SAPI25", "SAP6", "Bubblegum Vinyl"),
    "GM3": ("GM3 GLPA CER1 CER16 CER160 BGLC BGAL ANE5 ANE5A ANE5AC", "GLC GAL NMC CER", "Peach Pearl"),
}
SAME_RESIDUE = {"HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
                "CYSG": "CYS", "CYSP": "CYS", "CYS2": "CYS", "CYX": "CYS"}  # one residue under different force-field names
PROTEIN = "(protein or resname HSD HSE HSP CYSG CYSP CYS2 AIB LEM KTZ KRT KSM)"
WATER, IONS = "resname TIP3 SOL W", "resname SOD CLA ION"
# Camera orientations as VMD rotate matrices: side = membrane plane horizontal with +z up; top = looking down the normal.
VIEWS = {"side": [[1, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0], [0, 0, 0, 1]],
         "top": [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]}
CG_BEAD_RADIUS_A = 2.35  # half the Martini regular-bead sigma (0.47 nm)
SLAB_HALF_A = 12.0  # half-thickness of the lipid cross-section drawn in side views


def check_palette() -> None:
    """Refuse to render when two things would share a colour or a colour id."""
    named = {**PALETTE, **EXTRA_COLOURS}
    codes = [c.upper() for c in named.values()] + [c for _, c in ELEMENT_COLOURS.values() if c != "#FFFFFF"]
    ids = list(COLOUR_IDS.values()) + [i for i, _ in ELEMENT_COLOURS.values()]
    used = [LIPID_TYPES[k][2] for k in LIPID_TYPES] + CHAIN_COLOURS
    if len(set(codes)) != len(codes) or len(set(ids)) != len(ids) or sorted(used) != sorted(PALETTE):
        raise SystemExit("palette error: every lipid species and protein chain needs its own colour and VMD colour id "
                         f"({len(set(codes))}/{len(codes)} distinct colours, {len(set(ids))}/{len(ids)} distinct ids)")


class Scene:
    """One VMD session whose every command is recorded so the image can be replayed from a Tcl script."""

    def __init__(self):
        """Start with an empty command log."""
        from vmd import evaltcl
        self._evaltcl, self.commands, self.representations = evaltcl, [], []
        self.chains, self.reference_chains = {}, None  # chains per loaded molecule; sequences that fix the chain colours

    def tcl(self, command: str) -> str:
        """Run one Tcl command in VMD and remember it."""
        self.commands.append(command)
        return self._evaltcl(command)

    def install_palette(self) -> None:
        """Give every palette, extra and element colour its own VMD colour id with an explicit RGB value."""
        check_palette()
        for name, colour_id in COLOUR_IDS.items():
            r, g, b = rgb(name)
            self.tcl(f"color change rgb {colour_id} {r} {g} {b}")
        for colour_id, code in ELEMENT_COLOURS.values():
            r, g, b = rgb_of(code)
            self.tcl(f"color change rgb {colour_id} {r} {g} {b}")

    def create_materials(self) -> None:
        """Create (or reset) the explicit materials of MATERIALS."""
        existing = self._evaltcl("material list").split()
        for name, settings in MATERIALS.items():
            if name not in existing:
                self.tcl(f"material add {name} copy AOChalky")
            for key, value in settings.items():
                self.tcl(f"material change {key} {name} {value}")

    def configure(self, width: int, height: int, mode: RenderMode) -> None:
        """Apply the shared settings: palette, materials, white background, projection, sky light, shadows, sampling."""
        self.install_palette()
        self.create_materials()
        projection = "Perspective" if RENDER["dof"] else RENDER["projection"]
        for command in (f"display resize {width} {height}", f"color Display Background {ELEMENT_COLOURS['H'][0]}",
                        f"display projection {projection}", f"display depthcue {RENDER['depthcue']}",
                        f"axes location {RENDER['axes']}", f"display ambientocclusion {RENDER['ambient_occlusion']}",
                        f"display aoambient {RENDER['ao_ambient']}", f"display aodirect {RENDER['ao_direct']}",
                        f"display shadows {RENDER['shadows']}", f"display dof {'on' if RENDER['dof'] else 'off'}",
                        f"display dof_fnumber {RENDER['dof_fnumber']}", "display dof_focaldist 2.0",
                        f"render aasamples Tachyon {mode.aasamples}", f"render aosamples Tachyon {mode.aosamples}"):
            self.tcl(command)
        for k, (position, _intensity) in enumerate(LIGHTS.values()):  # directions only; intensities go in the scene file
            self.tcl(f"light {k} on")
            self.tcl(f"light {k} pos {{{position[0]} {position[1]} {position[2]}}}")
        for k in range(len(LIGHTS), 4):
            self.tcl(f"light {k} off")

    def load(self, path: Path) -> int:
        """Load a structure with no default representation and return its molecule id."""
        molid = int(self.tcl(f"mol new {{{path}}} type {path.suffix[1:]} waitfor all"))
        self.tcl(f"mol delrep 0 {molid}")
        return molid

    def add(self, molid: int, selection: str, style: str, colour: str, material: str) -> int:
        """Add one representation in a named colour ("element" = by atom name) and material; returns its atom count."""
        atoms = int(self.tcl(f"[atomselect {molid} {{{selection}}}] num"))
        if not atoms:
            return 0
        colouring = "Name" if colour == "element" else f"ColorID {COLOUR_IDS[colour]}"
        for command in (f"mol representation {style}", f"mol color {colouring}",
                        f"mol selection {{{selection}}}", f"mol material {material}", f"mol addrep {molid}"):
            self.tcl(command)
        self.representations.append({"molecule": molid, "selection": selection, "style": style, "colour": colour,
                                     "material": material, "atoms": atoms})
        return atoms

    def camera(self, centre: tuple, view: str, scale: float) -> dict:
        """Set every molecule's centre, rotation and scale explicitly (never an automatic fit) and return the matrices."""
        fmt = lambda m: "{" + " ".join("{" + " ".join(f"{v:.8g}" for v in row) + "}" for row in m) + "}"
        centre_m = [[1, 0, 0, -centre[0]], [0, 1, 0, -centre[1]], [0, 0, 1, -centre[2]], [0, 0, 0, 1]]
        scale_m = [[scale, 0, 0, 0], [0, scale, 0, 0], [0, 0, scale, 0], [0, 0, 0, 1]]
        identity = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
        for molid in self.tcl("molinfo list").split():
            self.tcl(f"molinfo {molid} set {{center_matrix rotate_matrix scale_matrix global_matrix}} "
                     f"{{{fmt(centre_m)} {fmt(VIEWS[view])} {fmt(scale_m)} {fmt(identity)}}}")
        return {"center_matrix": centre_m, "rotate_matrix": VIEWS[view], "scale_matrix": scale_m, "global_matrix": identity}

    def clear(self) -> None:
        """Delete every molecule and forget the recorded commands of the finished image."""
        for molid in self.tcl("molinfo list").split():
            self.tcl(f"mol delete {molid}")
        self.commands, self.representations, self.chains = [], [], {}


def find_tachyon(requested: str | None) -> str:
    """Return the Tachyon ray tracer to use, or fail with what to install."""
    path = requested or os.environ.get("MEMBRANEFORGER_TACHYON") or shutil.which("tachyon")
    if not path or not Path(path).is_file():
        raise SystemExit("tachyon not found; pass --tachyon or set MEMBRANEFORGER_TACHYON (conda-forge package 'tachyon')")
    return path


def detect_renderer_capabilities(tachyon: str) -> dict:
    """Ask the installed VMD and Tachyon what they can do, print a short report, and fail if rendering is impossible."""
    import vmd
    from vmd import evaltcl, render
    usage = subprocess.run([tachyon], text=True, capture_output=True)
    usage = usage.stdout + usage.stderr
    formats = sorted({line.split()[1] for line in usage.splitlines() if line.strip().startswith("-format")})
    option = lambda flag: flag in usage
    version = next((line.split("Version")[1].strip() for line in usage.splitlines() if "Version" in line), "unknown")
    caps = {
        "vmd_python": getattr(vmd, "__version__", "unknown"), "vmd": evaltcl("vmdinfo version"),
        "renderers": render.listall(), "tachyon_external": "Tachyon" in render.listall(),
        "tachyon_internal": "TachyonInternal" in render.listall(),
        "accelerated_backends": [r for r in render.listall() if "OptiX" in r or "OSPRay" in r],
        "tachyon_binary": tachyon, "tachyon_version": version, "output_formats": formats,
        "ppm48": "PPM48" in formats, "psd48": "PSD48" in formats,
        "aa_samples": option("-aasamples"), "ao_samples_scene_file": True, "ao_command_line": option("-skylight_samples"),
        "shadows": option("-fullshade"), "transparency_modes": [m for m in ("-trans_orig", "-trans_raster3d", "-trans_vmd") if option(m)],
        "transparent_shadow_filter": option("-shadow_filter_off"),
        "light_intensity": "scene file only (VMD has no per-light intensity)", "area_lights": False,
        "depth_of_field": "perspective projection only (VMD writes Perspective_DoF; orthographic export ignores it)",
        "materials": evaltcl("material list").split(), "colour_ids": int(evaltcl("colorinfo num")),
    }
    missing = [flag for flag in RENDER["tachyon_flags"] + ("-aasamples", "-res", "-numthreads") if not option(flag)]
    if not caps["tachyon_external"] or missing or not ({"PPM48", "PNG"} & set(formats)):
        raise SystemExit(f"cannot render: VMD Tachyon export {caps['tachyon_external']}, tachyon options missing {missing}, "
                         f"output formats {formats}")
    print(f"renderer capabilities: vmd-python {caps['vmd_python']}, VMD {caps['vmd']}, external Tachyon {version} "
          f"(TachyonInternal {'yes' if caps['tachyon_internal'] else 'no'}, accelerated {caps['accelerated_backends'] or 'none'})\n"
          f"  formats {' '.join(formats)}; master {'PPM48 (16 bit per channel)' if caps['ppm48'] else 'PNG (8 bit, no deeper format)'}\n"
          f"  AA samples: command line; AO samples, AO ambient/direct, shadows, projection: VMD scene file; "
          f"light intensity and soft key light: scene file rewrite\n"
          f"  depth of field: {caps['depth_of_field']}", flush=True)
    return caps


def protein_chains(scene: Scene, molid: int) -> list:
    """Protein chains in file order as {"residues", "sequence"}, split where the backbone C-N bond is missing."""
    # VMD's own "fragment" follows guessed bonds, which in a .gro can fuse two chains that touch; the peptide bond
    # between consecutive residues is the reliable test.
    from vmd import atomsel
    if molid in scene.chains:
        return scene.chains[molid]
    backbone = atomsel(f"{PROTEIN} and name N CA C", molid=molid)
    cell = np.array([float(v) for v in scene._evaltcl(f"molinfo {molid} get {{a b c}}").split()])
    atoms = {}
    for residue, name, resname, x, y, z in zip(backbone.residue, backbone.name, backbone.resname, backbone.x, backbone.y, backbone.z):
        atoms.setdefault(residue, {"resname": resname})[name] = np.array([x, y, z])
    chains, previous = [], None
    for residue in sorted(atoms):
        joined = False
        if previous is not None and "C" in atoms[previous] and "N" in atoms[residue]:
            step = atoms[residue]["N"] - atoms[previous]["C"]
            if cell.min() > 0:
                step -= cell * np.round(step / cell)  # a chain written across the periodic boundary is still one chain
            joined = np.linalg.norm(step) < 3.0  # a strained input bond (2.4 A) is still a bond; a missing residue leaves 3.8 A or more
        if not joined:
            chains.append({"residues": [], "sequence": [], "backbone": []})
        chains[-1]["residues"].append(int(residue))
        chains[-1]["sequence"].append(SAME_RESIDUE.get(atoms[residue]["resname"], atoms[residue]["resname"]))
        chains[-1]["backbone"].append(atoms[residue])
        previous = residue
    scene.chains[molid] = chains
    return chains


def chain_colours(scene: Scene, chains: list) -> list:
    """Palette colour of each chain: the colour of the reference chain with the same sequence, else file order."""
    import difflib
    if len(chains) > len(CHAIN_COLOURS) or not chains:
        raise SystemExit(f"{len(chains)} protein chains found; the palette has colours for 1 to {len(CHAIN_COLOURS)}")
    if not scene.reference_chains:
        return CHAIN_COLOURS[:len(chains)]
    scores = sorted(((difflib.SequenceMatcher(None, chain["sequence"], ref, autojunk=False).ratio(), i, j)
                     for i, chain in enumerate(chains) for j, ref in enumerate(scene.reference_chains)), reverse=True)
    chosen = {}
    for score, i, j in scores:  # best matches first; each reference chain gives its colour to one chain
        if i not in chosen and j not in chosen.values():
            if score < 0.8:
                raise SystemExit(f"protein chain {i + 1} matches no chain of the reference structure (best {score:.2f}); "
                                 "its colour would be a guess")
            chosen[i] = j
    return [CHAIN_COLOURS[chosen[i]] for i in range(len(chains))]


def assign_secondary_structure(scene: Scene, molid: int) -> dict:
    """Mark helices and strands from backbone dihedrals, because this VMD build has no STRIDE to do it."""
    def dihedral(p0, p1, p2, p3):
        """Signed dihedral angles (degrees) for arrays of four points."""
        b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
        b1 = b1 / np.linalg.norm(b1, axis=1)[:, None]
        v, w = b0 - (b0 * b1).sum(1)[:, None] * b1, b2 - (b2 * b1).sum(1)[:, None] * b1
        return np.degrees(np.arctan2((np.cross(b1, v) * w).sum(1), (v * w).sum(1)))

    found = {"H": 0, "E": 0}
    for chain in protein_chains(scene, molid):
        complete = [k for k, atoms in enumerate(chain["backbone"]) if all(name in atoms for name in ("N", "CA", "C"))]
        if len(complete) != len(chain["residues"]) or len(complete) < 4:
            continue  # a residue without a full backbone: leave this chain as coil rather than guess
        N, CA, C = (np.array([atoms[name] for atoms in chain["backbone"]]) for name in ("N", "CA", "C"))
        residue = np.array(chain["residues"])
        phi, psi = np.full(len(CA), np.nan), np.full(len(CA), np.nan)
        phi[1:] = dihedral(C[:-1], N[1:], CA[1:], C[1:])
        psi[:-1] = dihedral(N[:-1], CA[:-1], C[:-1], N[1:])
        helix = (phi > -100) & (phi < -30) & (psi > -80) & (psi < -5)
        strand = (phi > -180) & (phi < -45) & ((psi > 90) | (psi < -150))
        for code, flags, shortest in (("H", helix, 4), ("E", strand, 3)):
            chosen, start = [], None
            for i, flag in enumerate(list(flags) + [False]):
                if flag and start is None:
                    start = i
                elif start is not None and not flag:
                    if i - start >= shortest:
                        chosen += residue[start:i].tolist()
                    start = None
            if chosen:
                scene.tcl(f"[atomselect {molid} {{residue {' '.join(map(str, chosen))}}}] set structure {code}")
            found[code] += len(chosen)
    return found


def draw_protein(scene: Scene, molid: int, mode: RenderMode) -> None:
    """Cartoon for each protein chain in its own palette colour; a chain keeps its colour in every stage."""
    found = assign_secondary_structure(scene, molid)
    chains = protein_chains(scene, molid)
    if not found["H"] + found["E"]:
        raise SystemExit("no helix or strand was recognised in the protein; the cartoon would be a bare tube")
    for chain, colour in zip(chains, chain_colours(scene, chains)):
        scene.add(molid, f"{PROTEIN} and residue {chain['residues'][0]} to {chain['residues'][-1]}",
                  f"NewCartoon 0.38 {mode.cartoon_resolution} 4.5 0", colour, "MFProtein")


def draw_sticks(scene: Scene, molid: int, selection: str, carbon: str, radius: float, material: str) -> int:
    """Sticks coloured by element: the whole skeleton in the molecule's carbon colour, non-carbon atoms in element colours."""
    # VMD draws a bond only when both atoms are in the representation and cannot take a per-molecule carbon colour in
    # its by-element mode, so the skeleton is drawn once in the carbon colour and the heteroatoms again on top.
    atoms = scene.add(molid, f"({selection}) and noh", f"Licorice {radius} 16 16", carbon, material)
    scene.add(molid, f'({selection}) and noh and not name "C.*"', f"Licorice {radius + 0.04} 16 16", "element", material)
    return atoms


def draw_lipids(scene: Scene, molid: int, within: str, universe: str, radius: float = 0.28) -> None:
    """Every heavy atom of the selected lipids as sticks coloured by element, carbon in the lipid species' palette colour."""
    all_names = " ".join(names for names, _, _ in LIPID_TYPES.values())
    unknown = scene._evaltcl(f"lsort -unique [[atomselect {molid} {{({universe}) and not resname {all_names}}}] get resname]")
    if unknown.strip():
        raise SystemExit(f"unknown lipid species {unknown.split()}: add them to LIPID_TYPES with their own colour")
    drawn = sum(draw_sticks(scene, molid, f"resname {names} and ({within})", colour, radius, "MFLipid")
                for names, _cg, colour in LIPID_TYPES.values())
    if not drawn:
        raise SystemExit(f"no lipid atoms selected by: {within}")
    scene.add(molid, f"resname {all_names} and name P and ({within})", "VDW 0.55 24", "element", "MFLipid")


def draw_cg(scene: Scene, molid: int, within: str, material: str, protein: bool = True) -> None:
    """Martini beads as spheres: lipids coloured by the same lipid types, protein backbone beads in one colour."""
    # A CG file has no bonds, so VMD's "residue" is a single bead; whole molecules are selected by residue number.
    within = within.replace("same residue as", "same resid as")
    scene.tcl(f"[atomselect {molid} all] set radius {CG_BEAD_RADIUS_A}")
    size = 0.8 if protein else 0.45  # full-size beads on their own, small translucent ones under all-atom lipids
    for _aa, names, colour in LIPID_TYPES.values():
        scene.add(molid, f"resname {names} and not name V and ({within})", f"VDW {size} 20", colour, material)
    if protein and not scene.add(molid, "name BB", f"VDW {size} 20", "CG protein bead", material):
        raise SystemExit("no protein BB beads in the coarse-grained file")


def draw_box(scene: Scene, molid: int) -> None:
    """Draw the periodic cell as twelve thin dark cylinders."""
    a, b, c = (float(v) for v in scene._evaltcl(f"molinfo {molid} get {{a b c}}").split())
    if min(a, b, c) <= 0:
        return
    corners = [(x, y, z) for x in (0, a) for y in (0, b) for z in (0, c)]
    scene.tcl(f"graphics {molid} color {COLOUR_IDS['box']}")
    scene.tcl(f"graphics {molid} material Diffuse")
    for i, p in enumerate(corners):
        for q in corners[i + 1:]:
            if sum(abs(u - v) > 1e-6 for u, v in zip(p, q)) == 1:
                scene.tcl(f"graphics {molid} cylinder {{{p[0]} {p[1]} {p[2]}}} {{{q[0]} {q[1]} {q[2]}}} radius 0.35 resolution 12")


def draw_complex_whole(scene: Scene, molid: int, solute: str) -> int:
    """Shift each protein chain and ligand residue by whole box vectors next to the first chain; returns how many moved."""
    # mdrun may write a molecule in a neighbouring periodic image. That is the same structure, but a picture of it
    # shows a chain floating a box length away, so every image draws the complex in one piece.
    a, b, c = (float(v) for v in scene._evaltcl(f"molinfo {molid} get {{a b c}}").split())
    chains = [f"residue {c['residues'][0]} to {c['residues'][-1]}" for c in protein_chains(scene, molid)]
    if min(a, b, c) <= 0 or not chains:
        return 0
    centre = lambda selection: np.array([float(v) for v in scene._evaltcl(f"measure center [atomselect {molid} {{{selection}}}]").split()])
    reference, cell, moved = centre(chains[0]), np.array([a, b, c]), 0
    others = scene._evaltcl(f"lsort -unique -integer [[atomselect {molid} {{({solute}) and not {PROTEIN}}}] get residue]").split()
    for selection in chains[1:] + [f"residue {r}" for r in others]:
        shift = -cell * np.round((centre(selection) - reference) / cell)
        if np.abs(shift).max() > 1e-6:
            scene.tcl(f"[atomselect {molid} {{{selection}}}] moveby {{{shift[0]:.4f} {shift[1]:.4f} {shift[2]:.4f}}}")
            moved += 1
    return moved


def membrane_centre(scene: Scene, molid: int) -> tuple | None:
    """Centre of the cell in XY and of the lipid phosphates (or PO4 beads) in z; None when there is no membrane."""
    a, b, _c = (float(v) for v in scene._evaltcl(f"molinfo {molid} get {{a b c}}").split())
    lipid_names = " ".join(f"{aa} {cg}" for aa, cg, _ in LIPID_TYPES.values())
    z = scene._evaltcl(f"[atomselect {molid} {{resname {lipid_names} and name P PO4}}] get z").split()
    if not z or a <= 0:
        return None
    return a / 2, b / 2, float(np.mean([float(v) for v in z]))


def solute_centre(scene: Scene, molid: int) -> tuple:
    """Geometric centre of the protein CA atoms (or BB beads)."""
    xyz = scene._evaltcl(f"measure center [atomselect {molid} {{name CA BB}}]").split()
    return tuple(float(v) for v in xyz)


def tune_scene_file(scene_file: Path, mode: RenderMode, width: int, height: int) -> list:
    """Replace VMD's lights in the Tachyon scene file with the studio lights of LIGHTS; returns the lights written."""
    # Tachyon 0.99 has only point-like directional lights, whose shadows have a hard edge. Spreading the key light
    # over a small cone of weaker lights gives the edge a penumbra; total key intensity is unchanged.
    def unit(v):
        n = math.sqrt(sum(c * c for c in v))
        return [c / n for c in v]

    lights = []
    for k, (name, (position, intensity)) in enumerate(LIGHTS.items()):
        p, count = unit(position), (mode.key_lights if k == 0 else 1)
        a = unit([p[1], -p[0], 0.0] if abs(p[2]) < 0.9 else [1.0, 0.0, 0.0])
        b = [p[1] * a[2] - p[2] * a[1], p[2] * a[0] - p[0] * a[2], p[0] * a[1] - p[1] * a[0]]
        spread = math.tan(math.radians(RENDER["key_light_spread_deg"])) if count > 1 else 0.0
        for i in range(count):
            t = 2 * math.pi * i / count
            towards = unit([p[j] + spread * (math.cos(t) * a[j] + math.sin(t) * b[j]) for j in range(3)])
            lights.append({"light": name, "from": [round(c, 5) for c in towards], "intensity": round(intensity / count, 5)})
    rewritten, header, found = scene_file.with_suffix(".tuned.dat"), True, 0
    with open(scene_file) as source, open(rewritten, "w") as target:
        for line in source:
            if header and line.startswith("Directional_Light"):
                found += 1
                continue
            if header and RENDER["crop"] and line.startswith("  Zoom"):  # orthographic: pixels per unit = zoom * height
                line = f"  Zoom {float(line.split()[1]) * height / RENDER['crop'][2]:.8g}\n"
            if header and RENDER["crop"] and line.startswith("  Center"):
                unit_px = 0.5 * height  # VMD writes Zoom 0.5
                x, y = (RENDER["crop"][0] - width / 2) / unit_px, -(RENDER["crop"][1] - height / 2) / unit_px
                line = f"  Center  {x:.8g} {y:.8g} -2\n"
            if header and line.startswith("Background"):
                for light in lights:
                    x, y, z = (-c for c in light["from"])  # Tachyon wants the direction the light travels
                    i = light["intensity"]
                    target.write(f"Directional_Light Direction {x:.5f} {y:.5f} {z:.5f} Color {i} {i} {i}\n")
                header = False
            target.write(line)
    if header or not found:
        raise SystemExit(f"{scene_file.name}: not the Tachyon scene layout this code knows (lights {found}, Background missing {header})")
    rewritten.replace(scene_file)
    return lights


def read_ppm48(path: Path) -> np.ndarray:
    """Read Tachyon's 48-bit PPM as a (height, width, 3) uint16 array."""
    data = path.read_bytes()
    magic, size, maxval, pixels = data.split(b"\n", 3)
    width, height = (int(v) for v in size.split())
    if magic != b"P6" or int(maxval) != 65535 or len(pixels) != width * height * 6:
        raise SystemExit(f"{path.name}: not a 48-bit PPM ({magic!r}, maxval {maxval!r}, {len(pixels)} bytes)")
    return np.frombuffer(pixels, dtype=">u2").reshape(height, width, 3)


def write_png16(path: Path, pixels: np.ndarray, dpi: int) -> None:
    """Write a (height, width, 3) uint16 array as a 16-bit RGB PNG (PIL cannot), lossless, with the print density."""
    height, width, _ = pixels.shape
    raw = pixels.astype(">u2").view(np.uint8).reshape(height, width * 6)
    filtered = np.empty((height, width * 6 + 1), dtype=np.uint8)
    filtered[:, 0] = 1  # PNG "Sub" filter: each byte minus the same byte of the pixel to its left
    filtered[:, 1:7] = raw[:, :6]
    filtered[:, 7:] = raw[:, 6:] - raw[:, :-6]
    chunk = lambda kind, body: struct.pack(">I", len(body)) + kind + body + struct.pack(">I", zlib.crc32(kind + body))
    per_metre = round(dpi / 0.0254)
    path.write_bytes(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 16, 2, 0, 0, 0))
                     + chunk(b"pHYs", struct.pack(">IIB", per_metre, per_metre, 1))
                     + chunk(b"IDAT", zlib.compress(filtered.tobytes(), 6)) + chunk(b"IEND", b""))


def ray_trace(scene: Scene, png: Path, tachyon: str, width: int, height: int, mode: RenderMode, threads: int,
              caps: dict | None = None) -> dict:
    """Export the VMD scene, set its lights, ray-trace it, and write the 8-bit PNG (plus the 16-bit master it came from)."""
    from PIL import Image
    from vmd import display, render
    scene_file = png.with_suffix(".dat")
    display.update()  # headless VMD only builds geometry on an explicit update
    render.render("Tachyon", str(scene_file))
    lights = tune_scene_file(scene_file, mode, width, height)
    if RENDER["crop"]:
        width = height = RENDER["crop"][2]
    deep = bool(caps and caps["ppm48"])  # a real 16-bit render when Tachyon offers one, otherwise its native 8-bit PNG
    raw = png.with_suffix(".ppm") if deep else png
    command = [tachyon, str(scene_file), "-format", "PPM48" if deep else "PNG", "-o", str(raw), "-res", str(width), str(height),
               "-aasamples", str(mode.aasamples), "-numthreads", str(threads), *RENDER["tachyon_flags"]]
    started = time.time()
    result = subprocess.run(command, text=True, capture_output=True)
    scene_file.unlink(missing_ok=True)
    if result.returncode or not raw.is_file() or raw.stat().st_size == 0:
        raise SystemExit(f"tachyon failed for {png.name}: {(result.stdout + result.stderr)[-400:]}")
    record = {"lights": lights, "tachyon_command": " ".join(command), "seconds": round(time.time() - started, 1),
              "source_bit_depth": 16 if deep else 8, "master": None}
    if deep:
        pixels = read_ppm48(raw)
        raw.unlink()
        if mode.master:
            master = png.parent / "masters" / (png.stem + ".master16.png")
            master.parent.mkdir(exist_ok=True)
            write_png16(master, pixels, RENDER["dpi"])
            record["master"] = str(master)
        eight = ((pixels.astype(np.uint32) + 128) // 257).astype(np.uint8)  # exact 16 -> 8 bit rounding, no tone change
        Image.fromarray(eight, "RGB").save(png, format="PNG", dpi=(RENDER["dpi"], RENDER["dpi"]), optimize=True)
    else:
        with Image.open(png) as image:
            image.load()
            image.save(png, format="PNG", dpi=(RENDER["dpi"], RENDER["dpi"]), optimize=True)  # lossless; adds density only
    scene.commands.append("# ray traced with: " + " ".join(command))
    return record


def verify_png(png: Path, width: int, height: int) -> dict:
    """Decode the image and check size, content, framing and tone (no clipped colours, no crushed blacks)."""
    from PIL import Image
    with Image.open(png) as image:
        image.load()
        pixels = np.asarray(image.convert("RGB"))
        fmt, size = image.format, image.size
    drawn = np.argwhere((pixels < 250).any(axis=2))
    if fmt != "PNG" or size != (width, height) or not len(drawn):
        raise SystemExit(f"{png.name}: format {fmt}, size {size}, expected PNG {width}x{height} with content")
    (top, left), (bottom, right) = drawn.min(axis=0), drawn.max(axis=0)
    clipped = bool(top == 0 or left == 0 or bottom == height - 1 or right == width - 1)
    corner = pixels[:8, :8].reshape(-1, 3).mean(axis=0)
    molecule = pixels[(pixels < 250).any(axis=2)]
    blown = float(((molecule == 255).any(axis=1)).mean())  # a saturated channel inside a coloured pixel
    crushed = float((molecule.max(axis=1) < 16).mean())
    return {"format": fmt, "width": size[0], "height": size[1], "bytes": png.stat().st_size,
            "content_fraction": [round(float((right - left + 1) / width), 3), round(float((bottom - top + 1) / height), 3)],
            "touches_edge": clipped, "background_rgb": [int(v) for v in corner],
            "clipped_highlight_fraction": round(blown, 5), "crushed_black_fraction": round(crushed, 5),
            "median_molecule_rgb": [int(v) for v in np.median(molecule, axis=0)]}


def frame_scale(scene: Scene, png_dir: Path, tachyon: str, centre: tuple, view: str, width: int, height: int, threads: int) -> float:
    """Measure, on a small test render, the scale at which the drawn scene fills RENDER['fill'] of the image."""
    from PIL import Image
    probe = png_dir / f"_frame_probe_{os.getpid()}.png"
    small = (max(int(width / 6), 200), max(int(height / 6), 200))
    scene.tcl("display ambientocclusion off")
    scene.tcl("display shadows off")
    trial = 0.02
    for _ in range(12):  # shrink until the whole scene is inside the test image, then measure it
        scene.camera(centre, view, trial)
        ray_trace(scene, probe, tachyon, small[0], small[1], RenderMode("probe", max(small), 1, 1, 12, 1, False), threads)
        with Image.open(probe) as image:
            pixels = np.asarray(image.convert("RGB"))
        drawn = np.argwhere((pixels < 250).any(axis=2))
        inside = len(drawn) and drawn[:, 0].min() > 1 and drawn[:, 1].min() > 1 \
            and drawn[:, 0].max() < small[1] - 2 and drawn[:, 1].max() < small[0] - 2
        if inside:
            break
        trial *= 0.6
    else:
        raise SystemExit("could not frame the scene: it never fitted inside the test image")
    probe.unlink()
    scene.tcl(f"display ambientocclusion {RENDER['ambient_occlusion']}")
    scene.tcl(f"display shadows {RENDER['shadows']}")
    # extent measured about the image centre so the scene stays centred on `centre`
    half_h = np.abs(drawn[:, 0] - small[1] / 2).max() / (small[1] / 2)
    half_w = np.abs(drawn[:, 1] - small[0] / 2).max() / (small[0] / 2)
    return trial * RENDER["fill"] / max(half_h, half_w)


def stage_scenes(run: Path, all_atom: Path | None, coarse_grain: Path | None, solute_atoms: int) -> list[dict]:
    """Describe every image: which file, which view, and what is drawn."""
    slab = lambda yc: f"same residue as (abs(y - {yc:.2f}) < {SLAB_HALF_A})"
    solute = f"index < {solute_atoms}"
    scenes = [
        {"stage": "input_aa", "file": all_atom, "view": "side", "frame": "own", "draw": ["protein", "ligands"],
         "shows": "all-atom input complex in its own coordinate frame"},
        {"stage": "input_cg", "file": coarse_grain, "view": "side", "frame": "membrane", "draw": ["cg_slab", "box"],
         "shows": "Martini 3 input: protein BB beads and a 2.4 nm lipid cross-section, beads coloured by lipid type"},
        {"stage": "aligned_complex", "file": run_path(run, "membrane.pdb"), "view": "side", "frame": "membrane",
         "draw": ["protein", "ligands_pdb", "cg_protein_overlay"],
         "shows": "all-atom complex after the rigid fit, over the CG backbone beads it was fitted to (lipids hidden)"},
        {"stage": "backmapped_membrane", "file": run_path(run, "membrane.pdb"), "view": "side", "frame": "membrane",
         "draw": ["protein", "lipid_slab"], "shows": "backmapped all-atom lipids (2.4 nm cross-section) around the placed complex"},
        {"stage": "backmapped_membrane_top", "file": run_path(run, "membrane.pdb"), "view": "top", "frame": "membrane",
         "draw": ["protein", "lipids_all", "box"], "shows": "all backmapped lipids seen along the membrane normal"},
        {"stage": "backmapping_overlay", "file": run_path(run, "membrane.pdb"), "view": "side", "frame": "patch",
         "draw": ["lipid_patch", "cg_patch_overlay"],
         "shows": "close-up of a lipid patch: all-atom lipids with the Martini beads they were built from (translucent)"},
        {"stage": "boxed", "file": run_path(run, "boxed.gro"), "view": "side", "frame": "membrane", "draw": ["protein", "lipid_slab", "box"],
         "shows": "topology-ordered system in the rebuilt box"},
        {"stage": "solvated", "file": run_path(run, "solv.gro"), "view": "side", "frame": "membrane",
         "draw": ["protein", "lipid_slab", "water_slab", "box"],
         "shows": "water added (2.4 nm cross-section): the lipid core holds no water except cavity waters within 1.0 nm of the protein"},
        {"stage": "ionized", "file": run_path(run, "solv_ions.gro"), "view": "side", "frame": "membrane",
         "draw": ["protein", "lipid_slab", "ions", "box"], "shows": "0.15 M NaCl plus neutralizing ions (all ions shown)"},
        {"stage": "pre_em", "file": run_path(run, "solv_ions.gro"), "view": "side", "frame": "membrane", "draw": ["protein", "lipid_slab"],
         "shows": "lipid cross-section before energy minimization"},
        {"stage": "post_em", "file": run / "em.gro", "view": "side", "frame": "membrane", "draw": ["protein", "lipid_slab"],
         "shows": "the same cross-section after energy minimization"},
        {"stage": "post_em_patch", "file": run / "em.gro", "view": "side", "frame": "patch", "draw": ["lipid_patch"],
         "shows": "close-up of the same lipid patch after energy minimization"},
        {"stage": "final_validated_top", "file": run / "em.gro", "view": "top", "frame": "membrane",
         "draw": ["protein", "lipids_all", "box"], "shows": "validated system seen along the membrane normal"},
        {"stage": "final_validated", "file": run / "em.gro", "view": "side", "frame": "membrane",
         "draw": ["protein", "lipid_slab", "ions", "box"], "shows": "validated, energy-minimized system"},
    ]
    for item in scenes:
        item["slab"], item["solute"] = slab, solute
    return [s for s in scenes if s["file"] is not None and Path(s["file"]).is_file()]


def compose(scene: Scene, item: dict, coarse_grain: Path | None, shift_z: float, mode: RenderMode) -> tuple:
    """Load the structures of one image, draw what the stage asks for, and return the centre to look at."""
    molid = scene.load(Path(item["file"]))
    if Path(item["file"]).suffix == ".gro" and Path(item["file"]) != coarse_grain:
        draw_complex_whole(scene, molid, item["solute"])
    membrane = membrane_centre(scene, molid)
    centre = membrane if item["frame"] != "own" and membrane else solute_centre(scene, molid)
    slab = item["slab"](centre[1])
    # a lipid patch away from the protein: the corner region of the cell
    patch = (f"same residue as (name P PO4 ROH O3 and x > {0.04 * 2 * centre[0]:.2f} and x < {0.36 * 2 * centre[0]:.2f} "
             f"and abs(y - {0.2 * 2 * centre[1]:.2f}) < 9)")
    # everything that is neither solute, water nor ion must be a lipid species the palette knows
    pdb = Path(item["file"]).suffix == ".pdb"
    universe = "segname MEMB" if pdb else f"not ({item['solute']}) and not {WATER} and not {IONS}"
    for what in item["draw"]:
        if what == "protein":
            draw_protein(scene, molid, mode)
        elif what == "ligands":
            draw_sticks(scene, molid, f"not {PROTEIN}", "ligand carbon", 0.45, "MFStick")
        elif what == "ligands_pdb":
            draw_sticks(scene, molid, f"not {PROTEIN} and not segname MEMB", "ligand carbon", 0.45, "MFStick")
        elif what == "lipid_slab":
            draw_lipids(scene, molid, slab, universe)
        elif what == "lipids_all":
            draw_lipids(scene, molid, "all", universe)
        elif what == "lipid_patch":
            draw_lipids(scene, molid, patch, universe, radius=0.22)
        elif what == "cg_slab":
            draw_cg(scene, molid, slab, "MFBead")
        elif what in ("cg_protein_overlay", "cg_patch_overlay") and coarse_grain:
            cg = scene.load(coarse_grain)
            scene.tcl(f"[atomselect {cg} all] moveby {{0 0 {shift_z:.3f}}}")
            if what == "cg_protein_overlay":
                scene.tcl(f"[atomselect {cg} all] set radius {CG_BEAD_RADIUS_A}")
                scene.add(cg, "name BB", "VDW 0.6 20", "CG protein bead", "MFGhostBead")
            else:
                draw_cg(scene, cg, patch, "MFGhostBead", protein=False)
        elif what == "water_slab":
            scene.add(molid, f"name OH2 and abs(y - {centre[1]:.2f}) < {SLAB_HALF_A}", "VDW 0.35 12", "water", "MFBead")
        elif what == "ions":
            scene.add(molid, "resname SOD", "VDW 0.6 20", "sodium", "MFStick")
            scene.add(molid, "resname CLA", "VDW 0.6 20", "chloride", "MFStick")
        elif what == "box":
            draw_box(scene, molid)
    if item["frame"] == "patch":
        centre = (0.2 * 2 * centre[0], 0.2 * 2 * centre[1], centre[2])
    return centre


MODE_OVERRIDES = {}  # sampling changes requested with --set MODE.<field>=<n>


def apply_overrides(settings: list) -> dict:
    """Change central settings from the command line (tuning runs), e.g. RENDER.ao_ambient=0.8 LIGHTS.fill=0.5."""
    applied = {}
    for setting in settings or []:
        target, value = setting.split("=", 1)
        table, *keys = target.split(".")
        if table == "RENDER" and len(keys) == 1 and keys[0] in RENDER:
            RENDER[keys[0]] = type(RENDER[keys[0]])(float(value)) if not isinstance(RENDER[keys[0]], str) else value
        elif table == "LIGHTS" and len(keys) == 1 and keys[0] in LIGHTS:
            LIGHTS[keys[0]] = (LIGHTS[keys[0]][0], float(value))
        elif table == "MODE" and len(keys) == 1 and keys[0] in ("aasamples", "aosamples", "cartoon_resolution", "key_lights"):
            MODE_OVERRIDES[keys[0]] = int(value)
        elif table == "MATERIALS" and len(keys) == 2 and keys[1] in MATERIALS.get(keys[0], {}):
            MATERIALS[keys[0]][keys[1]] = float(value)
        else:
            raise SystemExit(f"unknown setting {target}; use RENDER.<key>, MODE.<field>, LIGHTS.<name> (intensity) or MATERIALS.<name>.<key>")
        applied[target] = value
    return applied


def source_commit() -> str | None:
    """Git commit of this source tree, or None when it is not a git checkout."""
    result = subprocess.run(["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "HEAD"], text=True, capture_output=True)
    return result.stdout.strip() if result.returncode == 0 else None


def render_run(run: Path, out: Path, tachyon: str, mode: RenderMode, all_atom: Path | None, coarse_grain: Path | None,
               threads: int, only: list | None = None, index: int | None = None, camera_from: Path | None = None,
               overrides: dict | None = None) -> list:
    """Render every stage of one build with a camera measured once per frame type, and return the manifest rows."""
    from vmd import evaltcl
    caps = detect_renderer_capabilities(tachyon)
    manifest = json.loads((run / "run_manifest.json").read_text()) if (run / "run_manifest.json").is_file() else {}
    solute_atoms = (manifest.get("topology") or {}).get("index_groups", {}).get("Protein_LIG", 0)
    out.mkdir(parents=True, exist_ok=True)
    scene, rows, scales = Scene(), [], {}
    everything = stage_scenes(run, all_atom, coarse_grain, solute_atoms)
    # the reference for each camera is the last (largest) stage that uses it, so no earlier stage can be clipped;
    # it is taken from the full stage list so a partial re-render keeps the same camera and file numbers
    reference = {(i["view"], i["frame"]): i for i in everything}
    if run_path(run, "membrane.pdb").is_file():  # its chain order fixes which chain gets which colour in every stage
        scene.reference_chains = [c["sequence"] for c in protein_chains(scene, scene.load(run_path(run, "membrane.pdb")))]
        scene.clear()
    for number, item in enumerate(everything, 1):
        if (only and item["stage"] not in only) or (index is not None and number - 1 != index):
            continue
        key = (item["view"], item["frame"])
        aspect = {"side": 0.75, "top": 1.0}[item["view"]] if item["frame"] != "patch" else 1.25  # width / height
        width, height = (int(mode.long_axis * aspect), mode.long_axis) if aspect <= 1 else (mode.long_axis, int(mode.long_axis / aspect))
        scene.configure(width // 6, height // 6, mode)
        earlier = camera_from / f"{number:03d}_{item['stage']}.render.json" if camera_from else None
        if earlier is not None:  # the camera of an earlier render, so the two images can be compared pixel for pixel
            if not earlier.is_file():
                raise SystemExit(f"--camera-from: {earlier} not found")
            scales[key] = json.loads(earlier.read_text())["scale"]
        elif key not in scales:
            centre = compose(scene, reference[key], coarse_grain, 0.0, mode)
            scales[key] = frame_scale(scene, out, tachyon, centre, item["view"], width, height, threads)
            scene.clear()
            scene.configure(width // 6, height // 6, mode)
        centre = compose(scene, item, coarse_grain, 0.0, mode)
        camera = scene.camera(centre, item["view"], scales[key])
        png = out / f"{number:03d}_{item['stage']}.png"
        traced = ray_trace(scene, png, tachyon, width, height, mode, threads, caps)
        checked = verify_png(png, *((RENDER["crop"][2],) * 2 if RENDER["crop"] else (width, height)))
        script = png.with_suffix(".vmd.tcl")
        script.write_text("# VMD scene that reproduces " + png.name + "\n" + "\n".join(scene.commands) + "\n")
        clipping = {k: float(evaltcl(f"display get {k}")) for k in ("nearclip", "farclip")}
        rows.append({"image": png.name, "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "category": "pipeline",
                     "stage_or_function": item["stage"], "shows": item["shows"], "view": item["view"],
                     "command_or_action": f"VMD {evaltcl('vmdinfo version')} scene -> Tachyon {mode.name} "
                                          f"{width}x{height}, aa {mode.aasamples}, ao {mode.aosamples}",
                     "source_input": str(item["file"]), "frame": 0, "output_artifact": str(png), "scene_script": script.name,
                     "master_16bit": traced["master"], "source_bit_depth": traced["source_bit_depth"], "output_format": "PNG (lossless RGB)",
                     "backend": {"renderer": RENDER["renderer"], "vmd": caps["vmd"], "vmd_python": caps["vmd_python"],
                                 "tachyon": caps["tachyon_version"], "tachyon_binary": caps["tachyon_binary"],
                                 "tachyon_command": traced["tachyon_command"], "render_seconds": traced["seconds"]},
                     "camera": camera, "centre_A": [round(v, 3) for v in centre], "scale": scales[key],
                     "camera_reused_from": str(camera_from) if camera_from else None, "clipping_planes": clipping,
                     "render_settings": {**RENDER, **asdict(mode)}, "lights": traced["lights"], "materials": MATERIALS,
                     "palette": {name: {"hex": code, "colour_id": COLOUR_IDS[name]} for name, code in {**PALETTE, **EXTRA_COLOURS}.items()},
                     "palette_assignment": {"lipids": {k: v[2] for k, v in LIPID_TYPES.items()}, "protein_chains_in_file_order": CHAIN_COLOURS,
                                            "elements": {e: {"colour_id": i, "hex": c} for e, (i, c) in ELEMENT_COLOURS.items()}},
                     "representations": list(scene.representations), "overrides": overrides or {}, "git_commit": source_commit(),
                     "png": checked, "validation_result": manifest.get("status", "unknown"),
                     "sha256": hashlib.sha256(png.read_bytes()).hexdigest()})
        print(f"rendered {png.name} {width}x{height} in {traced['seconds']} s fill {checked['content_fraction']} "
              f"edge {checked['touches_edge']} clipped {checked['clipped_highlight_fraction']} "
              f"crushed {checked['crushed_black_fraction']}", flush=True)
        if checked["touches_edge"] and not RENDER["crop"]:
            raise SystemExit(f"{png.name}: the drawn scene touches the image border")
        scene.clear()
    return rows


def main(argv: list | None = None) -> int:
    """Render one build directory (or one stage of it) and write a .render.json record next to each image."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="build output directory")
    parser.add_argument("--out", type=Path, required=True, help="directory for the PNG files")
    parser.add_argument("--mode", choices=sorted(MODES), default="preview")
    parser.add_argument("--all-atom", type=Path)
    parser.add_argument("--coarse-grain", type=Path)
    parser.add_argument("--tachyon")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--only", nargs="*", help="stage names to render (default: all)")
    parser.add_argument("--index", type=int, help="render only the stage with this position (0-based) in the full list; "
                                                  "for one Slurm array task per image")
    parser.add_argument("--camera-from", type=Path, help="directory of an earlier render whose cameras are reused exactly")
    parser.add_argument("--long-axis", type=int, help="override the mode's pixel size (test renders)")
    parser.add_argument("--crop", help="X,Y,SIZE: render only the SIZE-pixel square centred on pixel X,Y of the full image, "
                                       "at full pixel scale (test crops for tuning)")
    parser.add_argument("--dof", action="store_true", help="hero mode only: perspective projection with subtle depth of field")
    parser.add_argument("--set", nargs="*", metavar="KEY=VALUE", help="override central settings (tuning runs), "
                        "e.g. RENDER.ao_ambient=0.8 LIGHTS.fill=0.5 MATERIALS.MFLipid.ambient=0.3")
    args = parser.parse_args(argv)
    overrides = apply_overrides(args.set)
    mode = MODES[args.mode]
    if args.long_axis or MODE_OVERRIDES:
        mode = RenderMode(**{**asdict(mode), **MODE_OVERRIDES, **({"long_axis": args.long_axis} if args.long_axis else {})})
    if args.crop:
        if not args.camera_from:
            raise SystemExit("--crop needs --camera-from: a crop is a window of an image whose camera already exists")
        RENDER["crop"] = tuple(int(v) for v in args.crop.split(","))
    if args.dof:
        if args.mode != "hero":
            raise SystemExit("--dof is only allowed in hero mode; scientific figures stay sharp and orthographic")
        RENDER["dof"] = True
    rows = render_run(args.run.resolve(), args.out.resolve(), find_tachyon(args.tachyon), mode,
                      args.all_atom.resolve() if args.all_atom else None,
                      args.coarse_grain.resolve() if args.coarse_grain else None, args.threads, args.only, args.index,
                      args.camera_from.resolve() if args.camera_from else None, overrides)
    for row in rows:  # one record per image, so array tasks never write the same file
        (args.out / (Path(row["image"]).stem + ".render.json")).write_text(json.dumps(row, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
