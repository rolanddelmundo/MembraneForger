#!/usr/bin/env python3
"""Rewrite membrane lipid residues to match CellMembrane_Start templates exactly."""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import xml.etree.ElementTree as ET


SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent
DEFAULT_TEMPLATE_DIR = BASE_DIR.parent / "CellMembrane_Start"

# GM3 in membrane outputs uses unique names; CellMembrane_Start/GM3.pdb uses duplicated
# names grouped by template order. This list defines the unique source atom order that
# maps positionally onto the Start template order.
GM3_SOURCE_ORDER_FOR_START = [
    "N", "C", "O",
    "C1", "C7", "C14", "C1F", "C1S", "O1",
    "C2", "C8", "C15", "C2F", "C2S", "O2", "O9",
    "C3", "C10", "C16", "C3F", "C3S", "O3", "O10", "O34",
    "C4", "C11", "C17", "C4F", "C4S", "O4", "O11", "O17",
    "C5", "C12", "C20", "C5F", "C5S", "O5", "O12",
    "C6", "C13", "C18", "C6F", "C6S", "O6", "O13", "O18",
    "C19", "C7F", "C7S", "O19",
    "C21", "C8F", "C8S", "O21",
    "C22", "C9F", "C9S", "O22",
    "C10F", "C10S",
    "C11F", "C11S", "O14",
    "C12F", "C12S", "O15",
    "C13F", "C13S",
    "C14F", "C14S",
    "C15F", "C15S",
    "C16F", "C16S",
    "C17S",
    "C18S",
    "CT", "NF", "OF",
    "HN", "HNF",
    "HO2", "HO9",
    "HO3", "HO34",
    "HO11", "HO17",
    "HO6", "HO13", "HO19", "HO21", "HO22",
    "HT1", "HT2", "HT3",
    "H1", "H7",
    "H2", "H8",
    "H3", "H10",
    "H4", "H11", "H17",
    "H5", "H12", "H20",
    "H18", "H19", "H21", "H161", "H162",
    "H2F", "H2G", "H1S", "H1T",
    "H3F", "H3G", "H2S",
    "H61", "H131", "H4F", "H62", "H132", "H4G", "H3S",
    "H5F", "H5G", "H4S",
    "H6F", "H6G", "H5S",
    "H221", "H7F", "H222", "H7G", "H6S", "H6T",
    "H8F", "H8G", "H7S", "H7T",
    "H9F", "H9G", "H8S", "H8T",
    "H9S", "H9T",
    "H10F", "H10G",
    "H11F", "H11G", "H10S", "H10T",
    "H12F", "H12G", "H11S", "H11T",
    "H13F", "H13G", "H12S", "H12T",
    "H14F", "H14G", "H13S", "H13T",
    "H15F", "H15G", "H14S", "H14T",
    "H16F", "H16G", "H16H", "H15S", "H15T",
    "H16S", "H16T",
    "H17S", "H17T",
    "H18S", "H18T", "H18U",
]

TEMPLATE_RESNAME_ALIASES = {
    "CHL": "CHL1",
    "DOP": "DOPS",
    "POP": "POPS",
    "SAP": "SAPI",
    "SAP6": "SAPI",
    "SAPIS": "SAPI",
    "SAPI2": "SAPI",
    "SAPI24": "SAPI",
    "SAPI25": "SAPI",
    "TIP": "TIP3",
}

SAPI_NAME_NORMALIZER = {
    "HP42": "HP52",
}

CURRENT_NAME_NORMALIZERS = {
    "SAPI": SAPI_NAME_NORMALIZER,
    "SAPI2": SAPI_NAME_NORMALIZER,
    "SAPI24": SAPI_NAME_NORMALIZER,
    "SAPI25": SAPI_NAME_NORMALIZER,
    "SAP6": SAPI_NAME_NORMALIZER,
    "SAPIS": SAPI_NAME_NORMALIZER,
}

# Only membrane residues should be rewritten against CellMembrane_Start / lipid .itp files.
# Protein/cap residues are validated later against PRO*.itp in backmap_pipeline_2.py and
# must pass through unchanged here.
MEMBRANE_TEMPLATE_RESNAMES = {
    "CHL1",
    "DOPC",
    "DOPE",
    "DOPS",
    "GLPA",
    "GM3",
    "POPC",
    "POPE",
    "POPS",
    "PSM",
    "SAPI",
}

GLPA_MAX_BOND_LENGTH_A = 2.0
GLPA_MAX_ABS_BOND_DEV_A = 0.25

# This plan was derived once from a GM3.xml <-> GLPA.itp graph isomorphism using
# atom types, partial charges, and bond connectivity. It is fixed here so the
# output is reproducible across Python/networkx versions while validation still
# checks the final GLPA geometry against GLPA.itp and forcefield.itp.
GM3_TO_GLPA_ITP_PLAN = [
    ("C1S", "C1S"), ("H1T", "H1S"), ("H1S", "H1T"), ("NF", "NF"), ("HNF", "HNF"),
    ("C2S", "C2S"), ("H2S", "H2S"), ("C3S", "C3S"), ("H3S", "H3S"), ("O34", "O3"),
    ("HO34", "HO3"), ("C4S", "C4S"), ("H4S", "H4S"), ("C5S", "C5S"), ("H5S", "H5S"),
    ("C6S", "C6S"), ("H6S", "H6S"), ("H6T", "H6T"), ("C7S", "C7S"), ("H7T", "H7S"),
    ("H7S", "H7T"), ("C8S", "C8S"), ("H8T", "H8S"), ("H8S", "H8T"), ("C9S", "C9S"),
    ("H9T", "H9S"), ("H9S", "H9T"), ("C10S", "C10S"), ("H10T", "H10S"), ("H10S", "H10T"),
    ("C11S", "C11S"), ("H11T", "H11S"), ("H11S", "H11T"), ("C12S", "C12S"), ("H12T", "H12S"),
    ("H12S", "H12T"), ("C13S", "C13S"), ("H13S", "H13S"), ("H13T", "H13T"), ("C14S", "C14S"),
    ("H14S", "H14S"), ("H14T", "H14T"), ("C15S", "C15S"), ("H15S", "H15S"), ("H15T", "H15T"),
    ("C16S", "C16S"), ("H16T", "H16S"), ("H16S", "H16T"), ("C17S", "C17S"), ("H17T", "H17S"),
    ("H17S", "H17T"), ("C18S", "C18S"), ("H18S", "H18S"), ("H18U", "H18T"), ("H18T", "H18U"),
    ("C1F", "C1F"), ("OF", "OF"), ("C2F", "C2F"), ("H2G", "H2F"), ("H2F", "H2G"),
    ("C3F", "C3F"), ("H3F", "H3F"), ("H3G", "H3G"), ("C4F", "C4F"), ("H4G", "H4F"),
    ("H4F", "H4G"), ("C5F", "C5F"), ("H5F", "H5F"), ("H5G", "H5G"), ("C6F", "C6F"),
    ("H6G", "H6F"), ("H6F", "H6G"), ("C7F", "C7F"), ("H7F", "H7F"), ("H7G", "H7G"),
    ("C8F", "C8F"), ("H8G", "H8F"), ("H8F", "H8G"), ("C9F", "C9F"), ("H9G", "H9F"),
    ("H9F", "H9G"), ("C10F", "C10F"), ("H10G", "H10F"), ("H10F", "H10G"), ("C11F", "C11F"),
    ("H11F", "H11F"), ("H11G", "H11G"), ("C12F", "C12F"), ("H12F", "H12F"), ("H12G", "H12G"),
    ("C13F", "C13F"), ("H13F", "H13F"), ("H13G", "H13G"), ("C14F", "C14F"), ("H14F", "H14F"),
    ("H14G", "H14G"), ("C15F", "C15F"), ("H15G", "H15F"), ("H15F", "H15G"), ("C16F", "C16F"),
    ("H16G", "H16F"), ("H16H", "H16G"), ("H16F", "H16H"),
    ("C1", "C1"), ("H1", "H1"), ("O1", "O1"), ("C5", "C5"), ("H5", "H5"), ("O5", "O5"),
    ("C2", "C2"), ("H2", "H2"), ("O2", "O2"), ("HO2", "HO2"), ("C3", "C3"), ("H3", "H3"),
    ("O3", "O3"), ("HO3", "HO3"), ("C4", "C4"), ("H4", "H4"), ("O4", "O4"), ("C6", "C6"),
    ("H62", "H61"), ("H61", "H62"), ("O6", "O6"), ("HO6", "HO6"),
    ("C7", "C1"), ("H7", "H1"), ("C12", "C5"), ("H12", "H5"), ("O12", "O5"),
    ("C8", "C2"), ("H8", "H2"), ("O9", "O2"), ("HO9", "HO2"), ("C10", "C3"), ("H10", "H3"),
    ("O10", "O3"), ("C11", "C4"), ("H11", "H4"), ("O11", "O4"), ("HO11", "HO4"),
    ("C13", "C6"), ("H132", "H61"), ("H131", "H62"), ("O13", "O6"), ("HO13", "HO6"),
    ("C14", "C1"), ("O15", "O11"), ("O14", "O12"), ("C15", "C2"), ("C18", "C6"), ("H18", "H6"),
    ("O18", "O6"), ("C16", "C3"), ("H161", "H31"), ("H162", "H32"), ("C17", "C4"), ("H17", "H4"),
    ("O17", "O4"), ("HO17", "HO4"), ("C20", "C5"), ("H20", "H5"), ("N", "N"), ("HN", "HN"),
    ("C", "C"), ("O", "O"), ("CT", "CT"), ("HT1", "HT1"), ("HT2", "HT2"), ("HT3", "HT3"),
    ("C19", "C7"), ("H19", "H7"), ("O19", "O7"), ("HO19", "HO7"), ("C21", "C8"), ("H21", "H8"),
    ("O21", "O8"), ("HO21", "HO8"), ("C22", "C9"), ("H222", "H91"), ("H221", "H92"),
    ("O22", "O9"), ("HO22", "HO9"),
]


@dataclass(frozen=True)
class ItpAtom:
    index: int
    atom_type: str
    residu: str
    name: str
    charge: float


@dataclass(frozen=True)
class GlpaValidationSummary:
    topology_path: str
    residue_count: int
    bond_count: int
    mean_abs_bond_dev_A: float
    max_abs_bond_dev_A: float
    max_bond_length_A: float


def parse_itp_atom_order(path: Path) -> list[str]:
    return [atom.name for atom in parse_itp_atoms(path)]


def parse_itp_atoms(path: Path) -> list[ItpAtom]:
    atoms: list[ItpAtom] = []
    in_atoms = False
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith(";"):
                continue
            if stripped.startswith("["):
                in_atoms = stripped.lower() == "[ atoms ]"
                continue
            if in_atoms:
                parts = stripped.split()
                if len(parts) >= 7:
                    atoms.append(
                        ItpAtom(
                            index=int(parts[0]),
                            atom_type=parts[1],
                            residu=parts[3],
                            name=parts[4],
                            charge=float(parts[6]),
                        )
                    )
    if not atoms:
        raise ValueError(f"No atoms were found in {path}")
    return atoms


def parse_itp_bonds(path: Path) -> list[tuple[int, int]]:
    bonds: list[tuple[int, int]] = []
    in_bonds = False
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith(";"):
                continue
            if stripped.startswith("["):
                in_bonds = stripped.lower() == "[ bonds ]"
                continue
            if in_bonds:
                parts = stripped.split()
                if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                    bonds.append((int(parts[0]), int(parts[1])))
    if not bonds:
        raise ValueError(f"No bonds were found in {path}")
    return bonds


def parse_forcefield_bond_lengths(path: Path) -> dict[tuple[str, str], float]:
    bond_lengths: dict[tuple[str, str], float] = {}
    seen: set[Path] = set()

    def parse_one(current: Path) -> None:
        resolved = current.resolve()
        if resolved in seen:
            return
        seen.add(resolved)
        in_bondtypes = False
        with current.open() as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith(";"):
                    continue
                if stripped.startswith("#include"):
                    parts = stripped.split('"')
                    if len(parts) >= 2:
                        included = current.parent / parts[1]
                        if included.exists():
                            parse_one(included)
                    continue
                if stripped.startswith("["):
                    in_bondtypes = stripped.lower() == "[ bondtypes ]"
                    continue
                if not in_bondtypes:
                    continue
                parts = stripped.split()
                if len(parts) < 5 or parts[2] != "1":
                    continue
                # GROMACS stores b0 in nm; convert to angstroms for validation.
                bond_lengths[(parts[0], parts[1])] = float(parts[3]) * 10.0

    parse_one(path)
    if not bond_lengths:
        raise ValueError(f"No [ bondtypes ] were found in {path}")
    return bond_lengths


def parse_gm3_xml_atoms_and_bonds(path: Path) -> tuple[list[tuple[str, str, float]], list[tuple[str, str]]]:
    root = ET.parse(path).getroot()
    for residue in root.iter():
        if residue.tag.lower().endswith("residue") and residue.attrib.get("name") == "GM3":
            atoms: list[tuple[str, str, float]] = []
            bonds: list[tuple[str, str]] = []
            for child in residue:
                tag = child.tag.lower()
                if tag.endswith("atom"):
                    atoms.append(
                        (
                            child.attrib["name"],
                            child.attrib.get("type", ""),
                            float(child.attrib.get("charge", "0.0")),
                        )
                    )
                elif tag.endswith("bond"):
                    atom_1 = child.attrib.get("atomName1") or child.attrib.get("from")
                    atom_2 = child.attrib.get("atomName2") or child.attrib.get("to")
                    if atom_1 and atom_2:
                        bonds.append((atom_1, atom_2))
            if not atoms or not bonds:
                raise ValueError(f"GM3 residue in {path} is missing atoms or bonds")
            return atoms, bonds
    raise ValueError(f"GM3 residue was not found in {path}")


def _toppar_dir_key(toppar_dir: Path | None) -> str | None:
    return str(toppar_dir.resolve()) if toppar_dir else None


def _toppar_file_candidates(filename: str, toppar_dir: Path | None = None) -> list[Path]:
    candidates = []
    if toppar_dir is not None:
        candidates.append(toppar_dir / filename)
    candidates.extend([
        BASE_DIR / "charmmguitops" / "GPR_charmm_0207" / "gromacs" / "toppar" / filename,
        BASE_DIR / "charmmguitops" / "KOR_charmm_0207" / "gromacs" / "toppar" / filename,
        BASE_DIR / "toppar" / filename,
        BASE_DIR / "KOR" / "toppar" / filename,
        BASE_DIR / "GPR" / "toppar" / filename,
    ])
    return candidates


def find_residue_itp(resname: str, toppar_dir: Path | None = None) -> Path:
    filename = f"{resname}.itp"
    for candidate in _toppar_file_candidates(filename, toppar_dir):
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not locate {filename} for topology-order normalization")


def find_glpa_itp(toppar_dir: Path | None = None) -> Path:
    return find_residue_itp("GLPA", toppar_dir)


def find_forcefield_itp(toppar_dir: Path | None = None) -> Path:
    candidates = _toppar_file_candidates("forcefield.itp", toppar_dir) + [
        BASE_DIR / "resources" / "external" / "charmm36.ff" / "forcefield.itp",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Could not locate forcefield.itp for GLPA bond validation")


def find_gm3_xml(gm3_xml: Path | None = None) -> Path:
    candidates = ([gm3_xml] if gm3_xml is not None else []) + [
        BASE_DIR / "resources" / "templates" / "GM3.xml",
        BASE_DIR / "GM3.xml",
        BASE_DIR.parent / "GM3.xml",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Could not locate GM3.xml for GM3->GLPA topology mapping")


@lru_cache(maxsize=None)
def get_glpa_itp_atoms_and_bonds(toppar_dir_key: str | None = None) -> tuple[list[ItpAtom], list[tuple[int, int]], Path]:
    path = find_glpa_itp(Path(toppar_dir_key) if toppar_dir_key else None)
    return parse_itp_atoms(path), parse_itp_bonds(path), path


@lru_cache(maxsize=None)
def get_itp_atom_order(resname: str, toppar_dir_key: str | None = None) -> tuple[str, ...]:
    path = find_residue_itp(resname, Path(toppar_dir_key) if toppar_dir_key else None)
    return tuple(parse_itp_atom_order(path))


@lru_cache(maxsize=None)
def get_glpa_bond_lengths_A(toppar_dir_key: str | None = None) -> tuple[dict[tuple[str, str], float], Path]:
    path = find_forcefield_itp(Path(toppar_dir_key) if toppar_dir_key else None)
    return parse_forcefield_bond_lengths(path), path


@lru_cache(maxsize=None)
def get_gm3_xml_atoms_and_bonds(gm3_xml_key: str | None = None) -> tuple[list[tuple[str, str, float]], list[tuple[str, str]], Path]:
    path = find_gm3_xml(Path(gm3_xml_key) if gm3_xml_key else None)
    atoms, bonds = parse_gm3_xml_atoms_and_bonds(path)
    return atoms, bonds, path


@lru_cache(maxsize=None)
def get_gm3_to_glpa_itp_plan(
    toppar_dir_key: str | None = None,
    gm3_xml_key: str | None = None,
) -> tuple[tuple[str, str], ...]:
    glpa_atoms, _, _ = get_glpa_itp_atoms_and_bonds(toppar_dir_key)
    gm3_xml_atoms, _, _ = get_gm3_xml_atoms_and_bonds(gm3_xml_key)
    target_names = [atom.name for atom in glpa_atoms]
    source_names = [name for name, _, _ in gm3_xml_atoms]
    if len(GM3_TO_GLPA_ITP_PLAN) != len(target_names):
        raise ValueError("GM3_TO_GLPA_ITP_PLAN length does not match GLPA.itp atom count")
    if [target_name for _, target_name in GM3_TO_GLPA_ITP_PLAN] != target_names:
        raise ValueError("GM3_TO_GLPA_ITP_PLAN target order does not match GLPA.itp")
    if Counter(source_name for source_name, _ in GM3_TO_GLPA_ITP_PLAN) != Counter(source_names):
        raise ValueError("GM3_TO_GLPA_ITP_PLAN source inventory does not match GM3.xml")
    return tuple(GM3_TO_GLPA_ITP_PLAN)


@lru_cache(maxsize=None)
def get_glpa_itp_source_order(toppar_dir_key: str | None = None) -> tuple[str, ...]:
    return tuple(atom.name for atom in get_glpa_itp_atoms_and_bonds(toppar_dir_key)[0])


@dataclass
class GroAtomLine:
    resid_field: str
    resname_field: str
    atomname: str
    tail: str


def canonical_resname(resname: str) -> str:
    return TEMPLATE_RESNAME_ALIASES.get(resname, resname)


def topology_target_resname(raw_resname: str, canonical_name: str, gm3_target: str) -> str:
    if gm3_target != "itp":
        return canonical_name
    if raw_resname in {"GM3", "GLPA"}:
        return "GLPA"
    if canonical_name in {"SAPI", "SAPI25"}:
        return "SAPI"
    return canonical_name


def topology_order_resname(final_resname: str) -> str:
    if final_resname == "SAPI":
        return "SAPI25"
    return final_resname


def supported_template_names(
    raw_resname: str,
    canonical_name: str,
    template_orders: dict[str, list[str]],
) -> list[str] | None:
    if raw_resname not in MEMBRANE_TEMPLATE_RESNAMES and canonical_name not in MEMBRANE_TEMPLATE_RESNAMES:
        return None
    if raw_resname in {"GM3", "GLPA"}:
        return template_orders.get("GM3")
    return template_orders.get(canonical_name)


def parse_template_atom_orders(template_dir: Path) -> dict[str, list[str]]:
    orders: dict[str, list[str]] = {}
    for path in sorted(template_dir.glob("*.pdb")):
        names = [
            line[12:16].strip()
            for line in path.read_text().splitlines()
            if line.startswith(("ATOM", "HETATM"))
        ]
        if names:
            orders[path.stem] = names
    if not orders:
        raise ValueError(f"No template PDB files were found in {template_dir}")
    return orders


def format_pdb_atom_name(atom_name: str) -> str:
    if len(atom_name) >= 4:
        return atom_name[:4]
    if atom_name and atom_name[0].isdigit():
        return f"{atom_name:<4}"
    return f"{atom_name:>4}"


def format_pdb_resname(resname: str) -> str:
    return f"{resname:>4}"[:4]


def rename_pdb_line(line: str, atom_name: str, resname: str | None = None) -> str:
    renamed = f"{line[:12]}{format_pdb_atom_name(atom_name)}{line[16:]}"
    if resname is None:
        return renamed
    return f"{renamed[:17]}{format_pdb_resname(resname)}{renamed[21:]}"


def normalize_current_names(resname: str, names: list[str]) -> list[str]:
    normalizer = CURRENT_NAME_NORMALIZERS.get(resname, {})
    return [normalizer.get(name, name) for name in names]


def build_reorder_lookup_pdb(current_block: list[str], current_names: list[str]) -> dict[str, deque[str]]:
    lookup: dict[str, deque[str]] = defaultdict(deque)
    for line, name in zip(current_block, current_names):
        lookup[name].append(line)
    return lookup


def build_reorder_lookup_gro(atoms: list[GroAtomLine], current_names: list[str]) -> dict[str, deque[GroAtomLine]]:
    lookup: dict[str, deque[GroAtomLine]] = defaultdict(deque)
    for atom, name in zip(atoms, current_names):
        lookup[name].append(atom)
    return lookup


def reorder_names_for_template(
    resname: str,
    current_names: list[str],
    target_names: list[str],
    gm3_target: str = "start",
    glpa_itp_source_order: tuple[str, ...] | None = None,
    gm3_to_glpa_itp_plan: tuple[tuple[str, str], ...] | None = None,
    gm3_xml_key: str | None = None,
) -> list[tuple[str, str]]:
    if resname in {"GM3", "GLPA"}:
        glpa_order = glpa_itp_source_order or get_glpa_itp_source_order()
        gm3_plan = gm3_to_glpa_itp_plan or get_gm3_to_glpa_itp_plan()
        if current_names == target_names:
            return list(zip(current_names, target_names))
        if gm3_target == "itp":
            if current_names == list(glpa_order):
                return list(zip(current_names, target_names))
            gm3_xml_atoms, _, _ = get_gm3_xml_atoms_and_bonds(gm3_xml_key)
            gm3_xml_names = [name for name, _, _ in gm3_xml_atoms]
            if Counter(current_names) == Counter(gm3_xml_names):
                return list(gm3_plan)
            raise ValueError(
                "GM3 residue does not match GLPA.itp order or the GM3.xml atom inventory; "
                "refusing to remap coordinates onto GLPA atom names."
            )
        if Counter(current_names) != Counter(GM3_SOURCE_ORDER_FOR_START):
            raise ValueError("GM3 residue does not match a supported source atom naming/order")
        if len(target_names) != len(GM3_SOURCE_ORDER_FOR_START):
            raise ValueError("GM3 template atom count does not match the expected source order")
        return list(zip(GM3_SOURCE_ORDER_FOR_START, target_names))

    normalized = normalize_current_names(resname, current_names)
    if Counter(normalized) != Counter(target_names):
        raise ValueError(
            f"{resname} residue does not match the target atom inventory"
        )
    by_name = {name: original for name, original in zip(normalized, current_names)}
    return [(by_name[name], name) for name in target_names]


def process_pdb(
    input_path: Path,
    output_path: Path,
    template_orders: dict[str, list[str]],
    gm3_target: str = "start",
    toppar_dir_key: str | None = None,
    glpa_itp_source_order: tuple[str, ...] | None = None,
    gm3_to_glpa_itp_plan: tuple[tuple[str, str], ...] | None = None,
    gm3_xml_key: str | None = None,
) -> int:
    lines = input_path.read_text().splitlines(keepends=True)
    output_lines: list[str] = []
    current_block: list[str] = []
    current_key: tuple[str, str, str, str] | None = None
    converted = 0

    def flush_block() -> None:
        nonlocal converted
        if not current_block:
            return

        raw_resname = current_key[3].strip() if current_key else ""
        resname = canonical_resname(raw_resname)
        template_names = supported_template_names(raw_resname, resname, template_orders)
        if template_names is None:
            output_lines.extend(current_block)
            return

        final_resname = topology_target_resname(raw_resname, resname, gm3_target)
        if raw_resname in {"GM3", "GLPA"} and gm3_target == "itp":
            target_names = list(glpa_itp_source_order or get_glpa_itp_source_order(toppar_dir_key))
        elif gm3_target == "itp":
            target_names = list(get_itp_atom_order(topology_order_resname(final_resname), toppar_dir_key))
        else:
            target_names = template_names

        current_names = [line[12:16].strip() for line in current_block]
        if (
            raw_resname in {"GM3", "GLPA"}
            and gm3_target == "itp"
            and current_names == template_names
        ):
            renaming_plan = list(zip(current_names, target_names))
        else:
            renaming_plan = reorder_names_for_template(
                raw_resname,
                current_names,
                target_names,
                gm3_target=gm3_target,
                glpa_itp_source_order=glpa_itp_source_order,
                gm3_to_glpa_itp_plan=gm3_to_glpa_itp_plan,
                gm3_xml_key=gm3_xml_key,
            )
        by_name = build_reorder_lookup_pdb(current_block, current_names)
        reordered_lines = [
            rename_pdb_line(by_name[source_name].popleft(), target_name, final_resname)
            for source_name, target_name in renaming_plan
        ]
        output_lines.extend(reordered_lines)
        converted += 1

    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            key = (line[21], line[22:26], line[26], line[17:21])
            if current_key is None or key == current_key:
                current_block.append(line)
                current_key = key
            else:
                flush_block()
                current_block = [line]
                current_key = key
        else:
            flush_block()
            current_block = []
            current_key = None
            output_lines.append(line)
    flush_block()

    serial = 1
    renumbered_lines: list[str] = []
    for line in output_lines:
        if line.startswith(("ATOM", "HETATM")):
            renumbered_lines.append(f"{line[:6]}{serial % 100000:5d}{line[11:]}")
            serial += 1
        elif line.startswith("TER"):
            renumbered_lines.append(f"{line[:6]}{serial % 100000:5d}{line[11:]}")
            serial += 1
        else:
            renumbered_lines.append(line)

    output_path.write_text("".join(renumbered_lines))
    return converted


def parse_gro_blocks(atom_lines: list[str]) -> list[tuple[str, list[GroAtomLine]]]:
    blocks: list[tuple[str, list[GroAtomLine]]] = []
    current_key: tuple[str, str] | None = None
    current_atoms: list[GroAtomLine] = []
    current_resname = ""

    for raw in atom_lines:
        if len(raw) < 20:
            continue
        resid_field = raw[0:5]
        resname_field = raw[5:10]
        atomname = raw[10:15].strip()
        tail = raw[20:]
        key = (resid_field, resname_field)
        if current_key is None or key == current_key:
            current_key = key
        else:
            blocks.append((current_resname, current_atoms))
            current_atoms = []
            current_key = key
        current_resname = resname_field.strip()
        current_atoms.append(GroAtomLine(resid_field, resname_field, atomname, tail))

    if current_atoms:
        blocks.append((current_resname, current_atoms))
    return blocks


def rebuild_gro(title: str, blocks: list[tuple[str, list[GroAtomLine]]], box_line: str) -> str:
    total_atoms = sum(len(atoms) for _, atoms in blocks)
    lines = [title, f"{total_atoms:5d}\n"]
    atom_number = 1
    for _, atoms in blocks:
        for atom in atoms:
            lines.append(
                f"{atom.resid_field}{atom.resname_field}{atom.atomname:>5}{atom_number % 100000:5d}{atom.tail}"
            )
            atom_number += 1
    lines.append(box_line)
    return "".join(lines)


def process_gro(
    input_path: Path,
    output_path: Path,
    template_orders: dict[str, list[str]],
    gm3_target: str = "start",
    toppar_dir_key: str | None = None,
    glpa_itp_source_order: tuple[str, ...] | None = None,
    gm3_to_glpa_itp_plan: tuple[tuple[str, str], ...] | None = None,
    gm3_xml_key: str | None = None,
) -> int:
    with input_path.open() as handle:
        title = handle.readline()
        atom_count = int(handle.readline().strip())
        atom_lines = [handle.readline() for _ in range(atom_count)]
        box_line = handle.readline()

    blocks = parse_gro_blocks(atom_lines)
    converted = 0
    new_blocks: list[tuple[str, list[GroAtomLine]]] = []

    for raw_resname, atoms in blocks:
        resname = canonical_resname(raw_resname)
        template_names = supported_template_names(raw_resname, resname, template_orders)
        if template_names is None:
            new_blocks.append((raw_resname, atoms))
            continue

        final_resname = topology_target_resname(raw_resname, resname, gm3_target)
        if raw_resname in {"GM3", "GLPA"} and gm3_target == "itp":
            target_names = list(glpa_itp_source_order or get_glpa_itp_source_order(toppar_dir_key))
        elif gm3_target == "itp":
            target_names = list(get_itp_atom_order(topology_order_resname(final_resname), toppar_dir_key))
        else:
            target_names = template_names

        current_names = [atom.atomname for atom in atoms]
        if (
            raw_resname in {"GM3", "GLPA"}
            and gm3_target == "itp"
            and current_names == template_names
        ):
            renaming_plan = list(zip(current_names, target_names))
        else:
            renaming_plan = reorder_names_for_template(
                raw_resname,
                current_names,
                target_names,
                gm3_target=gm3_target,
                glpa_itp_source_order=glpa_itp_source_order,
                gm3_to_glpa_itp_plan=gm3_to_glpa_itp_plan,
                gm3_xml_key=gm3_xml_key,
            )
        by_name = build_reorder_lookup_gro(atoms, current_names)
        reordered_atoms: list[GroAtomLine] = []
        target_resname_field = f"{final_resname:>5}"[:5]
        for source_name, target_name in renaming_plan:
            atom = by_name[source_name].popleft()
            reordered_atoms.append(
                GroAtomLine(
                    resid_field=atom.resid_field,
                    resname_field=target_resname_field,
                    atomname=target_name,
                    tail=atom.tail,
                )
            )
        new_blocks.append((final_resname, reordered_atoms))
        converted += 1

    output_path.write_text(rebuild_gro(title, new_blocks, box_line))
    return converted


def apply_start_templates(
    pdb_path: Path | None,
    gro_path: Path | None = None,
    template_dir: Path = DEFAULT_TEMPLATE_DIR,
    gm3_target: str = "start",
    toppar_dir: Path | None = None,
    gm3_xml: Path | None = None,
) -> dict[str, int]:
    template_orders = parse_template_atom_orders(template_dir)
    toppar_dir_key = _toppar_dir_key(toppar_dir)
    glpa_itp_source_order = get_glpa_itp_source_order(toppar_dir_key) if gm3_target == "itp" else None
    gm3_xml_key = str(gm3_xml.resolve()) if gm3_xml is not None else None
    gm3_to_glpa_itp_plan = (
        get_gm3_to_glpa_itp_plan(toppar_dir_key, gm3_xml_key)
        if gm3_target == "itp"
        else None
    )
    converted = {"pdb": 0, "gro": 0}
    if pdb_path and pdb_path.exists():
        converted["pdb"] = process_pdb(
            pdb_path,
            pdb_path,
            template_orders,
            gm3_target=gm3_target,
            toppar_dir_key=toppar_dir_key,
            glpa_itp_source_order=glpa_itp_source_order,
            gm3_to_glpa_itp_plan=gm3_to_glpa_itp_plan,
            gm3_xml_key=gm3_xml_key,
        )
    if gro_path and gro_path.exists():
        converted["gro"] = process_gro(
            gro_path,
            gro_path,
            template_orders,
            gm3_target=gm3_target,
            toppar_dir_key=toppar_dir_key,
            glpa_itp_source_order=glpa_itp_source_order,
            gm3_to_glpa_itp_plan=gm3_to_glpa_itp_plan,
            gm3_xml_key=gm3_xml_key,
        )
    return converted


def iter_pdb_residue_blocks(path: Path, resname: str) -> list[list[str]]:
    blocks: list[list[str]] = []
    current_block: list[str] = []
    current_key: tuple[str, str, str, str] | None = None
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                if current_block:
                    blocks.append(current_block)
                    current_block = []
                    current_key = None
                continue
            key = (line[21], line[22:26], line[26], line[17:21])
            if line[17:21].strip() != resname:
                if current_block:
                    blocks.append(current_block)
                    current_block = []
                    current_key = None
                continue
            if current_key is None or key == current_key:
                current_block.append(line.rstrip("\n"))
                current_key = key
            else:
                blocks.append(current_block)
                current_block = [line.rstrip("\n")]
                current_key = key
    if current_block:
        blocks.append(current_block)
    return blocks


def iter_gro_residue_blocks(path: Path, resname: str) -> list[list[GroAtomLine]]:
    with path.open() as handle:
        _ = handle.readline()
        atom_count = int(handle.readline().strip())
        atom_lines = [handle.readline() for _ in range(atom_count)]
    return [atoms for block_resname, atoms in parse_gro_blocks(atom_lines) if block_resname.strip() == resname]


def validate_glpa_blocks(
    blocks: list[tuple[list[str], list[tuple[float, float, float]]]],
    toppar_dir: Path | None = None,
) -> GlpaValidationSummary:
    toppar_dir_key = _toppar_dir_key(toppar_dir)
    glpa_atoms, glpa_bonds, topology_path = get_glpa_itp_atoms_and_bonds(toppar_dir_key)
    target_names = [atom.name for atom in glpa_atoms]
    bond_lengths_A, _ = get_glpa_bond_lengths_A(toppar_dir_key)

    residue_count = 0
    bond_count = 0
    abs_deviations: list[float] = []
    bond_lengths: list[float] = []

    for names, positions in blocks:
        residue_count += 1
        if names != target_names:
            raise ValueError("GLPA residue atom order does not match GLPA.itp exactly")
        if len(positions) != len(glpa_atoms):
            raise ValueError("GLPA residue atom count does not match GLPA.itp")
        for atom_i, atom_j in glpa_bonds:
            left = glpa_atoms[atom_i - 1]
            right = glpa_atoms[atom_j - 1]
            ref_b0_A = bond_lengths_A.get((left.atom_type, right.atom_type), bond_lengths_A.get((right.atom_type, left.atom_type)))
            if ref_b0_A is None:
                raise ValueError(
                    f"Bond type for GLPA atoms {left.name}-{right.name} "
                    f"({left.atom_type}-{right.atom_type}) was not found in forcefield.itp"
                )
            x1, y1, z1 = positions[atom_i - 1]
            x2, y2, z2 = positions[atom_j - 1]
            bond_length_A = ((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) ** 0.5
            abs_deviation = abs(bond_length_A - ref_b0_A)
            abs_deviations.append(abs_deviation)
            bond_lengths.append(bond_length_A)
            bond_count += 1

    if residue_count and (
        max(bond_lengths) > GLPA_MAX_BOND_LENGTH_A
        or max(abs_deviations) > GLPA_MAX_ABS_BOND_DEV_A
    ):
        raise ValueError(
            "GLPA bond geometry validation failed: "
            f"max bond length = {max(bond_lengths):.3f} Å, "
            f"max |bond-b0| = {max(abs_deviations):.3f} Å"
        )

    return GlpaValidationSummary(
        topology_path=str(topology_path),
        residue_count=residue_count,
        bond_count=bond_count,
        mean_abs_bond_dev_A=(sum(abs_deviations) / len(abs_deviations)) if abs_deviations else 0.0,
        max_abs_bond_dev_A=max(abs_deviations) if abs_deviations else 0.0,
        max_bond_length_A=max(bond_lengths) if bond_lengths else 0.0,
    )


def validate_glpa_outputs(
    pdb_path: Path | None,
    gro_path: Path | None = None,
    toppar_dir: Path | None = None,
) -> dict[str, GlpaValidationSummary]:
    results: dict[str, GlpaValidationSummary] = {}
    if pdb_path and pdb_path.exists():
        pdb_blocks = [
            (
                [line[12:16].strip() for line in block],
                [
                    (
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    )
                    for line in block
                ],
            )
            for block in iter_pdb_residue_blocks(pdb_path, "GLPA")
        ]
        results["pdb"] = validate_glpa_blocks(pdb_blocks, toppar_dir=toppar_dir)
    if gro_path and gro_path.exists():
        gro_blocks = [
            (
                [atom.atomname for atom in block],
                [
                    (
                        float(atom.tail[0:8]) * 10.0,
                        float(atom.tail[8:16]) * 10.0,
                        float(atom.tail[16:24]) * 10.0,
                    )
                    for atom in block
                ],
            )
            for block in iter_gro_residue_blocks(gro_path, "GLPA")
        ]
        results["gro"] = validate_glpa_blocks(gro_blocks, toppar_dir=toppar_dir)
    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, help="Path to membrane.pdb")
    parser.add_argument("--gro", type=Path, help="Path to membrane.gro")
    parser.add_argument(
        "--templates",
        type=Path,
        default=DEFAULT_TEMPLATE_DIR,
        help=f"Directory of CellMembrane_Start template PDBs (default: {DEFAULT_TEMPLATE_DIR})",
    )
    parser.add_argument(
        "--gm3-target",
        choices=["start", "itp"],
        default="start",
        help="Use CellMembrane_Start naming/order or topology (.itp) atom order for supported membrane residues.",
    )
    parser.add_argument("--toppar-dir", type=Path, help="Optional toppar directory containing GLPA.itp and forcefield.itp.")
    args = parser.parse_args()

    result = apply_start_templates(
        args.pdb,
        args.gro,
        args.templates.resolve(),
        gm3_target=args.gm3_target,
        toppar_dir=args.toppar_dir.resolve() if args.toppar_dir else None,
    )
    print(f"Converted lipid residues -> pdb: {result['pdb']}, gro: {result['gro']}")
