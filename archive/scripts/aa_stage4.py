"""Repository-local all-atom topology, slicing, solvation, and EM helpers."""
from __future__ import annotations

from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Any, Callable
import math
import shutil
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from membraneforger.chemistry import component_for


class Stage4Error(RuntimeError):
    pass


PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HSD", "HSE", "HSP", "ACE", "NME", "NMA", "CT3", "NTER", "CTER",
}
CHARMM_EXPLICIT_TERMINAL_CAP_RESNAMES = {"ACE", "NME", "NMA", "CT3"}
WATER_IONS = {"W", "WF", "SOL", "TIP3", "NA", "CL", "SOD", "CLA", "ION"}
GLPA_COMPONENT_RESNAMES = {"CER1", "CER160", "BGLC", "BGAL", "ANE5", "ANE5AC"}
GLPA_COMPONENT_SEQUENCE = ("CER1", "BGLC", "BGAL", "ANE5")
MOLECULE_ALIASES = {
    "CHOL": "CHL1", "CHL": "CHL1", "CLR": "CHL1", "CHL1": "CHL1",
    "GM3": "GLPA", "GLPA": "GLPA",
    "SAPI": "SAPI25", "SAPI24": "SAPI25", "SAPI25": "SAPI25", "SAP6": "SAPI25",
    "TIP3": "TIP3", "SOL": "TIP3", "SOD": "SOD", "NA": "SOD", "CLA": "CLA", "CL": "CLA",
}
MEMBRANE_MOLECULE_ORDER = ("CHL1", "DOPC", "DOPE", "DOPS", "GLPA", "POPC", "POPE", "POPS", "PSM", "SAPI25")
SOLVENT_ION_MOLECULE_ORDER = ("TIP3", "SOD", "CLA")


def _pdb_atoms(path: Path) -> tuple[list[str], list[dict[str, Any]], list[str]]:
    prefix: list[str] = []
    suffix: list[str] = []
    atoms: list[dict[str, Any]] = []
    seen_atoms = False
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.startswith(("ATOM  ", "HETATM")):
            seen_atoms = True
            atoms.append({
                "record": raw[:6], "name": raw[12:16].strip(), "resname": raw[17:21].strip(),
                "chain": raw[21:22].strip() or " ", "resid": int(raw[22:26]), "icode": raw[26:27].strip(),
                "x": float(raw[30:38]), "y": float(raw[38:46]), "z": float(raw[46:54]),
                "element": raw[76:78].strip() if len(raw) >= 78 else "",
            })
        elif not seen_atoms:
            prefix.append(raw)
        elif raw not in {"END", "TER"}:
            suffix.append(raw)
    if not atoms:
        raise Stage4Error(f"no ATOM/HETATM records found in {path}")
    return prefix, atoms, suffix


def _box_from_prefix(prefix: list[str]) -> np.ndarray:
    for line in prefix:
        if line.startswith("CRYST1"):
            return np.array([float(line[6:15]), float(line[15:24]), float(line[24:33])], dtype=float)
    raise Stage4Error("Stage 4 input PDB has no CRYST1 box record")


def _write_pdb(path: Path, prefix: list[str], atoms: list[dict[str, Any]], box_a: np.ndarray) -> None:
    out = [line for line in prefix if not line.startswith("CRYST1")]
    out.append(f"CRYST1{box_a[0]:9.3f}{box_a[1]:9.3f}{box_a[2]:9.3f}  90.00  90.00  90.00 P 1           1")
    previous_chain = None
    for serial, atom in enumerate(atoms, start=1):
        if previous_chain is not None and atom["chain"] != previous_chain:
            out.append("TER")
        previous_chain = atom["chain"]
        element = atom["element"] or atom["name"][:1]
        serial_field = serial % 100000
        out.append(
            f"{atom['record']:<6}{serial_field:5d} {atom['name']:>4} {atom['resname']:>4}{atom['chain'][:1]:1}"
            f"{int(atom['resid']):4d}{atom['icode'][:1]:1}   {atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}          {element:>2}"
        )
    out.extend(["TER", "END"])
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def strip_inherited_solvent(input_pdb: Path, output_pdb: Path) -> dict[str, int]:
    prefix, atoms, _ = _pdb_atoms(input_pdb)
    box_a = _box_from_prefix(prefix)
    kept = [atom for atom in atoms if atom["resname"] not in WATER_IONS]
    if not kept:
        raise Stage4Error("all Stage 3 atoms were classified as solvent or ions")
    _write_pdb(output_pdb, prefix, kept, box_a)
    return {"removed_atoms": len(atoms) - len(kept), "kept_atoms": len(kept)}


def explicit_charmm_terminal_cap_residues(path: Path) -> list[str]:
    _, atoms, _ = _pdb_atoms(path)
    seen: set[tuple[str, int, str, str]] = set()
    labels: list[str] = []
    for atom in atoms:
        resname = atom["resname"]
        if resname not in CHARMM_EXPLICIT_TERMINAL_CAP_RESNAMES:
            continue
        key = (atom["chain"], int(atom["resid"]), atom["icode"], resname)
        if key in seen:
            continue
        seen.add(key)
        icode = atom["icode"] or ""
        labels.append(f"{atom['chain']}:{resname}{int(atom['resid'])}{icode}")
    return labels


def reject_explicit_charmm_terminal_cap_residues(path: Path, context: str) -> None:
    labels = explicit_charmm_terminal_cap_residues(path)
    if labels:
        raise Stage4Error(
            f"{context} contains explicit terminal cap residues {labels[:20]}. "
            "Stage 4 uses CHARMM terminal patches instead: N-terminal ACE and C-terminal CT3. "
            "Do not add standalone ACE/NME residues for the CHARMM36 workflow; preserve or resolve this chemistry explicitly before Stage 4."
        )


def align_bilayer_to_z(input_pdb: Path, output_pdb: Path) -> dict[str, float]:
    prefix, atoms, _ = _pdb_atoms(input_pdb)
    box_a = _box_from_prefix(prefix)
    phosphates = np.array([[atom["x"], atom["y"], atom["z"]] for atom in atoms if atom["name"] == "P"])
    if len(phosphates) < 3:
        raise Stage4Error("cannot align bilayer: fewer than three phosphate atoms")
    _, _, vt = np.linalg.svd(phosphates - phosphates.mean(axis=0), full_matrices=False)
    normal = vt[-1]
    if normal[2] < 0.0:
        normal = -normal
    normal /= np.linalg.norm(normal)
    target = np.array([0.0, 0.0, 1.0])
    cosine = float(np.clip(np.dot(normal, target), -1.0, 1.0))
    angle = math.degrees(math.acos(cosine))
    axis = np.cross(normal, target)
    axis_norm = np.linalg.norm(axis)
    if axis_norm > 1e-10:
        axis /= axis_norm
        skew = np.array([[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]])
        rotation = np.eye(3) + skew * axis_norm + skew @ skew * (1.0 - cosine)
        center = np.mean([[atom["x"], atom["y"], atom["z"]] for atom in atoms], axis=0)
        for atom in atoms:
            xyz = rotation @ (np.array([atom["x"], atom["y"], atom["z"]]) - center) + center
            atom["x"], atom["y"], atom["z"] = (float(xyz[0]), float(xyz[1]), float(xyz[2]))
    _write_pdb(output_pdb, prefix, atoms, box_a)
    return {"initial_normal_z_cosine": cosine, "initial_angle_deg": angle}


def _residue_blocks(atoms: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    blocks: list[list[dict[str, Any]]] = []
    current: list[dict[str, Any]] = []
    key: tuple[str, int, str, str] | None = None
    for atom in atoms:
        atom_key = (atom["chain"], atom["resid"], atom["icode"], atom["resname"])
        if key is None or key == atom_key:
            current.append(atom)
        else:
            blocks.append(current)
            current = [atom]
        key = atom_key
    if current:
        blocks.append(current)
    return blocks


def _component_key(resname: str) -> str:
    if resname.startswith("CER1"):
        return "CER1"
    if resname.startswith("ANE5"):
        return "ANE5"
    return resname


def _molecule_blocks(atoms: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    residue_blocks = _residue_blocks(atoms)
    molecule_blocks: list[list[dict[str, Any]]] = []
    index = 0
    while index < len(residue_blocks):
        block = residue_blocks[index]
        if block[0]["resname"] in GLPA_COMPONENT_RESNAMES:
            component_blocks = residue_blocks[index:index + len(GLPA_COMPONENT_SEQUENCE)]
            sequence = tuple(_component_key(item[0]["resname"]) for item in component_blocks)
            if sequence != GLPA_COMPONENT_SEQUENCE:
                raise Stage4Error(f"incomplete or out-of-order GLPA component sequence near residue {block[0]['resid']}")
            molecule_blocks.append([atom for item in component_blocks for atom in item])
            index += len(GLPA_COMPONENT_SEQUENCE)
            continue
        molecule_blocks.append(block)
        index += 1
    return molecule_blocks


def _gro_atoms(path: Path, chain: str = "X") -> tuple[list[dict[str, Any]], np.ndarray]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        raise Stage4Error(f"invalid GRO file: {path}")
    try:
        natoms = int(lines[1].strip())
    except ValueError as exc:
        raise Stage4Error(f"invalid GRO atom count in {path}") from exc
    if len(lines) != natoms + 3:
        raise Stage4Error(f"GRO atom count does not match file length: {path}")
    atoms: list[dict[str, Any]] = []
    for raw in lines[2:2 + natoms]:
        atoms.append({
            "record": "ATOM",
            "name": raw[10:15].strip(),
            "resname": raw[5:10].strip(),
            "chain": chain[:1] or "X",
            "resid": int(raw[0:5]),
            "icode": "",
            "x": float(raw[20:28]) * 10.0,
            "y": float(raw[28:36]) * 10.0,
            "z": float(raw[36:44]) * 10.0,
            "element": raw[10:15].strip()[:1],
        })
    box = np.array([float(value) * 10.0 for value in lines[-1].split()[:3]], dtype=float)
    return atoms, box


def _molecule_for_block(block: list[dict[str, Any]]) -> str | None:
    resname = block[0]["resname"]
    if resname in WATER_IONS:
        return None
    if resname in GLPA_COMPONENT_RESNAMES:
        return "GLPA"
    if resname in PROTEIN_RESNAMES:
        return "PROTEIN"
    return MOLECULE_ALIASES.get(resname, resname)


def _fit_atoms_to_reference(mobile_atoms: list[dict[str, Any]], reference_atoms: list[dict[str, Any]]) -> dict[str, float]:
    reference_by_key: dict[tuple[str, int, str, str], deque[dict[str, Any]]] = defaultdict(deque)
    for atom in reference_atoms:
        if atom["name"].startswith("H"):
            continue
        reference_by_key[(atom["chain"], atom["resid"], atom["icode"], atom["name"])].append(atom)

    mobile_points: list[list[float]] = []
    reference_points: list[list[float]] = []
    for atom in mobile_atoms:
        if atom["name"].startswith("H"):
            continue
        matches = reference_by_key.get((atom["chain"], atom["resid"], atom["icode"], atom["name"]))
        if matches:
            ref = matches.popleft()
            mobile_points.append([atom["x"], atom["y"], atom["z"]])
            reference_points.append([ref["x"], ref["y"], ref["z"]])

    if len(mobile_points) < 3:
        raise Stage4Error("cannot fit pdb2gmx protein to Stage 3 protein: fewer than three common heavy atoms")
    mobile = np.array(mobile_points, dtype=float)
    reference = np.array(reference_points, dtype=float)
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    centered_mobile = mobile - mobile_center
    centered_reference = reference - reference_center
    u, _, vt = np.linalg.svd(centered_mobile.T @ centered_reference)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0.0:
        vt[-1, :] *= -1.0
        rotation = vt.T @ u.T
    before = float(np.sqrt(np.mean(np.sum((mobile - reference) ** 2, axis=1))))
    after_points = centered_mobile @ rotation + reference_center
    after = float(np.sqrt(np.mean(np.sum((after_points - reference) ** 2, axis=1))))
    for atom in mobile_atoms:
        xyz = (np.array([atom["x"], atom["y"], atom["z"]], dtype=float) - mobile_center) @ rotation + reference_center
        atom["x"], atom["y"], atom["z"] = (float(xyz[0]), float(xyz[1]), float(xyz[2]))
    return {"fit_atoms": len(mobile_points), "rmsd_before_a": before, "rmsd_after_a": after}


def _reuse_reference_protein_coordinates(
    topology_atoms: list[dict[str, Any]],
    reference_atoms: list[dict[str, Any]],
) -> tuple[dict[str, int], list[int]]:
    reference_by_key: dict[tuple[str, int, str, str], deque[dict[str, Any]]] = defaultdict(deque)
    for atom in reference_atoms:
        reference_by_key[(atom["chain"], atom["resid"], atom["icode"], atom["name"])].append(atom)
    aliases = {"H1": "HT1", "H2": "HT2", "H3": "HT3"}
    reused = 0
    fallback = 0
    fallback_indices: list[int] = []
    for atom_index, atom in enumerate(topology_atoms, start=1):
        keys = [(atom["chain"], atom["resid"], atom["icode"], atom["name"])]
        if atom["name"] in aliases:
            keys.append((atom["chain"], atom["resid"], atom["icode"], aliases[atom["name"]]))
        ref = None
        for key in keys:
            matches = reference_by_key.get(key)
            if matches:
                ref = matches.popleft()
                break
        if ref is None:
            fallback += 1
            fallback_indices.append(atom_index)
            continue
        atom["x"], atom["y"], atom["z"] = ref["x"], ref["y"], ref["z"]
        reused += 1
    return {"reused_stage3_protein_coordinates": reused, "fitted_fallback_protein_coordinates": fallback}, fallback_indices


def _assign_generated_protein_identity(generated_atoms: list[dict[str, Any]], reference_atoms: list[dict[str, Any]]) -> None:
    chains_by_resid: dict[int, set[str]] = defaultdict(set)
    residue_identity: dict[int, tuple[str, str]] = {}
    for atom in reference_atoms:
        chains_by_resid[atom["resid"]].add(atom["chain"])
        residue_identity.setdefault(atom["resid"], (atom["chain"], atom["icode"]))
    ambiguous = sorted(resid for resid, chains in chains_by_resid.items() if len(chains) > 1)
    if ambiguous:
        raise Stage4Error(
            "pdb2gmx GRO lacks chain identifiers and Stage 3 protein has duplicate residue numbers across chains: "
            f"{ambiguous[:20]}; provide uniquely numbered chains or an explicit multichain mapping"
        )
    for atom in generated_atoms:
        identity = residue_identity.get(atom["resid"])
        if identity:
            atom["chain"], atom["icode"] = identity


def _protein_hydrogen_parents(protein_top: Path) -> dict[int, int]:
    atom_names: dict[int, str] = {}
    bonds: list[tuple[int, int]] = []
    section: str | None = None
    for raw in protein_top.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line.strip("[] ").lower()
            continue
        parts = line.split()
        if section == "atoms" and len(parts) >= 5 and parts[0].isdigit():
            atom_names[int(parts[0])] = parts[4]
        elif section == "bonds" and len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
            bonds.append((int(parts[0]), int(parts[1])))
    parents: dict[int, int] = {}
    for left, right in bonds:
        left_name = atom_names.get(left, "")
        right_name = atom_names.get(right, "")
        if left_name.startswith("H") and right_name and not right_name.startswith("H"):
            parents[left] = right
        elif right_name.startswith("H") and left_name and not left_name.startswith("H"):
            parents[right] = left
    return parents


def _place_fallback_hydrogens_from_parents(
    topology_atoms: list[dict[str, Any]],
    template_atoms: list[dict[str, Any]],
    parent_by_atom: dict[int, int],
    fallback_indices: list[int],
) -> dict[str, int]:
    placed = 0
    left_as_fit = 0
    for atom_index in fallback_indices:
        atom = topology_atoms[atom_index - 1]
        parent_index = parent_by_atom.get(atom_index)
        if not atom["name"].startswith("H") or parent_index is None:
            left_as_fit += 1
            continue
        template_atom = template_atoms[atom_index - 1]
        template_parent = template_atoms[parent_index - 1]
        parent = topology_atoms[parent_index - 1]
        vector = np.array([
            template_atom["x"] - template_parent["x"],
            template_atom["y"] - template_parent["y"],
            template_atom["z"] - template_parent["z"],
        ], dtype=float)
        length = float(np.linalg.norm(vector))
        if not np.isfinite(length) or length <= 0.0 or length > 2.0:
            left_as_fit += 1
            continue
        atom["x"] = parent["x"] + float(vector[0])
        atom["y"] = parent["y"] + float(vector[1])
        atom["z"] = parent["z"] + float(vector[2])
        placed += 1
    return {"parent_placed_fallback_hydrogens": placed, "unplaced_fallback_protein_coordinates": left_as_fit}


def replace_protein_with_pdb2gmx_coordinates(
    input_pdb: Path,
    protein_gro: Path,
    output_pdb: Path,
    membrane_toppar: Path | None = None,
    protein_top: Path | None = None,
) -> dict[str, int]:
    prefix, atoms, _ = _pdb_atoms(input_pdb)
    box_a = _box_from_prefix(prefix)
    blocks = _molecule_blocks(atoms)
    reference_protein = [atom for block in blocks if _molecule_for_block(block) == "PROTEIN" for atom in block]
    protein_chains = [block[0]["chain"] for block in blocks if _molecule_for_block(block) == "PROTEIN"]
    chain = protein_chains[0] if protein_chains else "X"
    generated_protein, _ = _gro_atoms(protein_gro, chain=chain)
    generated_caps = sorted({atom["resname"] for atom in generated_protein if atom["resname"] in CHARMM_EXPLICIT_TERMINAL_CAP_RESNAMES})
    if generated_caps:
        raise Stage4Error(
            f"pdb2gmx generated explicit terminal cap residues {generated_caps}; "
            "Stage 4 requires CHARMM ACE/CT3 terminal patches on terminal amino-acid residues, not standalone cap residues"
        )
    _assign_generated_protein_identity(generated_protein, reference_protein)
    generated_template = [atom.copy() for atom in generated_protein]
    fit_report = _fit_atoms_to_reference(generated_protein, reference_protein)
    coordinate_report, fallback_indices = _reuse_reference_protein_coordinates(generated_protein, reference_protein)
    protein_top = protein_top or protein_gro.with_name("protein.top")
    parent_by_atom = _protein_hydrogen_parents(protein_top) if protein_top.exists() else {}
    fallback_report = _place_fallback_hydrogens_from_parents(
        generated_protein,
        generated_template,
        parent_by_atom,
        fallback_indices,
    )
    membrane_blocks = [block for block in blocks if _molecule_for_block(block) != "PROTEIN"]
    if membrane_toppar is not None:
        itp_by_molecule = {_moleculetype(path): path for path in membrane_toppar.glob("*.itp")}
        grouped: dict[str, list[list[dict[str, Any]]]] = defaultdict(list)
        for block in membrane_blocks:
            molecule = _molecule_for_block(block)
            if molecule is None:
                continue
            if molecule not in itp_by_molecule:
                raise Stage4Error(f"missing all-atom ITP for molecule: {molecule}")
            grouped[molecule].append(_reorder_block_to_names(block, _itp_atom_order(itp_by_molecule[molecule]), molecule))
        membrane_atoms = [atom for molecule in sorted(grouped) for block in grouped[molecule] for atom in block]
    else:
        membrane_atoms = [atom for block in membrane_blocks for atom in block]
    _write_pdb(output_pdb, prefix, [*generated_protein, *membrane_atoms], box_a)
    return {
        "replaced_stage3_protein_atoms": len(atoms) - len(membrane_atoms),
        "pdb2gmx_protein_atoms": len(generated_protein),
        "retained_membrane_atoms": len(membrane_atoms),
        **fit_report,
        **coordinate_report,
        **fallback_report,
    }


def slice_membrane(input_pdb: Path, output_pdb: Path, buffer_nm: float, z_slab_padding_nm: float = 1.5) -> dict[str, Any]:
    if not np.isfinite(z_slab_padding_nm) or z_slab_padding_nm <= 0.0:
        raise Stage4Error("Stage 4 z_slab_padding_nm must be positive")
    prefix, atoms, _ = _pdb_atoms(input_pdb)
    box_a = _box_from_prefix(prefix)
    blocks = _molecule_blocks(atoms)
    protein = [atom for block in blocks if _molecule_for_block(block) == "PROTEIN" for atom in block]
    if not protein:
        raise Stage4Error("cannot slice Stage 4 system: no protein atoms found")
    protein_xyz = np.array([[atom["x"], atom["y"], atom["z"]] for atom in protein])
    buffer_a = buffer_nm * 10.0
    lower = np.array([
        max(0.0, float(protein_xyz[:, 0].min() - buffer_a)),
        max(0.0, float(protein_xyz[:, 1].min() - buffer_a)),
        0.0,
    ])
    upper = np.array([
        min(float(box_a[0]), float(protein_xyz[:, 0].max() + buffer_a)),
        min(float(box_a[1]), float(protein_xyz[:, 1].max() + buffer_a)),
        float(box_a[2]),
    ])
    crop_box = upper - lower
    if np.any(crop_box <= 0.0):
        raise Stage4Error("invalid Stage 4 crop dimensions")

    kept_blocks: list[tuple[str, list[dict[str, Any]]]] = []
    counts: Counter[str] = Counter()
    protein_chains: set[str] = set()
    removed: Counter[str] = Counter()
    for block in blocks:
        molecule = _molecule_for_block(block)
        if molecule is None:
            removed[block[0]["resname"]] += 1
            continue
        keep = molecule == "PROTEIN"
        if not keep:
            ref = next((atom for atom in block if atom["name"].startswith("P")), block[0])
            ref_xyz = np.array([ref["x"], ref["y"], ref["z"]])
            block_xyz = np.array([[atom["x"], atom["y"], atom["z"]] for atom in block])
            ref_inside = bool(np.all(ref_xyz[:2] >= lower[:2]) and np.all(ref_xyz[:2] <= upper[:2]))
            block_inside = bool(np.all(block_xyz[:, :2] >= lower[:2]) and np.all(block_xyz[:, :2] <= upper[:2]))
            keep = ref_inside and block_inside
        if keep:
            kept_blocks.append((molecule, block))
            if molecule == "PROTEIN":
                protein_chains.add(block[0]["chain"])
            else:
                counts[molecule] += 1
        else:
            removed[molecule] += 1
    kept_source = [atom for _, block in kept_blocks for atom in block]
    kept_xyz = np.array([[atom["x"], atom["y"], atom["z"]] for atom in kept_source])
    z_padding_a = z_slab_padding_nm * 10.0
    lower[2] = float(kept_xyz[:, 2].min() - z_padding_a)
    upper[2] = float(kept_xyz[:, 2].max() + z_padding_a)
    crop_box = upper - lower
    if np.any(crop_box <= 0.0):
        raise Stage4Error("invalid Stage 4 crop dimensions after molecule-boundary filtering")
    kept: list[dict[str, Any]] = []
    for atom in kept_source:
        shifted = atom.copy()
        shifted["x"] -= lower[0]
        shifted["y"] -= lower[1]
        shifted["z"] -= lower[2]
        kept.append(shifted)
    _write_pdb(output_pdb, prefix, kept, crop_box)
    shifted_z = np.array([atom["z"] for atom in kept], dtype=float)
    return {
        "input_atoms": len(atoms), "output_atoms": len(kept), "crop_box_nm": [round(float(x / 10.0), 6) for x in crop_box],
        "buffer_nm": buffer_nm, "z_slab_padding_nm": z_slab_padding_nm,
        "z_lower_clearance_nm": round(float(shifted_z.min() / 10.0), 6),
        "z_upper_clearance_nm": round(float((crop_box[2] - shifted_z.max()) / 10.0), 6),
        "original_box_z_nm": round(float(box_a[2] / 10.0), 6),
        "kept_molecules": dict(sorted({**counts, "PROTEIN": len(protein_chains)}.items())),
        "removed_molecules": dict(sorted(removed.items())),
    }


def _moleculetype(itp: Path) -> str:
    in_section = False
    for raw in itp.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.split(";", 1)[0].strip()
        if line.lower() == "[ moleculetype ]":
            in_section = True
            continue
        if in_section and line and not line.startswith("["):
            return line.split()[0]
    raise Stage4Error(f"cannot determine moleculetype from {itp}")


def _itp_atom_order(itp: Path) -> list[str]:
    names: list[str] = []
    in_atoms = False
    for raw in itp.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.lower() == "[ atoms ]":
            in_atoms = True
            continue
        if in_atoms and line.startswith("["):
            break
        if in_atoms:
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                names.append(parts[4])
    if not names:
        raise Stage4Error(f"cannot determine atom order from {itp}")
    return names


def _reorder_block_to_names(block: list[dict[str, Any]], target_names: list[str], molecule: str) -> list[dict[str, Any]]:
    current_names = [atom["name"] for atom in block]
    if Counter(current_names) != Counter(target_names):
        raise Stage4Error(f"{molecule} atom inventory does not match topology atom order")
    by_name: dict[str, deque[dict[str, Any]]] = defaultdict(deque)
    for atom in block:
        by_name[atom["name"]].append(atom)
    return [by_name[name].popleft() for name in target_names]


def _protein_topology_body(protein_top: Path) -> list[str]:
    lines = protein_top.read_text(encoding="utf-8", errors="replace").splitlines()
    start = next((i for i, raw in enumerate(lines) if raw.strip().lower() == "[ moleculetype ]"), None)
    if start is None:
        return []
    stop_markers = ("; include water topology", "; include topology for ions")
    body: list[str] = []
    for raw in lines[start:]:
        stripped = raw.strip().lower()
        if stripped == "[ system ]" or stripped in stop_markers:
            break
        if stripped.startswith("#include") and (
            "charmm36.ff" in stripped
            or "tip3p.itp" in stripped
            or "ions.itp" in stripped
        ):
            continue
        body.append(raw)
    return body


def write_topology(
    output_top: Path,
    protein_top: Path,
    protein_itps: list[Path],
    membrane_toppar: Path,
    molecule_counts: dict[str, int],
) -> None:
    itp_by_molecule = {_moleculetype(path): path for path in membrane_toppar.glob("*.itp")}
    missing = sorted(set(molecule_counts) - {"PROTEIN"} - set(itp_by_molecule))
    if missing:
        raise Stage4Error(f"missing all-atom ITPs for molecules: {missing}")
    protein_molecules = []
    in_molecules = False
    for raw in protein_top.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.strip().lower() == "[ molecules ]":
            in_molecules = True
            continue
        if in_molecules:
            line = raw.split(";", 1)[0].strip()
            if line and not line.startswith("["):
                protein_molecules.append(line.split()[0])
    if molecule_counts.get("PROTEIN", 0) != len(protein_molecules):
        raise Stage4Error(
            f"protein molecule count mismatch: finalized PDB has {molecule_counts.get('PROTEIN', 0)}, "
            f"pdb2gmx generated {len(protein_molecules)}"
        )
    lines = ['#include "charmm36.ff/forcefield.itp"']
    lines.extend(f'#include "{path.name}"' for path in protein_itps)
    if not protein_itps:
        protein_body = _protein_topology_body(protein_top)
        if not protein_body:
            raise Stage4Error(f"pdb2gmx generated no protein topology body in {protein_top}")
        lines.extend(["", *protein_body])
    lines.extend(f'#include "toppar/{itp_by_molecule[name].name}"' for name in sorted(itp_by_molecule))
    lines.extend(["", "[ system ]", "MembraneForger Stage 4", "", "[ molecules ]"])
    lines.extend(f"{name:<20} 1" for name in protein_molecules)
    ordered_names = [
        *MEMBRANE_MOLECULE_ORDER,
        *sorted(set(molecule_counts) - {"PROTEIN", *MEMBRANE_MOLECULE_ORDER, *SOLVENT_ION_MOLECULE_ORDER}),
        *SOLVENT_ION_MOLECULE_ORDER,
    ]
    for name in ordered_names:
        count = molecule_counts.get(name, 0)
        if count:
            lines.append(f"{name:<20} {count}")
    output_top.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _ensure_charmm_r2b_aliases(charmm36: Path) -> None:
    r2b = charmm36 / "merged.r2b"
    if not r2b.exists():
        return
    aliases = {
        "ACE": "ACE",
        "CT3": "CT3",
        "HSD": "HSD",
        "HSE": "HSE",
        "HSP": "HSP",
    }
    rows = r2b.read_text(encoding="utf-8", errors="replace").splitlines()
    existing = {row.split()[0] for row in rows if row.strip() and not row.lstrip().startswith(";")}
    additions = [f"{source}\t{target}" for source, target in aliases.items() if source not in existing]
    if additions:
        r2b.write_text("\n".join([*rows, *additions]) + "\n", encoding="utf-8")


def copy_stage4_resources(charmm36: Path, membrane_toppar: Path, work: Path) -> None:
    if not (charmm36 / "forcefield.itp").exists():
        raise Stage4Error(f"CHARMM36 force field is incomplete: {charmm36}")
    if not membrane_toppar.exists():
        raise Stage4Error(f"membrane topology directory is missing: {membrane_toppar}")
    shutil.copytree(charmm36, work / "charmm36.ff")
    _ensure_charmm_r2b_aliases(work / "charmm36.ff")
    shutil.copytree(membrane_toppar, work / "toppar")


def molecule_counts_from_pdb(path: Path) -> dict[str, int]:
    _, atoms, _ = _pdb_atoms(path)
    counts: Counter[str] = Counter()
    protein_chains: set[str] = set()
    for block in _molecule_blocks(atoms):
        molecule = _molecule_for_block(block)
        if molecule == "PROTEIN":
            protein_chains.add(block[0]["chain"])
        elif molecule is not None:
            counts[molecule] += 1
    counts["PROTEIN"] = len(protein_chains)
    return dict(counts)


def membrane_z_bounds_nm(path: Path) -> tuple[float, float]:
    _, atoms, _ = _pdb_atoms(path)
    membrane_z = [
        atom["z"] / 10.0
        for block in _molecule_blocks(atoms)
        if _molecule_for_block(block) not in {None, "PROTEIN"}
        for atom in block
    ]
    if not membrane_z:
        raise Stage4Error(f"cannot determine membrane z bounds from {path}")
    return min(membrane_z), max(membrane_z)


def phosphate_headgroup_core_z_bounds_nm(path: Path) -> tuple[float, float]:
    _, atoms, _ = _pdb_atoms(path)
    headgroup_z = sorted(
        atom["z"] / 10.0
        for block in _molecule_blocks(atoms)
        if _molecule_for_block(block) not in {None, "PROTEIN"}
        for atom in block
        if _is_representative_headgroup_atom(_molecule_for_block(block), atom["name"])
    )
    if len(headgroup_z) < 4:
        raise Stage4Error(f"cannot determine representative headgroup z bounds from {path}")
    gaps = [(headgroup_z[index + 1] - headgroup_z[index], index) for index in range(len(headgroup_z) - 1)]
    gap, split_index = max(gaps, key=lambda item: item[0])
    if gap <= 0.5:
        raise Stage4Error(f"cannot separate representative headgroup leaflets in {path}")
    return headgroup_z[split_index], headgroup_z[split_index + 1]


def _is_representative_headgroup_atom(molecule: str | None, atom_name: str) -> bool:
    if molecule is None:
        return False
    component = component_for(molecule)
    if component and component.headgroup_beads and atom_name in component.headgroup_beads:
        return True
    return atom_name in {"P", "PO4", "PO3", "ROH"}


def normalize_solvated_gro(
    input_gro: Path,
    output_gro: Path,
    excluded_water_z_nm: tuple[float, float] | None = None,
) -> dict[str, int | float | None]:
    lines = input_gro.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        raise Stage4Error(f"invalid solvated GRO file: {input_gro}")
    natoms = int(lines[1].strip())
    if len(lines) != natoms + 3:
        raise Stage4Error(f"solvated GRO atom count does not match file length: {input_gro}")
    rename = {"OW": "OH2", "HW1": "H1", "HW2": "H2"}
    output_atoms: list[str] = []
    waters = 0
    removed_waters = 0
    index = 0
    atom_lines = lines[2:2 + natoms]
    while index < len(atom_lines):
        raw = atom_lines[index]
        resname = raw[5:10].strip()
        if resname != "SOL":
            output_atoms.append(raw)
            index += 1
            continue
        resid = raw[0:5]
        water_block: list[str] = []
        while index < len(atom_lines) and atom_lines[index][0:5] == resid and atom_lines[index][5:10].strip() == "SOL":
            water_block.append(atom_lines[index])
            index += 1
        oxygen = next((item for item in water_block if item[10:15].strip() == "OW"), None)
        if oxygen is None:
            raise Stage4Error(f"solvated water residue {resid.strip()} has no OW atom")
        oxygen_z = float(oxygen[36:44])
        if excluded_water_z_nm is not None and excluded_water_z_nm[0] <= oxygen_z <= excluded_water_z_nm[1]:
            removed_waters += 1
            continue
        waters += 1
        for item in water_block:
            name = item[10:15].strip()
            output_atoms.append(f"{item[:5]}{'TIP3':>5}{rename.get(name, name):>5}{item[15:]}")
    output = [lines[0], f"{len(output_atoms):5d}", *output_atoms]
    output.append(lines[-1])
    output_gro.write_text("\n".join(output) + "\n", encoding="utf-8")
    lower = excluded_water_z_nm[0] if excluded_water_z_nm is not None else None
    upper = excluded_water_z_nm[1] if excluded_water_z_nm is not None else None
    return {
        "waters": waters,
        "removed_waters_in_membrane_z": removed_waters,
        "excluded_z_min_nm": lower,
        "excluded_z_max_nm": upper,
    }


def filter_solvent_ions_by_z(
    input_gro: Path,
    output_gro: Path,
    excluded_z_nm: tuple[float, float],
) -> dict[str, int | float]:
    lines = input_gro.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        raise Stage4Error(f"invalid GRO file: {input_gro}")
    natoms = int(lines[1].strip())
    if len(lines) != natoms + 3:
        raise Stage4Error(f"GRO atom count does not match file length: {input_gro}")
    output_records: list[list[str]] = []
    removed_waters = 0
    removed_sod = 0
    removed_cla = 0
    removed_other_ions = 0
    index = 0
    atom_lines = lines[2:2 + natoms]
    water_resnames = {"TIP3", "SOL", "W", "WF"}
    ion_resnames = {"SOD", "CLA", "NA", "CL", "ION"}
    oxygen_names = {"OH2", "OW", "O"}
    outside_sod: list[tuple[float, int]] = []
    outside_cla: list[tuple[float, int]] = []
    while index < len(atom_lines):
        raw = atom_lines[index]
        resname = raw[5:10].strip()
        if resname in water_resnames:
            resid = raw[0:5]
            water_block: list[str] = []
            while index < len(atom_lines) and atom_lines[index][0:5] == resid and atom_lines[index][5:10].strip() in water_resnames:
                water_block.append(atom_lines[index])
                index += 1
            oxygen = next((item for item in water_block if item[10:15].strip() in oxygen_names), water_block[0])
            oxygen_z = float(oxygen[36:44])
            if excluded_z_nm[0] <= oxygen_z <= excluded_z_nm[1]:
                removed_waters += 1
                continue
            output_records.append(water_block)
            continue
        if resname in ion_resnames:
            ion_z = float(raw[36:44])
            if excluded_z_nm[0] <= ion_z <= excluded_z_nm[1]:
                if resname in {"SOD", "NA"}:
                    removed_sod += 1
                elif resname in {"CLA", "CL"}:
                    removed_cla += 1
                else:
                    removed_other_ions += 1
                index += 1
                continue
            record_index = len(output_records)
            output_records.append([raw])
            distance_to_slab = min(abs(ion_z - excluded_z_nm[0]), abs(ion_z - excluded_z_nm[1]))
            if resname in {"SOD", "NA"}:
                outside_sod.append((distance_to_slab, record_index))
            elif resname in {"CLA", "CL"}:
                outside_cla.append((distance_to_slab, record_index))
            index += 1
            continue
        output_records.append([raw])
        index += 1
    balanced_sod = 0
    balanced_cla = 0
    drop_records: set[int] = set()
    if removed_sod > removed_cla:
        need = removed_sod - removed_cla
        candidates = sorted(outside_cla)[:need]
        if len(candidates) != need:
            raise Stage4Error("cannot preserve charge after membrane-z ion cleanup: insufficient outside CLA counterions")
        drop_records.update(record_index for _, record_index in candidates)
        balanced_cla = need
    elif removed_cla > removed_sod:
        need = removed_cla - removed_sod
        candidates = sorted(outside_sod)[:need]
        if len(candidates) != need:
            raise Stage4Error("cannot preserve charge after membrane-z ion cleanup: insufficient outside SOD counterions")
        drop_records.update(record_index for _, record_index in candidates)
        balanced_sod = need
    output_atoms = [line for record_index, record in enumerate(output_records) if record_index not in drop_records for line in record]
    output = [lines[0], f"{len(output_atoms):5d}", *output_atoms, lines[-1]]
    output_gro.write_text("\n".join(output) + "\n", encoding="utf-8")
    return {
        "removed_waters_in_membrane_z": removed_waters,
        "removed_sod_in_membrane_z": removed_sod,
        "removed_cla_in_membrane_z": removed_cla,
        "removed_other_ions_in_membrane_z": removed_other_ions,
        "removed_sod_for_charge_balance": balanced_sod,
        "removed_cla_for_charge_balance": balanced_cla,
        "excluded_z_min_nm": excluded_z_nm[0],
        "excluded_z_max_nm": excluded_z_nm[1],
    }


def solvent_ion_counts_from_gro(path: Path) -> dict[str, int]:
    atoms, _ = _gro_atoms(path)
    counts: Counter[str] = Counter()
    for atom in atoms:
        resname = atom["resname"]
        name = atom["name"]
        if resname in {"TIP3", "SOL", "W", "WF"} and name in {"OH2", "OW", "O"}:
            counts["TIP3"] += 1
        elif resname in {"SOD", "NA"}:
            counts["SOD"] += 1
        elif resname in {"CLA", "CL"}:
            counts["CLA"] += 1
    return dict(counts)


def hard_clash_report_gro(path: Path, cutoff_nm: float = 0.06) -> dict[str, Any]:
    atoms, box = _gro_atoms(path)
    box_nm = box / 10.0
    heavy: list[tuple[int, str, str, np.ndarray]] = []
    for index, atom in enumerate(atoms, start=1):
        name = atom["name"]
        if name.startswith("H"):
            continue
        xyz = np.array([atom["x"], atom["y"], atom["z"]], dtype=float) / 10.0
        heavy.append((index, atom["resname"], name, xyz))
    cells: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    cell_size = cutoff_nm
    dims = np.maximum(np.floor(box_nm / cell_size).astype(int), 1)
    for heavy_index, (_, _, _, xyz) in enumerate(heavy):
        wrapped = np.mod(xyz, box_nm)
        cell = tuple(np.floor(wrapped / box_nm * dims).astype(int))
        cells[cell].append(heavy_index)
    clashes = 0
    min_dist = float("inf")
    examples: list[dict[str, Any]] = []
    offsets = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
    for cell, occupants in cells.items():
        neighbor_indices: list[int] = []
        for offset in offsets:
            neighbor = tuple((cell[axis] + offset[axis]) % int(dims[axis]) for axis in range(3))
            neighbor_indices.extend(cells.get(neighbor, []))
        for left in occupants:
            left_atom = heavy[left]
            for right in neighbor_indices:
                if right <= left:
                    continue
                right_atom = heavy[right]
                delta = right_atom[3] - left_atom[3]
                delta -= box_nm * np.round(delta / box_nm)
                dist = float(np.linalg.norm(delta))
                min_dist = min(min_dist, dist)
                if dist < cutoff_nm:
                    clashes += 1
                    if len(examples) < 10:
                        examples.append({
                            "distance_nm": round(dist, 5),
                            "left": [left_atom[0], left_atom[1], left_atom[2]],
                            "right": [right_atom[0], right_atom[1], right_atom[2]],
                        })
    return {
        "cutoff_nm": cutoff_nm,
        "heavy_atoms_checked": len(heavy),
        "minimum_heavy_atom_distance_nm": None if not np.isfinite(min_dist) else round(min_dist, 5),
        "hard_clashes": clashes,
        "examples": examples,
    }


def write_index_ini(gro_path: Path, output_ndx: Path) -> dict[str, int]:
    atoms, _ = _gro_atoms(gro_path)
    groups: dict[str, list[int]] = {
        "System": [],
        "Protein": [],
        "Membrane": [],
        "Water": [],
        "Ions": [],
        "Water_and_ions": [],
    }
    for index, atom in enumerate(atoms, start=1):
        resname = atom["resname"]
        groups["System"].append(index)
        if resname in PROTEIN_RESNAMES:
            groups["Protein"].append(index)
        elif resname in {"TIP3", "SOL", "W", "WF"}:
            groups["Water"].append(index)
            groups["Water_and_ions"].append(index)
        elif resname in {"SOD", "CLA", "NA", "CL", "ION"}:
            groups["Ions"].append(index)
            groups["Water_and_ions"].append(index)
        else:
            groups["Membrane"].append(index)

    lines: list[str] = []
    for name, indices in groups.items():
        lines.append(f"[ {name} ]")
        for start in range(0, len(indices), 15):
            lines.append(" ".join(str(value) for value in indices[start:start + 15]))
        lines.append("")
    output_ndx.write_text("\n".join(lines), encoding="utf-8")
    return {name: len(indices) for name, indices in groups.items()}
