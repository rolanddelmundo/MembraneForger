from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import aa_stage4


def test_pdb_round_trip_preserves_four_character_glpa(tmp_path: Path) -> None:
    path = tmp_path / "glpa.pdb"
    atoms = [{
        "record": "HETATM", "name": "C1S", "resname": "GLPA", "chain": "A", "resid": 10,
        "icode": "", "x": 1.0, "y": 2.0, "z": 3.0, "element": "C",
    }]
    aa_stage4._write_pdb(path, [], atoms, np.array([50.0, 60.0, 70.0]))
    _, parsed, _ = aa_stage4._pdb_atoms(path)
    assert parsed[0]["resname"] == "GLPA"
    assert aa_stage4.molecule_counts_from_pdb(path) == {"GLPA": 1, "PROTEIN": 0}


def test_pdb_writer_keeps_columns_after_atom_serial_overflow(tmp_path: Path) -> None:
    path = tmp_path / "large.pdb"
    atoms = [
        {
            "record": "ATOM", "name": "CA", "resname": "POPE", "chain": "S", "resid": 1091,
            "icode": "", "x": 1.0, "y": 2.0, "z": 3.0, "element": "C",
        }
        for _ in range(100001)
    ]
    aa_stage4._write_pdb(path, [], atoms, np.array([50.0, 60.0, 70.0]))
    _, parsed, _ = aa_stage4._pdb_atoms(path)
    assert len(parsed) == 100001
    assert parsed[-1]["resname"] == "POPE"
    assert parsed[-1]["chain"] == "S"
    assert parsed[-1]["resid"] == 1091


def test_write_topology_embeds_inline_protein_without_builtin_includes(tmp_path: Path) -> None:
    protein_top = tmp_path / "protein.top"
    protein_top.write_text(
        "\n".join([
            '#include "./charmm36.ff/forcefield.itp"',
            "",
            "[ moleculetype ]",
            "Protein_chain_X 3",
            "",
            "[ atoms ]",
            "1 CT1 1 ALA CA 1 0.0 12.011",
            "",
            "; Include water topology",
            '#include "./charmm36.ff/tip3p.itp"',
            '#include "./charmm36.ff/ions.itp"',
            "",
            "[ system ]",
            "protein",
            "",
            "[ molecules ]",
            "Protein_chain_X 1",
            "",
        ]),
        encoding="utf-8",
    )
    toppar = tmp_path / "toppar"
    toppar.mkdir()
    (toppar / "POPE.itp").write_text("[ moleculetype ]\nPOPE 2\n", encoding="utf-8")
    output = tmp_path / "topol.top"
    aa_stage4.write_topology(output, protein_top, [], toppar, {"PROTEIN": 1, "POPE": 2})
    text = output.read_text(encoding="utf-8")
    assert "[ moleculetype ]\nProtein_chain_X 3" in text
    assert '#include "tip3p.itp"' not in text
    assert '#include "ions.itp"' not in text
    assert '#include "toppar/POPE.itp"' in text


def test_replace_protein_with_pdb2gmx_coordinates_keeps_membrane(tmp_path: Path) -> None:
    pdb = tmp_path / "dry.pdb"
    aa_stage4._write_pdb(
        pdb,
        [],
        [
            {
                "record": "ATOM", "name": "CA", "resname": "ALA", "chain": "A", "resid": 1,
                "icode": "", "x": 1.0, "y": 2.0, "z": 3.0, "element": "C",
            },
            {
                "record": "ATOM", "name": "N", "resname": "ALA", "chain": "A", "resid": 1,
                "icode": "", "x": 0.5, "y": 2.0, "z": 3.0, "element": "N",
            },
            {
                "record": "ATOM", "name": "C", "resname": "ALA", "chain": "A", "resid": 1,
                "icode": "", "x": 1.5, "y": 2.0, "z": 3.0, "element": "C",
            },
            {
                "record": "HETATM", "name": "P", "resname": "POPE", "chain": "B", "resid": 10,
                "icode": "", "x": 4.0, "y": 5.0, "z": 6.0, "element": "P",
            },
        ],
        np.array([50.0, 60.0, 70.0]),
    )
    gro = tmp_path / "protein.gro"
    gro.write_text(
        "\n".join([
            "protein",
            "3",
            "    1ALA      N    1   0.700   0.800   0.900",
            "    1ALA     CA    2   0.750   0.800   0.900",
            "    1ALA      C    3   0.800   0.800   0.900",
            "   5.00000   6.00000   7.00000",
        ]) + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "matched.pdb"
    report = aa_stage4.replace_protein_with_pdb2gmx_coordinates(pdb, gro, out)
    _, atoms, _ = aa_stage4._pdb_atoms(out)
    assert report["replaced_stage3_protein_atoms"] == 3
    assert report["pdb2gmx_protein_atoms"] == 3
    assert report["retained_membrane_atoms"] == 1
    assert report["fit_atoms"] == 3
    assert report["rmsd_after_a"] < 1e-6
    assert [atom["name"] for atom in atoms] == ["N", "CA", "C", "P"]
    assert atoms[-1]["resname"] == "POPE"


def test_replace_protein_rejects_ambiguous_multichain_duplicate_resids(tmp_path: Path) -> None:
    pdb = tmp_path / "dry.pdb"
    atoms = []
    for chain in ("A", "B"):
        for name, x in (("N", 0.0), ("CA", 1.0), ("C", 2.0)):
            atoms.append({
                "record": "ATOM", "name": name, "resname": "ALA", "chain": chain, "resid": 1,
                "icode": "", "x": x, "y": 2.0 if chain == "A" else 5.0, "z": 3.0, "element": name[0],
            })
    aa_stage4._write_pdb(pdb, [], atoms, np.array([50.0, 60.0, 70.0]))
    gro = tmp_path / "protein.gro"
    gro.write_text(
        "\n".join([
            "protein",
            "6",
            "    1ALA      N    1   0.000   0.000   0.000",
            "    1ALA     CA    2   0.100   0.000   0.000",
            "    1ALA      C    3   0.200   0.000   0.000",
            "    1ALA      N    4   0.000   0.300   0.000",
            "    1ALA     CA    5   0.100   0.300   0.000",
            "    1ALA      C    6   0.200   0.300   0.000",
            "   5.00000   6.00000   7.00000",
        ])
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(aa_stage4.Stage4Error, match="duplicate residue numbers across chains"):
        aa_stage4.replace_protein_with_pdb2gmx_coordinates(pdb, gro, tmp_path / "matched.pdb")


def test_slice_membrane_crops_z_to_retained_complex_with_padding(tmp_path: Path) -> None:
    source = tmp_path / "aligned.pdb"
    aa_stage4._write_pdb(
        source,
        [],
        [
            {
                "record": "ATOM", "name": "CA", "resname": "ALA", "chain": "A", "resid": 1,
                "icode": "", "x": 50.0, "y": 50.0, "z": 10.0, "element": "C",
            },
            {
                "record": "HETATM", "name": "P", "resname": "POPE", "chain": "B", "resid": 10,
                "icode": "", "x": 55.0, "y": 55.0, "z": 90.0, "element": "P",
            },
            {
                "record": "HETATM", "name": "P", "resname": "DOPC", "chain": "C", "resid": 11,
                "icode": "", "x": 95.0, "y": 95.0, "z": 50.0, "element": "P",
            },
        ],
        np.array([100.0, 100.0, 120.0]),
    )
    out = tmp_path / "sliced.pdb"
    report = aa_stage4.slice_membrane(source, out, 1.0)
    prefix, atoms, _ = aa_stage4._pdb_atoms(out)
    box = aa_stage4._box_from_prefix(prefix)
    assert box[2] == 110.0
    assert report["z_slab_padding_nm"] == 1.5
    assert report["z_lower_clearance_nm"] == 1.5
    assert report["z_upper_clearance_nm"] == 1.5
    assert report["kept_molecules"]["POPE"] == 1
    assert report["removed_molecules"]["DOPC"] == 1
    assert {atom["resname"] for atom in atoms} == {"ALA", "POPE"}
    assert min(atom["z"] for atom in atoms) == 15.0


def test_slice_keeps_glpa_components_as_complete_molecule(tmp_path: Path) -> None:
    source = tmp_path / "components.pdb"
    atoms = [{
        "record": "ATOM", "name": "CA", "resname": "ALA", "chain": "A", "resid": 1,
        "icode": "", "x": 50.0, "y": 50.0, "z": 60.0, "element": "C",
    }]
    for resid, resname in enumerate(["CER1", "BGLC", "BGAL", "ANE5"], start=2):
        atoms.append({
            "record": "HETATM", "name": "C1", "resname": resname, "chain": "B", "resid": resid,
            "icode": "", "x": 52.0, "y": 52.0, "z": 60.0, "element": "C",
        })
    for resid, resname in enumerate(["CER1", "BGLC", "BGAL", "ANE5"], start=6):
        atoms.append({
            "record": "HETATM", "name": "C1", "resname": resname, "chain": "C", "resid": resid,
            "icode": "", "x": 95.0, "y": 95.0, "z": 60.0, "element": "C",
        })
    aa_stage4._write_pdb(source, [], atoms, np.array([100.0, 100.0, 120.0]))
    out = tmp_path / "sliced.pdb"
    report = aa_stage4.slice_membrane(source, out, 1.0)
    assert report["kept_molecules"]["GLPA"] == 1
    assert report["removed_molecules"]["GLPA"] == 1
    assert aa_stage4.molecule_counts_from_pdb(out) == {"GLPA": 1, "PROTEIN": 1}


def test_write_index_ini_groups_final_gro(tmp_path: Path) -> None:
    gro = tmp_path / "ionized.gro"
    gro.write_text(
        "\n".join([
            "system",
            "5",
            "    1ALA     CA    1   0.100   0.100   0.100",
            "    2POPE     P    2   0.200   0.200   0.200",
            "    3TIP3   OH2    3   0.300   0.300   0.300",
            "    4SOD    SOD    4   0.400   0.400   0.400",
            "    5CLA    CLA    5   0.500   0.500   0.500",
            "   1.00000   1.00000   1.00000",
        ]) + "\n",
        encoding="utf-8",
    )
    ndx = tmp_path / "index_ini.ndx"
    report = aa_stage4.write_index_ini(gro, ndx)
    text = ndx.read_text(encoding="utf-8")
    assert report == {"System": 5, "Protein": 1, "Membrane": 1, "Water": 1, "Ions": 2, "Water_and_ions": 3}
    assert "[ Protein ]" in text
    assert "[ Membrane ]" in text
    assert "[ Water_and_ions ]" in text


def test_normalize_solvated_gro_removes_water_inside_membrane_z(tmp_path: Path) -> None:
    gro = tmp_path / "raw.gro"
    gro.write_text(
        "\n".join([
            "solvated",
            "7",
            "    1POPE     P    1   0.100   0.100   0.500",
            "    2SOL     OW    2   0.200   0.200   0.500",
            "    2SOL    HW1    3   0.210   0.200   0.500",
            "    2SOL    HW2    4   0.200   0.210   0.500",
            "    3SOL     OW    5   0.300   0.300   0.800",
            "    3SOL    HW1    6   0.310   0.300   0.800",
            "    3SOL    HW2    7   0.300   0.310   0.800",
            "   1.00000   1.00000   1.00000",
        ]) + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "filtered.gro"
    report = aa_stage4.normalize_solvated_gro(gro, out, (0.4, 0.6))
    lines = out.read_text(encoding="utf-8").splitlines()
    assert report["waters"] == 1
    assert report["removed_waters_in_membrane_z"] == 1
    assert int(lines[1]) == 4
    assert any(raw[5:10].strip() == "TIP3" and raw[10:15].strip() == "OH2" for raw in lines[2:-1])
    water_oxygen_z = [
        float(raw[36:44])
        for raw in lines[2:-1]
        if raw[5:10].strip() == "TIP3" and raw[10:15].strip() == "OH2"
    ]
    assert water_oxygen_z == [0.8]


def test_phosphate_headgroup_core_z_bounds_use_leaflet_gap(tmp_path: Path) -> None:
    path = tmp_path / "phosphates.pdb"
    atoms = []
    for resid, z in enumerate([10.0, 10.4, 10.8, 25.0, 25.3, 25.7], start=1):
        atoms.append({
            "record": "HETATM", "name": "P", "resname": "POPC", "chain": "M", "resid": resid,
            "icode": "", "x": float(resid), "y": 1.0, "z": z, "element": "P",
        })
    aa_stage4._write_pdb(path, [], atoms, np.array([50.0, 50.0, 250.0]))
    assert aa_stage4.phosphate_headgroup_core_z_bounds_nm(path) == (1.08, 2.5)


def test_headgroup_core_z_bounds_do_not_require_phosphorus_atoms(tmp_path: Path) -> None:
    path = tmp_path / "chol.pdb"
    atoms = []
    for resid, z in enumerate([10.0, 10.4, 10.8, 25.0, 25.3, 25.7], start=1):
        atoms.append({
            "record": "HETATM", "name": "ROH", "resname": "CHOL", "chain": "M", "resid": resid,
            "icode": "", "x": float(resid), "y": 1.0, "z": z, "element": "O",
        })
    aa_stage4._write_pdb(path, [], atoms, np.array([50.0, 50.0, 250.0]))
    assert aa_stage4.phosphate_headgroup_core_z_bounds_nm(path) == (1.08, 2.5)


def test_charmm_histidine_is_protein_like_and_explicit_caps_are_flagged(tmp_path: Path) -> None:
    path = tmp_path / "caps.pdb"
    aa_stage4._write_pdb(
        path,
        [],
        [
            {
                "record": "ATOM", "name": "CH3", "resname": "ACE", "chain": "A", "resid": 0,
                "icode": "", "x": 1.0, "y": 1.0, "z": 1.0, "element": "C",
            },
            {
                "record": "ATOM", "name": "CA", "resname": "HSD", "chain": "A", "resid": 1,
                "icode": "", "x": 2.0, "y": 2.0, "z": 2.0, "element": "C",
            },
            {
                "record": "ATOM", "name": "CH3", "resname": "CT3", "chain": "A", "resid": 2,
                "icode": "", "x": 3.0, "y": 3.0, "z": 3.0, "element": "C",
            },
        ],
        np.array([50.0, 50.0, 50.0]),
    )
    assert aa_stage4.molecule_counts_from_pdb(path) == {"PROTEIN": 1}
    assert aa_stage4.explicit_charmm_terminal_cap_residues(path) == ["A:ACE0", "A:CT32"]


def test_charmm_r2b_aliases_include_ct3_terminal_cap(tmp_path: Path) -> None:
    ff = tmp_path / "charmm36.ff"
    ff.mkdir()
    r2b = ff / "merged.r2b"
    r2b.write_text("; aliases\nHISD\tHSD\n", encoding="utf-8")
    aa_stage4._ensure_charmm_r2b_aliases(ff)
    text = r2b.read_text(encoding="utf-8")
    assert "ACE\tACE" in text
    assert "CT3\tCT3" in text
    assert "CT3\tNME" not in text
    assert "HSD\tHSD" in text


def test_filter_solvent_ions_by_z_preserves_charge_balance(tmp_path: Path) -> None:
    gro = tmp_path / "raw.gro"
    gro.write_text(
        "\n".join([
            "system",
            "7",
            "    1TIP3   OH2    1   0.100   0.100   0.500",
            "    1TIP3    H1    2   0.110   0.100   0.500",
            "    1TIP3    H2    3   0.100   0.110   0.500",
            "    2SOD    SOD    4   0.200   0.200   0.500",
            "    3CLA    CLA    5   0.300   0.300   0.300",
            "    4CLA    CLA    6   0.400   0.400   0.800",
            "    5SOD    SOD    7   0.500   0.500   0.800",
            "   1.00000   1.00000   1.00000",
        ]) + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "filtered.gro"
    report = aa_stage4.filter_solvent_ions_by_z(gro, out, (0.4, 0.6))
    text = out.read_text(encoding="utf-8")
    assert report["removed_waters_in_membrane_z"] == 1
    assert report["removed_sod_in_membrane_z"] == 1
    assert report["removed_cla_for_charge_balance"] == 1
    assert int(text.splitlines()[1]) == 2
    assert aa_stage4.solvent_ion_counts_from_gro(out) == {"CLA": 1, "SOD": 1}
