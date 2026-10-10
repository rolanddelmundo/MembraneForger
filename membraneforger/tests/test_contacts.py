"""The force-field contact validator and bond-crossing detector, and the gates built on them."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from .common import failure, mf

# A small CHARMM-like force field: real CHARMM36 sigma/epsilon of the types used, plus one NBFIX row.
FORCEFIELD = """[ defaults ]
1 2 yes 1.0 1.0

[ atomtypes ]
 CTL2  6 12.011000 0.000 A 0.358141284692 0.23430
 HAL2  1  1.008000 0.000 A 0.238760856462 0.11715
   OT  8 15.999400 0.000 A 0.315057422683 0.63639
   HT  1  1.008000 0.000 A 0.0400013524445 0.19246
  SOD 11 22.989770 0.000 A 0.242992625373 0.19623
  CLA 17 35.450000 0.000 A 0.404468018036 0.62760

[ nonbond_params ]
  SOD  OT  1  0.288700000000  0.313883680000
"""
# LIP: a four-carbon chain with one hydrogen, nrexcl 3 (C1..C4 is a 1-4 pair and excluded, as in CHARMM).
ITPS = {
    "LIP": """[ moleculetype ]
LIP 3
[ atoms ]
1 CTL2 1 LIP C1 1 0.0 12.011
2 HAL2 1 LIP H1 2 0.0 1.008
3 CTL2 1 LIP C2 3 0.0 12.011
4 CTL2 1 LIP C3 4 0.0 12.011
5 CTL2 1 LIP C4 5 0.0 12.011
[ bonds ]
1 2
1 3
3 4
4 5
""",
    "PROT": """[ moleculetype ]
PROT 3
[ atoms ]
1 CTL2 1 ALA CA 1 0.0 12.011
2 HAL2 1 ALA HA 2 0.0 1.008
[ bonds ]
1 2
""",
    "TIP3": """[ moleculetype ]
TIP3 2
[ atoms ]
1 OT 1 TIP3 OH2 1 -0.834 15.9994
2 HT 1 TIP3 H1 1 0.417 1.008
3 HT 1 TIP3 H2 1 0.417 1.008
[ settles ]
1 1 0.09572 0.15139
""",
    "SOD": "[ moleculetype ]\nSOD 1\n[ atoms ]\n1 SOD 1 SOD SOD 1 1.0 22.98977\n",
    "CLA": "[ moleculetype ]\nCLA 1\n[ atoms ]\n1 CLA 1 CLA CLA 1 -1.0 35.45\n",
}
NAMES = {"LIP": ["C1", "H1", "C2", "C3", "C4"], "PROT": ["CA", "HA"], "TIP3": ["OH2", "H1", "H2"], "SOD": ["SOD"], "CLA": ["CLA"]}
RESNAME = {"PROT": "ALA"}
BOX_A = (40.0, 40.0, 40.0)


def chain(x, y, z, axis=(1.0, 0.0, 0.0), step=1.53):
    """LIP coordinates (A): C1 at (x, y, z), carbons every `step` along `axis`, H1 1.1 A off C1 in z."""
    a = np.array(axis, float)
    c = [np.array([x, y, z]) + k * step * a for k in range(4)]
    return [c[0], c[0] + np.array([0.0, 0.0, 1.1]), c[1], c[2], c[3]]


def water(x, y, z, towards=(1.0, 0.0, 0.0)):
    """TIP3 coordinates (A): O at (x, y, z), one H 0.96 A along `towards`, the other off in y."""
    o, t = np.array([x, y, z]), np.array(towards, float)
    return [o, o + 0.9572 * t / np.linalg.norm(t), o + np.array([0.0, 0.9572, 0.0])]


def build(out: Path, molecules: list[tuple], structure: str = "boxed.gro", box_a=BOX_A) -> Path:
    """Write topol.top, the itps and one .gro with the given (moltype, resid, coordinates in A) molecules, in order."""
    (out / "toppar").mkdir(exist_ok=True)
    (out / "toppar" / "forcefield.itp").write_text(FORCEFIELD)
    for name, text in ITPS.items():
        (out / "toppar" / f"{name}.itp").write_text(text)
    counts, atoms = [], []
    for mol, resid, xyz in molecules:
        if not counts or counts[-1][0] != mol:
            counts.append([mol, 0])
        counts[-1][1] += 1
        for name, p in zip(NAMES[mol], xyz):
            atoms.append({"resid": resid, "resname": RESNAME.get(mol, mol), "atom": name, "x": p[0] / 10.0, "y": p[1] / 10.0, "z": p[2] / 10.0})
    (out / "topol.top").write_text('#include "toppar/forcefield.itp"\n' + "".join(f'#include "toppar/{n}.itp"\n' for n in ITPS)
                                   + "[ system ]\ntest\n[ molecules ]\n" + "".join(f"{m} {n}\n" for m, n in counts))
    mf.write_gro(atoms, [v / 10.0 for v in box_a], out / structure, "test")
    return out / structure


def clean_system() -> list[tuple]:
    """Protein, two packed lipids, a folded lipid, hydrogen-bonded waters and a salt pair: every contact legitimate."""
    folded = chain(30.0, 30.0, 30.0)
    folded[4] = folded[0] + np.array([2.0, 0.0, 0.0])  # C4 2.0 A from C1: a 1-4 pair, excluded by nrexcl 3
    return [("PROT", 1, [np.array([5.0, 5.0, 5.0]), np.array([5.0, 5.0, 6.1])]),
            ("LIP", 2, chain(10.0, 10.0, 10.0)), ("LIP", 3, chain(10.0, 14.5, 10.0)),            # 4.5 A apart
            ("LIP", 4, chain(10.0, 7.4, 10.0)),                                                    # 2.6 A beside lipid 2: packed
            ("LIP", 5, folded),
            ("TIP3", 6, water(20.0, 20.0, 20.0)), ("TIP3", 7, water(22.8, 20.0, 20.0, towards=(-1, 0, 0))),  # O..H 1.84 A, O..O 2.8 A
            ("TIP3", 8, water(30.0, 30.0, 31.1 + 1.3 + 0.9572, towards=(0, 0, -1))),                 # water H 1.3 A from lipid 5's H1
            ("SOD", 9, [np.array([30.0, 10.0, 10.0])]), ("TIP3", 10, water(32.3, 10.0, 10.0)),     # Na+..O 2.3 A
            ("CLA", 11, [np.array([30.0, 10.0, 13.0])])]                                           # Cl-..Na+ 3.0 A


class ContactValidator(unittest.TestCase):
    """Pairs are judged by their own Lennard-Jones Rmin, exclusions by the topology, distances under PBC."""

    def setUp(self):
        self.out = Path(tempfile.mkdtemp())

    def report(self, molecules, structure="boxed.gro", solute_atoms=2, **settings):
        build(self.out, molecules, structure)
        return mf.check_contacts(self.out, structure, solute_atoms, mf.Settings(**settings))

    def test_a_clean_system_has_no_findings(self):
        report = self.report(clean_system())
        self.assertEqual((report["contacts"], report["crossings"]), ([], []))
        self.assertEqual(report["atoms"], 2 + 4 * 5 + 4 * 3 + 2)
        self.assertEqual(report["contacts_by_class"], {"solute": 0, "membrane": 0, "solvent": 0})

    def test_the_limit_is_a_fraction_of_the_pair_rmin(self):
        tables = mf.system_tables(build(self.out, clean_system()).parent / "topol.top")
        c, h = tables["type_names"].index("CTL2"), tables["type_names"].index("HAL2")
        self.assertAlmostEqual(tables["rmin"][c, c], 2 ** (1 / 6) * 0.358141284692, places=9)          # 4.02 A
        self.assertAlmostEqual(tables["rmin"][c, h], 2 ** (1 / 6) * 0.5 * (0.358141284692 + 0.238760856462), places=9)
        self.assertAlmostEqual(tables["rmin"][tables["type_names"].index("SOD"), tables["type_names"].index("OT")],
                               2 ** (1 / 6) * 0.2887, places=9)                                         # the NBFIX row, not the rule
        self.assertEqual(len(tables["bonds"]), 4 * 4 + 1 + 4 * 2)                                       # settles count as two bonds

    def test_lipid_on_lipid_and_lipid_on_solute(self):
        system = clean_system()
        system[2] = ("LIP", 3, chain(10.0, 12.3, 10.0))             # C..C 2.3 A < 0.6 * 4.02 = 2.41 A
        system[4] = ("LIP", 5, chain(5.0, 7.2, 5.0))                 # C1 2.2 A from the protein CA
        report = self.report(system)
        pairs = {tuple(sorted(a.split("(")[0] for a in c["atoms"])) for c in report["contacts"]}
        self.assertIn(("LIP2:C1", "LIP3:C1"), pairs)
        self.assertIn(("ALA1:CA", "LIP5:C1"), pairs)
        self.assertEqual({c["class"] for c in report["contacts"]}, {"membrane"})
        worst = report["contacts"][0]
        self.assertEqual(worst["class"], "membrane")
        self.assertAlmostEqual(worst["limit_A"], 10 * 0.6 * 2 ** (1 / 6) * 0.358141284692, places=2)
        self.assertEqual(report["crossings"], [])
        self.assertIn("inside their LJ contact limit", mf.describe_contacts(report))

    def test_exclusions_follow_the_topology_not_the_residue_number(self):
        system = clean_system()
        system[1] = ("LIP", 7, chain(10.0, 10.0, 10.0))              # two different molecules share resid 7...
        system[2] = ("LIP", 7, chain(10.0, 12.0, 10.0))              # ...and are 2.0 A apart: a real pair
        system[3] = ("LIP", 123456, chain(10.0, 17.0, 10.0))         # a residue number beyond the .gro field
        report = self.report(system)
        self.assertEqual(len(report["contacts"]), 4)                 # C1..C1, C2..C2, C3..C3, C4..C4 of the two resid-7 lipids
        self.assertEqual({tuple(c["molecules"]) for c in report["contacts"]}, {(2, 3)})
        # the folded lipid's C1..C4 at 2.0 A is a 1-4 pair of one molecule: excluded, not reported
        self.assertFalse(any("LIP5" in a for c in report["contacts"] for a in c["atoms"]))

    def test_water_and_ion_overlaps(self):
        system = clean_system()
        system[5] = ("TIP3", 6, water(10.0, 10.0, 7.8))              # O 2.2 A from lipid 2's C1 (limit 2.27 A)
        system[6] = ("TIP3", 7, water(10.0, 10.0, 11.1 + 1.5))       # O 1.5 A from lipid 2's H1 (limit 1.87 A)
        system[8] = ("SOD", 9, [np.array([30.0, 10.0, 10.0])])
        system[9] = ("TIP3", 10, water(31.9, 10.0, 10.0))            # Na+..O 1.90 A, inside the NBFIX limit 1.94 A: solvent on solvent, not judged
        system[10] = ("CLA", 11, [np.array([14.59, 7.4, 12.9])])     # Cl- 2.9 A above lipid 4's C4 (limit 0.6 * 4.28 = 2.57 A): fine
        system.append(("TIP3", 12, water(14.59, 10.0, 12.4)))        # O 2.4 A above lipid 2's C4 (limit 2.27 A): fine; 2.65 A from the
        #                                                              Cl- (0.66 of their 4.04 A Rmin, a genion neighbour): not judged
        system.append(("CLA", 13, [np.array([30.0, 10.0, 13.0])]))   # 3.0 A from the Na+: a salt pair, not judged
        report = self.report(system)
        found = {tuple(sorted(a.split("(")[0] for a in c["atoms"])): c for c in report["contacts"]}
        self.assertEqual(set(found), {("LIP2:C1", "TIP36:OH2"), ("LIP2:H1", "TIP37:OH2")})
        self.assertEqual({c["class"] for c in found.values()}, {"solvent"})
        system[10] = ("CLA", 11, [np.array([14.59, 7.4, 12.2])])     # Cl- 2.2 A above lipid 4's C4: inside its limit
        found = {tuple(sorted(a.split("(")[0] for a in c["atoms"])) for c in self.report(system)["contacts"]}
        self.assertIn(("CLA11:CLA", "LIP4:C4"), found)
        self.assertEqual(len(found), 3)

    def test_overlapping_marks_the_water_atoms_that_solvation_must_drop(self):
        lj = mf.lj_table([Path(build(self.out, clean_system()).parent / "toppar" / "forcefield.itp")])
        lipid = np.array(chain(10.0, 10.0, 10.0)) / 10.0
        waters = [water(10.0, 10.0, 7.8), water(10.0, 10.0, 11.1 + 1.5), water(10.0, 10.0, 11.1 + 1.3 + 0.9572, towards=(0, 0, -1)),
                  water(20.0, 20.0, 20.0)]
        xyz = np.vstack(waters) / 10.0
        bad, worst = mf.overlapping(xyz, ["OT", "HT", "HT"] * 4, lipid, ["CTL2", "HAL2", "CTL2", "CTL2", "CTL2"], np.array(BOX_A) / 10.0, lj, 0.6)
        self.assertEqual(bad.reshape(-1, 3).any(axis=1).tolist(), [True, True, False, False])   # a water H 1.3 A from a lipid H has no LJ to speak of
        self.assertEqual(worst[0] // 3, 1)                                                        # the deepest pair is water 2's O on H1
        self.assertEqual(worst[1], 1)
        self.assertAlmostEqual(worst[2], 0.15, places=3)

    def test_periodic_images(self):
        system = clean_system()
        system[1] = ("LIP", 2, chain(0.2, 10.0, 10.0, axis=(0, 1, 0)))
        system[2] = ("LIP", 3, chain(BOX_A[0] - 1.8, 10.0, 10.0, axis=(0, 1, 0)))   # 2.0 A apart across x
        system[3] = ("LIP", 4, chain(BOX_A[0] - 3.5, 25.0, 10.0, axis=(0, 1, 0)))   # 3.5 A from lipid 5 across x: fine
        system[4] = ("LIP", 5, chain(0.0, 25.0, 10.0, axis=(0, 1, 0)))
        report = self.report(system)
        self.assertEqual({tuple(c["molecules"]) for c in report["contacts"]}, {(2, 3)})
        self.assertAlmostEqual(report["contacts"][0]["distance_A"], 2.0, places=2)

    def test_crossing_chains(self):
        system = clean_system()
        # lipid 3 runs along y and passes 0.8 A above the middle of lipid 2's C2-C3 bond
        x_mid = 10.0 + 1.53 * 1.5
        system[2] = ("LIP", 3, chain(x_mid, 10.0 - 1.53 * 1.5, 10.8, axis=(0, 1, 0)))
        system[3] = ("LIP", 4, chain(10.0, 20.0, 10.0))
        report = self.report(system)
        self.assertEqual(len(report["crossings"]), 1)
        crossing = report["crossings"][0]
        self.assertAlmostEqual(crossing["distance_A"], 0.8, places=2)
        self.assertEqual(sorted(crossing["molecules"]), [2, 3])
        self.assertEqual(crossing["class"], "membrane")
        self.assertIn("bond crossing", mf.describe_contacts(report))
        # the same two chains 2.6 A apart are packed, not crossed: no crossing, no contact
        system[2] = ("LIP", 3, chain(x_mid, 10.0 - 1.53 * 1.5, 12.6, axis=(0, 1, 0)))
        report = self.report(system)
        self.assertEqual((report["contacts"], report["crossings"]), ([], []))

    def test_crossing_across_the_periodic_boundary(self):
        system = clean_system()
        system[1] = ("LIP", 2, chain(BOX_A[0] - 2.3, 10.0, 10.0))                        # C2-C3 straddles x = 0
        system[2] = ("LIP", 3, chain(0.0, 10.0 - 1.53 * 1.5, 10.8, axis=(0, 1, 0)))     # crosses it 0.8 A above
        system[3] = ("LIP", 4, chain(10.0, 20.0, 10.0))
        report = self.report(system)
        self.assertEqual([sorted(c["molecules"]) for c in report["crossings"]], [[2, 3]])

    def test_a_bad_fraction_is_refused(self):
        self.assertIn("contact_rmin_fraction", failure(mf.check_fraction, 0.9))
        self.assertIsNone(failure(mf.check_fraction, 0.6))


class ContactGates(unittest.TestCase):
    """What the pipeline does with the findings: retry a backmap, stop before EM, refuse an EM result, name the culprit."""

    def setUp(self):
        self.out = Path(tempfile.mkdtemp())

    def crossed(self) -> list[tuple]:
        system = clean_system()
        x_mid = 10.0 + 1.53 * 1.5
        system[2] = ("LIP", 3, chain(x_mid, 10.0 - 1.53 * 1.5, 10.8, axis=(0, 1, 0)))
        system[3] = ("LIP", 4, chain(10.0, 20.0, 10.0))
        return system

    def staged(self, report):
        return {"rings": {"pierced": []}, "solute_atoms": 2, "contact": {"distance_nm": 0.3, "atoms": ["a", "b"]}, "overlaps": report}

    def test_backmap_verdict_retries_on_a_membrane_overlap_and_tolerates_the_input_itself(self):
        build(self.out, self.crossed())
        report = mf.check_contacts(self.out, "boxed.gro", 2)
        verdict = mf.backmap_verdict(self.staged(report), {}, mf.Settings(), last=False)
        self.assertFalse(verdict["accept"])
        self.assertEqual(verdict["stage"], "contacts")
        self.assertIn("bond crossing", verdict["problem"])
        self.assertEqual(verdict["counts"]["membrane_bond_crossings"], 1)
        self.assertGreater(verdict["counts"]["membrane_contact_violations"], 0)
        self.assertIn("after 5 backmapping attempt(s)", failure(mf.pipeline.refuse_backmap, verdict["problem"], 5))
        # a contact inside the all-atom input is not the membrane's: another seed cannot change it
        solute = {"contacts": [{"class": "solute", "atoms": ["ALA1:CA(C, atom 1)", "ALA2:CB(C, atom 3)"], "distance_A": 2.0,
                                "limit_A": 2.4, "rmin_A": 4.0, "molecules": [1, 1]}], "crossings": []}
        self.assertTrue(mf.backmap_verdict(self.staged(solute), {}, mf.Settings(), last=False)["accept"])
        self.assertTrue(mf.backmap_verdict({**self.staged(solute), "overlaps": None}, {}, mf.Settings(), last=False)["accept"])

    def test_refuse_overlaps_before_em(self):
        build(self.out, self.crossed(), "solv_ions.gro")
        report = mf.check_contacts(self.out, "solv_ions.gro", 2)
        message = failure(mf.refuse_overlaps, report, ("membrane", "solvent"), "before EM")
        self.assertIn("before EM", message)
        self.assertIn("minimization cannot be trusted", message)
        clean = mf.check_contacts(self.out, "solv_ions.gro", 2) if build(self.out, clean_system(), "solv_ions.gro") else None
        self.assertIsNone(failure(mf.refuse_overlaps, clean, ("membrane", "solvent"), "before EM"))
        # a solute-internal contact passes this gate; a crossing of any class never does
        solute = dict(clean, contacts=[{"class": "solute", "atoms": ["a", "b"], "distance_A": 2.0, "limit_A": 2.4, "rmin_A": 4.0,
                                        "molecules": [1, 1]}])
        self.assertIsNone(failure(mf.refuse_overlaps, solute, ("membrane", "solvent"), "before EM"))
        crossed = dict(clean, crossings=[{"class": "solute", "bonds": [["a", "b"], ["c", "d"]], "distance_A": 0.5, "molecules": [1, 1]}])
        self.assertIn("bond crossing", failure(mf.refuse_overlaps, crossed, ("membrane", "solvent"), "before EM"))

    def validate(self, system, fmax, atom=12):
        build(self.out, system, "solv_ions.gro")
        build(self.out, system, "em.unverified.gro")
        build(self.out, system, "boxed.gro")
        (self.out / "em.log").write_text(f"Potential Energy  = -1.0e+05\nMaximum force     =  {fmax:.5e} on atom {atom}\n"
                                         "Steepest Descents converged to Fmax < 500 in 10 steps\n")
        names = [a["atom"] for a in mf.read_gro(self.out / "em.unverified.gro")[0]]
        index = {"System": len(names), "Protein_LIG": 2, "MEMB": 20, "SOL_ION": len(names) - 22}
        return failure(mf.validate_em, {"out": self.out, "inputs": {}}, {"names": names, "charge": 0.0}, index)

    def test_em_may_not_leave_an_overlap_or_a_crossing(self):
        message = self.validate(self.crossed(), 100.0)
        self.assertIn("EM reached Fmax 100.0 but left", message)
        self.assertIn("bond crossing", message)

    def test_a_clean_minimized_system_passes(self):
        self.assertIsNone(self.validate(clean_system(), 100.0))

    def test_an_fmax_failure_names_the_culprit_with_elements_and_its_overlaps(self):
        message = self.validate(self.crossed(), 20878.1, atom=3 + 5 + 2)  # lipid 3's C2
        self.assertIn("EM final Fmax 20878.1 >= 500 kJ/mol/nm", message)
        self.assertIn("largest force on atom 10 LIP3:C2 (C)", message)
        self.assertRegex(message, r"LIP2:C\d\(C\) \d\.\d\d A")
        self.assertIn("force-field overlaps of LIP3 after EM:", message)
        self.assertIn("bond crossing", message)
        self.assertIn("force-field overlaps of LIP3 before EM:", message)


if __name__ == "__main__":
    unittest.main()
