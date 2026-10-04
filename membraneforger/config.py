#//=============================================================
#// MembraneForger - (c) 2026 Roland Del Mundo
#// CHARMM36/GROMACS constants and tunable acceptance criteria.
#//=============================================================
"""Constants shared by every stage, and the Settings that carry the tunable acceptance criteria."""
from dataclasses import dataclass

__all__ = ['ORIENTATION', 'ORIENT_CHAINS', 'NTERM_SIDE', 'PDB_ID', 'PPM_MEMBRANE', 'BOX',
           'AMINO', 'SOLVENT', 'N_CAP', 'C_CAP', 'RENAME_RESIDUE', 'RENAME_ATOM', 'RENAME_MOLECULE', 'TERMINUS_MENU',
           'DEFAULT_LIGANDS', 'NSTEPS', 'DISULFIDE_MAX_A', 'DISULFIDE_OK_A', 'SLAB_Z_PAD_NM', 'MIN_Z_PAD_TOTAL_NM',
           'SALT_M', 'ION_RMIN_NM', 'GENION_ATTEMPTS', 'WATER_CLASH_NM', 'WATER_PROTECT_NM', 'LIPID_SCAN_A',
           'LIPID_DELETE_A', 'GENERATED', 'INDEX_GROUPS', 'FMAX_TARGET', 'H_BOND_RANGE', 'GLPA_MAX_BOND_A',
           'LYS_BACKBONE', 'ONE_LETTER', 'GM3_XML_TO_GLPA', 'LIPIDATED', 'CYSG_HDB', 'EM_MDP', 'Settings']

# ============================================================
# MEMBRANE ORIENTATION
# ============================================================
# Select the chain that actually spans or associates with the membrane. MembraneForger determines its membrane
# orientation (from the exact OPM entry of the structure when a PDB ID is known and matches, otherwise with a local
# PPM 3.0 run on the anchor's own coordinates) and moves the complete complex as ONE rigid object. Command-line
# options (--orientation, --orient-chain(s), --nterm-side, --pdb-id, --ppm-membrane) override these defaults.
ORIENTATION = "auto"       # auto | ppm | opm | none      (none: use the input coordinates as given)
ORIENT_CHAINS = ""         # e.g. "R", "A", or "A,B"      (empty: the only protein chain, else fail and ask)
NTERM_SIDE = "auto"        # auto | in | out              (side of the N terminus of the first anchor chain)
PDB_ID = ""                # optional exact PDB ID         (else the HEADER record of the input, else none)
PPM_MEMBRANE = ""          # advanced: PPM 3.0 membrane code, "" = undefined flat bilayer (e.g. "PMm", "GnI")

# ============================================================
# BOX
# ============================================================
# auto: x and y come from the membrane cell; z is sized around the bilayer midplane so that the complex and the
# membrane are covered by SLAB_Z_PAD_NM of water above and below. A user box is opt-in: "x,y,z" in nm, where x and
# y must equal the membrane cell (the membrane is periodic in it) and z must be at least the automatic minimum.
BOX = "auto"               # auto | "x,y,z" (nm)

AMINO = {"ALA", "ARG", "ASN", "ASP", "CYS", "CYSG", "CYSP", "GLN", "GLU", "GLY", "HIS", "HSD", "HSE", "HSP", "ILE",
         "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "AIB", "LEM", "KTZ", "KRT", "KSM"}


SOLVENT = {"TIP3", "SOD", "CLA", "POT"}


N_CAP = ["CAY", "HY1", "HY2", "HY3", "CY", "OY"]


C_CAP = ["NT", "HNT", "CAT", "HT1", "HT2", "HT3"]


RENAME_RESIDUE = {"HIS": "HSD", "HID": "HSD", "HIE": "HSE", "HIP": "HSP"}


RENAME_ATOM = {"ILE": {"CD1": "CD"}}


RENAME_MOLECULE = {"SAPI": "SAPI25", "GM3": "GLPA"}  # structure resname -> toppar ITP


TERMINUS_MENU = {"NTER": "NH3+", "ACE": "None", "CTER": "COO-", "CT2": "CT2", "CT3": "None"}


DEFAULT_LIGANDS = frozenset({"GDP", "MG", "V6G"})  # non-protein residues treated as ligands when no all-atom input says so
NSTEPS = 50000  # default maximum EM steps
DISULFIDE_MAX_A, DISULFIDE_OK_A = 3.0, (1.8, 2.2)  # SG-SG below 3.0 A is a disulfide; outside 1.8-2.2 A it is rebuilt
SLAB_Z_PAD_NM, MIN_Z_PAD_TOTAL_NM = 1.5, 2.0  # water above/below the solute, and the least total z padding
SALT_M, ION_RMIN_NM, GENION_ATTEMPTS = 0.15, 0.60, 3
WATER_CLASH_NM, WATER_PROTECT_NM = 0.18, 0.45
LIPID_SCAN_A, LIPID_DELETE_A = 1.0, 0.10  # EM survived 0.15 A (6WHC_MORF) but not 0.05 A (7RA3_MRTR)
# Every file a build writes into the output directory; removed at the start so nothing stale survives a rerun.
GENERATED = ("oriented.pdb", "orientation_report.json", "membrane.pdb", "aa_cg_mapping.tsv", "prot-memb.pdb",
             "topol.top", "topol.pre_genion.top", "boxed.gro",
             "solv_raw.gro", "solv.gro", "solv_ions.gro", "index_ini.ndx", "genion.ndx", "ions.mdp", "ions.tpr",
             "ions_mdout.mdp", "em.mdp", "mdout.mdp", "em.tpr", "em.log", "em.edr", "em.trr", "em.gro",
             "em.unverified.gro", "toppar", "audit.json", "run_manifest.json", "ring_piercing.json", "work")


INDEX_GROUPS = ("System", "Protein_LIG", "MEMB", "SOL_ION")


FMAX_TARGET = 500.0


H_BOND_RANGE = (0.60, 1.35)


GLPA_MAX_BOND_A = 2.0


LYS_BACKBONE = {"N": ("NH1", -.47), "HN": ("H", .31), "HA": ("HB1", .09), "C": ("C", .51), "O": ("O", -.51)}


ONE_LETTER = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
              "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
              "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "HSD": "H", "HSE": "H", "HSP": "H", "CYSG": "C",
              "CYSP": "C", "KTZ": "K", "KRT": "K", "KSM": "K"}  # anything else aligns as X


# Step-4 GM3 (GM3.xml names) -> GLPA.itp names, in GLPA.itp order. Verbatim from finalize-templates
# (GM3_TO_GLPA_ITP_PLAN); several H names swap, so never match GM3 by name.
GM3_XML_TO_GLPA = [
    ("C1S", "C1S"), ("H1T", "H1S"), ("H1S", "H1T"), ("NF", "NF"), ("HNF", "HNF"), ("C2S", "C2S"),
    ("H2S", "H2S"), ("C3S", "C3S"), ("H3S", "H3S"), ("O34", "O3"), ("HO34", "HO3"), ("C4S", "C4S"),
    ("H4S", "H4S"), ("C5S", "C5S"), ("H5S", "H5S"), ("C6S", "C6S"), ("H6S", "H6S"), ("H6T", "H6T"),
    ("C7S", "C7S"), ("H7T", "H7S"), ("H7S", "H7T"), ("C8S", "C8S"), ("H8T", "H8S"), ("H8S", "H8T"),
    ("C9S", "C9S"), ("H9T", "H9S"), ("H9S", "H9T"), ("C10S", "C10S"), ("H10T", "H10S"), ("H10S", "H10T"),
    ("C11S", "C11S"), ("H11T", "H11S"), ("H11S", "H11T"), ("C12S", "C12S"), ("H12T", "H12S"), ("H12S", "H12T"),
    ("C13S", "C13S"), ("H13S", "H13S"), ("H13T", "H13T"), ("C14S", "C14S"), ("H14S", "H14S"), ("H14T", "H14T"),
    ("C15S", "C15S"), ("H15S", "H15S"), ("H15T", "H15T"), ("C16S", "C16S"), ("H16T", "H16S"), ("H16S", "H16T"),
    ("C17S", "C17S"), ("H17T", "H17S"), ("H17S", "H17T"), ("C18S", "C18S"), ("H18S", "H18S"), ("H18U", "H18T"),
    ("H18T", "H18U"), ("C1F", "C1F"), ("OF", "OF"), ("C2F", "C2F"), ("H2G", "H2F"), ("H2F", "H2G"),
    ("C3F", "C3F"), ("H3F", "H3F"), ("H3G", "H3G"), ("C4F", "C4F"), ("H4G", "H4F"), ("H4F", "H4G"),
    ("C5F", "C5F"), ("H5F", "H5F"), ("H5G", "H5G"), ("C6F", "C6F"), ("H6G", "H6F"), ("H6F", "H6G"),
    ("C7F", "C7F"), ("H7F", "H7F"), ("H7G", "H7G"), ("C8F", "C8F"), ("H8G", "H8F"), ("H8F", "H8G"),
    ("C9F", "C9F"), ("H9G", "H9F"), ("H9F", "H9G"), ("C10F", "C10F"), ("H10G", "H10F"), ("H10F", "H10G"),
    ("C11F", "C11F"), ("H11F", "H11F"), ("H11G", "H11G"), ("C12F", "C12F"), ("H12F", "H12F"), ("H12G", "H12G"),
    ("C13F", "C13F"), ("H13F", "H13F"), ("H13G", "H13G"), ("C14F", "C14F"), ("H14F", "H14F"), ("H14G", "H14G"),
    ("C15F", "C15F"), ("H15G", "H15F"), ("H15F", "H15G"), ("C16F", "C16F"), ("H16G", "H16F"), ("H16H", "H16G"),
    ("H16F", "H16H"), ("C1", "C1"), ("H1", "H1"), ("O1", "O1"), ("C5", "C5"), ("H5", "H5"), ("O5", "O5"),
    ("C2", "C2"), ("H2", "H2"), ("O2", "O2"), ("HO2", "HO2"), ("C3", "C3"), ("H3", "H3"), ("O3", "O3"),
    ("HO3", "HO3"), ("C4", "C4"), ("H4", "H4"), ("O4", "O4"), ("C6", "C6"), ("H62", "H61"), ("H61", "H62"),
    ("O6", "O6"), ("HO6", "HO6"), ("C7", "C1"), ("H7", "H1"), ("C12", "C5"), ("H12", "H5"), ("O12", "O5"),
    ("C8", "C2"), ("H8", "H2"), ("O9", "O2"), ("HO9", "HO2"), ("C10", "C3"), ("H10", "H3"), ("O10", "O3"),
    ("C11", "C4"), ("H11", "H4"), ("O11", "O4"), ("HO11", "HO4"), ("C13", "C6"), ("H132", "H61"),
    ("H131", "H62"), ("O13", "O6"), ("HO13", "HO6"), ("C14", "C1"), ("O15", "O11"), ("O14", "O12"),
    ("C15", "C2"), ("C18", "C6"), ("H18", "H6"), ("O18", "O6"), ("C16", "C3"), ("H161", "H31"),
    ("H162", "H32"), ("C17", "C4"), ("H17", "H4"), ("O17", "O4"), ("HO17", "HO4"), ("C20", "C5"),
    ("H20", "H5"), ("N", "N"), ("HN", "HN"), ("C", "C"), ("O", "O"), ("CT", "CT"), ("HT1", "HT1"),
    ("HT2", "HT2"), ("HT3", "HT3"), ("C19", "C7"), ("H19", "H7"), ("O19", "O7"), ("HO19", "HO7"),
    ("C21", "C8"), ("H21", "H8"), ("O21", "O8"), ("HO21", "HO8"), ("C22", "C9"), ("H222", "H91"),
    ("H221", "H92"), ("O22", "O9"), ("HO22", "HO9"),
]


# CGenFF lipidated lysines -> peptide RTP entries; tail=None maps the remaining CGenFF atoms by name.
LIPIDATED = {
    "KTZ": dict(package="MTZP_cgenff_model_pH7_explicitH_gromacs", top="MTZP_cgenff_model_pH7_explicitH_gmx.top",
                count=139, unmapped=(137, 138, 139), tail=range(19, 135), offset=1, junction="C19", ca_charge=0.04,
                mapping={1: "O", 2: "C", 3: "CA", 4: "HA", 5: "CB", 6: "HB1", 7: "HB2", 8: "CG", 9: "HG1", 10: "HG2",
                         11: "CD", 12: "HD1", 13: "HD2", 14: "CE", 15: "HE1", 16: "HE2", 17: "NZ", 18: "HZ1",
                         135: "N", 136: "HN"}),
    "KRT": dict(package="MRTR_cgenff_model_pH7_explicitH_gromacs", top="MRTR_cgenff_model_pH7_explicitH_gmx.top",
                count=118, unmapped=(4, 5, 118), tail=range(23, 118), offset=-3, junction="C23", ca_charge=0.04,
                mapping={1: "C", 2: "N", 3: "HN", 6: "CA", 7: "HA", 8: "CB", 9: "HB1", 10: "HB2", 11: "CG",
                         12: "HG1", 13: "HG2", 14: "CD", 15: "HD1", 16: "HD2", 17: "CE", 18: "HE1", 19: "HE2",
                         20: "NZ", 21: "HZ1", 22: "O"}),
    "KSM": dict(package="MSEM_cgenff_model_Lys26_explicitH_gromacs", top="semaglutide_Lys26_CGenFF_ready_gmx.top",
                count=142, unmapped=(1, 2, 3, 63, 64, 65, 66, 67, 139, 140, 141, 142), tail=None, junction="C11",
                ca_charge=0.084,
                mapping={4: "N", 68: "HN", 5: "CA", 69: "HA", 6: "CB", 70: "HB1", 71: "HB2", 7: "CG", 72: "HG1",
                         73: "HG2", 8: "CD", 74: "HD1", 75: "HD2", 9: "CE", 76: "HE1", 77: "HE2", 10: "NZ",
                         78: "HZ1", 61: "C", 62: "O"}),
}


CYSG_HDB = """CYSG		37
1	1	HN	N	-C	CA
1	5	HA	CA	N	C	CB
1	6	HB1	CB	SG	CA
1	6	HB2	CB	CA	SG
1	6	H1A	C1	SG	C2
1	6	H1B	C1	C2	SG
1	6	H2A	C2	C3	C1
1	6	H4A	C4	C3	C2
1	6	H4B	C4	C3	H4A
1	6	H4C	C4	C3	H4B
1	6	H5A	C5	C3	C6
1	6	H5B	C5	C6	C3
1	6	H6A	C6	C5	C7
1	6	H6B	C6	C7	C5
1	6	H7A	C7	C8	C6
1	6	H9A	C9	C8	C7
1	6	H9B	C9	C8	H9A
1	6	H9C	C9	C8	H9B
1	6	H10A	C10	C8	C11
1	6	H10B	C10	C11	C8
1	6	H11A	C11	C10	C12
1	6	H11B	C11	C12	C10
1	6	H12A	C12	C13	C11
1	6	H14A	C14	C13	C12
1	6	H14B	C14	C13	H14A
1	6	H14C	C14	C13	H14B
1	6	H15A	C15	C13	C16
1	6	H15B	C15	C16	C13
1	6	H16A	C16	C15	C17
1	6	H16B	C16	C17	C15
1	6	H17A	C17	C18	C16
1	6	H19A	C19	C18	C17
1	6	H19B	C19	C18	H19A
1	6	H19C	C19	C18	H19B
1	4	H20A	C20	C18	C19
1	4	H20B	C20	C18	H20A
1	4	H20C	C20	C18	H20B
"""


EM_MDP = """integrator      = steep
emtol           = {emtol}
emstep          = 0.005
nsteps          = {nsteps}
cutoff-scheme   = Verlet
nstlist         = 20
rlist           = 1.2
coulombtype     = {coulombtype}
rcoulomb        = 1.2
vdwtype         = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw            = 1.2
pme-order       = 4
fourierspacing  = 0.12
pbc             = xyz
"""


@dataclass(frozen=True)
class Settings:
    """Tunable acceptance criteria; every default is documented in docs/tutorial section 8."""
    # Largest core RMSD (A) between all-atom backbone centres and CG BB beads after the rigid fit. Martini 3
    # elastic-network proteins stay within 1.2-2.6 A of their all-atom parent over 10 us in the repository
    # examples, so 5 A separates "same fold, thermal drift" from "different conformation or wrong chain".
    fit_max_core_rmsd_a: float = 5.0
    # Smallest fraction of backbone pairs the trimmed fit must retain. Below one half the "core" would be a
    # minority of the complex and the transform would no longer describe the assembly as a whole.
    fit_min_core_fraction: float = 0.5
    # The trimmed fit drops a pair when its deviation exceeds max(fit_trim_factor x median, fit_trim_floor_a). By
    # construction that keeps at least half of all pairs, so fit_min_core_fraction only guards against a changed
    # trimming rule; misplaced chains are caught by fit_max_chain_median_a.
    # A matched chain may differ from its CG counterpart (a peptide or helix that shifted during the CG run keeps
    # its all-atom pose), but not be somewhere else: its median backbone deviation in the complex fit must stay
    # below this. Observed in the repository examples: up to 7.9 A for genuine matches; a misplaced chain is tens of A.
    fit_max_chain_median_a: float = 10.0
    fit_trim_factor: float = 2.5
    fit_trim_floor_a: float = 2.0
    # An accepted sequence match must cover at least this fraction of the shorter of the two sequences.
    fit_min_coverage: float = 0.5
    # Chains whose match to the same CG segment is within this fraction of the best are interchangeable; geometry
    # then decides, and two assignments closer than fit_dead_heat_a (A, RMSD) are refused as ambiguous.
    fit_rival_identity: float = 0.95
    fit_dead_heat_a: float = 0.1
    # A chain maps onto a CG segment only above this sequence identity over at least this many aligned residues.
    # Unrelated sequences align at 5-10% identity; engineered peptides/G proteins in the examples align at >= 51%.
    fit_min_identity: float = 0.30
    fit_min_pairs: int = 5
    # Consecutive BB beads further apart than this start a new CG segment. Bonded BB beads sit at 0.35-0.40 nm
    # and one skipped residue (e.g. a non-standard residue the CG model omits) leaves about 0.55-0.70 nm.
    cg_chain_break_nm: float = 1.0
    # mstool minimizer iterations, as in the validated step-4 backmaps.
    backmap_em_steps: int = 20000
    # How many backmapping seeds to try when a heavy-atom bond threads a ring (stochastic atom placement).
    backmap_attempts: int = 5
    # A membrane atom (hydrogens included) starting closer than this to a solute atom makes steepest descent diverge:
    # EM survived 0.15 A and failed at 0.05-0.08 A in the repository runs. mstool's relaxation uses a soft-core
    # repulsion that is finite and force-free at zero distance, so such overlaps can survive it.
    min_start_contact_nm: float = 0.012
    # Independent-audit limits: protein CA RMSD across EM (nm) and the closest lipid-solute heavy-atom pair (nm).
    audit_max_ca_rmsd_nm: float = 0.15
    audit_min_contact_nm: float = 0.10
    # Membrane orientation. An OPM/OPRLM entry is accepted as the orientation reference only when the anchor chain
    # matches it by sequence (identity over the aligned residues) AND its membrane-embedded backbone (reference
    # residues with CA inside the hydrophobic slab) superposes within orient_opm_max_core_rmsd_a after trimming at
    # orient_opm_trim_floor_a; a 1 A core deviation over a 30 A bundle changes the normal by about 2 degrees, which
    # is below PPM's own tilt uncertainty. Beyond that the structure has moved and PPM is run on its real coordinates.
    orient_min_identity: float = 0.95
    orient_opm_max_core_rmsd_a: float = 1.0
    orient_opm_min_core_fraction: float = 0.7
    orient_opm_trim_floor_a: float = 1.0
    orient_opm_min_slab_residues: int = 20
    # PPM only repositions the submitted atoms (it writes 3 decimals), so its output must be a rigid copy of the
    # anchor to within rounding; anything larger means the wrong atoms were matched or the output is not the input's.
    orient_ppm_max_fit_rmsd_a: float = 0.05
    orient_ppm_min_matched_fraction: float = 0.9
    orient_ppm_timeout_s: int = 3600
    # CG registration after orientation. The CG membrane has a cavity shaped around the CG protein's own pose; a
    # rigid complex rotated by theta about the midplane moves the ends of a 17 A half-height bundle by 17 sin(theta)
    # A at the membrane surfaces, which at 20 degrees is 6 A, about one lipid diameter: beyond that the backmapped
    # membrane cannot be expected to accommodate the oriented complex. The depth offset between the CG pose and the
    # midplane placement is limited likewise (0.5 nm, about one phosphate-plane width).
    register_max_tilt_deg: float = 20.0
    register_max_depth_offset_nm: float = 0.5
    # The lateral registration is fitted on the anchor residues inside the hydrophobic slab (what must sit in the
    # lipid cavity) when at least this many are matched; otherwise on every matched anchor residue.
    register_min_embedded_pairs: int = 20
    # Water above and below the complex/membrane in the automatic box (nm).
    box_z_pad_nm: float = SLAB_Z_PAD_NM
