"""Configuration information for PDB2PQR."""

# PDB2PQR version number.
VERSION = "3.0"

# How to format PDB2PQR title in output
TITLE_FORMAT_STRING = "PDB2PQR v{version} - biomolecular structure conversion software"

# Citation strings for PDB2PQR
CITATIONS = [("Please cite:  Jurrus E, et al.  Improvements to the APBS biomolecular "
              "solvation software suite.  Protein Sci 27 112-128 (2018)."),
             ("Please cite:  Dolinsky TJ, et al.  PDB2PQR: expanding and upgrading "
              "automated preparation of biomolecular structures for molecular simulations. "
              "Nucleic Acids Res 35 W522-W525 (2007).")]

# Standard force field names
FORCE_FIELDS = ["amber", "charmm", "parse", "tyl06", "peoepb", "swanson"]

# Standard amino acid names
AA_NAMES = ["ALA", "ARG", "ASH", "ASN", "ASP", "CYS", "CYM", "GLN", "GLU", "GLH",
            "GLY", "HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "ILE", "LEU",
            "LYS", "LYN", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "TYM",
            "VAL"]

# Standard nucleic acid names
NA_NAMES = ["A", "A5", "A3", "C", "C5", "C3", "G", "G5", "G3", "T", "T5", "T3",
            "U", "U5", "U3", "RA", "RG", "RC", "RU", "DA", "DG", "DC", "DT"]

# Standard backbone atom names
BACKBONE = ["N", "CA", "C", "O", "O2", "HA", "HN", "H", "tN"]

# A small number used by some math routines.
SMALL_NUMBER = 1.0e-7

# A number of unknown origin used in dihedral angle calculations
DIHEDRAL_WTF = 57.2958

# The start of warning strings to be filtered.
FILTER_WARNINGS = ["Skipped atom during water optimization",
                   "The best donorH was not picked",
                   "Multiple occupancies found"]

# The number of times one of the warning strings should be printed before
# supressing further output.
FILTER_WARNINGS_LIMIT = 20

# Expected location for topology definition file
TOPOLOGY_DEF_PATH = "dat/TOPOLOGY.xml"

# Expected location for amino acid topology definition file
AA_DEF_PATH = "dat/AA.xml"

# Expected location for nucleic acid topology definition file
NA_DEF_PATH = "dat/NA.xml"

# Expected location for hydrogens topology definition file
HYD_DEF_PATH = "dat/HYDROGENS.xml"

# Expected location for topology patch definition file
PATCH_DEF_PATH = "dat/PATCHES.xml"

# Number of angle steps to scan when debumping
DEBUMP_ANGLE_STEPS = 72

# Size of debumping step
DEBUMP_ANGLE_STEP_SIZE = float(360 // DEBUMP_ANGLE_STEPS)

# Debump angle test
DEBUMP_ANGLE_TEST_COUNT = 10

# Size of cells used for neighbor lookups
CELL_SIZE = 2

# Debumping hydrogen size
BUMP_HYDROGEN_SIZE = 0.5

# Debumping heavy atom size
BUMP_HEAVY_SIZE = 1.0

# Disulfide bond distance limit
BONDED_SS_LIMIT = 2.5

# Peptide bond distance limit
PEPTIDE_DIST = 1.7

# Limit on fraction of molecule missing before giving up on repairs
REPAIR_LIMIT = 0.1

# Cutoff for A - D - H(D) hydrogen bond angle
ANGLE_CUTOFF = 20.0

# Cutoff for H(D) to A hydrogen bond distance
DIST_CUTOFF = 3.3
