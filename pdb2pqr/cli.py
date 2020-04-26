"""Command line option parser for PDB2PQR"""
import logging
import argparse
from . import __version__
from . import aa


_LOGGER = logging.getLogger(__name__)
FIELD_NAMES = ('amber', 'charmm', 'parse', 'tyl06', 'peoepb', 'swanson')


def build_parser():
    """Build an argument parser.

    Return:
        ArgumentParser() object
    """

    desc = "PDB2PQR {version}.  Wields awesome powers to turn PDBs into PQRs."
    desc = desc.format(version=__version__)
    p = argparse.ArgumentParser(description=desc,
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("input_pdb",
                   help="Input PDB path or ID (to be retrieved from RCSB database")
    p.add_argument("output_pqr", help="Output PQR path")
    p.add_argument("--log-level", help="Logging level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    g1 = p.add_argument_group(title="Mandatory options",
                             description="One of the following options must be used")
    g1.add_argument("--ff", choices=[x.upper() for x in FIELD_NAMES],
                   help="The forcefield to use.")
    g1.add_argument("--userff",
                   help="The user-created forcefield file to use. Requires --usernames and overrides --ff")
    g1.add_argument("--clean", action='store_true', default=False,
                   help="Do no optimization, atom addition, or parameter assignment, just return the original PDB file in aligned format. Overrides --ff and --userff")
    g2 = p.add_argument_group(title="General options")
    g2.add_argument('--nodebump', dest='debump', action='store_false',
                    default=True, help='Do not perform the debumping operation')
    g2.add_argument('--noopt', dest='opt', action='store_false', default=True,
                    help='Do not perform hydrogen optimization')
    g2.add_argument('--chain', action='store_true', default=False,
                    help='Keep the chain ID in the output PQR file')
    g2.add_argument('--assign-only', action='store_true', default=False,
                    help='Only assign charges and radii - do not add atoms, debump, or optimize.')
    g2.add_argument('--ffout', choices=[x.upper() for x in FIELD_NAMES],
                    help='Instead of using the standard canonical naming scheme for residue and atom names, use the names from the given forcefield')
    g2.add_argument('--usernames', 
                    help='The user created names file to use. Required if using --userff')
    g2.add_argument('--apbs-input', action='store_true', default=False,
                    help='Create a template APBS input file based on the generated PQR file.  Also creates a Python pickle for using these parameters in other programs.')
    g2.add_argument('--ligand',
                    help='Calculate the parameters for the specified MOL2-format ligand at the path specified by this option.  PDB2PKA must be compiled.')
    g2.add_argument('--whitespace', action='store_true', default=False,
                    help='Insert whitespaces between atom name and residue name, between x and y, and between y and z.')
    g2.add_argument('--typemap', action='store_true', default=False,
                    help='Create Typemap output.')
    g2.add_argument('--neutraln', action='store_true', default=False,
                    help='Make the N-terminus of this protein neutral (default is charged). Requires PARSE force field.')
    g2.add_argument('--neutralc', action='store_true', default=False,
                    help='Make the C-terminus of this protein neutral (default is charged). Requires PARSE force field.')
    g2.add_argument('--drop-water', action='store_true', default=False,
                    help='Drop waters (%s) before processing protein.' % aa.WAT.water_residue_names)

    g2.add_argument('--include-header', action='store_true', default=False,
                    help='Include pdb header in pqr file. WARNING: The resulting PQR file will not work with APBS versions prior to 1.5')
    g3 = p.add_argument_group(title="pKa options",
                              description="Options for titration calculations")
    g3.add_argument('--titration-state-method', dest="pka_method", 
                    choices=('propka', 'pdb2pka'),
                    help='Method used to calculate titration states. If a titration state method is selected, titratable residue charge states will be set by the pH value supplied by --with_ph')
    g3.add_argument('--with-ph', dest='ph', type=float, action='store', default=7.0,
                    help='pH values to use when applying the results of the selected pH calculation method.')
    g4 = p.add_argument_group(title="PDB2PKA method options")
    g4.add_argument('--pdb2pka-out', default='pdb2pka_output',
                    help='Output directory for PDB2PKA results.')
    g4.add_argument('--pdb2pka-resume', action="store_true", default=False,
                    help='Resume run from state saved in output directory.')
    g4.add_argument('--pdie', default=8.0,
                    help='Protein dielectric constant.')
    g4.add_argument('--sdie', default=80.0,
                    help='Solvent dielectric constant.')
    g4.add_argument('--pairene',  default=1.0,
                    help='Cutoff energy in kT for pairwise pKa interaction energies.')
    g5 = p.add_argument_group(title="PROPKA method options")
    g5.add_argument("--propka-reference", default="neutral", choices=('neutral','low-pH'),
                    help="Setting which reference to use for stability calculations. See PROPKA 3.0 documentation.")
    return p