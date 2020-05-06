"""Perform functions related to _main_ execution of PDB2PQR.

This module is intended for functions that directly touch arguments provided at
the invocation of PDB2PQR.  It was created to avoid cluttering the __init__.py
file.
"""
import logging
import argparse
from pathlib import Path
from . import run
from . import utilities
from .io import get_molecule, test_dat_file, dump_apbs, DuplicateFilter
from .definitions import Definition
from .config import VERSION, TITLE_FORMAT_STRING, CITATIONS, FORCE_FIELDS


_LOGGER = logging.getLogger(__name__)
_LOGGER.addFilter(DuplicateFilter())


def build_parser():
    """Build an argument parser.

    Return:
        ArgumentParser() object
    """

    desc = TITLE_FORMAT_STRING.format(version=VERSION)
    pars = argparse.ArgumentParser(description=desc,
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pars.add_argument("input_path",
                      help="Input PDB path or ID (to be retrieved from RCSB database")
    pars.add_argument("output_pqr", help="Output PQR path")
    pars.add_argument("--log-level", help="Logging level", default="INFO",
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    grp1 = pars.add_argument_group(title="Mandatory options",
                                   description="One of the following options must be used")
    grp1.add_argument("--ff", choices=[x.upper() for x in FORCE_FIELDS],
                      help="The forcefield to use.")
    grp1.add_argument("--userff",
                      help=("The user-created forcefield file to use. Requires "
                            "--usernames and overrides --ff"))
    grp1.add_argument("--clean", action='store_true', default=False,
                      help=("Do no optimization, atom addition, or parameter "
                            "assignment, just return the original PDB file in "
                            "aligned format. Overrides --ff and --userff"))
    grp2 = pars.add_argument_group(title="General options")
    grp2.add_argument('--nodebump', dest='debump', action='store_false',
                      default=True, help='Do not perform the debumping operation')
    grp2.add_argument('--noopt', dest='opt', action='store_false', default=True,
                      help='Do not perform hydrogen optimization')
    grp2.add_argument('--chain', action='store_true', default=False,
                      help='Keep the chain ID in the output PQR file')
    grp2.add_argument('--assign-only', action='store_true', default=False,
                      help=("Only assign charges and radii - do not add atoms, "
                            "debump, or optimize."))
    grp2.add_argument('--ffout', choices=[x.upper() for x in FORCE_FIELDS],
                      help=('Instead of using the standard canonical naming '
                            'scheme for residue and atom names, use the names '
                            'from the given forcefield'))
    grp2.add_argument('--usernames',
                      help=('The user-created names file to use. Required if '
                            'using --userff'))
    grp2.add_argument('--apbs-input', action='store_true', default=False,
                      help=('Create a template APBS input file based on the '
                            'generated PQR file.'))
    grp2.add_argument('--ligand',
                      help=('Calculate the parameters for the specified '
                            'MOL2-format ligand at the path specified by this '
                            'option.  PDB2PKA must be compiled.'))
    grp2.add_argument('--whitespace', action='store_true', default=False,
                      help=('Insert whitespaces between atom name and residue '
                            'name, between x and y, and between y and z.'))
    grp2.add_argument('--typemap', action='store_true', default=False,
                      help='Create Typemap output.')
    grp2.add_argument('--neutraln', action='store_true', default=False,
                      help=('Make the N-terminus of this protein neutral '
                            '(default is charged). Requires PARSE force field.'))
    grp2.add_argument('--neutralc', action='store_true', default=False,
                      help=('Make the C-terminus of this protein neutral '
                            '(default is charged). Requires PARSE force field.'))
    grp2.add_argument('--drop-water', action='store_true', default=False,
                      help='Drop waters before processing protein.')
    grp2.add_argument('--include-header', action='store_true', default=False,
                      help=('Include pdb header in pqr file. WARNING: The '
                            'resulting PQR file will not work with APBS versions '
                            'prior to 1.5'))
    grp3 = pars.add_argument_group(title="pKa options",
                                   description="Options for titration calculations")
    grp3.add_argument('--titration-state-method', dest="pka_method",
                      choices=('propka', 'pdb2pka'),
                      help=('Method used to calculate titration states. If a '
                            'titration state method is selected, titratable '
                            'residue charge states will be set by the pH value '
                            'supplied by --with_ph'))
    grp3.add_argument('--with-ph', dest='ph', type=float, action='store',
                      default=7.0,
                      help=('pH values to use when applying the results of the '
                            'selected pH calculation method.'))
    grp4 = pars.add_argument_group(title="PDB2PKA method options")
    grp4.add_argument('--pdb2pka-out', default='pdb2pka_output',
                      help='Output directory for PDB2PKA results.')
    grp4.add_argument('--pdb2pka-resume', action="store_true", default=False,
                      help='Resume run from state saved in output directory.')
    grp4.add_argument('--pdie', default=8.0,
                      help='Protein dielectric constant.')
    grp4.add_argument('--sdie', default=80.0,
                      help='Solvent dielectric constant.')
    grp4.add_argument('--pairene', default=1.0,
                      help='Cutoff energy in kT for pairwise pKa interaction energies.')
    grp5 = pars.add_argument_group(title="PROPKA method options")
    grp5.add_argument("--propka-reference", default="neutral",
                      choices=('neutral', 'low-pH'),
                      help=("Setting which reference to use for stability "
                            "calculations. See PROPKA 3.0 documentation."))
    return pars


def print_splash_screen(args):
    """Print argument overview and citation information.

    Args:
        args:  argparse namespace
    """
    _LOGGER.debug("Args:  %s", args)
    _LOGGER.info("%s", TITLE_FORMAT_STRING.format(version=VERSION))
    for citation in CITATIONS:
        _LOGGER.info(citation)


def check_files(args):
    """Check for other necessary files.

    Args:
        args:  argparse namespace
    Raises:
        FileNotFoundError:  necessary files not found
        RuntimeError:  input argument or file parsing problems
    """
    if args.usernames is not None:
        usernames = Path(args.usernames)
        if not usernames.is_file():
            error = "User-provided names file does not exist: %s" % usernames
            raise FileNotFoundError(error)

    if args.userff is not None:
        userff = Path(args.userff)
        if not userff.is_file():
            error = "User-provided forcefield file does not exist: %s" % userff
            raise FileNotFoundError(error)
        if args.usernames is None:
            raise RuntimeError('--usernames must be specified if using --userff')
    elif args.ff is not None:
        if test_dat_file(args.ff) == "":
            raise RuntimeError("Unable to load parameter file for forcefield %s" % args.ff)

    if args.ligand is not None:
        ligand = Path(args.ligand)
        if not ligand.is_file():
            error = "Unable to find ligand file: %s" % ligand
            raise FileNotFoundError(error)


def check_options(args):
    """Sanity check options.

    Args:
        args:  argparse namespace
    Raises:
        RuntimeError:  silly option combinations were encountered.
    """
    if (args.ph < 0) or (args.ph > 14):
        raise RuntimeError(("Specified pH (%s) is outside the range [1, 14] "
                            "of this program") % args.ph)

    if args.neutraln and (args.ff is None or args.ff.lower() != 'parse'):
        raise RuntimeError('--neutraln option only works with PARSE forcefield!')

    if args.neutralc and (args.ff is None or args.ff.lower() != 'parse'):
        raise RuntimeError('--neutralc option only works with PARSE forcefield!')


def print_pqr(args, pqr_lines, header_lines, missing_lines, is_cif):
    """Print output to specified file

    TODO - move this to another module (utilities)

    Args:
        args:  argparse namespace
        pqr_lines:  output lines (records)
        header_lines:  header lines
        missing_lines:  lines describing missing atoms (should go in header)
        is_cif:  flag indicating CIF-format
    """
    with open(args.output_pqr, "wt") as outfile:
        # Adding whitespaces if --whitespace is in the options
        if header_lines:
            _LOGGER.warning("Ignoring %d header lines in output.", len(header_lines))
        if missing_lines:
            _LOGGER.warning("Ignoring %d missing lines in output.", len(missing_lines))
        for line in pqr_lines:
            if args.whitespace:
                if line[0:4] == 'ATOM':
                    newline = line[0:6] + ' ' + line[6:16] + ' ' + \
                        line[16:38] + ' ' + line[38:46] + ' ' + line[46:]
                    outfile.write(newline)
                elif line[0:6] == 'HETATM':
                    newline = line[0:6] + ' ' + line[6:16] + ' ' + \
                        line[16:38] + ' ' + line[38:46] + ' ' + line[46:]
                    outfile.write(newline)
                elif line[0:3] == "TER" and is_cif:
                    pass
            else:
                if line[0:3] == "TER" and is_cif:
                    pass
                else:
                    outfile.write(line)
        if is_cif:
            outfile.write("#\n")


def transform_arguments(args):
    """Transform arguments with logic not provided by argparse.

    TODO - I wish this could be done with argparse.

    Args:
        args:  argparse namespace
    Returns:
        argparse namespace
    """
    if args.assign_only or args.clean:
        args.debump = False
        args.opt = False
    return args


def main(args):
    """Main driver for running program from the command line.

    Validate inputs, launch PDB2PQR, handle output.

    Args:
        args:  argument namespace object (e.g., as returned by argparse).
    """
    logging.basicConfig(level=getattr(logging, args.log_level))
    args = transform_arguments(args)
    print_splash_screen(args)
    check_files(args)
    check_options(args)

    pdblist, is_cif = get_molecule(args.input_path)
    results = run.run_pdb2pqr(pdblist, args, is_cif)

    print_pqr(args=args, pqr_lines=results["lines"], header_lines=results["header"],
              missing_lines=results["missed_ligands"], is_cif=is_cif)

    if args.apbs_input:
        dump_apbs(args.output_pqr)
