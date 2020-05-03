"""PDB2PQR

This package takes a PDB file as input and performs optimizations before
yielding a new PDB-style file as output.

For more information, see http://www.poissonboltzmann.org/
"""
__version__ = "3.0"
from sys import version_info
assert version_info >= (3, 5)
import logging
from pathlib import Path
from . import run
from . import pdb, cif, utilities, structures
from .propka import lib as propka_lib
from . import inputgen, psize


TITLE_TEXT = "PDB2PQR v{version} - biomolecular structure conversion software"
TITLE_TEXT = TITLE_TEXT.format(version=__version__)
CITE_TEXTS = [
    ("Please cite:  Jurrus E, et al.  Improvements to the APBS biomolecular "
     "solvation software suite.  Protein Sci 27 112-128 (2018)."),
    ("Please cite:  Dolinsky TJ, et al.  PDB2PQR: expanding and upgrading "
     "automated preparation of biomolecular structures for molecular "
     "simulations.  Nucleic Acids Res 35 W522-W525 (2007).")
]


_LOGGER = logging.getLogger(__name__)
logging.captureWarnings(True)



def main(args):
    """Main driver for running program from the command line.

    Args:
        args:  argument namespace object (e.g., as returned by argparse).
    """
    logging.basicConfig(level=getattr(logging, args.log_level))

    _LOGGER.debug("Args:  %s", args)
    _LOGGER.info(TITLE_TEXT)
    for citation in CITE_TEXTS:
        _LOGGER.info(citation)

    if args.assign_only or args.clean:
        args.debump = False
        args.opt = False

    if not args.clean:
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
        elif utilities.test_dat_file(args.ff) == "":
            raise RuntimeError("Unable to load parameter file for forcefield %s" % args.ff)
        if (args.ph < 0) or (args.ph > 14):
            raise RuntimeError(("Specified pH (%s) is outside the range [1, 14] "
                                "of this program") % args.ph)

    if args.ligand is not None:
        ligand = Path(args.ligand)
        if not ligand.is_file():
            error = "Unable to find ligand file: %s" % ligand
            raise FileNotFoundError(error)

    # TODO - it appears none of the following code is actually used
    # if args.pka_method == 'propka':
    #     ph_calc_options, _ = propka_lib.loadOptions('--quiet')
    # elif args.pka_method == 'pdb2pka':
    #     if args.ff.lower() != 'parse':
    #         raise RuntimeError('PDB2PKA requires the PARSE force field.')
    #     ph_calc_options = {'output_dir': args.output_pqr,
    #                        'clean_output': not args.pdb2pka_resume,
    #                        'pdie': args.pdie,
    #                        'sdie': args.sdie,
    #                        'pairene': args.pairene}
    # else:
    #     ph_calc_options = None

    if args.neutraln and (args.ff is None or args.ff.lower() != 'parse'):
        raise RuntimeError('--neutraln option only works with PARSE forcefield!')

    if args.neutralc and (args.ff is None or args.ff.lower() != 'parse'):
        raise RuntimeError('--neutralc option only works with PARSE forcefield!')


    path = Path(args.input_pdb)
    pdb_file = utilities.get_pdb_file(args.input_pdb)

    args.is_cif = False
    if path.suffix.lower() == "cif":
        pdblist, errlist = cif.read_cif(pdb_file)
        args.is_cif = True
    else:
        pdblist, errlist = pdb.read_pdb(pdb_file)

    if len(pdblist) == 0 and len(errlist) == 0:
        raise RuntimeError("Unable to find file %s!" % path)

    if len(errlist) != 0:
        if args.is_cif:
            _LOGGER.warning("Warning: %s is a non-standard CIF file.\n", path)
        else:
            _LOGGER.warning("Warning: %s is a non-standard PDB file.\n", path)
        _LOGGER.error(errlist)

    args.outname = args.output_pqr

    results_dict = run.run_pdb2pqr(pdblist, args)
    _ = results_dict["header"]
    lines = results_dict["lines"]
    _ = results_dict["missed_ligands"]

    # Print the PQR file
    # TODO - move this to another module.
    with open(args.output_pqr, "wt") as outfile:
        # Adding whitespaces if --whitespace is in the options
        for line in lines:
            if args.whitespace:
                if line[0:4] == 'ATOM':
                    newline = line[0:6] + ' ' + line[6:16] + ' ' + \
                        line[16:38] + ' ' + line[38:46] + ' ' + line[46:]
                    outfile.write(newline)
                elif line[0:6] == 'HETATM':
                    newline = line[0:6] + ' ' + line[6:16] + ' ' + \
                        line[16:38] + ' ' + line[38:46] + ' ' + line[46:]
                    outfile.write(newline)
                elif line[0:3] == "TER" and args.is_cif:
                    pass
            else:
                if line[0:3] == "TER" and args.is_cif:
                    pass
                else:
                    outfile.write(line)
        if args.is_cif:
            outfile.write("#\n")

    if args.apbs_input:
        method = "mg-auto"
        size = psize.Psize()
        size.parse_input(args.output_pqr)
        size.run_pize(args.output_pqr)
        input_ = inputgen.Input(args.output_pqr, size, method, 0, potdx=True)
        input_.print_input_files()
        input_.dump_pickle()


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    logging.captureWarnings(True)
