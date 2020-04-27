"""PDB2PQR

This package takes a PDB file as input and performs optimizations before
yielding a new PDB-style file as output.

For more information, see http://www.poissonboltzmann.org/
"""
__version__ = "3.0"


import logging
from pathlib import Path
from . import run
from . import pdb, cif, utilities, structures
from .errors import PDB2PQRError
from .propka import lib as propka_lib
from . import extensions
from .pdb2pka.ligandclean import ligff
from . import inputgen, psize


TITLE_TEXT = "PDB2PQR v{version} - biomolecular structure conversion software"
TITLE_TEXT = TITLE_TEXT.format(version=__version__)
CITE_TEXTS = [
    "Please cite:  Jurrus E, et al.  Improvements to the APBS biomolecular solvation software suite.  Protein Sci 27 112-128 (2018).  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5734301/",
    "Please cite:  Dolinsky TJ, et al.  PDB2PQR: expanding and upgrading automated preparation of biomolecular structures for molecular simulations.  Nucleic Acids Res 35 W522-W525 (2007).  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1933214/"
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
            # TODO - it makes me sad to open a file without a close() statement
            user_names_file = open(args.usernames, 'rt', encoding="utf-8")
        if args.userff is not None:
            # TODO - it makes me sad to open a file without a close() statement
            user_ff_file = open(args.userff, "rt", encoding="utf-8")
            if args.usernames is None:
                raise RuntimeError('--usernames must be specified if using --userff')
        if utilities.getFFfile(args.ff) == "":
            raise RuntimeError("Unable to load parameter file for forcefield %s" % args.ff)
        if (args.ph < 0) or (args.ph > 14):
            raise RuntimeError("Specified pH (%s) is outside the range [1, 14] of this program" % args.ph)
    
    ph_calc_options = None

    if args.pka_method == 'propka':
        ph_calc_options, _ = propka_lib.loadOptions('--quiet')
    elif args.pka_method == 'pdb2pka':
        if args.ff.lower() != 'parse':
            raise RuntimeError('PDB2PKA requires the PARSE force field.')
        ph_calc_options = {'output_dir': args.output_pqr,
                          'clean_output': not args.pdb2pka_resume,
                          'pdie': args.pdie,
                          'sdie': args.sdie,
                          'pairene': args.pairene}

    if args.ligand is not None:
        try:
            # TODO - it makes me sad to open a file without a close() statement
            ligand_file = open(args.ligand, 'rt', encoding="utf-8")
        except IOError:
            raise RuntimeError('Unable to find ligand file %s!' % args.ligand)

    if args.neutraln and (args.ff is None or args.ff.lower() != 'parse'):
        raise RuntimeError('--neutraln option only works with PARSE forcefield!')

    if args.neutralc and (args.ff is None or args.ff.lower() != 'parse'):
        raise RuntimeError('--neutralc option only works with PARSE forcefield!')


    path = Path(args.input_pdb)
    pdbFile = utilities.getPDBFile(args.input_pdb)

    args.isCIF = False
    if path.suffix.lower() == "cif":
        pdblist, errlist = cif.readCIF(pdbFile)
        args.isCIF = True
    else:
        pdblist, errlist = pdb.readPDB(pdbFile)

    if len(pdblist) == 0 and len(errlist) == 0:
        raise RuntimeError("Unable to find file %s!" % path)

    if len(errlist) != 0:
        if(args.isCIF):
            _LOGGER.warn("Warning: %s is a non-standard CIF file.\n", path)
        else:
            _LOGGER.warn("Warning: %s is a non-standard PDB file.\n", path)
        _LOGGER.error(errlist)

    args.outname = args.output_pqr

    # In case no extensions were specified or no extensions exist.
    # TODO - there are no command line options for extensions so I'm not sure what this does
    if not hasattr(args, 'active_extensions'):
        args.active_extensions = []
    elif args.active_extensions is None:
        args.active_extensions = []
    extensionOpts = args

    try:
        results_dict = run.runPDB2PQR(pdblist, args)
        header = results_dict["header"]
        lines = results_dict["lines"]
        missedligands = results_dict["missed_ligands"]
    except PDB2PQRError as error:
        _LOGGER.error(error)
        raise PDB2PQRError(error)

    # Print the PQR file
    # TODO - move this to another function... this function is already way too long.
    outfile = open(args.output_pqr,"w")
    outfile.write(header)
    # Adding whitespaces if --whitespace is in the options
    for line in lines:
        if args.whitespace:
            if line[0:4] == 'ATOM':
                newline = line[0:6] + ' ' + line[6:16] + ' ' + line[16:38] + ' ' + line[38:46] + ' ' + line[46:]
                outfile.write(newline)
            elif line[0:6] == 'HETATM':
                newline = line[0:6] + ' ' + line[6:16] + ' ' + line[16:38] + ' ' + line[38:46] + ' ' + line[46:]
                outfile.write(newline)
            elif line[0:3] == "TER" and args.isCIF:
                pass
        else:
            if line[0:3] == "TER" and args.isCIF:
                pass
            else:
                outfile.write(line)
    if(args.isCIF):
        outfile.write("#\n")
    outfile.close()

    if args.apbs_input:
        method = "mg-auto"
        size = psize.Psize()
        size.parseInput(args.output_pqr)
        size.runPsize(args.output_pqr)
        #async = 0 # No async files here!
        input = inputgen.Input(args.output_pqr, size, method, 0, potdx=True)
        input.printInputFiles()
        input.dumpPickle()


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    logging.captureWarnings(True)
