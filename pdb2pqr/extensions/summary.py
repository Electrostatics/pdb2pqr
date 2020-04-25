"""Summary extension

Print protein summary file. 
Currently prints a list of all residue in protein.

Author:  Kyle Monson
"""
import logging


_LOGGER = logging.getLogger(__name__)


def usage():
    """
    Returns usage text for summary.
    """
    return 'Print protein summary information to {output-path}.summary.'


def create_summary_output(routines, outfile):
    """
    Output the interaction energy between each possible residue pair.
    """
    _LOGGER.debug("Printing protein summary...")
    
    output = extensions.extOutputHelper(routines, outfile)
    
    output.write(routines.protein.getSummary() + '\n')
    
    for residue in routines.protein.getResidues():
        output.write(str(residue)+'\n')


def run_extension(routines, outroot, options):
    outname = outroot + ".summary"
    with open(outname, "w") as outfile:
        create_summary_output(routines, outfile)