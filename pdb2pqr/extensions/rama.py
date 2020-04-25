"""
    Ramachandran extension

    Print both the phi and psi angles to standard out.  See the individual
    functions for more info.

    Author:  Mike Bradley and Todd Dolinsky
"""


import logging
from src.utilities import getDihedral
import extensions


_LOGGER = logging.getLogger(__name__)


def addExtensionOptions(extensionGroup):
    """
        Add options to set output type.
    """
    extensionGroup.parser.set_defaults(rama_output='rama')
    extensionGroup.add_option('--phi_only', dest='rama_output', action='store_const', const = 'phi',
                              help='Only include phi angles in output. '+
                                   'Rename output file {output-path}.phi')
    
    extensionGroup.add_option('--psi_only', dest='rama_output', action='store_const', const = 'psi',
                              help='Only include psi angles in output. '+
                                   'Rename output file {output-path}.psi')


def usage():
    return 'Print the per-residue phi and psi angles to {output-path}.rama for Ramachandran plots'

def create_rama_output(routines, outfile, outputtype='rama'):

    _LOGGER.debug("\nPrinting %s angles for each residue..." % (outputtype if outputtype != 'rama' else 'phi and psi'))
    verboseHeader = "Residue        %s"  % (outputtype.capitalize() if outputtype != 'rama' else 'Phi          Psi')
    _LOGGER.debug(verboseHeader)
    _LOGGER.debug('-' * len(verboseHeader))
    
    # Initialize some variables

    protein = routines.protein

    for residue in protein.getResidues():
        if residue.hasAtom("N"): 
            ncoords = residue.getAtom("N").getCoords()
        else: 
            continue

        if residue.hasAtom("CA"): 
            cacoords = residue.getAtom("CA").getCoords()
        else: 
            continue

        if residue.hasAtom("C"): 
            ccoords = residue.getAtom("C").getCoords()
        else: 
            continue

        try:
            if residue.peptideN != None:
                pepncoords = residue.peptideN.getCoords()
            else: 
                continue

            if residue.peptideC != None:
                pepccoords = residue.peptideC.getCoords()
            else: 
                continue
        except AttributeError: # Non amino acids
            continue

        _LOGGER.debug(str(residue))
        
        if outputtype in ('rama', 'phi'):
            phi = getDihedral(pepccoords, ncoords, cacoords, ccoords)
            _LOGGER.debug("\t%.4f" % phi)
            
        if outputtype in ('rama', 'psi'):
            psi = getDihedral(ncoords, cacoords, ccoords, pepncoords)
            _LOGGER.debug("\t%.4f" % psi)
            

    
def run_extension(routines, outroot, options):
    """
        Print the list of phi and psi angles for use in a Ramachandran plot.

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
            options:   options object 
    """
    outputType = options.rama_output
    outname = outroot + '.' + outputType
    with open(outname, "w") as outfile:
        create_rama_output(routines, outfile, outputtype=outputType)
    
