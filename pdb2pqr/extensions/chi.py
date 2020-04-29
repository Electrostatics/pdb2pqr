"""Chi extension

Print the backbone chi angle for each residue in the structure. Chi angle is
determined by the coordinates of the N, CA, CB (if available), and CG/OG/SG
atoms (if available).

Author:  Todd Dolinsky
"""
import logging
from ..utilities import getDihedral


_LOGGER = logging.getLogger(__name__)


def usage():
    return 'Print the per-residue backbone chi angle to {output-path}.chi'


def run_extension(routines, outroot, options):
    """
        Print the list of psi angles

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
            options:   options object 
    """

    outname = outroot + ".chi"
    outfile = open(outname, "w")

    _LOGGER.debug("Printing chi angles for each residue...")
    _LOGGER.debug("Residue     chi")
    _LOGGER.debug("----------------")
    
    # Initialize some variables

    protein = routines.protein

    for residue in protein.get_residues():
        if residue.has_atom("N"): 
            ncoords = residue.get_atom("N").getCoords()
        else: 
            continue

        if residue.has_atom("CA"):
            cacoords = residue.get_atom("CA").getCoords()
        else: 
            continue

        if residue.has_atom("CB"): 
            cbcoords = residue.get_atom("CB").getCoords()
        else: 
            continue

        if residue.has_atom("CG"): 
            gcoords = residue.get_atom("CG").getCoords()
        elif residue.has_atom("OG"): 
            gcoords = residue.get_atom("OG").getCoords()
        elif residue.has_atom("SG"): 
            gcoords = residue.get_atom("SG").getCoords()
        else: 
            continue

        chi = getDihedral(ncoords, cacoords, cbcoords, gcoords)
        _LOGGER.debug("%s\t%.4f" % (residue, chi))
        outfile.write("%s\t%.4f\n" % (residue, chi))
        
    outfile.close()