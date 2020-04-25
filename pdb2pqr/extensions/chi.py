"""Chi extension

Print the backbone chi angle for each residue in the structure. Chi angle is
determined by the coordinates of the N, CA, CB (if available), and CG/OG/SG
atoms (if available).

Author:  Todd Dolinsky
"""


import logging
from src.utilities import getDihedral


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

    for residue in protein.getResidues():
        if residue.hasAtom("N"): 
            ncoords = residue.getAtom("N").getCoords()
        else: 
            continue

        if residue.hasAtom("CA"):
            cacoords = residue.getAtom("CA").getCoords()
        else: 
            continue

        if residue.hasAtom("CB"): 
            cbcoords = residue.getAtom("CB").getCoords()
        else: 
            continue

        if residue.hasAtom("CG"): 
            gcoords = residue.getAtom("CG").getCoords()
        elif residue.hasAtom("OG"): 
            gcoords = residue.getAtom("OG").getCoords()
        elif residue.hasAtom("SG"): 
            gcoords = residue.getAtom("SG").getCoords()
        else: 
            continue

        chi = getDihedral(ncoords, cacoords, cbcoords, gcoords)
        _LOGGER.debug("%s\t%.4f" % (residue, chi))
        outfile.write("%s\t%.4f\n" % (residue, chi))
        
    outfile.close()
