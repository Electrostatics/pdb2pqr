"""Chi extension

Print the backbone chi angle for each residue in the structure. Chi angle is
determined by the coordinates of the N, CA, CB (if available), and CG/OG/SG
atoms (if available).

Author:  Todd Dolinsky
"""
import logging
from ..utilities import dihedral


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

    for residue in protein.residues:
        if residue.has_atom("N"): 
            ncoords = residue.get_atom("N").coords
        else: 
            continue

        if residue.has_atom("CA"):
            cacoords = residue.get_atom("CA").coords
        else: 
            continue

        if residue.has_atom("CB"): 
            cbcoords = residue.get_atom("CB").coords
        else: 
            continue

        if residue.has_atom("CG"): 
            gcoords = residue.get_atom("CG").coords
        elif residue.has_atom("OG"): 
            gcoords = residue.get_atom("OG").coords
        elif residue.has_atom("SG"): 
            gcoords = residue.get_atom("SG").coords
        else: 
            continue

        chi = dihedral(ncoords, cacoords, cbcoords, gcoords)
        _LOGGER.debug("%s\t%.4f" % (residue, chi))
        outfile.write("%s\t%.4f\n" % (residue, chi))
        
    outfile.close()
