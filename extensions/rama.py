"""
    Ramachandran extension

    Print both the phi and psi angles to standard out.  See the individual
    functions for more info.

    Author:  Mike Bradley and Todd Dolinsky
"""

__date__ = "17 February 2006"
__author__ = "Mike Bradley, Todd Dolinsky"
  
from src.utilities import getDihedral
import extensions

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

    routines.write("\nPrinting %s angles for each residue...\n" % (outputtype if outputtype != 'rama' else 'phi and psi'))
    verboseHeader = "Residue        %s\n"  % (outputtype.capitalize() if outputtype != 'rama' else 'Phi          Psi')
    routines.write(verboseHeader)
    routines.write('-' * len(verboseHeader) + '\n')
    
    output = extensions.extOutputHelper(routines, outfile)

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

        output.write(str(residue))
        
        if outputtype in ('rama', 'phi'):
            phi = getDihedral(pepccoords, ncoords, cacoords, ccoords)
            output.write("\t%.4f" % phi)
            
        if outputtype in ('rama', 'psi'):
            psi = getDihedral(ncoords, cacoords, ccoords, pepncoords)
            output.write("\t%.4f" % psi)
            
        output.write('\n')

    routines.write("\n")
    
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
    
