"""Hbond extension

Find all hydrogen bonds as determined by the cutoffs specified. Uses PDB2PQR to
determine donors and acceptors, and displays all available bonds to file. 

The original bonding parameters were an angle of 20.0, distance of 3.30, and
using the old method for calculating distance.

The original parameters for WHAT-IF output was an angle of 90.0, distance of
3.30, and using the old method for calculating distance. 

Authors:  Todd Dolinsky, Michael J Bradley, Julie Mitchell, and Kyle Monson
"""


# NOTE: This extension edited and updated on 05 August 2008 by Michael J Bradley
# and again on 17 October by Kyle Monson.


# NOTE: The current defaults for hbonds used below were utilized in the following study:
#Bradley MJ, Chivers PT, Baker NA. Molecular dynamics simulation of the 
#Escherichia coli NikR protein: Equilibrium conformational fluctuations reveal 
#inter-domain allosteric communication pathways. Journal of Molecular Biology, 378, 
#1155-1173, 2008.  http://dx.doi.org/10.1016/j.jmb.2008.03.010


import logging
from ..utilities import distance, getAngle
from ..routines import Cells
from math import cos


# TODO - this extension used to write to a file as well as stdout.  This logger may need to be modified to recreate the writing to a file.
_LOGGER = logging.getLogger(__name__)


ANGLE_CUTOFF = 30.0       # A - D - H(D) angle
DIST_CUTOFF = 3.4         # D to A distance


def addExtensionOptions(extensionGroup):
    """
        Add options to set output type, angle cutoff, distance cutoff, and distance calculating method.
    """
    extensionGroup.add_option('--whatif', dest='whatif', action='store_true', default=False,
                              help='Change hbond output to WHAT-IF format.')
    
    extensionGroup.add_option('--angle_cutoff', dest='angle_cutoff',  type="float", action='store', default=ANGLE_CUTOFF,
                              help='Angle cutoff to use when creating hbond data (default %s)' % ANGLE_CUTOFF)
    
    extensionGroup.add_option('--distance_cutoff', dest='distance_cutoff',  type="float", action='store', default=DIST_CUTOFF,
                              help='Distance cutoff to use when creating hbond data (default %s)' % DIST_CUTOFF)
    
    extensionGroup.add_option('--old_distance_method', dest='old_distance_method', action='store_true', default=False,
                              help='Use distance from donor hydrogen to acceptor to calculate distance used with --distance_cutoff.')


def usage():
    return 'Print a list of hydrogen bonds to {output-path}.hbond'


#TODO: replace this with a ''.format() call.
def _residueString(residue, name):
    return '%4d %-4s (%4d  ) %s     %-4s' % \
        (residue.res_seq, residue.name, residue.res_seq, residue.chain_id, name)


def create_hbond_output(routines, outfile, whatif=False, 
                                           angleCutoff=ANGLE_CUTOFF, 
                                           distanceCutoff=DIST_CUTOFF, 
                                           oldDistanceMethod=False):

    _LOGGER.debug("Printing hydrogen bond list...")
    
    cellsize = int(distanceCutoff + 1.0 + 1.0) 
    protein = routines.protein
    routines.set_donors_acceptors()
    routines.cells = Cells(cellsize)
    routines.cells.assign_cells(protein)

    for donor in protein.get_atoms():

        # Grab the list of donors
        if not donor.hdonor: 
            continue
        donorhs = []
        for bond in donor.bonds:
            if bond.isHydrogen(): 
                donorhs.append(bond)
        if donorhs == []: 
            continue

        # For each donor, grab all acceptors
        closeatoms = routines.cells.get_near_cells(donor)
        for acc in closeatoms:
            if not acc.hacceptor: 
                continue
            if donor.residue == acc.residue: 
                continue
            
            #TODO: do we need to do this for plain hbond stuff?
            if whatif and (donor.residue.chain_id == acc.residue.chain_id): 
                continue
            
            # Do new style distance check
            if not oldDistanceMethod:
                dist = distance(donor.getCoords(), acc.getCoords())
                if dist > distanceCutoff:
                    continue
            
            for donorh in donorhs:

                # Do old style distance check       
                if oldDistanceMethod:
                    dist = distance(donorh.getCoords(), acc.getCoords())
                    if dist > distanceCutoff: 
                        continue
                    
                # Do angle check
                angle = getAngle(acc.getCoords(), donor.getCoords(), donorh.getCoords())
                if angle > angleCutoff: 
                    continue
                
                if whatif:
                    if (donor.temp_factor > 60.0): 
                        continue
                    if (acc.temp_factor > 60.0): 
                        continue
                    
                    thisBstring='B' if donor.isBackbone() else 'S'
                    thatBstring='B' if acc.isBackbone() else 'S'

                    score= (1.7/dist) * cos(angle * 3.142 / 180.0)
                    _LOGGER.debug(_residueString(donor.residue, donor.name))
                    _LOGGER.debug('-> ')
                    _LOGGER.debug(_residueString(acc.residue, acc.name))
                    _LOGGER.debug('Sym=   1 Val= %6.3lf  DA=%6.2f  DHA=%6.2f (%s-%s)\n' % 
                                 (score, dist, angle, thisBstring, thatBstring))
#                    outfile.write("%4d %-4s (%4d  ) %s     %-4s-> %4d %-4s (%4d  ) %s     %-4sSym=   1 Val= %6.3lf  DA=%6.2f  DHA=%6.2f (%s-%s)\n" % \
#                      (donor.residue.res_seq,donor.residue.name,donor.residue.res_seq, donor.residue.chain_id,donor.name,acc.residue.res_seq,acc.residue.name,acc.residue.res_seq, acc.residue.chain_id,acc.name, score, dist, angle, thisBstring, thatBstring)) 

                else:
                    s = "Donor: %s %s\tAcceptor: %s %s\tdist: %.2f\tAngle: %.2f\n" % \
                        (donor.residue, donor.name, acc.residue, acc.name, dist, angle)
                    _LOGGER.debug(s) 
                    

def run_extension(routines, outroot, options):
    """
        Print a list of hydrogen bonds.

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
            options:   options object 
    """
    
    outname = outroot + ".hbond"
    with open(outname, "w") as outfile:
        create_hbond_output(routines, outfile, whatif=options.whatif, 
                            angleCutoff=options.angle_cutoff,
                            distanceCutoff=options.distance_cutoff,
                            oldDistanceMethod=options.old_distance_method)

