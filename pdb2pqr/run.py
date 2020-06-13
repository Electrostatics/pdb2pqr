"""Routines for running the code with a given set of options and PDB files."""
import logging
import time
import tempfile
from pathlib import Path
from . import debump
from . import hydrogens
from . import forcefield
from . import aa
from . import na
from . import input_output as io


_LOGGER = logging.getLogger(__name__)


def run_pdb2pka(ph, force_field, pdb_list, ligand, pdb2pka_params):
    """Run PDB2PKA"""
    # TODO - we are not ready to deal with PDB2PKA yet
    raise NotImplementedError("TODO - fix and re-enable PDB2PKA")
    # if force_field.lower() != 'parse':
    #     PDB2PKAError('PDB2PKA can only be run with the PARSE force field.')

    # _LOGGER.info("Running PDB2PKA and applying at pH %.2f... ", ph)
    # init_params = pdb2pka_params.copy()
    # init_params.pop('pairene')
    # init_params.pop('clean_output')
    # results = pka.pre_init(original_pdb_list=pdb_list, ff=force_field, ligand=ligand,
    #                        **init_params)
    # TODO - this is a messed-up variable unpacking:
    # output_dir, protein, routines, forcefield, apbs_setup, \
    #     ligand_titratable_groups, maps, sd = results
    # mypkaRoutines = pka_routines.pKaRoutines(protein, routines, forcefield,
    #                                          apbs_setup, output_dir, maps, sd,
    #                                          restart=pdb2pka_params.get('clean_output'),
    #                                          pairene=pdb2pka_params.get('pairene'))

    # _LOGGER.info('Doing full pKa calculation')
    # mypkaRoutines.runpKa()
    # pdb2pka_warnings = mypkaRoutines.warnings[:]
    # _LOGGER.warning(pdb2pka_warnings)

    # residue_ph = {}
    # for pka_residue_tuple, calc_ph in mypkaRoutines.ph_at_0_5.items():
    #     tit_type, chain_id, number_str = pka_residue_tuple
    #     if tit_type == 'NTR':
    #         tit_type = 'N+'
    #     elif tit_type == 'CTR':
    #         tit_type = 'C-'

    #     key = ' '.join([tit_type, number_str, chain_id])
    #     residue_ph[key] = calc_ph
    # pformat(residue_ph)
    # protein.apply_pka_values(ff, ph, residue_ph)
    # _LOGGER.debug('Finished running PDB2PKA.')

def run_pdb2pqr(pdblist, my_protein, my_definition, options, is_cif):
    """Run the PDB2PQR Suite"""
    raise DeprecationWarning("TODO - This function is deprecated")
#     """Run the PDB2PQR Suite
#     Args:
#         pdblist: The list of objects that was read from the PDB file given as
#                  input (list)
#         my_protein: Protein object
#         options: The name of the forcefield (string)
#         is_cif:  Boolean indicating whether input is CIF

#     Returns
#         A dictionary with the following elements:
#         * header:  The PQR file header (string)
#         * lines:  The PQR file atoms (list)
#         * missed_ligands:  A list of ligand residue names whose charges could
#                            not be assigned (ligand)
#         * protein:  The protein object
#     """
#     pkaname = ""
#     lines = []
#     ligand = None
#     atomcount = 0
#     output_pqr = Path(options.output_pqr)
#     outroot = output_pqr.stem

#     if options.pka_method == 'propka':
#         pkaname = Path(outroot + ".propka")
#         if pkaname.is_file():
#             _LOGGER.warning("PROPKA file already exists: %s", pkaname)

#     start = time.time()

#     ligsuccess = 0
#     if options.ligand is not None:
#         # If this is independent, we can assign charges and radii here
#         for residue in my_protein.residues:
#             if isinstance(residue, aa.LIG):
#                 templist = []
#                 ligand.make_up2date(residue)
#                 for atom in residue.atoms:
#                     atom.ffcharge = ligand.ligand_props[atom.name]["charge"]
#                     atom.radius = ligand.ligand_props[atom.name]["radius"]
#                     if atom in misslist:
#                         misslist.pop(misslist.index(atom))
#                         templist.append(atom)

#                 charge = residue.charge
#                 if abs(charge - int(charge)) > 0.001:
#                     # Ligand parameterization failed
#                     _LOGGER.warning(("WARNING: PDB2PQR could not successfully "
#                                      "parameterize the desired ligand; it has "
#                                      "been left out of the PQR file."))

#                     # remove the ligand
#                     my_protein.residues.remove(residue)
#                     for my_chain in my_protein.chains:
#                         if residue in my_chain.residues:
#                             my_chain.residues.remove(residue)
#                 else:
#                     ligsuccess = 1
#                     # Mark these atoms as hits
#                     hitlist = hitlist + templist

#     # Temporary fix; if ligand was successful, pull all ligands from misslist
#     if ligsuccess:
#         templist = misslist[:]
#         for atom in templist:
#             if isinstance(atom.residue, (aa.Amino, na.Nucleic)):
#                 continue
#             misslist.remove(atom)


#     # Determine if any of the atoms in misslist were ligands
#     missedligandresidues = []
#     for atom in misslist:
#         if isinstance(atom.residue, (aa.Amino, na.Nucleic)):
#             continue
#         if atom.res_name not in missedligandresidues:
#             missedligandresidues.append(atom.res_name)

#     _LOGGER.debug("Total time taken: %.2f seconds", (time.time() - start))
#     result_dict = {}
#     result_dict["header"] = header
#     result_dict["lines"] = lines
#     result_dict["missed_ligands"] = missedligandresidues
#     result_dict["protein"] = my_protein
#     return result_dict
