"""Routines for running the code with a given set of options and PDB files."""
import logging
import time
import io
from pathlib import Path
from . import definitions
from . import protein
from . import routines
from . import hydrogens
from . import forcefield
from . import aa
from . import na
from . import pdb
from .io import print_pqr_header_cif, print_pqr_header


_LOGGER = logging.getLogger(__name__)


def run_pdb2pqr(pdblist, my_definition, options, is_cif):
    """Run the PDB2PQR Suite

    Args:
        pdblist: The list of objects that was read from the PDB file given as
                 input (list)
        options: The name of the forcefield (string)
        is_cif:  Boolean indicating whether input is CIF

    Returns
        A dictionary with the following elements:
        * header:  The PQR file header (string)
        * lines:  The PQR file atoms (list)
        * missed_ligands:  A list of ligand residue names whose charges could
                           not be assigned (ligand)
        * protein:  The protein object
    """
    pkaname = ""
    lines = []
    ligand = None
    atomcount = 0
    output_pqr = Path(options.output_pqr)
    outroot = output_pqr.stem

    if options.pka_method == 'propka':
        pkaname = Path(outroot + ".propka")
        if pkaname.is_file():
            _LOGGER.warning("PROPKA file already exists: %s", pkaname)

    start = time.time()
    _LOGGER.info("Beginning PDB2PQR...")


    if options.drop_water:
        pdblist = routines.drop_water(pdblist)

    # Check for the presence of a ligand!  This code is taken from pdb2pka/pka.py
    if options.ligand is not None:
        raise NotImplementedError("Ligand functionality is temporarily disabled.")
        with open(options.ligand, "rt", encoding="utf-8") as ligand_file:
            my_protein, my_definition, ligand = ligff.initialize(my_definition,
                                                                 ligand_file,
                                                                 pdblist)
        for atom in my_protein.atoms:
            if atom.type == "ATOM":
                atomcount += 1
    else:
        my_protein = protein.Protein(pdblist, my_definition)

    _LOGGER.info("Created protein object:")
    _LOGGER.info("  Number of residues in protein: %s", len(my_protein.residues))
    _LOGGER.info("  Number of atoms in protein   : %s", len(my_protein.atoms))

    my_routines = routines.Routines(my_protein)
    for residue in my_protein.residues:
        multoccupancy = 0
        for atom in residue.atoms:
            if atom.alt_loc != "":
                multoccupancy = 1
                txt = "Warning: multiple occupancies found: %s in %s." % (atom.name,
                                                                          residue)
                _LOGGER.warning(txt)
        if multoccupancy == 1:
            _LOGGER.warning(("WARNING: multiple occupancies found in %s at least "
                             "one of the instances is being ignored."), residue)

    my_routines.set_termini(options.neutraln, options.neutralc)
    my_routines.update_bonds()

    if options.clean:
        header = ""
        lines = my_protein.print_atoms(my_protein.atoms, options.chain)
        _LOGGER.debug("Total time taken: %.2f seconds", (time.time() - start))

        #Be sure to include None for missed ligand residues
        results_dict = {"header": header, "lines": lines,
                        "missed_ligands": None, "protein": my_protein}
        return results_dict

    #remove any future need to convert to lower case
    if options.userff is not None:
        force_field = options.userff.lower()
    elif options.ff is not None:
        force_field = options.ff.lower()
    if options.ffout is not None:
        ffout = options.ffout.lower()

    if not options.assign_only:
        # It is OK to process ligands with no ATOM records in the pdb
        if atomcount == 0 and ligand is not None:
            pass
        else:
            my_routines.find_missing_heavy()
        my_routines.update_ss_bridges()

        if options.debump:
            my_routines.debump_protein()

        # TODO - both PROPKA and PDB2PKA are messed up
        if options.pka_method == 'propka':
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
            my_routines.run_propka(options.ph, options.ff, options={})
        elif options.pka_method == 'pdb2pka':
            raise NotImplementedError("PROPKA is broken.")
            # my_routines.run_pdb2pka(options.ph, options.ff, pdblist, ligand, ph_calc_options)

        my_routines.add_hydrogens()
        my_hydrogen_routines = hydrogens.HydrogenRoutines(my_routines)

        if options.debump:
            my_routines.debump_protein()

        if options.opt:
            my_hydrogen_routines.set_optimizeable_hydrogens()
            my_routines.hold_residues(None)
            my_hydrogen_routines.initialize_full_optimization()
            my_hydrogen_routines.optimize_hydrogens()
        else:
            my_hydrogen_routines.initialize_wat_optimization()
            my_hydrogen_routines.optimize_hydrogens()

        # Special for GLH/ASH, since both conformations were added
        my_hydrogen_routines.cleanup()

    else:  # Special case for HIS if using assign-only
        for residue in my_protein.residues:
            if isinstance(residue, aa.HIS):
                my_routines.apply_patch("HIP", residue)

    my_routines.set_states()
    my_forcefield = forcefield.Forcefield(force_field, my_definition, options.userff,
                                          options.usernames)
    hitlist, misslist = my_routines.apply_force_field(my_forcefield)

    ligsuccess = 0
    if options.ligand is not None:
        # If this is independent, we can assign charges and radii here
        for residue in my_protein.residues:
            if isinstance(residue, aa.LIG):
                templist = []
                ligand.make_up2date(residue)
                for atom in residue.atoms:
                    atom.ffcharge = ligand.ligand_props[atom.name]["charge"]
                    atom.radius = ligand.ligand_props[atom.name]["radius"]
                    if atom in misslist:
                        misslist.pop(misslist.index(atom))
                        templist.append(atom)

                charge = residue.charge
                if abs(charge - int(charge)) > 0.001:
                    # Ligand parameterization failed
                    _LOGGER.warning(("WARNING: PDB2PQR could not successfully "
                                     "parameterize the desired ligand; it has "
                                     "been left out of the PQR file."))

                    # remove the ligand
                    my_protein.residues.remove(residue)
                    for my_chain in my_protein.chains:
                        if residue in my_chain.residues:
                            my_chain.residues.remove(residue)
                else:
                    ligsuccess = 1
                    # Mark these atoms as hits
                    hitlist = hitlist + templist

    # Temporary fix; if ligand was successful, pull all ligands from misslist
    if ligsuccess:
        templist = misslist[:]
        for atom in templist:
            if isinstance(atom.residue, (aa.Amino, na.Nucleic)):
                continue
            misslist.remove(atom)

    # Create the Typemap
    if options.typemap:
        typemapname = "%s-typemap.html" % outroot
        my_protein.create_html_typemap(my_definition, typemapname)

    # Grab the protein charge
    reslist, charge = my_protein.charge

    # If we want a different naming scheme, use that
    if options.ffout is not None:
        scheme = ffout
        userff = None # Currently not supported
        if scheme != force_field:
            my_name_scheme = forcefield.Forcefield(scheme, my_definition, userff)
        else:
            my_name_scheme = my_forcefield
        my_routines.apply_name_scheme(my_name_scheme)

    if is_cif:
        header = print_pqr_header_cif(misslist, reslist, charge, force_field,
                                      options.pka_method, options.ph, options.ffout,
                                      include_old_header=options.include_header)
    else:
        header = print_pqr_header(pdblist, misslist, reslist, charge, force_field,
                                  options.pka_method, options.ph, options.ffout,
                                  include_old_header=options.include_header)

    lines = my_protein.print_atoms(hitlist, options.chain)

    # Determine if any of the atoms in misslist were ligands
    missedligandresidues = []
    for atom in misslist:
        if isinstance(atom.residue, (aa.Amino, na.Nucleic)):
            continue
        if atom.res_name not in missedligandresidues:
            missedligandresidues.append(atom.res_name)

    _LOGGER.debug("Total time taken: %.2f seconds", (time.time() - start))
    result_dict = {}
    result_dict["header"] = header
    result_dict["lines"] = lines
    result_dict["missed_ligands"] = missedligandresidues
    result_dict["protein"] = my_protein
    return result_dict
