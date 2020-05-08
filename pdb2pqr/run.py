"""Routines for running the code with a given set of options and PDB files."""
import logging
import time
from pathlib import Path
from . import routines
from . import hydrogens
from . import forcefield
from . import aa
from . import na
from . import io


_LOGGER = logging.getLogger(__name__)


def run_propka_31(protein, pka_options):
    """Run PROPKA 3.1 on the current protein, setting protonation states to
    the correct values. pH is set in pka_options

    Parameters
        pka_options: Options for propKa 3.1, including pH

    Returns
        pka_molecule: pKa's internal molecule object (including pKa's, etc)
        not_found:    dict of residues found in pka_molecule but not in PDB2PQR (with pKa)
    """
    # See https://github.com/jensengroup/propka-3.1/blob/master/scripts/propka31.py

    ph = pka_options.ph
    _LOGGER.info("Running propka 3.1 at pH %.2f... ", ph)

    # Initialize some variables
    pkadic = {}

    # Reorder the atoms in each residue to start with N - TONI is this necessary?
    for residue in protein.residues:
        residue.reorder()

    # TONI Make a string with all non-hydrogen atoms. Previously it was removing the "element"
    # column and hydrogens. This does not seem to be necessary in propKa 3.1 .
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb") as h_free_file:
        for atom in protein.atoms:
            if not atom.is_hydrogen:
                atomtxt = atom.get_pdb_string()
                h_free_file.write(atomtxt + '\n')

        # Run PropKa 3.1 -------------
        # Creating protein object. Annoyingly, at this stage propka generates a
        # *.propka_input file in PWD and does not delete it (irrespective of the original
        # .pdb location)
        pka_molecule = propka.molecular_container.Molecular_container(h_free_file.name,
                                                                        pka_options)

    # calculating pKa values for ionizable residues -
    pka_molecule.calculate_pka()

    ##  pka_molecule.write_pka()
    for grp in pka_molecule.conformations['AVR'].groups:
        key = str.strip('%s %s %s' % (grp.residue_type, grp.atom.resNumb, grp.atom.chain_id))
        pkadic[key] = grp.pka_value

    protein.pka_protein = pka_molecule
    return pkadic


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


def run_propka(protein, ph, force_field, options, version=30):
    """Run PROPKA on the current protein, setting protonation states to the correct values

    Parameters
        ph:  The desired pH of the system
        force_field:  The forcefield name to be used
        outname: The name of the PQR outfile
        options: Options to propka
        version: may be 30 or 31 (uses external propka 3.1)
    """
    _LOGGER.info("Running PROPKA v%d and applying at pH %.2f... ", version, ph)
    pkadic = self.run_propka_31(options)

    if len(pkadic) == 0:
        raise ValueError("PROPKA returned empty results!")

    # Now apply each pka to the appropriate residue
    protein.apply_pka_values(force_field, ph, pkadic)
    _LOGGER.debug("Done running PROPKA")


def run_pdb2pqr(pdblist, my_protein, my_definition, options, is_cif):
    """Run the PDB2PQR Suite

    Args:
        pdblist: The list of objects that was read from the PDB file given as
                 input (list)
        my_protein: Protein object
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

    my_routines = routines.Debump(my_protein)

    if not options.assign_only:

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

        my_protein.add_hydrogens()
        my_hydrogen_routines = hydrogens.HydrogenRoutines(my_routines)

        if options.debump:
            my_routines.debump_protein()

        if options.opt:
            my_hydrogen_routines.set_optimizeable_hydrogens()
            my_protein.hold_residues(None)
            my_hydrogen_routines.initialize_full_optimization()
            my_hydrogen_routines.optimize_hydrogens()
        else:
            my_hydrogen_routines.initialize_wat_optimization()
            my_hydrogen_routines.optimize_hydrogens()

        # Special for GLH/ASH, since both conformations were added
        my_hydrogen_routines.cleanup()

    my_protein.set_states()
    my_forcefield = forcefield.Forcefield(options.ff, my_definition, options.userff,
                                          options.usernames)
    hitlist, misslist = my_protein.apply_force_field(my_forcefield)

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
        scheme = options.ffout
        userff = None # Currently not supported
        if scheme != options.ff:
            my_name_scheme = forcefield.Forcefield(scheme, my_definition, userff)
        else:
            my_name_scheme = my_forcefield
        my_protein.apply_name_scheme(my_name_scheme)

    if is_cif:
        header = io.print_pqr_header_cif(misslist, reslist, charge, options.ff,
                                      options.pka_method, options.ph, options.ffout,
                                      include_old_header=options.include_header)
    else:
        header = io.print_pqr_header(pdblist, misslist, reslist, charge, options.ff,
                                  options.pka_method, options.ph, options.ffout,
                                  include_old_header=options.include_header)

    lines = io.print_protein_atoms(hitlist, options.chain)

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
