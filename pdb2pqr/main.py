"""Perform functions related to _main_ execution of PDB2PQR.

This module is intended for functions that directly touch arguments provided at
the invocation of PDB2PQR.
It was created to avoid cluttering the __init__.py file.

.. codeauthor:: Nathan Baker (et al.)
"""

import argparse
import logging
import sys
from collections import OrderedDict
from collections.abc import Sequence
from io import StringIO
from os import PathLike
from pathlib import Path

import propka.input as pk_in
import propka.lib
import propka.output as pk_out
from propka.molecular_container import MolecularContainer
from propka.parameters import Parameters

from . import aa, debump, forcefield, hydrogens, io
from . import biomolecule as biomol
from .config import (
    CITATIONS,
    FORCE_FIELDS,
    IGNORED_PROPKA_OPTIONS,
    REPAIR_LIMIT,
    TITLE_STR,
    VERSION,
)
from .ligand.mol2 import Mol2Molecule
from .utilities import noninteger_charge

_LOGGER = logging.getLogger(f"PDB2PQR{VERSION}")


def build_main_parser():
    """Build an argument parser.

    .. todo::
        Need separate argparse groups for PDB2PKA and PROPKA.
        These exist but need real options.

    :returns:  argument parser
    :rtype:  argparse.ArgumentParser
    """
    desc = TITLE_STR
    pars = argparse.ArgumentParser(
        prog="pdb2pqr",
        description=desc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        conflict_handler="resolve",
    )
    pars.add_argument(
        "input_path",
        help="Input PDB path or ID (to be retrieved from RCSB database",
    )
    pars.add_argument("output_pqr", help="Output PQR path")
    pars.add_argument(
        "--log-level",
        help="Logging level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    grp1 = pars.add_argument_group(
        title="Mandatory options",
        description="One of the following options must be used",
    )
    grp1.add_argument(
        "--ff",
        choices=[ff.upper() for ff in FORCE_FIELDS],
        default="PARSE",
        help="The forcefield to use.",
    )
    grp1.add_argument(
        "--userff",
        help=(
            "The user-created forcefield file to use. Requires "
            "--usernames and overrides --ff"
        ),
    )
    grp1.add_argument(
        "--clean",
        action="store_true",
        default=False,
        help=(
            "Do no optimization, atom addition, or parameter assignment, "
            "just return the original PDB file in aligned format. Overrides "
            "--ff and --userff"
        ),
    )
    grp2 = pars.add_argument_group(title="General options")
    grp2.add_argument(
        "--nodebump",
        dest="debump",
        action="store_false",
        default=True,
        help="Do not perform the debumping operation",
    )
    grp2.add_argument(
        "--noopt",
        dest="opt",
        action="store_false",
        default=True,
        help="Do not perform hydrogen optimization",
    )
    grp2.add_argument(
        "--keep-chain",
        action="store_true",
        default=False,
        help="Keep the chain ID in the output PQR file",
    )
    grp2.add_argument(
        "--assign-only",
        action="store_true",
        default=False,
        help=(
            "Only assign charges and radii - do not add atoms, "
            "debump, or optimize."
        ),
    )
    grp2.add_argument(
        "--ffout",
        choices=[ff.upper() for ff in FORCE_FIELDS],
        help=(
            "Instead of using the standard canonical naming scheme for "
            "residue and atom names, use the names from the given forcefield"
        ),
    )
    grp2.add_argument(
        "--usernames",
        help=(
            "The user-created names file to use. Required if using --userff"
        ),
    )
    grp2.add_argument(
        "--apbs-input",
        help=(
            "Create a template APBS input file based on the generated PQR "
            "file at the specified location."
        ),
    )
    grp2.add_argument(
        "--pdb-output",
        default=None,
        help=(
            "Create a PDB file based on input. This will be missing charges "
            "and radii"
        ),
    )
    grp2.add_argument(
        "--ligand",
        help=(
            "Calculate the parameters for a single MOL2-format ligand at the "
            "path specified by this option.  The MOL2 ligand name should "
            "match only one ligand in the PDB file."
        ),
    )
    grp2.add_argument(
        "--whitespace",
        action="store_true",
        default=False,
        help=(
            "Insert whitespaces between atom name and residue name, between x "
            "and y, and between y and z."
        ),
    )
    grp2.add_argument(
        "--neutraln",
        action="store_true",
        default=False,
        help=(
            "Make the N-terminus of a protein neutral (default is "
            "charged). Requires PARSE force field."
        ),
    )
    grp2.add_argument(
        "--neutralc",
        action="store_true",
        default=False,
        help=(
            "Make the C-terminus of a protein neutral (default is "
            "charged). Requires PARSE force field."
        ),
    )
    grp2.add_argument(
        "--drop-water",
        action="store_true",
        default=False,
        help="Drop waters before processing biomolecule.",
    )
    grp2.add_argument(
        "--include-header",
        action="store_true",
        default=False,
        help=(
            "Include pdb header in pqr file. WARNING: The resulting PQR file "
            "will not work with APBS versions prior to 1.5"
        ),
    )
    grp3 = pars.add_argument_group(
        title="pKa options", description="Options for titration calculations"
    )
    grp3.add_argument(
        "--titration-state-method",
        dest="pka_method",
        choices=(["propka"]),
        help=(
            "Method used to calculate titration states. If a titration state "
            "method is selected, titratable residue charge states will be set "
            "by the pH value supplied by --with_ph"
        ),
    )
    grp3.add_argument(
        "--with-ph",
        dest="ph",
        type=float,
        action="store",
        default=7.0,
        help=(
            "pH values to use when applying the results of the selected pH "
            "calculation method."
        ),
    )
    pars: argparse.ArgumentParser = propka.lib.build_parser(pars)

    # Override version flag set by PROPKA
    pars.add_argument(
        "--version", action="version", version=f"%(prog)s {VERSION}"
    )

    return pars


def print_splash_screen(args):
    """Print argument overview and citation information.

    :param args:  command-line arguments
    :type args:  argparse.Namespace
    """
    _LOGGER.debug(f"Args:  {args}")
    _LOGGER.info(TITLE_STR)
    for citation in CITATIONS:
        _LOGGER.info(citation)


def check_files(args):
    """Check for other necessary files.

    :param args:  command-line arguments
    :type args:  argparse.Namespace
    :raises FileNotFoundError:  necessary files not found
    :raises RuntimeError:  input argument or file parsing problems
    """
    if args.usernames is not None:
        usernames = Path(args.usernames)
        if not usernames.is_file():
            error = f"User-provided names file does not exist: {usernames}"
            raise FileNotFoundError(error)
    if args.userff is not None:
        userff = Path(args.userff)
        if not userff.is_file():
            error = f"User-provided forcefield file does not exist: {userff}"
            raise FileNotFoundError(error)
        if args.usernames is None:
            err = "--usernames must be specified if using --userff"
            raise RuntimeError(err)
    elif args.ff is not None:
        io.test_dat_file(args.ff)
    if args.ligand is not None:
        ligand = Path(args.ligand)
        if not ligand.is_file():
            error = f"Unable to find ligand file: {ligand}"
            raise FileNotFoundError(error)


def check_options(args):
    """Sanity check options.

    :param args:  command-line arguments
    :type args:  argparse.Namespace
    :raises RuntimeError:  silly option combinations were encountered.
    """
    for option, new_value in IGNORED_PROPKA_OPTIONS.items():
        if option in args:
            value = getattr(args, option)
            if value:
                _LOGGER.warning(
                    f"PROPKA option '{option}' {getattr(args, option)} is not "
                    f"processed correctly by PDB2PQR. Ignoring."
                )
                setattr(args, option, new_value)
    if (args.ph < 0) or (args.ph > 14):
        err = (
            f"Specified pH ({args.ph}) is outside the range "
            "[1, 14] of this program"
        )
        raise RuntimeError(err)
    if args.neutraln and (args.ff is None or args.ff.lower() != "parse"):
        err = "--neutraln option only works with PARSE forcefield!"
        raise RuntimeError(err)
    if args.neutralc and (args.ff is None or args.ff.lower() != "parse"):
        err = "--neutralc option only works with PARSE forcefield!"
        raise RuntimeError(err)


def print_pqr(args, pqr_lines, header_lines, missing_lines, is_cif):
    """Print PQR-format output to specified file

    .. todo::  Move this to another module (io)

    :param argparse.Namespace args:  command-line arguments
    :param [str] pqr_lines:  output lines (records)
    :param [str] header_lines:  header lines
    :param [str] missing_lines:  lines describing missing atoms (should go
        in header)
    :param bool is_cif:  flag indicating CIF format
    """
    with open(args.output_pqr, "w") as outfile:
        # Adding whitespaces if --whitespace is in the options
        if header_lines:
            _LOGGER.warning(
                f"Ignoring {len(header_lines)} header lines in output."
            )
        if missing_lines:
            _LOGGER.warning(
                f"Ignoring {len(missing_lines)} missing lines in output."
            )
        for line in pqr_lines:
            if args.whitespace:
                if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                    newline = (
                        line[0:6]
                        + " "
                        + line[6:16]
                        + " "
                        + line[16:38]
                        + " "
                        + line[38:46]
                        + " "
                        + line[46:]
                    )
                    outfile.write(newline)
            elif line[0:3] != "TER" or not is_cif:
                outfile.write(line)
        if is_cif:
            outfile.write("#\n")


def print_pdb(args, pdb_lines, header_lines, missing_lines, is_cif):
    """Print PDB-format output to specified file

    .. todo::  Move this to another module (io)

    :param argparse.Namespace args:  command-line arguments
    :param [str]] pdb_lines:  output lines (records)
    :param [str] header_lines:  header lines
    :param [str] missing_lines:  lines describing missing atoms (should go in
        header)
    :param bool is_cif:  flag indicating CIF format
    """
    with open(args.pdb_output, "w") as outfile:
        # Adding whitespaces if --whitespace is in the options
        if header_lines:
            _LOGGER.warning(
                f"Ignoring {len(header_lines)} header lines in output."
            )
        if missing_lines:
            _LOGGER.warning(
                f"Ignoring {len(missing_lines)} missing lines in output."
            )
        for line in pdb_lines:
            if line[0:3] != "TER" or not is_cif:
                outfile.write(line)


def transform_arguments(args):
    """Transform arguments with logic not provided by argparse.

    .. todo::  I wish this could be done with argparse.

    :param args:  command-line arguments
    :type args:  argparse.Namespace
    :return:  modified arguments
    :rtype:  argparse.Namespace
    """
    if args.assign_only or args.clean:
        args.debump = False
        args.opt = False
    if args.userff is not None:
        args.userff = args.userff
    elif args.ff is not None:
        args.ff = args.ff.lower()
    if args.ffout is not None:
        args.ffout = args.ffout.lower()
    return args


def setup_molecule(pdblist, definition, ligand_path):
    """Set up the molecular system.

    :param pdblist:  list of PDB records
    :type pdblist:  list
    :param definition:  topology definition
    :type definition:  Definition
    :param ligand_path:  path to ligand (may be None)
    :type ligand_path:  str
    :return: (biomolecule object, definition object--revised if ligand was
        parsed, ligand object--may be None)
    :rtype: (Biomolecule, Definition, Ligand)
    """
    if ligand_path is not None:
        ligand = Mol2Molecule()
        with open(ligand_path, encoding="utf-8") as ligand_file:
            ligand.read(ligand_file)
    else:
        ligand = None
    biomolecule = biomol.Biomolecule(pdblist, definition)
    _LOGGER.info(
        f"Created biomolecule object with {len(biomolecule.residues)} "
        f"residues and {len(biomolecule.atoms)} atoms."
    )
    for residue in biomolecule.residues:
        multoccupancy = False
        for atom in residue.atoms:
            if atom.alt_loc != "":
                multoccupancy = True
                txt = f"Multiple occupancies found: {atom.name} in {residue}."
                _LOGGER.warning(txt)
        if multoccupancy:
            err = (
                f"Multiple occupancies found in {residue}. At least one of "
                "the instances is being ignored."
            )
            _LOGGER.warning(err)
    return biomolecule, definition, ligand


def is_repairable(biomolecule, has_ligand):
    """Determine if the biomolecule can be (or needs to be) repaired.

    :param biomolecule:  biomolecule object
    :type biomolecule:  biomol.Biomolecule
    :param has_ligand:  does the system contain a ligand?
    :type has_ligand:  bool
    :return:  indication of whether biomolecule can be repaired
    :rtype:  bool
    :raises ValueError: if there are insufficient heavy atoms or a significant
        part of the biomolecule is missing
    """
    num_heavy = biomolecule.num_heavy
    num_missing = biomolecule.num_missing_heavy
    if num_heavy == 0:
        if not has_ligand:
            err = (
                "No biomolecule heavy atoms found and no ligand present. "
                "Unable to proceed.  You may also see this message if "
                "PDB2PQR does not have parameters for any residue in your "
                "biomolecule."
            )
            raise ValueError(err)
        else:
            err = (
                "No heavy atoms found but a ligand is present. Proceeding "
                "with caution."
            )
            _LOGGER.warning(err)
            return False
    if num_missing == 0:
        _LOGGER.info("This biomolecule is clean.  No repair needed.")
        return False
    miss_frac = float(num_missing) / float(num_heavy)
    if miss_frac > REPAIR_LIMIT:
        error = f"This PDB file is missing too many ({num_missing} out of "
        error += f"{num_heavy:d}, {miss_frac:g}) "
        error += "heavy atoms to accurately repair the file."
        error += f"The current repair limit is set at {REPAIR_LIMIT:g}. "
        error += "You may also see this message if PDB2PQR does not have "
        error += "parameters for enough residues in your biomolecule."
        _LOGGER.error(error)
        return False
    return True


def drop_water(pdblist):
    """Drop waters from a list of PDB records.

    .. todo:: this module is already too long but this function fits better
        here. Other possible place would be utilities.

    :param pdb_list:  list of PDB records as returned by io.get_molecule
    :type pdb_list:  [str]
    :return:  new list of PDB records with waters removed.
    :rtype:  [str]
    """
    pdblist_new = []
    for record in pdblist:
        record_type = record.record_type()
        if (
            record_type in ["HETATM", "ATOM", "SIGATM", "SEQADV"]
            and record.res_name in aa.WAT.water_residue_names
        ):
            continue
        pdblist_new.append(record)
    return pdblist_new


def run_propka(args, biomolecule):
    """Run a PROPKA calculation.

    :param args:  command-line arguments
    :type args:  argparse.Namespace
    :param biomolecule:  biomolecule object
    :type biomolecule:  Biomolecule
    :return:  (DataFrame-convertible table of assigned pKa values,
               pKa information from PROPKA)
    :rtype:  (list, str)
    """
    lines = io.print_biomolecule_atoms(
        atomlist=biomolecule.atoms, chainflag=args.keep_chain, pdbfile=True
    )

    with StringIO() as fpdb:
        fpdb.writelines(lines)
        parameters = pk_in.read_parameter_file(args.parameters, Parameters())
        molecule = MolecularContainer(parameters, args)
        # needs a mock name with .pdb extension to work with stream data, hence the "input.pdb"
        molecule = pk_in.read_molecule_file("input.pdb", molecule, fpdb)

    molecule.calculate_pka()

    # Extract pKa information from PROPKA
    # write_pka(
    #    self,
    #    self.version.parameters,
    #    filename=filename,
    #    conformation="AVR",
    #    reference=reference,
    # )
    lines = []
    # lines.append(f"{pk_out.get_propka_header()}")
    # lines.append(f"{pk_out.get_references_header()}")
    # lines.append(
    #     pk_out.get_determinant_section(
    #         molecule, 'AVR', molecule.version.parameters
    #     )
    # )
    # lines.append(
    #     pk_out.get_summary_section(
    #         molecule, "AVR", molecule.version.parameters
    #     )
    # )
    # lines.append(pk_out.get_the_line())
    lines.append(
        pk_out.get_folding_profile_section(
            molecule,
            conformation="AVR",
            reference="neutral",
            window=[0.0, 14.0, 1.0],
        )
    )
    lines.append(
        pk_out.get_charge_profile_section(molecule, conformation="AVR")
    )
    lines.append(pk_out.get_the_line())
    pka_str = "\n".join(lines)

    # Summarize in pKas in DataFrame for later use
    conformation = molecule.conformations["AVR"]
    rows = []
    for group in conformation.groups:
        row_dict = OrderedDict()
        atom = group.atom
        row_dict["res_num"] = atom.res_num
        row_dict["ins_code"] = atom.icode
        row_dict["res_name"] = atom.res_name
        row_dict["chain_id"] = atom.chain_id
        row_dict["group_label"] = group.label
        row_dict["group_type"] = getattr(group, "type", None)
        row_dict["pKa"] = group.pka_value
        row_dict["model_pKa"] = group.model_pka
        row_dict["buried"] = group.buried
        if group.coupled_titrating_group:
            row_dict["coupled_group"] = group.coupled_titrating_group.label
        else:
            row_dict["coupled_group"] = None
        rows.append(row_dict)
    return rows, pka_str


def non_trivial(args, biomolecule, ligand, definition, is_cif):
    """Perform a non-trivial PDB2PQR run.

    .. todo::
       These routines should be generalized to biomolecules; none of them are
       specific to biomolecules.

    :param args:  command-line arguments
    :type args:  argparse.Namespace
    :param biomolecule:  biomolecule
    :type biomolecule:  Biomolecule
    :param ligand:  ligand object or None
    :type ligand:  Mol2Molecule
    :param definition:  topology definition
    :type definition:  Definition
    :param is_cif:  indicates whether file is CIF format
    :type is_cif:  bool
    :raises ValueError:  for missing atoms that prevent debumping
    :return:  dictionary with results
    :rtype:  dict
    """
    _LOGGER.info("Loading forcefield.")
    forcefield_ = forcefield.Forcefield(
        args.ff, definition, args.userff, args.usernames
    )
    _LOGGER.info("Loading hydrogen topology definitions.")
    hydrogen_handler = hydrogens.create_handler()
    debumper = debump.Debump(biomolecule)
    pka_df = None
    if args.assign_only:
        # TODO - I don't understand why HIS needs to be set to HIP for
        # assign-only
        biomolecule.set_hip()
    else:
        if is_repairable(biomolecule, args.ligand is not None):
            _LOGGER.info(
                f"Attempting to repair {biomolecule.num_missing_heavy:d} "
                "missing atoms in biomolecule."
            )
            biomolecule.repair_heavy()
        _LOGGER.info("Updating disulfide bridges.")
        biomolecule.update_ss_bridges()
        if args.debump:
            _LOGGER.info("Debumping biomolecule.")
            try:
                debumper.debump_biomolecule()
            except ValueError as err:
                err = f"Unable to debump biomolecule. {err}"
                raise ValueError(err)
        if args.pka_method == "propka":
            _LOGGER.info("Assigning titration states with PROPKA.")
            biomolecule.remove_hydrogens()
            pka_df, pka_str = run_propka(args, biomolecule)
            _LOGGER.info(f"PROPKA information:\n{pka_str}")
            biomolecule.apply_pka_values(
                forcefield_.name,
                args.ph,
                {
                    f"{row['res_name']} {row['res_num']} {row['chain_id']}": row[
                        "pKa"
                    ]
                    for row in pka_df
                    if row["group_label"].startswith(row["res_name"])
                },
            )

        _LOGGER.info("Adding hydrogens to biomolecule.")
        biomolecule.add_hydrogens()
        if args.debump:
            _LOGGER.info("Debumping biomolecule (again).")
            debumper.debump_biomolecule()
        _LOGGER.info("Optimizing hydrogen bonds")
        hydrogen_routines = hydrogens.HydrogenRoutines(
            debumper, hydrogen_handler
        )
        if args.opt:
            hydrogen_routines.set_optimizeable_hydrogens()
            biomolecule.hold_residues(None)
            hydrogen_routines.initialize_full_optimization()
        else:
            hydrogen_routines.initialize_wat_optimization()
        hydrogen_routines.optimize_hydrogens()
        hydrogen_routines.cleanup()
    _LOGGER.info("Applying force field to biomolecule states.")
    biomolecule.set_states()
    matched_atoms, missing_atoms = biomolecule.apply_force_field(forcefield_)
    if args.ligand is not None:
        _LOGGER.info("Processing ligand.")
        _LOGGER.warning("Using ZAP9 forcefield for ligand radii.")
        ligand.assign_parameters()
        lig_atoms = []
        for residue in biomolecule.residues:
            tot_charge = 0
            for pdb_atom in residue.atoms:
                # Only check residues with HETATM
                if pdb_atom.type == "ATOM":
                    break
                try:
                    mol2_atom = ligand.atoms[pdb_atom.name]
                    pdb_atom.radius = mol2_atom.radius
                    pdb_atom.ffcharge = mol2_atom.charge
                    tot_charge += mol2_atom.charge
                    lig_atoms.append(pdb_atom)
                except KeyError:
                    err = (
                        f"Can't find HETATM {residue.name} {residue.res_seq} "
                        f"{pdb_atom.name} in MOL2 file"
                    )
                    _LOGGER.warning(err)
                    missing_atoms.append(pdb_atom)
        matched_atoms += lig_atoms
    total_charge = 0
    for residue in biomolecule.residues:
        charge = residue.charge
        charge_err = noninteger_charge(charge)
        if charge_err:
            _LOGGER.warning(
                f"Residue {residue} has non-integer charge:  {charge_err}"
            )
        total_charge += charge
    charge_err = noninteger_charge(total_charge)
    if charge_err:
        raise ValueError(charge_err)
    if args.ffout is not None:
        _LOGGER.info(f"Applying custom naming scheme ({args.ffout}).")
        if args.ffout != args.ff:
            name_scheme = forcefield.Forcefield(args.ffout, definition, None)
        else:
            name_scheme = forcefield_
        biomolecule.apply_name_scheme(name_scheme)
    _LOGGER.info("Regenerating headers.")
    reslist, charge = biomolecule.charge
    if is_cif:
        header = io.print_pqr_header_cif(
            missing_atoms,
            reslist,
            charge,
            args.ff,
            args.pka_method,
            args.ph,
            args.ffout,
            include_old_header=args.include_header,
        )
    else:
        header = io.print_pqr_header(
            biomolecule.pdblist,
            missing_atoms,
            reslist,
            charge,
            args.ff,
            args.pka_method,
            args.ph,
            args.ffout,
            include_old_header=args.include_header,
        )
    _LOGGER.info("Regenerating PDB lines.")
    lines = io.print_biomolecule_atoms(matched_atoms, args.keep_chain)
    return {
        "lines": lines,
        "header": header,
        "missed_residues": missing_atoms,
        "pka_df": pka_df,
    }


def main_driver(args: argparse.Namespace):
    """Main driver for running program from the command line.

    Validate inputs, launch PDB2PQR, handle output.

    :param args:  command-line arguments
    :type args:  argparse.Namespace
    """
    _LOGGER.debug(f"Invoked with arguments: {args}")
    print_splash_screen(args)
    _LOGGER.info("Checking and transforming input arguments.")
    args = transform_arguments(args)
    check_files(args)
    check_options(args)
    _LOGGER.info("Loading topology files.")
    definition = io.get_definitions()
    _LOGGER.info(f"Loading molecule: {args.input_path}")
    pdblist, is_cif = io.get_molecule(args.input_path)
    if args.drop_water:
        _LOGGER.info("Dropping water from structure.")
        pdblist = drop_water(pdblist)
    _LOGGER.info("Setting up molecule.")
    biomolecule, definition, ligand = setup_molecule(
        pdblist, definition, args.ligand
    )
    _LOGGER.info("Setting termini states for biomolecule chains.")
    biomolecule.set_termini(neutraln=args.neutraln, neutralc=args.neutralc)
    biomolecule.update_bonds()
    if args.clean:
        _LOGGER.info(
            "Arguments specified cleaning only; skipping remaining steps."
        )
        results = {
            "header": "",
            "missed_residues": None,
            "biomolecule": biomolecule,
            "lines": io.print_biomolecule_atoms(
                biomolecule.atoms, args.keep_chain
            ),
            "pka_df": None,
        }
    else:
        try:
            results = non_trivial(
                args=args,
                biomolecule=biomolecule,
                ligand=ligand,
                definition=definition,
                is_cif=is_cif,
            )
        except ValueError as err:
            _LOGGER.critical(err)
            _LOGGER.critical("Giving up.")
            raise RuntimeError from err
    print_pqr(
        args=args,
        pqr_lines=results["lines"],
        header_lines=results["header"],
        missing_lines=results["missed_residues"],
        is_cif=is_cif,
    )
    if args.pdb_output:
        print_pdb(
            args=args,
            pdb_lines=io.print_biomolecule_atoms(
                biomolecule.atoms, chainflag=args.keep_chain, pdbfile=True
            ),
            header_lines=results["header"],
            missing_lines=results["missed_residues"],
            is_cif=is_cif,
        )
    if args.apbs_input:
        io.dump_apbs(args.output_pqr, args.apbs_input)
    return results["missed_residues"], results["pka_df"], biomolecule


def main():
    """Hook for command-line usage."""
    parser = build_main_parser()
    args = parser.parse_args()
    io.setup_logger(args.output_pqr, args.log_level)
    if main_driver(args) == 1:
        sys.exit(1)


def run_pdb2pqr(args: Sequence[str | PathLike]):
    """Run PDB2PQR with a list of arguments.

    Logger is not set up so that it can be called multiple times.

    :param args:  list of command-line arguments
    :type args:  list
    :return:  results of PDB2PQR run
    :rtype:  tuple
    """
    args_strlist = [str(arg) for arg in args]
    parser = build_main_parser()
    args_parsed = parser.parse_args(args_strlist)
    return main_driver(args_parsed)


def dx_to_cube():
    """Convert DX file format to Cube file format.

    The OpenDX file format is defined at
    <https://www.idvbook.com/wp-content/uploads/2010/12/opendx.pdf` and the
    Cube file format is defined at
    <https://docs.chemaxon.com/display/Gaussian_Cube_format.html>.

    .. todo:: This function should be moved into the APBS code base.
    """
    desc = f"{TITLE_STR}\ndx2cube: converting OpenDX-format files to "
    desc += "Gaussian Cube format since (at least) 2015"
    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("dx_input", help="name of the dx_input file")
    parser.add_argument("pqr_input", help="name of the pqr_input file")
    parser.add_argument("output", help="name of the output file")
    parser.add_argument(
        "--log-level",
        help="set logging level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    args = parser.parse_args()
    log_level = getattr(logging, args.log_level)
    logging.basicConfig(level=log_level)
    _LOGGER.debug(f"Got arguments: {args}", args)
    _LOGGER.info(f"Reading PQR from {args.pqr_input}...")
    with open(args.pqr_input) as pqr_file:
        atom_list = io.read_pqr(pqr_file)
    _LOGGER.info(f"Reading DX from {args.dx_input}...")
    with open(args.dx_input) as dx_file:
        dx_dict = io.read_dx(dx_file)
    _LOGGER.info(f"Writing Cube to {args.output}...")
    with open(args.output, "w") as cube_file:
        io.write_cube(cube_file, dx_dict, atom_list)
