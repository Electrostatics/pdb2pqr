"""Implements the PEOE method.

The PEOE method is described in:  Paul Czodrowski, Ingo Dramburg,
Christoph A. Sotriffer  Gerhard Klebe. Development, validation, and
application of adapted PEOE charges to estimate pKa values of functional
groups in protein-ligand complexes. Proteins, 65, 424-437, 2006.
https://doi.org/10.1002/prot.21110
"""

import logging
from math import isclose

_LOGGER = logging.getLogger(__name__)


# The terms of the third-order polynomial fit for the electronegativity.
# See https://doi.org/10.1002/prot.21110 for more information.
# NOTE - this data has no meaning outside of this module; do not move.
POLY_TERMS = {
    "H": (7.17, 6.24, -0.56, 12.85),
    "C.3": (7.98, 9.18, 1.88, 19.04),
    "C.CAT": (7.98, 9.18, 1.88, 19.04),
    "C.2": (8.79 + 0.5, 9.32, 1.51, 19.62),
    "C.AR": (7.98 + 0.55, 9.18, 1.88, 19.04),
    "C.1": (10.39, 9.45, 0.73, 20.57),
    "N.3": (11.54 + 6.0, 10.28, 1.36, 28.00),
    "N.4": (11.54 + 6.0, 10.28, 1.36, 28.00),
    "N.AR": (12.87 - 1.29, 11.15, 0.85, 24.87),
    "N.2": (12.87, 11.15, 0.85, 24.87),
    "N.PL3": (12.87 + 0.5, 11.15, 0.85, 24.87),
    "N.AM": (12.87 + 3.5, 11.15, 0.85, 24.87),
    "N.1": (15.68, 11.70, -0.27, 27.11),
    "O.OH": (14.18 + 0.8, 12.92, 1.39, 28.49),
    "O.3": (14.18 - 3.1, 12.92, 1.39, 28.49),
    "O.2": (14.18, 12.92, 1.39, 28.49),
    "O.CO2": (15.25, 13.79, 0.47, 31.33),
    "F": (12.36, 13.85, 2.31, 30.82),
    "CL": (9.38 + 1.0, 9.69, 1.35, 22.04),
    "BR": (10.08 + 0.8, 8.47, 1.16, 19.71),
    "I": (9.90 + 1.0, 7.96, 0.96, 18.82),
    "S.3": (10.13 + 0.5, 9.13, 1.38, 20.65),
    "S.2": (10.13 + 0.5, 9.13, 1.38, 20.65),
    "S.O2": (10.13 + 0.5, 9.13, 1.38, 20.65),
    "P.3": (10.13 + 0.5, 9.13, 1.38, 20.65),
}
# Maximum (absolute) value of charge after which contribution to polynomial
# is capped
MAX_CHARGE = 1.1
DEFAULT_H_ELECTRONEG = 20.02
DEFAULT_H_CHARGE = 1.0
# These next values are from the "Adaptation of the PEOE Procedure" section of
# https://doi.org/10.1002/prot.21110.
DAMPING_FACTOR = 0.778
SCALING_FACTOR = 1.56
NUM_CYCLES = 6


def electronegativity(charge, poly_terms, atom_type):
    """Calculate the electronegativity.

    Calculation is based on a third-order polynomial in the atomic charge as
    described in Equation 2 of https://doi.org/10.1002/prot.21110.

    :param charge:  charge of atom
    :type charge:  float
    :param poly_terms:  polynomial terms ordered from 0th- to 3rd-order
    :param atom_type:  string with atom type
    :type atom_type:  str
    :return: electronegativity value
    :rtype:  float
    :raises IndexError:  if incorrect number of poly_terms given
    """
    chi = None
    if abs(charge) > MAX_CHARGE:
        charge = -1.0 * MAX_CHARGE if charge < 0 else MAX_CHARGE
    if (atom_type == "H") and isclose(charge, DEFAULT_H_CHARGE):
        chi = DEFAULT_H_ELECTRONEG
    else:
        if len(poly_terms) == 4:
            chi = (
                poly_terms[0]
                + poly_terms[1] * charge
                + poly_terms[2] * charge * charge
                + poly_terms[3] * charge * charge * charge
            )
        elif len(poly_terms) == 3:
            chi = (
                poly_terms[0]
                + poly_terms[1] * charge
                + poly_terms[2] * charge * charge
            )
        else:
            err = f"Cannot parse length-{len(poly_terms):d} polynomial"
            raise IndexError(err)
    return chi


def assign_terms(atoms, term_dict):
    """Assign polynomial terms to each atom.

    :param atoms:  list of Mol2Atom atoms
    :type atoms:  list
    :param term_dict:  dictionary of polynomial terms
    :type term_dict:  dict
    :return:  modified list of atoms
    :rtype:  list
    """
    for atom in atoms:
        atom_type = atom.type.upper()
        if atom_type == "O.3":
            atom_type = "O.OH"
        try:
            atom.poly_terms = term_dict[atom_type]
        except KeyError as e:
            raise KeyError(
                f"Unable to find polynomial terms for atom type {atom_type}"
            ) from e
    return atoms


def equilibrate(
    atoms,
    damp=DAMPING_FACTOR,
    scale=SCALING_FACTOR,
    num_cycles=NUM_CYCLES,
    term_dict=POLY_TERMS,
):
    """Equilibrate the atomic charges.

    :param atoms:  list of Mol2Atom atoms to equilibrate
    :type atoms:  list
    :param damp:  damping factor for equilibration process
    :type damp:  float
    :param scale:  scaling factor for equilibration process
    :type scale:  float
    :param num_cycles:  number of PEOE cycles
    :type num_cycles:  int
    :param term_dict:  dictionary of polynomial terms
    :type term_dict:  dict
    :return: revised list of atoms
    :rtype:  list
    """
    atoms = assign_terms(atoms, term_dict)
    # Reset or accumulate charges
    abs_qges = 0.0
    for atom in atoms:
        if isclose(atom.charge, 0.0):
            atom.equil_formal_charge = 0.0
        else:
            # PEOE multiples all atoms by a scaling factor at the end to
            # account for increased polarizability.  The initial formal
            # charge needs to be reduced to account for this scaling.
            atom.equil_formal_charge = atom.charge * (1.0 / scale)
            abs_qges += abs(atom.charge)
        atom.charge = 0
    # A finite number of cycles is used to prevent complete equilibration of
    # the molecule.  I'm not sure why this is a good idea but people have
    # been doing it since the original 1978 Tetrahedron paper with Gasteiger
    # & Marsili
    for icycle in range(num_cycles):
        for atom1 in atoms:
            chi1 = electronegativity(
                atom1.charge, atom1.poly_terms, atom1.type
            )
            atom1.delta_charge = 0.0
            for atom2 in atom1.bonded_atoms:
                chi2 = electronegativity(
                    atom2.charge, atom2.poly_terms, atom2.type
                )
                chi_diff = chi2 - chi1
                if chi2 > chi1:
                    chi_norm = electronegativity(
                        +1, atom1.poly_terms, atom1.type
                    )
                else:
                    chi_norm = electronegativity(
                        +1, atom2.poly_terms, atom2.type
                    )
                # Damping is used in PEOE to accelerate convergence
                atom1.delta_charge += (chi_diff / chi_norm) * (
                    damp ** (icycle + 1)
                )
        for atom in atoms:
            if isclose(abs_qges, 0.0):
                atom.charge += atom.delta_charge
            else:
                atom.charge += (
                    atom.delta_charge
                    + (1.0 / num_cycles) * atom.equil_formal_charge
                )
    for atom in atoms:
        atom.charge = scale * atom.charge
    return atoms
