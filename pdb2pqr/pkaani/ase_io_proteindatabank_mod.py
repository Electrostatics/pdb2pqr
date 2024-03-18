# CODE REUSED FROM https://github.com/isayevlab/pKa-ANI
# PLEASE CREDIT THE Isayev Lab if you use pKa-ANI to assign titration states!
#
# CITATION:
# Gokcan, H.; Isayev, O. Prediction of Protein pKa with Representation Learning.
# Chemical Science, 2022, 13, 2462–2474. https://doi.org/10.1039/d1sc05610g.

"""Module to read and write atoms in PDB file format.

See::

    http://www.wwpdb.org/documentation/file-format

Note: The PDB format saves cell lengths and angles; hence the absolute
orientation is lost when saving.  Saving and loading a file will
conserve the scaled positions, not the absolute ones.
"""

import warnings

import numpy as np

from ase.atoms import Atoms
from ase.geometry import cellpar_to_cell
from ase.io.espresso import label_to_symbol
from ase.utils import writer


def read_atom_line(line_full):
    """
    Read atom line from pdb format
    HETATM    1  H14 ORTE    0       6.301   0.693   1.919  1.00  0.00        H
    """

    line = line_full.rstrip("\n")
    type_atm = line[0:6]
    if type_atm in ("ATOM  ", "HETATM"):

        name = line[12:16].strip()

        altloc = line[16]
        # HATICE
        # check pdb file format
        # columns 18-20 is residue name
        # column 22 chain identifier
        # resname = line[17:21]
        resname = line[17:20]
        chainid = line[21]
        # HATICE

        # chainid = line[21]        # Not used
        # resseq = int(line[22:26].split()[0])  # sequence identifier
        # in some cases there is an insertion mutation
        # and in some cases an insertion code is used
        # to keep the residue numbering shceme same as in WT
        # i.e. 1GOA.pdb GLN80 GLY80B TRP81
        # if we get 22:26 we dont get the insertion code
        # this results to have same residue number which causes errors
        resseq = line[22:27].split()[0]  # int(line[22:26].split()[0])
        # icode = line[26]          # insertion code, not used

        # atomic coordinates
        try:
            coord = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=np.float64,
            )
        except ValueError:
            raise ValueError("Invalid or missing coordinate(s)")

        # occupancy & B factor
        try:
            occupancy = float(line[54:60])
        except ValueError:
            occupancy = None  # Rather than arbitrary zero or one

        if occupancy is not None and occupancy < 0:
            warnings.warn("Negative occupancy in one or more atoms")

        try:
            bfactor = float(line[60:66])
        except ValueError:
            # The PDB use a default of zero if the data is missing
            bfactor = 0.0

        # segid = line[72:76] # not used
        symbol = line[76:78].strip().upper()

    else:
        raise ValueError("Only ATOM and HETATM supported")

    # HATICE
    # return symbol, name, altloc, resname, coord, occupancy, bfactor, resseq
    return (
        symbol,
        name,
        altloc,
        resname,
        coord,
        occupancy,
        bfactor,
        resseq,
        chainid,
        type_atm,
    )
    # HATICE


def read_proteindatabank(fileobj, index=-1, read_arrays=True):
    """Read PDB files."""
    images = []
    orig = np.identity(3)
    trans = np.zeros(3)
    occ = []
    bfactor = []
    residuenames = []
    residuenumbers = []
    atomtypes = []

    symbols = []
    positions = []
    chain_id = []  # HATICE
    type_atm = []
    cell = None
    pbc = None

    def build_atoms():
        atoms = Atoms(symbols=symbols, cell=cell, pbc=pbc, positions=positions)

        if not read_arrays:
            return atoms
        # START: HATICE
        # info = {'occupancy': occ,
        #        'bfactor': bfactor,
        #        'residuenames': residuenames,
        #        'atomtypes': atomtypes,
        #        'residuenumbers': residuenumbers}

        info = {
            "occupancy": occ,
            "bfactor": bfactor,
            "residuenames": residuenames,
            "atomtypes": atomtypes,
            "residuenumbers": residuenumbers,
            "chainid": chain_id,
            "type_atm": type_atm,
        }

        # END: HATICE
        for name, array in info.items():
            if len(array) == 0:
                pass
            elif len(array) != len(atoms):
                warnings.warn(
                    "Length of {} array, {}, "
                    "different from number of atoms {}".format(
                        name, len(array), len(atoms)
                    )
                )
            else:
                atoms.set_array(name, np.array(array))
        return atoms

    for line in fileobj.readlines():
        if line.startswith("CRYST1"):
            cellpar = [
                float(line[6:15]),  # a
                float(line[15:24]),  # b
                float(line[24:33]),  # c
                float(line[33:40]),  # alpha
                float(line[40:47]),  # beta
                float(line[47:54]),
            ]  # gamma
            cell = cellpar_to_cell(cellpar)
            pbc = True
        for c in range(3):
            if line.startswith("ORIGX" + "123"[c]):
                orig[c] = [
                    float(line[10:20]),
                    float(line[20:30]),
                    float(line[30:40]),
                ]
                trans[c] = float(line[45:55])

        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Atom name is arbitrary and does not necessarily
            # contain the element symbol.  The specification
            # requires the element symbol to be in columns 77+78.
            # Fall back to Atom name for files that do not follow
            # the spec, e.g. packmol.

            # line_info = symbol, name, altloc, resname, coord, occupancy,
            #             bfactor, resseq
            line_info = read_atom_line(line)

            try:
                symbol = label_to_symbol(line_info[0])
            except (KeyError, IndexError):
                symbol = label_to_symbol(line_info[1])

            position = np.dot(orig, line_info[4]) + trans
            atomtypes.append(line_info[1])
            residuenames.append(line_info[3])
            if line_info[5] is not None:
                occ.append(line_info[5])
            bfactor.append(line_info[6])
            residuenumbers.append(line_info[7])

            symbols.append(symbol)
            positions.append(position)

            # START: HATICE
            # print(line_info)
            # ('O', 'O', ' ', 'ASP', array([24.933, 46.981, 51.939]), 1.0, 0.0, 12, 'B')
            chain_id.append(line_info[8])
            type_atm.append(line_info[9])
            # END: HATICE
        if line.startswith("END"):
            # End of configuration reached
            # According to the latest PDB file format (v3.30),
            # this line should start with 'ENDMDL' (not 'END'),
            # but in this way PDB trajectories from e.g. CP2K
            # are supported (also VMD supports this format).
            atoms = build_atoms()
            images.append(atoms)
            occ = []
            bfactor = []
            residuenames = []
            atomtypes = []
            symbols = []
            positions = []
            cell = None
            pbc = None

    if len(images) == 0:
        atoms = build_atoms()
        images.append(atoms)
    return images[index]


@writer
def write_proteindatabank(fileobj, images, write_arrays=True):
    """Write images to PDB-file."""
    if hasattr(images, "get_positions"):
        images = [images]

    rotation = None
    if images[0].get_pbc().any():
        from ase.geometry import cell_to_cellpar, cellpar_to_cell

        currentcell = images[0].get_cell()
        cellpar = cell_to_cellpar(currentcell)
        exportedcell = cellpar_to_cell(cellpar)
        rotation = np.linalg.solve(currentcell, exportedcell)
        # ignoring Z-value, using P1 since we have all atoms defined explicitly
        format = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1\n"
        fileobj.write(
            format
            % (
                cellpar[0],
                cellpar[1],
                cellpar[2],
                cellpar[3],
                cellpar[4],
                cellpar[5],
            )
        )

    #     1234567 123 6789012345678901   89   67   456789012345678901234567 890
    format = (
        "ATOM  %5d %4s MOL     1    %8.3f%8.3f%8.3f%6.2f%6.2f"
        "          %2s  \n"
    )

    # RasMol complains if the atom index exceeds 100000. There might
    # be a limit of 5 digit numbers in this field.
    MAXNUM = 100000

    symbols = images[0].get_chemical_symbols()
    natoms = len(symbols)

    for n, atoms in enumerate(images):
        fileobj.write("MODEL     " + str(n + 1) + "\n")
        p = atoms.get_positions()
        occupancy = np.ones(len(atoms))
        bfactor = np.zeros(len(atoms))
        if write_arrays:
            if "occupancy" in atoms.arrays:
                occupancy = atoms.get_array("occupancy")
            if "bfactor" in atoms.arrays:
                bfactor = atoms.get_array("bfactor")
        if rotation is not None:
            p = p.dot(rotation)
        for a in range(natoms):
            x, y, z = p[a]
            occ = occupancy[a]
            bf = bfactor[a]
            fileobj.write(
                format
                % (
                    (a + 1) % MAXNUM,
                    symbols[a],
                    x,
                    y,
                    z,
                    occ,
                    bf,
                    symbols[a].upper(),
                )
            )
        fileobj.write("ENDMDL\n")
