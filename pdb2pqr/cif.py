"""CIF parsing methods.

This methods use the pdbx/cif parser provided by WWPDB
(http://mmcif.wwpdb.org/docs/sw-examples/python/html/index.html)

.. todo:  Why do we have this module when we have pdbx?
.. codeauthor::  Juan Brandi
"""

import logging
from datetime import datetime

import pdbx
from numpy import ceil, minimum

from . import pdb

_LOGGER = logging.getLogger(__name__)


def atom_site(block):
    """Handle ATOM_SITE block.

    Data items in the ATOM_SITE category record details about the atom sites
    in a macromolecular crystal structure, such as the positional coordinates,
    atomic displacement parameters, magnetic moments and directions.
    (Source: https://j.mp/2Zprx41)

    :param block:  PDBx data block
    :type block:  [str]
    :return:  (array of pdb.ATOM objects, array of things that weren't handled
        by parser)
    :rtype:  ([Atom], [str])
    """
    line = 0
    pdb_arr = []
    err_arr = []
    atoms = block.get_object("atom_site")
    num_model_arr = count_models(block)
    if len(num_model_arr) == 1:
        # TODO - this part of the conditional should be a separate function
        for i in range(atoms.row_count):
            if atoms.get_value("group_PDB", i) == "ATOM":
                try:
                    line = ""
                    # 1  - 6 RECORD NAME (ATOM)
                    line += atoms.get_value("group_PDB", i) + " " * (
                        6 - len(atoms.get_value("group_PDB", i))
                    )
                    # 7  - 11 ATOM SERIAL
                    line += " " * (
                        5 - len(str(atoms.get_value("id", i)))
                    ) + str(atoms.get_value("id", i))
                    # 12 - 13
                    line += "  "
                    # 14 - 16 ATOM NAME
                    line += atoms.get_value("label_atom_id", i) + " " * (
                        3 - len(atoms.get_value("label_atom_id", i))
                    )
                    # 17 ALT LOCATION
                    if atoms.get_value("label_alt_id", i) == ".":
                        line += " "
                    else:
                        atoms.get_value("label_alt_id", i)
                    # 18 - 20 RES NAME
                    line += " " * (
                        3 - len(atoms.get_value("label_comp_id", i))
                    ) + atoms.get_value("label_comp_id", i)
                    # 21
                    line += " "
                    # 22 CHAIN ID
                    line += " " * (
                        1 - len(atoms.get_value("label_asym_id", i))
                    ) + atoms.get_value("label_asym_id", i)
                    # 23 - 26 RES SEQ ID
                    line += " " * (
                        4 - len(str(atoms.get_value("auth_seq_id", i)))
                    ) + str(atoms.get_value("auth_seq_id", i))
                    # 27 - 30
                    line += " " * 3
                    # 31 - 38 X Coords
                    line += " " * (
                        8 - len(str(atoms.get_value("Cartn_x", i)))
                    ) + str(atoms.get_value("Cartn_x", i))
                    # 39 - 46 Y Coords
                    line += " " * (
                        8 - len(str(atoms.get_value("Cartn_y", i)))
                    ) + str(atoms.get_value("Cartn_y", i))
                    # 47 - 54 Z Coords
                    line += " " * (
                        8 - len(str(atoms.get_value("Cartn_z", i)))
                    ) + str(atoms.get_value("Cartn_z", i))
                    # 55 - 60 OCCUPANCY
                    line += " " * (
                        6 - len(str(atoms.get_value("occupancy", i)))
                    ) + str(atoms.get_value("occupancy", i))
                    # 61 - 66 TEMP FACTOR
                    line += " " * (
                        6 - len(str(atoms.get_value("B_iso_or_equiv", i)))
                    ) + str(atoms.get_value("B_iso_or_equiv", i))
                    # 67 - 76
                    line += " " * 10
                    # 77 - 78 ELEMENT SYMBOL
                    line += " " * (
                        2 - len(atoms.get_value("type_symbol", i))
                    ) + atoms.get_value("type_symbol", i)
                    # 79 - 80 CHARGE OF ATOM
                    if atoms.get_value("pdbx_formal_charge", i) == "?":
                        line += " " * 2
                    else:
                        atoms.get_value("pdbx_formal_charge", i)
                    pdb_arr.append(pdb.ATOM(line))
                except KeyError:
                    _LOGGER.error(f"atom_site: Error reading line: #{line}#\n")
            elif atoms.get_value("group_PDB", i) == "HETATM":
                try:
                    line = ""
                    # 1  - 6 RECORD NAME (HETATM)
                    line += atoms.get_value("group_PDB", i) + "" * (
                        6 - len(atoms.get_value("group_PDB", i))
                    )
                    # 7  - 11 ATOM SERIAL
                    line += " " * (
                        5 - len(str(atoms.get_value("id", i)))
                    ) + str(atoms.get_value("id", i))
                    # 12 - 13
                    line += "  "
                    # 14 - 16 ATOM NAME
                    line += atoms.get_value("label_atom_id", i) + " " * (
                        3 - len(atoms.get_value("label_atom_id", i))
                    )
                    # 17 ALT LOCATION
                    if atoms.get_value("label_alt_id", i) == ".":
                        line += " "
                    else:
                        atoms.get_value("label_alt_id", i)
                    # 18 - 20 RES NAME
                    line += " " * (
                        3 - len(atoms.get_value("label_comp_id", i))
                    ) + atoms.get_value("label_comp_id", i)
                    # 21
                    line += " "
                    # 22 CHAIN ID
                    line += " " * (
                        1 - len(atoms.get_value("label_asym_id", i))
                    ) + atoms.get_value("label_asym_id", i)
                    # 23 - 26 RES SEQ ID
                    line += " " * (
                        4 - len(str(atoms.get_value("auth_seq_id", i)))
                    ) + str(atoms.get_value("auth_seq_id", i))
                    # 27 - 30
                    line += " " * 3
                    # 31 - 38 X Coords
                    line += " " * (
                        8 - len(str(atoms.get_value("Cartn_x", i)))
                    ) + str(atoms.get_value("Cartn_x", i))
                    # 39 - 46 Y Coords
                    line += " " * (
                        8 - len(str(atoms.get_value("Cartn_y", i)))
                    ) + str(atoms.get_value("Cartn_y", i))
                    # 47 - 54 Z Coords
                    line += " " * (
                        8 - len(str(atoms.get_value("Cartn_z", i)))
                    ) + str(atoms.get_value("Cartn_z", i))
                    # 55 - 60 OCCUPANCY
                    line += " " * (
                        6 - len(str(atoms.get_value("occupancy", i)))
                    ) + str(atoms.get_value("occupancy", i))
                    # 61 - 66 TEMP FACTOR
                    line += " " * (
                        6 - len(str(atoms.get_value("B_iso_or_equiv", i)))
                    ) + str(atoms.get_value("B_iso_or_equiv", i))
                    # 67 - 76
                    line += " " * (10)
                    # 77 - 78 ELEMENT SYMBOL
                    line += " " * (
                        2 - len(atoms.get_value("type_symbol", i))
                    ) + atoms.get_value("type_symbol", i)
                    # 79 - 80 CHARGE OF ATOM
                    if atoms.get_value("pdbx_formal_charge", i) == "?":
                        line += " " * 2
                    else:
                        atoms.get_value("pdbx_formal_charge", i)
                    pdb_arr.append(pdb.HETATM(line))
                except KeyError:
                    _LOGGER.error(f"atom_site: Error reading line:\n{line}")
        return pdb_arr, err_arr
    # TODO - Given the return statement above, is this "else" ever reached?
    else:
        # TODO - this part of the conditional should be a separate function
        for j in num_model_arr:
            try:
                line = "MODEL "
                line += " " * 4
                line += " " * (4 - len(str(j))) + str(j)
                pdb_arr.append(pdb.MODEL(line))
            except ValueError:
                _LOGGER.error(f"atom_site: Error readline line:\n{line}")
                err_arr.append("MODEL")

            for i in range(atoms.row_count):
                if atoms.get_value("pdbx_PDB_model_num", i) == j:
                    if atoms.get_value("group_PDB", i) == "ATOM":
                        try:
                            line = ""
                            # 1  - 6 RECORD NAME (ATOM)
                            line += atoms.get_value("group_PDB", i) + " " * (
                                6 - len(atoms.get_value("group_PDB", i))
                            )
                            # 7  - 11 ATOM SERIAL
                            line += " " * (
                                5 - len(str(atoms.get_value("id", i)))
                            ) + str(atoms.get_value("id", i))
                            # 12 - 13
                            line += "  "
                            # 14 - 16 ATOM NAME
                            line += atoms.get_value(
                                "label_atom_id", i
                            ) + " " * (
                                3 - len(atoms.get_value("label_atom_id", i))
                            )
                            # 17 ALT LOCATION
                            if atoms.get_value("label_alt_id", i) == ".":
                                line += " "
                            else:
                                atoms.get_value("label_alt_id", i)
                            # 18 - 20 RES NAME
                            line += " " * (
                                3 - len(atoms.get_value("label_comp_id", i))
                            ) + atoms.get_value("label_comp_id", i)
                            # 21
                            line += " "
                            # 22 CHAIN ID
                            line += " " * (
                                1 - len(atoms.get_value("label_asym_id", i))
                            ) + atoms.get_value("label_asym_id", i)
                            # 23 - 26 RES SEQ ID
                            line += " " * (
                                4 - len(str(atoms.get_value("auth_seq_id", i)))
                            ) + str(atoms.get_value("auth_seq_id", i))
                            # 27 - 30
                            line += " " * 3
                            # 31 - 38 X Coords
                            line += " " * (
                                8 - len(str(atoms.get_value("Cartn_x", i)))
                            ) + str(atoms.get_value("Cartn_x", i))
                            # 39 - 46 Y Coords
                            line += " " * (
                                8 - len(str(atoms.get_value("Cartn_y", i)))
                            ) + str(atoms.get_value("Cartn_y", i))
                            # 47 - 54 Z Coords
                            line += " " * (
                                8 - len(str(atoms.get_value("Cartn_z", i)))
                            ) + str(atoms.get_value("Cartn_z", i))
                            # 55 - 60 OCCUPANCY
                            line += " " * (
                                6 - len(str(atoms.get_value("occupancy", i)))
                            ) + str(atoms.get_value("occupancy", i))
                            # 61 - 66 TEMP FACTOR
                            line += " " * (
                                6
                                - len(
                                    str(atoms.get_value("B_iso_or_equiv", i))
                                )
                            ) + str(atoms.get_value("B_iso_or_equiv", i))
                            # 67 - 76
                            line += " " * 10
                            # 77 - 78 ELEMENT SYMBOL
                            line += " " * (
                                2 - len(atoms.get_value("type_symbol", i))
                            ) + atoms.get_value("type_symbol", i)
                            # 79 - 80 CHARGE OF ATOM
                            if atoms.get_value("pdbx_formal_charge", i) == "?":
                                line += " " * 2
                            else:
                                atoms.get_value("pdbx_formal_charge", i)
                            pdb_arr.append(pdb.ATOM(line))
                        except KeyError:
                            _LOGGER.error(
                                f"atom_site: Error reading line:\n{line}"
                            )
                            err_arr.append("ATOM")
                    elif atoms.get_value("group_PDB", i) == "HETATM":
                        try:
                            line = ""
                            # 1  - 6 RECORD NAME (HETATM)
                            line += atoms.get_value("group_PDB", i) + "" * (
                                6 - len(atoms.get_value("group_PDB", i))
                            )
                            # 7  - 11 ATOM SERIAL
                            line += " " * (
                                5 - len(str(atoms.get_value("id", i)))
                            ) + str(atoms.get_value("id", i))
                            # 12 - 13
                            line += "  "
                            # 14 - 16 ATOM NAME
                            line += atoms.get_value(
                                "label_atom_id", i
                            ) + " " * (
                                3 - len(atoms.get_value("label_atom_id", i))
                            )
                            # 17      ALT LOCATION
                            if atoms.get_value("label_alt_id", i) == ".":
                                line += " "
                            else:
                                atoms.get_value("label_alt_id", i)
                            # 18 - 20 RES NAME
                            line += " " * (
                                3 - len(atoms.get_value("label_comp_id", i))
                            ) + atoms.get_value("label_comp_id", i)
                            # 21
                            line += " "
                            # 22      CHAIN ID
                            line += " " * (
                                1 - len(atoms.get_value("label_asym_id", i))
                            ) + atoms.get_value("label_asym_id", i)
                            # 23 - 26 RES SEQ ID
                            line += " " * (
                                4 - len(str(atoms.get_value("auth_seq_id", i)))
                            ) + str(atoms.get_value("auth_seq_id", i))
                            # 27 - 30
                            line += " " * 3
                            # 31 - 38 X Coords
                            line += " " * (
                                8 - len(str(atoms.get_value("Cartn_x", i)))
                            ) + str(atoms.get_value("Cartn_x", i))
                            # 39 - 46 Y Coords
                            line += " " * (
                                8 - len(str(atoms.get_value("Cartn_y", i)))
                            ) + str(atoms.get_value("Cartn_y", i))
                            # 47 - 54 Z Coords
                            line += " " * (
                                8 - len(str(atoms.get_value("Cartn_z", i)))
                            ) + str(atoms.get_value("Cartn_z", i))
                            # 55 - 60 OCCUPANCY
                            line += " " * (
                                6 - len(str(atoms.get_value("occupancy", i)))
                            ) + str(atoms.get_value("occupancy", i))
                            # 61 - 66 TEMP FACTOR
                            line += " " * (
                                6
                                - len(
                                    str(atoms.get_value("B_iso_or_equiv", i))
                                )
                            ) + str(atoms.get_value("B_iso_or_equiv", i))
                            # 67 - 76
                            line += " " * 10
                            # 77 - 78 ELEMENT SYMBOL
                            line += " " * (
                                2 - len(atoms.get_value("type_symbol", i))
                            ) + atoms.get_value("type_symbol", i)
                            # 79 - 80 CHARGE OF ATOM
                            if atoms.get_value("pdbx_formal_charge", i) == "?":
                                line += " " * 2
                            else:
                                atoms.get_value("pdbx_formal_charge", i)
                            pdb_arr.append(pdb.HETATM(line))
                        except KeyError:
                            _LOGGER.error(
                                f"atom_site: Error reading line:\n{line}"
                            )
                            err_arr.append("HETATOM")
            try:
                line = "ENDMDL"
                pdb_arr.append(pdb.ENDMDL(line))
            except KeyError:
                _LOGGER.error(f"atom_site: Error reading line:\n{line}")
                err_arr.append("ENDMDL")
        return pdb_arr, err_arr


def conect(block):
    """Handle CONECT block.

    Data items in the STRUCT_CONN category record details about the connections
    between portions of the structure.
    These can be hydrogen bonds, salt bridges, disulfide bridges and so on.

    The ``STRUCT_CONN_TYPE`` records define the criteria used to identify these
    connections.
    (Source: https://j.mp/3gPkJT5)

    :param block:  PDBx data block
    :type block:  [str]
    :return:  (array of pdb.conect objects, array of things that did not parse)
    :rtype:  ([pdb.CONECT], [str])
    """
    pdb_arr = []
    err_arr = []
    struct_conn = block.get_object("struct_conn")
    atoms = block.get_object("atom_site")
    if struct_conn is None or atoms is None:
        return pdb_arr, err_arr
    for index in range(struct_conn.row_count):
        atom_pair = []
        for partner in ["ptnr1_", "ptnr2_"]:
            # Retrieve all the information necessary to uniquely identify atom4
            atom_dict = {
                "auth_seq_id": struct_conn.get_value(
                    partner + "auth_seq_id", index
                ),
                "auth_comp_id": struct_conn.get_value(
                    partner + "auth_comp_id", index
                ),
                "auth_asym_id": struct_conn.get_value(
                    partner + "auth_asym_id", index
                ),
                "label_atom_id": struct_conn.get_value(
                    partner + "label_atom_id", index
                ),
            }
            for i in range(atoms.row_count):
                found = all(
                    atoms.get_value(key, i) == atom_dict[key]
                    for key in atom_dict
                )
                if found:
                    atom_pair.append(atoms.get_value("id", i))
        if len(atom_pair) == 2:
            pline = (
                "CONECT"
                + " " * (5 - len(str(atom_pair[0])))
                + str(atom_pair[0])
                + " " * (5 - len(str(atom_pair[1])))
                + str(atom_pair[1])
            )
            try:
                pdb_arr.append(pdb.CONECT(pline))
            except KeyError:
                _LOGGER.error(f"conect:   Error parsing line: \n{pline}")
                err_arr.append("conect")
    return pdb_arr, err_arr


def header(block):
    """Handle HEADER block.

    :param block:  PDBx data block
    :type block:  [str]
    :return:  (array of pdb.conect objects, array of things that did not parse)
    :rtype:  ([pdb.CONECT], [str])
    """
    header_arr = []
    header_err = []

    struct_obj = block.get_object("struct_keywords")
    database_obj = block.get_object("pdbx_database_status")
    entry_obj = block.get_object("entry")
    ridd = database_obj.get_value("recvd_initial_deposition_date")
    if len(ridd) > 9:
        ridd = datetime.strptime(ridd, "%Y-%m-%d").strftime("%d-%b-%y").upper()
    line = "HEADER"
    line += " " * 4
    line += struct_obj.get_value("pdbx_keywords") + " " * (
        40 - len(struct_obj.get_value("pdbx_keywords"))
    )
    line += " " * (9 - len(ridd)) + ridd
    line += " " * 3
    line += " " * (4 - len(entry_obj.get_value("id"))) + entry_obj.get_value(
        "id"
    )
    try:
        header_arr.append(pdb.HEADER(line))
    except KeyError:
        _LOGGER.error(f"header:   Error parsing line: #{line}#")
        header_err.append("header")
    return header_arr, header_err


def title(block):
    """Handle TITLE block.

    :param block:  PDBx data block
    :type block:  [str]
    :return:  (array of pdb.conect objects, array of things that did not parse)
    :rtype:  ([pdb.CONECT], [str])
    """
    title_arr = []
    title_err = []
    struct_obj = block.get_object("struct")
    title_string = struct_obj.get_value("title")
    title_chunk = int(ceil(len(title_string) / 70.0))
    for i in range(title_chunk):
        line = "TITLE  "
        line += " " * (2 - len(str(i + 1))) + str(i + 1) if i > 0 else "  "
        line += title_string[
            (i * 70) : minimum(len(title_string), (i + 1) * 70)
        ]
        try:
            title_arr.append(pdb.TITLE(line))
        except KeyError:
            _LOGGER.error(f"TITLE:    Error parsing line:\n{line}")
            title_err.append("title")
    return title_arr, title_err


def compnd(block):
    """Handle COMPND block.

    :param block:  PDBx data block
    :type block:  [str]
    :return:  (array of pdb.conect objects, array of things that did not parse)
    :rtype:  ([pdb.CONECT], [str])
    """
    compnd_arr = []
    compnd_err = []
    entity_obj = block.get_object("entity")
    cont = 1
    for i in range(entity_obj.row_count):
        line1 = "COMPND "
        line1 += " " * (3 - len(str(cont))) + str(cont) if cont > 1 else "   "
        line1 += "MOL_ID: " + str(entity_obj.get_value("id", i)) + ""
        try:
            compnd_arr.append(pdb.COMPND(line1))
        except KeyError:
            _LOGGER.error(f"compnd:    Error parsing line:\n{line1}")
            compnd_err.append("compnd")
        cont += 1
        line2 = "COMPND "
        line2 += " " * (3 - len(str(cont))) + str(cont) if cont > 1 else "   "
        line2 += (
            "MOLECULE: " + entity_obj.get_value("pdbx_description", i) + ""
        )
        try:
            compnd_arr.append(pdb.COMPND(line2))
        except KeyError:
            _LOGGER.error(f"compnd:    Error parsing line:\n{line2}")
            compnd_err.append("compnd")
        cont += 1
    return compnd_arr, compnd_err


def source(block):
    """Handle SOURCE block.

    :param block:  PDBx data block
    :type block:  [str]
    :return:  (array of pdb.conect objects, array of things that did not parse)
    :rtype:  ([pdb.CONECT], [str])
    """
    src_arr = []
    src_err = []
    src_obj = block.get_object("entity_src_gen")
    if src_obj is None:
        return src_arr, src_err
    cont = 1
    for i in range(src_obj.row_count):
        if src_obj.get_value("entity_id", 0) != "?":
            line = "SOURCE "
            line += (
                " " * (3 - len(str(cont))) + str(cont) if cont > 1 else "   "
            )
            line += "MOL_ID: " + str(src_obj.get_value("entity_id", i)) + ""
            cont += 1
            try:
                src_arr.append(pdb.SOURCE(line))
            except KeyError:
                _LOGGER.error(f"source:  Error parsing line:\n{line}")
                src_err.append("source")
        if src_obj.get_value("pdbx_gene_src_scientific_name", i) != "?":
            line = "SOURCE "
            line += (
                " " * (3 - len(str(cont))) + str(cont) if cont > 1 else "   "
            )
            line += (
                "ORGANISM_SCIENTIFIC: "
                + src_obj.get_value("pdbx_gene_src_scientific_name", i)
                + ""
            )
            cont += 1
            try:
                src_arr.append(pdb.SOURCE(line))
            except KeyError:
                _LOGGER.error(f"source:  Error parsing line:\n{line}")
                src_err.append("source")
        if src_obj.get_value("gene_src_common_name", i) != "?":
            line = "SOURCE "
            line += (
                " " * (3 - len(str(cont))) + str(cont) if cont > 1 else "   "
            )
            line += (
                "ORGANISM_COMMON: "
                + src_obj.get_value("gene_src_common_name", i)
                + ""
            )
            cont += 1
            try:
                src_arr.append(pdb.SOURCE(line))
            except KeyError:
                _LOGGER.error(f"source:  Error parsing line:\n{line}")
                src_err.append("source")
        if src_obj.get_value("pdbx_gene_src_ncbi_taxonomy_id", i) != "?":
            line = "SOURCE "
            line += (
                " " * (3 - len(str(cont))) + str(cont) if cont > 1 else "   "
            )
            line += (
                "ORGANISM_TAXID: "
                + src_obj.get_value("pdbx_gene_src_ncbi_taxonomy_id", i)
                + ""
            )
            cont += 1
            try:
                src_arr.append(pdb.SOURCE(line))
            except KeyError:
                _LOGGER.error(f"source:    Error parsing line:\n{line}")
                src_err.append("source")
    return src_arr, src_err


def keywds(block):
    """Handle KEYWDS block.

    :param block:  PDBx data block
    :type block:  [str]
    :return:  (array of pdb.conect objects, array of things that did not parse)
    :rtype:  ([pdb.CONECT], [str])
    """
    key_arr = []
    key_err = []
    key_obj = block.get_object("struct_keywords")
    key_string = key_obj.get_value("text")
    key_chunk = int(ceil(len(key_string) / 69.0))
    for i in range(key_chunk):
        line = "KEYWDS  "
        line += " " * (2 - len(str(i + 1))) + str(i + 1) if i > 0 else "  "
        line += key_string[(i * 69) : minimum(len(key_string), (i + 1) * 69)]
        try:
            key_arr.append(pdb.KEYWDS(line))
        except KeyError:
            _LOGGER.error(f"keywds:    Error parsing line:\n{line}")
            key_err.append("keywds")
    return key_arr, key_err


def expdata(block):
    """Handle EXPDTA block.

    :param block:  PDBx data block
    :type block:  [str]
    :return:  (array of pdb.conect objects, array of things that did not parse)
    :rtype:  ([pdb.CONECT], [str])
    """
    ex_arr = []
    ex_err = []
    ex_obj = block.get_object("exptl")
    line = "EXPDTA  "
    line += "  "
    line += ex_obj.get_value("method", 0)
    try:
        ex_arr.append(pdb.EXPDTA(line))
    except KeyError:
        _LOGGER.error(f"expdata:    Error parsing line:\n{line}\n")
        ex_err.append("expdata")
    return ex_arr, ex_err


def author(block):
    """Handle AUTHOR block.

    :param block:  PDBx data block
    :type block:  [str]
    :return:  (array of pdb.conect objects, array of things that did not parse)
    :rtype:  ([pdb.CONECT], [str])
    """
    aut_arr = []
    aut_err = []
    aut_obj = block.get_object("audit_author")
    for i in range(aut_obj.row_count):
        line = "AUTHOR  "
        line += "  " * (
            2 - len(str(aut_obj.get_value("pdbx_ordinal", i)))
        ) + str(aut_obj.get_value("pdbx_ordinal", i))
        line += aut_obj.get_value("name", i)
        try:
            aut_arr.append(pdb.AUTHOR(line))
        except KeyError:
            _LOGGER.error(f"author:  Error parsing line:\n{line}")
            aut_err.append("author")
    return aut_arr, aut_err


def cryst1(block):
    """Handle CRYST1 block.

    :param block:  PDBx data block
    :type block:  [str]
    :return:  (array of pdb.conect objects, array of things that did not parse)
    :rtype:  ([pdb.CONECT], [str])
    """
    cry_arr = []
    cry_err = []
    cry_obj = block.get_object("cell")
    sym_obj = block.get_object("symmetry")
    line = "CRYST1"
    line += " " * (
        9 - len(str(cry_obj.get_value("length_a", 0)))
    ) + cry_obj.get_value("length_a", 0)
    line += " " * (
        9 - len(str(cry_obj.get_value("length_b", 0)))
    ) + cry_obj.get_value("length_b", 0)
    line += " " * (
        9 - len(str(cry_obj.get_value("length_c", 0)))
    ) + cry_obj.get_value("length_c", 0)
    line += " " * (
        7 - len(str(cry_obj.get_value("angle_alpha", 0)))
    ) + cry_obj.get_value("angle_alpha", 0)
    line += " " * (
        7 - len(str(cry_obj.get_value("angle_beta", 0)))
    ) + cry_obj.get_value("angle_beta", 0)
    line += " " * (
        7 - len(str(cry_obj.get_value("angle_gamma", 0)))
    ) + cry_obj.get_value("angle_gamma", 0)
    line += " " * (
        11 - len(str(sym_obj.get_value("space_group_name_H-M", 0)))
    ) + sym_obj.get_value("space_group_name_H-M", 0)
    line += " " * (
        4 - len(str(cry_obj.get_value("Z_PDB", 0)))
    ) + cry_obj.get_value("Z_PDB", 0)
    try:
        cry_arr.append(pdb.CRYST1(line))
    except KeyError:
        _LOGGER.error(f"cif.cryst1:    Error parsing line:\n{line}")
        cry_err.append(cryst1)
    return cry_arr, cry_err


def scalen(block):
    """Handle SCALEn block.

    :param block:  PDBx data block
    :type block:  [str]
    :return:  (array of pdb.conect objects, array of things that did not parse)
    :rtype:  ([pdb.CONECT], [str])
    """
    sc_arr = []
    sc_err = []
    sc_obj = block.get_object("atom_sites")
    scale1 = ""
    scale1 += "SCALE1    "
    scale1 += " " * (
        10 - len(str(sc_obj.get_value("fract_transf_matrix[1][1]", 0)))
    ) + str(sc_obj.get_value("fract_transf_matrix[1][1]", 0))
    scale1 += " " * (
        10 - len(str(sc_obj.get_value("fract_transf_matrix[1][2]", 0)))
    ) + str(sc_obj.get_value("fract_transf_matrix[1][2]", 0))
    scale1 += " " * (
        10 - len(str(sc_obj.get_value("fract_transf_matrix[1][3]", 0)))
    ) + str(sc_obj.get_value("fract_transf_matrix[1][3]", 0))
    scale1 += "     "
    scale1 += " " * (
        10 - len(str(sc_obj.get_value("fract_transf_vector[1]", 0)))
    ) + str(sc_obj.get_value("fract_transf_vector[1]", 0))
    scale2 = ""
    scale2 += "SCALE2    "
    scale2 += " " * (
        10 - len(str(sc_obj.get_value("fract_transf_matrix[2][1]", 0)))
    ) + str(sc_obj.get_value("fract_transf_matrix[2][1]", 0))
    scale2 += " " * (
        10 - len(str(sc_obj.get_value("fract_transf_matrix[2][2]", 0)))
    ) + str(sc_obj.get_value("fract_transf_matrix[2][2]", 0))
    scale2 += " " * (
        10 - len(str(sc_obj.get_value("fract_transf_matrix[2][3]", 0)))
    ) + str(sc_obj.get_value("fract_transf_matrix[2][3]", 0))
    scale2 += "     "
    scale2 += " " * (
        10 - len(str(sc_obj.get_value("fract_transf_vector[2]", 0)))
    ) + str(sc_obj.get_value("fract_transf_vector[2]", 0))
    scale3 = ""
    scale3 += "SCALE3    "
    scale3 += " " * (
        10 - len(str(sc_obj.get_value("fract_transf_matrix[3][1]", 0)))
    ) + str(sc_obj.get_value("fract_transf_matrix[3][1]", 0))
    scale3 += " " * (
        10 - len(str(sc_obj.get_value("fract_transf_matrix[3][2]", 0)))
    ) + str(sc_obj.get_value("fract_transf_matrix[3][2]", 0))
    scale3 += " " * (
        10 - len(str(sc_obj.get_value("fract_transf_matrix[3][3]", 0)))
    ) + str(sc_obj.get_value("fract_transf_matrix[3][3]", 0))
    scale3 += "     "
    scale3 += " " * (
        10 - len(str(sc_obj.get_value("fract_transf_vector[3]", 0)))
    ) + str(sc_obj.get_value("fract_transf_vector[3]", 0))
    try:
        sc_arr.append(pdb.SCALE1(scale1))
    except KeyError:
        _LOGGER.error(f"cif.scalen:    Error parsing line:\n{scale1}")
        sc_err.append("SCALE1")
    try:
        sc_arr.append(pdb.SCALE2(scale2))
    except KeyError:
        _LOGGER.error(f"cif.scalen:    Error parsing line:\n{scale2}")
        sc_err.append("SCALE2")
    try:
        sc_arr.append(pdb.SCALE3(scale3))
    except KeyError:
        _LOGGER.error(f"cif.scalen:    Error parsing line:\n{scale3}")
        sc_err.append("SCALE3")
    return sc_arr, sc_err


def origxn(block):
    """Handle ORIGXn block.

    :param block:  PDBx data block
    :type block:  [str]
    :return:  (array of pdb.conect objects, array of things that did not parse)
    :rtype:  ([pdb.CONECT], [str])
    """
    or_arr = []
    or_err = []
    or_obj = block.get_object("database_PDB_matrix")
    orig1 = "ORIGX1    "
    orig1 += " " * (10 - len(str(or_obj.get_value("origx[1][1]", 0)))) + str(
        or_obj.get_value("origx[1][1]", 0)
    )
    orig1 += " " * (10 - len(str(or_obj.get_value("origx[1][2]", 0)))) + str(
        or_obj.get_value("origx[1][2]", 0)
    )
    orig1 += " " * (10 - len(str(or_obj.get_value("origx[1][3]", 0)))) + str(
        or_obj.get_value("origx[1][3]", 0)
    )
    orig1 += "     "
    orig1 += " " * (
        10 - len(str(or_obj.get_value("origx_vector[1]", 0)))
    ) + str(or_obj.get_value("origx_vector[1]", 0))
    orig2 = "ORIGX2    "
    orig2 += " " * (10 - len(str(or_obj.get_value("origx[2][1]", 0)))) + str(
        or_obj.get_value("origx[2][1]", 0)
    )
    orig2 += " " * (10 - len(str(or_obj.get_value("origx[2][2]", 0)))) + str(
        or_obj.get_value("origx[2][2]", 0)
    )
    orig2 += " " * (10 - len(str(or_obj.get_value("origx[2][3]", 0)))) + str(
        or_obj.get_value("origx[2][3]", 0)
    )
    orig2 += "     "
    orig2 += " " * (
        10 - len(str(or_obj.get_value("origx_vector[2]", 0)))
    ) + str(or_obj.get_value("origx_vector[2]", 0))
    orig3 = "ORIGX3    "
    orig3 += " " * (10 - len(str(or_obj.get_value("origx[3][1]", 0)))) + str(
        or_obj.get_value("origx[3][1]", 0)
    )
    orig3 += " " * (10 - len(str(or_obj.get_value("origx[3][2]", 0)))) + str(
        or_obj.get_value("origx[3][2]", 0)
    )
    orig3 += " " * (10 - len(str(or_obj.get_value("origx[3][3]", 0)))) + str(
        or_obj.get_value("origx[3][3]", 0)
    )
    orig3 += "     "
    orig3 += " " * (
        10 - len(str(or_obj.get_value("origx_vector[3]", 0)))
    ) + str(or_obj.get_value("origx_vector[3]", 0))
    try:
        or_arr.append(pdb.ORIGX1(orig1))
    except KeyError:
        _LOGGER.error(f"cif.origxn:  Error parsing line:\n{orig1}")
        or_err.append("ORIGX1")
    try:
        or_arr.append(pdb.ORIGX2(orig2))
    except KeyError:
        _LOGGER.error(f"cif.origxn:  Error parsing line:\n{orig2}")
        or_err.append("ORIGX2")
    try:
        or_arr.append(pdb.ORIGX3(orig3))
    except KeyError:
        _LOGGER.error(f"cif.origxn:  Error parsing line:\n{orig3}")
        or_err.append("ORIGX3")
    return or_arr, or_err


def cispep(block):
    """Handle CISPEP block.

    :param block:  PDBx data block
    :type block:  [str]
    :return:  (array of pdb.conect objects, array of things that did not parse)
    :rtype:  ([pdb.CONECT], [str])
    """
    cis_arr = []
    cis_err = []
    cis_obj = block.get_object("struct_mon_prot_cis")
    if cis_obj is None:
        return cis_arr, cis_err
    for i in range(cis_obj.row_count):
        line = "CISPEP "
        line += " " * (3 - len(str(cis_obj.get_value("pdbx_id", i)))) + str(
            cis_obj.get_value("pdbx_id", i)
        )
        line += " "
        line += " " * (
            3 - len(cis_obj.get_value("auth_comp_id", i))
        ) + cis_obj.get_value("auth_comp_id", i)
        line += " "
        line += cis_obj.get_value("auth_asym_id", i)
        line += " "
        line += " " * (
            4 - len(str(cis_obj.get_value("auth_seq_id", i)))
        ) + str(cis_obj.get_value("auth_seq_id", i))
        value = cis_obj.get_value("pdbx_PDB_ins_code", i)
        if value not in ["?", None]:
            line += value
        else:
            line += " "
        line += "   "
        line += " " * (
            3 - len(cis_obj.get_value("pdbx_auth_comp_id_2", i))
        ) + cis_obj.get_value("pdbx_auth_comp_id_2", i)
        line += " "
        line += cis_obj.get_value("pdbx_auth_asym_id_2", i)
        line += " "
        line += " " * (
            4 - len(str(cis_obj.get_value("pdbx_auth_seq_id_2", i)))
        ) + str(cis_obj.get_value("pdbx_auth_seq_id_2", i))
        value = cis_obj.get_value("pdbx_PDB_ins_code_2", i)
        if value not in ["?", None]:
            line += value
        else:
            line += " "
        line += " " * 7
        line += " " * (
            3 - len(str(cis_obj.get_value("pdbx_PDB_model_num", i)))
        ) + str(cis_obj.get_value("pdbx_PDB_model_num", i))
        line += " " * 7
        line += " " * (
            6 - len(str(cis_obj.get_value("pdbx_omega_angle", i)))
        ) + str(cis_obj.get_value("pdbx_omega_angle", i))
        try:
            cis_arr.append(pdb.CISPEP(line))
        except KeyError:
            _LOGGER.error(f"cif.cispep:    Erro parsing line:\n{line}")
            cis_err.append("cispep")
    return cis_arr, cis_err


def ssbond(block):
    """Handle SSBOND block.

    :param block:  PDBx data block
    :type block:  [str]
    :return:  (array of pdb.conect objects, array of things that did not parse)
    :rtype:  ([pdb.CONECT], [str])
    """
    ssb_arr = []
    ssb_err = []
    ssb_obj = block.get_object("struct_conn")
    if ssb_obj is None:
        return ssb_arr, ssb_err
    for i in range(ssb_obj.row_count):
        line = "SSBOND "
        line += " " * (3 - len(str(ssb_obj.get_value("id", i)[-1]))) + str(
            ssb_obj.get_value("id", i)[-1]
        )
        line += " "
        line += " " * (
            3 - len(ssb_obj.get_value("ptnr1_auth_comp_id", i))
        ) + ssb_obj.get_value("ptnr1_auth_comp_id", i)
        line += " "
        line += ssb_obj.get_value("ptnr1_auth_asym_id", i)
        line += " "
        line += " " * (
            4 - len(str(ssb_obj.get_value("ptnr1_auth_seq_id", i)))
        ) + str(ssb_obj.get_value("ptnr1_auth_seq_id", i))
        value = ssb_obj.get_value("pdbx_ptnr1_PDB_ins_code", i)
        if value not in ["?", None]:
            line += value
        else:
            line += " "
        line += " " * 3
        line += " " * (
            3 - len(ssb_obj.get_value("ptnr2_auth_comp_id", i))
        ) + ssb_obj.get_value("ptnr2_auth_comp_id", i)
        line += " "
        line += ssb_obj.get_value("ptnr2_auth_asym_id", i)
        line += " "
        line += " " * (
            4 - len(str(ssb_obj.get_value("ptnr2_auth_seq_id", i)))
        ) + str(ssb_obj.get_value("ptnr2_auth_seq_id", i))
        value = ssb_obj.get_value("pdbx_ptnr2_PDB_ins_code", i)
        if value not in ["?", None]:
            line += value
        else:
            line += " "
        line += " " * 23
        line += " " * (
            6 - len(ssb_obj.get_value("ptnr1_symmetry", i).replace("_", ""))
        ) + ssb_obj.get_value("ptnr1_symmetry", i).replace("_", "")
        line += " "
        line += " " * (
            6 - len(ssb_obj.get_value("ptnr2_symmetry", i).replace("_", ""))
        ) + ssb_obj.get_value("ptnr2_symmetry", i).replace("_", "")
        line += " "
        line += " " * (
            5 - len(str(ssb_obj.get_value("pdbx_dist_value", i)))
        ) + str(ssb_obj.get_value("pdbx_dist_value", i))
        try:
            ssb_arr.append(pdb.SSBOND(line))
        except KeyError:
            _LOGGER.error(f"cif.ssbond:    Error parsing line:\n{line}")
            ssb_err.append("ssbond")
    return ssb_arr, ssb_err


def count_models(block):
    """Count models in structure file block.

    :param block:  PDBx data block
    :type block:  [str]
    :return:  number of models in block
    :rtype:  int
    """
    atom_obj = block.get_object("atom_site")
    model_num = []
    for i in range(atom_obj.row_count):
        tmp = atom_obj.get_value("pdbx_PDB_model_num", i)
        if tmp not in model_num:
            model_num.append(tmp)
    return model_num


def read_cif(cif_file):
    """Parse CIF-format data into array of Atom objects.

    .. todo::  Manage several blocks of data.

    :param file:  open file-like object
    :type file:  file
    :return:  (a dictionary indexed by PDBx/CIF record names, a list of record
        names that couldn't be parsed)
    :rtype:  (dict, [str])
    """
    pdblist = []  # Array of parsed lines (as objects)
    errlist = []  # List of record names that couldn't be parsed.
    if cif_file is None:
        return pdblist, errlist
    pdbdata = pdbx.load(cif_file)
    if len(pdbdata) > 0:
        for block in pdbdata:
            head_pdb, head_err = header(block)
            title_pdb, title_err = title(block)
            cmpnd_pdb, cmpnd_err = compnd(block)
            src_pdb, src_err = source(block)
            key_pdb, key_err = keywds(block)
            ex_pdb, ex_err = expdata(block)
            aut_pdb, aut_err = author(block)
            ssb_pdb, ssb_err = ssbond(block)
            cis_pdb, cis_err = cispep(block)
            cry_pdb, cry_err = cryst1(block)
            or_pdb, or_err = origxn(block)
            sc_pdb, sc_err = scalen(block)
            ato_pdb, ato_err = atom_site(block)
            con_pdb, con_err = conect(block)
            pdblist = (
                head_pdb
                + title_pdb
                + cmpnd_pdb
                + src_pdb
                + key_pdb
                + ex_pdb
                + aut_pdb
                + ssb_pdb
                + cis_pdb
                + cry_pdb
                + or_pdb
                + sc_pdb
                + ato_pdb
                + con_pdb
            )
            errlist = (
                head_err
                + title_err
                + cmpnd_err
                + src_err
                + key_err
                + ex_err
                + aut_err
                + ssb_err
                + cis_err
                + cry_err
                + or_err
                + sc_err
                + ato_err
                + con_err
            )
    else:
        _LOGGER.error("Unknown error while reading CIF file.")

    return pdblist, errlist
