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


def atom_site(block: pdbx.containers.ContainerBase):
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
    if not isinstance(atoms, pdbx.containers.DataCategory):
        _LOGGER.error("atom_site: No lines were found\n")
        raise ValueError(
            "Did not find atom_site in cif file, terminal failure."
        )
    num_model_arr = count_models(block)
    if len(num_model_arr) == 1:
        # TODO - this part of the conditional should be a separate function
        for i in range(atoms.row_count):
            result = convert_cif_atom_site_to_pdb_line(atoms=atoms, row_index=i)
            if result is None:
                continue
            line, serial, chain, res_seq = result
            try:
                if atoms.get_value("group_PDB", i) == "ATOM":
                    record = pdb.ATOM(line)
                elif atoms.get_value("group_PDB", i) == "HETATM":
                    record = pdb.HETATM(line)
                else:
                    continue
                # Restore the true, untruncated values that do not fit the
                # fixed PDB columns.
                record.serial = serial
                record.chain_id = chain
                record.res_seq = res_seq
                pdb_arr.append(record)
            except KeyError:
                _LOGGER.error(f"atom_site: Error reading line: #{line}#\n")
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
            except KeyError:
                _LOGGER.error(f"atom_site: Error readline line:\n{line}")
                err_arr.append("MODEL")

            for i in range(atoms.row_count):
                if atoms.get_value("pdbx_PDB_model_num", i) == j:
                    result = convert_cif_atom_site_to_pdb_line(
                        atoms=atoms, row_index=i
                    )
                    if result is None:
                        continue
                    line, serial, chain, res_seq = result
                    try:
                        if atoms.get_value("group_PDB", i) == "ATOM":
                            record = pdb.ATOM(line)
                        elif atoms.get_value("group_PDB", i) == "HETATM":
                            record = pdb.HETATM(line)
                        else:
                            continue
                        # Restore the true, untruncated values that do not fit
                        # the fixed PDB columns.
                        record.serial = serial
                        record.chain_id = chain
                        record.res_seq = res_seq
                        pdb_arr.append(record)
                    except (KeyError, ValueError):
                        _LOGGER.error(
                            f"atom_site: Error reading line: #{line}#\n"
                        )
            try:
                line = "ENDMDL"
                pdb_arr.append(pdb.ENDMDL(line))
            except KeyError:
                _LOGGER.error(f"atom_site: Error reading line:\n{line}")
                err_arr.append("ENDMDL")
        return pdb_arr, err_arr


def convert_cif_atom_site_to_pdb_line(
    atoms: pdbx.containers.DataCategory,
    row_index: int,
) -> tuple[str, int, str, int] | None:
    """Converts cif data for one row into a pdb line.

    Extracts all the relevant fields from atoms, does basic formatting,
    and converts them into a pdb line.

    The atom serial, chain id, and residue sequence number do not always fit
    in the fixed PDB columns (serial > 99999, multi-character chain ids,
    res_seq outside [-999, 9999]).  For those, the line carries a clamped
    placeholder so it stays a valid 80-char record, and the true untruncated
    values are returned alongside it so the caller can restore them on the
    parsed :class:`pdb.ATOM`/:class:`pdb.HETATM` record.

    :param atoms: DataCategory containing all atoms.
    :type atoms: pdbx.containers.DataCategory
    :param row_index: Index corresponding to the cif line we want to convert to a pdb line
    :type row_index: int
    :return: tuple of (pdb line exactly 80 chars long, true serial, true chain
        id, true res_seq), or None if the row failed to process properly
    :rtype: tuple[str, int, str, int] | None
    """
    # Extract and cast values
    group = atoms.get_value("group_PDB", row_index)
    serial = int(atoms.get_value("id", row_index))
    name = atoms.get_value("label_atom_id", row_index)
    alt_id = atoms.get_value("label_alt_id", row_index)
    res_name = atoms.get_value("label_comp_id", row_index)
    chain = atoms.get_value("label_asym_id", row_index)
    res_seq = int(atoms.get_value("auth_seq_id", row_index))
    x = float(atoms.get_value("Cartn_x", row_index))
    y = float(atoms.get_value("Cartn_y", row_index))
    z = float(atoms.get_value("Cartn_z", row_index))
    occ_raw = atoms.get_value("occupancy", row_index)
    try:
        occ = float(occ_raw)
    except (TypeError, ValueError):
        _LOGGER.error(
            f"Invalid occupancy {occ_raw} at row {row_index}, setting to 1.0"
        )
        occ = 1.0
    temp_raw = atoms.get_value("B_iso_or_equiv", row_index)
    try:
        temp = float(temp_raw)
    except (TypeError, ValueError):
        _LOGGER.error(
            f"Invalid B-factor {temp_raw} at row {row_index}, setting to 0.0"
        )
        temp = 0.0
    element = atoms.get_value("type_symbol", row_index)

    # Handle the '?' or '.' cases for alt_id and charge
    alt_id = (
        alt_id if (alt_id is not None and alt_id not in (".", "?")) else " "
    )
    charge = (
        atoms.get_value("pdbx_formal_charge", row_index)
        if "pdbx_formal_charge" in atoms.attribute_list
        else "  "
    )
    if charge in ["?", None, "."]:
        charge = "  "

    # Format the Atom Name (4 chars, specifically aligned)
    # 1-2 char elements start at the 2nd position of the 4-char field
    if len(element) < 2 and len(name) < 4:
        name = f" {name:<3}"
    else:
        name = f"{name:<4}"

    # The PDB line is a fixed-column format that cannot represent serials
    # > 99999, multi-character chain ids, or residue sequence numbers outside
    # [-999, 9999].  We only need the line to parse cleanly for the remaining
    # fields (name, coordinates, element, ...); the true, untruncated values
    # for these three fields are returned alongside the line so the caller can
    # overwrite them on the parsed record.  See atom_site().
    line_serial = serial if 0 <= serial <= 99999 else 99999
    line_chain = chain[0] if chain else " "
    line_res_seq = res_seq if -999 <= res_seq <= 9999 else 9999

    # THE PDB LINE MAPPING
    line = (
        f"{group:<6}"  #    1-6 Record name
        f"{line_serial:>5}"  #   7-11 Atom serial number
        f" "  #             12 Space
        f"{name}"  #        13-16 Atom name
        f"{alt_id:1}"  #    17 Alternate location indicator
        f"{res_name:>3}"  # 18-20 Residue name
        f" "  #             21 Space
        f"{line_chain:1}"  #     22 Chain Id
        f"{line_res_seq:>4}"  #  23-26 Residue sequence number
        f" "  #             27 Code for insertion of residues
        f"   "  #           28-30 Spaces
        f"{x:8.3f}"  #      31-38 X-coordinates
        f"{y:8.3f}"  #      39-46 Y-coordinates
        f"{z:8.3f}"  #      47-54 Z-coordinates
        f"{occ:6.2f}"  #    55-60 Occupancy
        f"{temp:6.2f}"  #   61-66 Temperature factor
        f"          "  #    67-76 Spaces
        f"{element:>2}"  #  77-78 Element Symbol
        f"{charge:>2}"  #   79-80 Charge of the atom
    )
    if len(line) != 80:
        # This would happen if any of the values above are larger than the max space they are allocated, e.g.:
        # if x < -999.999
        # if name = LIG1
        _LOGGER.error(
            f"atom_site: pdb line extracted from cif record is not exactly 80 char long, dropping line:\n{line}"
        )
        return None
    return line, serial, chain, res_seq


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
    line = "HEADER"
    line += " " * 4
    if struct_obj is not None:
        line += struct_obj.get_value("pdbx_keywords") + " " * (
            40 - len(struct_obj.get_value("pdbx_keywords"))
        )
    else:
        line += " " * 40
    if database_obj is not None:
        ridd = database_obj.get_value("recvd_initial_deposition_date")
        if isinstance(ridd, str) and len(ridd) > 9:
            ridd = (
                datetime.strptime(ridd, "%Y-%m-%d")
                .strftime("%d-%b-%y")
                .upper()
            )
            line += " " * (9 - len(ridd)) + ridd
    else:
        line += " " * 9
    line += " " * 3
    if entry_obj is not None:
        line += " " * (
            4 - len(entry_obj.get_value("id"))
        ) + entry_obj.get_value("id")
    else:
        line += " " * 4
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


def source(block: pdbx.containers.ContainerBase):
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
        if ("entity_id" in src_obj.attribute_list) and (
            src_obj.get_value("entity_id", i) not in ["?", None]
        ):
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
        if ("pdbx_gene_src_scientific_name" in src_obj.attribute_list) and (
            src_obj.get_value("pdbx_gene_src_scientific_name", i)
            not in ["?", None]
        ):
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
        if ("gene_src_common_name" in src_obj.attribute_list) and (
            src_obj.get_value("gene_src_common_name", i) not in ["?", None]
        ):
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
        if ("pdbx_gene_src_ncbi_taxonomy_id" in src_obj.attribute_list) and (
            src_obj.get_value("pdbx_gene_src_ncbi_taxonomy_id", i)
            not in ["?", None]
        ):
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
    if key_obj is None:
        return key_arr, key_err
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
    if ex_obj is None:
        return ex_arr, ex_err
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
    if aut_obj is None:
        return aut_arr, aut_err
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
    if len(pdbdata) == 0:
        _LOGGER.error("Unknown error while reading CIF file.")
        return pdblist, errlist

    # Only ``atom_site`` is required to build the biomolecule.  The descriptive
    # header categories are best-effort: they feed the optional regenerated
    # output header and nothing in the calculation, so a missing or unusual
    # category must not abort the read.
    #
    # Connectivity/geometry records (SSBOND/CISPEP/CONECT/CRYST1/SCALE/ORIGX)
    # are intentionally NOT parsed.  pdb2pqr never consumes them -- disulfides
    # are detected by distance in :func:`Biomolecule.update_ss_bridges`, bonds
    # come from the reference topology, and none of these records are written
    # back out.  Their fixed-column reconstruction also crashes on large
    # assemblies (multi-character chains, residue numbers > 9999).
    header_parsers = (header, title, compnd, source, keywds, expdata, author)
    for block in pdbdata:
        for parser in header_parsers:
            try:
                recs, errs = parser(block)
            except Exception as exc:
                _LOGGER.warning(
                    f"Skipping CIF header category '{parser.__name__}': {exc}"
                )
                continue
            pdblist += recs
            errlist += errs
        atom_pdb, atom_err = atom_site(block)
        pdblist += atom_pdb
        errlist += atom_err

    return pdblist, errlist
