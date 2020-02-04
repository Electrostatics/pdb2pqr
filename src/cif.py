""" CIF parsing methods

    This methods use the pdbx/cif parser provided by WWPDB
    (http://mmcif.wwpdb.org/docs/sw-examples/python/html/index.html)

    ----------------------------

    PDB2PQR -- An automated pipeline for the setup, execution, and analysis of
    Poisson-Boltzmann electrostatics calculations

    Copyright (c) 2002-2011, Jens Erik Nielsen, University College Dublin;
    Nathan A. Baker, Battelle Memorial Institute, Developed at the Pacific
    Northwest National Laboratory, operated by Battelle Memorial Institute,
    Pacific Northwest Division for the U.S. Department Energy.;
    Paul Czodrowski & Gerhard Klebe, University of Marburg.

	All rights reserved.

	Redistribution and use in source and binary forms, with or without modification,
	are permitted provided that the following conditions are met:

		* Redistributions of source code must retain the above copyright notice,
		  this list of conditions and the following disclaimer.
		* Redistributions in binary form must reproduce the above copyright notice,
		  this list of conditions and the following disclaimer in the documentation
		  and/or other materials provided with the distribution.
        * Neither the names of University College Dublin, Battelle Memorial Institute,
          Pacific Northwest National Laboratory, US Department of Energy, or University
          of Marburg nor the names of its contributors may be used to endorse or promote
          products derived from this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
	ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
	IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
	INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
	BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
	OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
	OF THE POSSIBILITY OF SUCH DAMAGE.

    ----------------------------
"""

__date__ = "2 January 2020"
__author__ = "Juan Brandi"

import os
import sys
import pdb
from pdbx.reader.PdbxReader import PdbxReader;

def ATOM_SITE(block):
    """
    Data items in the ATOM_SITE category record details about
    the atom sites in a macromolecular crystal structure, such as
    the positional coordinates, atomic displacement parameters,
    magnetic moments and directions.
    (source: http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/atom_site.html)

    Parameters:
        block: Pdbx data block
    Returs:
        pdblist: array of pdb.ATOM objects;
        errlist: array of thigs that
    """
    pdb_arr = [];
    err_arr = [];

    atoms= block.getObj("atom_site");

    for i in range(atoms.getRowCount()):

        if(atoms.getValue("group_PDB", i) == "ATOM"):
            try:
                line = "";
                line += atoms.getValue("group_PDB", i) + " "*(6 - len(atoms.getValue("group_PDB", i)));                           # 1  - 6  RECORD NAME (ATOM)
                line += " "*(5 - len(str(atoms.getValue("id", i))))       + str(atoms.getValue("id", i));                         # 7  - 11 ATOM SERIAL
                line += "  ";                                                                                                     # 12 - 13
                line += atoms.getValue("label_atom_id", i) + " "*(3 - len(atoms.getValue("label_atom_id", i)));                   # 14 - 16 ATOM NAME
                line += " " if(atoms.getValue("label_alt_id", i)==".") else atoms.getValue("label_alt_id", i);                    # 17      ALT LOCATION
                line += " "*(3 - len(atoms.getValue("label_comp_id", i))) + atoms.getValue("label_comp_id", i);                   # 18 - 20 RES NAME
                line += " ";                                                                                                      # 21
                line += " "*(1 - len(atoms.getValue("label_asym_id", i))) + atoms.getValue("label_asym_id", i);                   # 22      CHAIN ID
                line += " "*(4 - len(str(atoms.getValue("auth_seq_id", i)))) +  str(atoms.getValue("auth_seq_id", i));            # 23 - 26 RES SEQ ID
                line += " "*3;                                                                                                    # 27 - 30
                line += " "*(8 - len(str(atoms.getValue("Cartn_x", i)))) + str(atoms.getValue("Cartn_x"));                        # 31 - 38 X Coords
                line += " "*(8 - len(str(atoms.getValue("Cartn_y", i)))) + str(atoms.getValue("Cartn_y"));                        # 39 - 46 Y Coords
                line += " "*(8 - len(str(atoms.getValue("Cartn_z", i)))) + str(atoms.getValue("Cartn_z"));                        # 47 - 54 Z Coords
                line += " "*(6 - len(str(atoms.getValue("occupancy", i)))) + str(atoms.getValue("occupancy", i));                 # 55 - 60 OCCUPANCY
                line += " "*(6 - len(str(atoms.getValue("B_iso_or_equiv", i)))) + str(atoms.getValue("B_iso_or_equiv", i));       # 61 - 66 TEMP FACTOR
                line += " "*(10);                                                                                                 # 67 - 76
                line += " "*(2 - len(atoms.getValue("type_symbol", i))) + atoms.getValue("type_symbol", i);                       # 77 - 78 ELEMENT SYMBOL
                line += " "*2 if(atoms.getValue("pdbx_formal_charge", i) == "?") else atoms.getValue("pdbx_formal_charge", i);    # 79 - 80 CHARGE OF ATOM
                pdb_arr.append(pdb.ATOM(line));
            except:
                print("cif.ATOM_SITE: Error reading line:\n%s" % line);
                
        elif(atoms.getValue("group_PDB", i) == "HETATM"):
            try:
                line = "";
                line += atoms.getValue("group_PDB", i) + ""*(6 - len(atoms.getValue("group_PDB", i)));                            # 1  - 6  RECORD NAME (HETATM)
                line += " "*(5 - len(str(atoms.getValue("id", i))))       + str(atoms.getValue("id", i));                         # 7  - 11 ATOM SERIAL
                line += "  ";                                                                                                     # 12 - 13
                line += atoms.getValue("label_atom_id", i) + " "*(3 - len(atoms.getValue("label_atom_id", i)));                   # 14 - 16 ATOM NAME
                line += " " if(atoms.getValue("label_alt_id", i)==".") else atoms.getValue("label_alt_id", i);                    # 17      ALT LOCATION
                line += " "*(3 - len(atoms.getValue("label_comp_id", i))) + atoms.getValue("label_comp_id", i);                   # 18 - 20 RES NAME
                line += " ";                                                                                                      # 21
                line += " "*(1 - len(atoms.getValue("label_asym_id", i))) + atoms.getValue("label_asym_id", i);                   # 22      CHAIN ID
                line += " "*(4 - len(str(atoms.getValue("auth_seq_id", i)))) +  str(atoms.getValue("auth_seq_id", i));            # 23 - 26 RES SEQ ID
                line += " "*3;                                                                                                    # 27 - 30
                line += " "*(8 - len(str(atoms.getValue("Cartn_x", i)))) + str(atoms.getValue("Cartn_x"));                        # 31 - 38 X Coords
                line += " "*(8 - len(str(atoms.getValue("Cartn_y", i)))) + str(atoms.getValue("Cartn_y"));                        # 39 - 46 Y Coords
                line += " "*(8 - len(str(atoms.getValue("Cartn_z", i)))) + str(atoms.getValue("Cartn_z"));                        # 47 - 54 Z Coords
                line += " "*(6 - len(str(atoms.getValue("occupancy", i)))) + str(atoms.getValue("occupancy", i));                 # 55 - 60 OCCUPANCY
                line += " "*(6 - len(str(atoms.getValue("B_iso_or_equiv", i)))) + str(atoms.getValue("B_iso_or_equiv", i));       # 61 - 66 TEMP FACTOR
                line += " "*(10);                                                                                                 # 67 - 76
                line += " "*(2 - len(atoms.getValue("type_symbol", i))) + atoms.getValue("type_symbol", i);                       # 77 - 78 ELEMENT SYMBOL
                line += " "*2 if(atoms.getValue("pdbx_formal_charge", i) == "?") else atoms.getValue("pdbx_formal_charge", i);    # 79 - 80 CHARGE OF ATOM
                pdb_arr.append(pdb.HETATM(line));
            except ValueError as e:
                print("cif.ATOM_SITE: Error reading line:\n%s" % line);

    return pdb_arr, err_arr;

def CONECT(block):
    """
        Data items in the STRUCT_CONN category record details about
        the connections between portions of the structure. These can be
        hydrogen bonds, salt bridges, disulfide bridges and so on.

        The STRUCT_CONN_TYPE records define the criteria used to
        identify these connections.
        (source: http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/struct_conn.html)

        Parameters:
            block: Pdbx data block
        Returs:
            pdblist: array of pdb.CONECT objects;
            errlist: array of thigs that
    """
    pdb_arr = [];
    err_arr = [];

    struct_conn = block.getObj("struct_conn");
    atoms = block.getObj("atom_site");

    for index in range(struct_conn.getRowCount()):
        atom_pair = [];
        for partner in ["ptnr1_", "ptnr2_"] :
            # Retrieve all the information necessary to uniquely identify the atom4
            atom = {"auth_seq_id" : struct_conn.getValue(partner + "auth_seq_id", index),
                "auth_comp_id" : struct_conn.getValue(partner + "auth_comp_id", index),
                "auth_asym_id" : struct_conn.getValue(partner + "auth_asym_id", index),
                "label_atom_id" : struct_conn.getValue(partner + "label_atom_id", index),
                }

            for i in range(atoms.getRowCount()):
                found = True;
                for key in atom.keys() :
                    if atoms.getValue(key, i) != atom[key] :
                        found = False

                if(found):
                    atom_pair.append(atoms.getValue("id", i));


        if(len(atom_pair) == 2):
            p1 = "CONECT" + " "*(5 - len(str(atom_pair[0]))) + str(atom_pair[0]) + \
                ' '*(5 - len(str(atom_pair[1]))) + str(atom_pair[1]);

            try:
                pdb_arr.append(pdb.CONECT(p1));
            except:
                sys.stderr.write("cif.CONECT:   Error parsing line: \n%s" % p1);
                err_arr.append("CONECT")

    return pdb_arr, err_arr;

def HEADER(block):

    header_arr = [];
    header_err = [];

    struct_obj   = block.getObj("struct_keywords");
    database_obj = block.getObj("pdbx_database_status");
    entry_obj    = block.getObj("entry");

    line = "HEADER";
    line += " "*4;
    line += struct_obj.getValue("pdbx_keywords") + " "*(40 - len(struct_obj.getValue("pdbx_keywords")));
    line += " "*(9 - len(database_obj.getValue("recvd_initial_deposition_date"))) + database_obj.getValue("recvd_initial_deposition_date");
    line += " "*(4 - len(entry_obj.getValue("id"))) + entry_obj.getValue("id");

    try:
        header_arr.append(pdb.HEADER(line));
    except:
        sys.stderr.write("cif.HEADER:   Error parsing line: \n%s" % line);
        header_err.append("HEADER");

    return header_arr, header_err;


def readCIF(file):
    """ Parse CIF-format data into array of Atom objects.
        Parameters:
            file: open file object
        Returns (dict, errlist):
            dict:       a dictionary indexed by PDBx/CIF record names.
            errlist:    a list of record names that couldn't be parsed.
    """

    pdblist = []; # Array of parsed lines (as objects)
    errlist = []; # List of record names that couldn't be parsed.

    if file is None:
        return pdblist, errlist;

    pdbdata = [];

    reader = PdbxReader(file);
    reader.read(pdbdata);

    # @TODO manage several blocks of data.
    if(len(pdbdata) > 0):

        for i in range(len(pdbdata)):
            block = pdbdata[i];
            head_pdb, head_err = HEADER(block);
            ato_pdb, ato_err = ATOM_SITE(block);
            con_pdb, con_err = CONECT(block);

            pdblist = ato_pdb + con_pdb;
            errlist = ato_err + con_err;

        return pdblist, errlist

    else:
        print("Unknown error while reading cif file.");
        return pdblist, errlist;
