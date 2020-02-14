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
from numpy import minimum, ceil
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
    
    numModel_arr = CountModels(block);
    
    if(len(numModel_arr) == 1):

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
                    line += " "*(8 - len(str(atoms.getValue("Cartn_x", i)))) + str(atoms.getValue("Cartn_x", i));                     # 31 - 38 X Coords
                    line += " "*(8 - len(str(atoms.getValue("Cartn_y", i)))) + str(atoms.getValue("Cartn_y", i));                     # 39 - 46 Y Coords
                    line += " "*(8 - len(str(atoms.getValue("Cartn_z", i)))) + str(atoms.getValue("Cartn_z", i));                     # 47 - 54 Z Coords
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
                    line += " "*(8 - len(str(atoms.getValue("Cartn_x", i)))) + str(atoms.getValue("Cartn_x", i));                     # 31 - 38 X Coords
                    line += " "*(8 - len(str(atoms.getValue("Cartn_y", i)))) + str(atoms.getValue("Cartn_y", i));                     # 39 - 46 Y Coords
                    line += " "*(8 - len(str(atoms.getValue("Cartn_z", i)))) + str(atoms.getValue("Cartn_z", i));                     # 47 - 54 Z Coords
                    line += " "*(6 - len(str(atoms.getValue("occupancy", i)))) + str(atoms.getValue("occupancy", i));                 # 55 - 60 OCCUPANCY
                    line += " "*(6 - len(str(atoms.getValue("B_iso_or_equiv", i)))) + str(atoms.getValue("B_iso_or_equiv", i));       # 61 - 66 TEMP FACTOR
                    line += " "*(10);                                                                                                 # 67 - 76
                    line += " "*(2 - len(atoms.getValue("type_symbol", i))) + atoms.getValue("type_symbol", i);                       # 77 - 78 ELEMENT SYMBOL
                    line += " "*2 if(atoms.getValue("pdbx_formal_charge", i) == "?") else atoms.getValue("pdbx_formal_charge", i);    # 79 - 80 CHARGE OF ATOM
                    pdb_arr.append(pdb.HETATM(line));
                     
                except:
                    print("cif.ATOM_SITE: Error reading line:\n%s" % line);
    
        return pdb_arr, err_arr;
    
    else:
        for j in numModel_arr:
            try:
                line  = "MODEL ";
                line += " "*4;
                line += " "*(4 - len(str(j))) + str(j);
                pdb_arr.append(pdb.MODEL(line));
            except:
                sys.stderr.write("cif.ATOM_SITE: Error readline line:\n%s\n" % line);
                err_arr.append("MODEL");
                
            for i in range(atoms.getRowCount()):
                if(atoms.getValue("pdbx_PDB_model_num", i) == j):
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
                            line += " "*(8 - len(str(atoms.getValue("Cartn_x", i)))) + str(atoms.getValue("Cartn_x", i));                     # 31 - 38 X Coords
                            line += " "*(8 - len(str(atoms.getValue("Cartn_y", i)))) + str(atoms.getValue("Cartn_y", i));                     # 39 - 46 Y Coords
                            line += " "*(8 - len(str(atoms.getValue("Cartn_z", i)))) + str(atoms.getValue("Cartn_z", i));                     # 47 - 54 Z Coords
                            line += " "*(6 - len(str(atoms.getValue("occupancy", i)))) + str(atoms.getValue("occupancy", i));                 # 55 - 60 OCCUPANCY
                            line += " "*(6 - len(str(atoms.getValue("B_iso_or_equiv", i)))) + str(atoms.getValue("B_iso_or_equiv", i));       # 61 - 66 TEMP FACTOR
                            line += " "*(10);                                                                                                 # 67 - 76
                            line += " "*(2 - len(atoms.getValue("type_symbol", i))) + atoms.getValue("type_symbol", i);                       # 77 - 78 ELEMENT SYMBOL
                            line += " "*2 if(atoms.getValue("pdbx_formal_charge", i) == "?") else atoms.getValue("pdbx_formal_charge", i);    # 79 - 80 CHARGE OF ATOM
                            pdb_arr.append(pdb.ATOM(line));
                        except:
                            print("cif.ATOM_SITE: Error reading line:\n%s" % line);
                            err_arr.append("ATOM");
                            
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
                            line += " "*(8 - len(str(atoms.getValue("Cartn_x", i)))) + str(atoms.getValue("Cartn_x", i));                     # 31 - 38 X Coords
                            line += " "*(8 - len(str(atoms.getValue("Cartn_y", i)))) + str(atoms.getValue("Cartn_y", i));                     # 39 - 46 Y Coords
                            line += " "*(8 - len(str(atoms.getValue("Cartn_z", i)))) + str(atoms.getValue("Cartn_z", i));                     # 47 - 54 Z Coords
                            line += " "*(6 - len(str(atoms.getValue("occupancy", i)))) + str(atoms.getValue("occupancy", i));                 # 55 - 60 OCCUPANCY
                            line += " "*(6 - len(str(atoms.getValue("B_iso_or_equiv", i)))) + str(atoms.getValue("B_iso_or_equiv", i));       # 61 - 66 TEMP FACTOR
                            line += " "*(10);                                                                                                 # 67 - 76
                            line += " "*(2 - len(atoms.getValue("type_symbol", i))) + atoms.getValue("type_symbol", i);                       # 77 - 78 ELEMENT SYMBOL
                            line += " "*2 if(atoms.getValue("pdbx_formal_charge", i) == "?") else atoms.getValue("pdbx_formal_charge", i);    # 79 - 80 CHARGE OF ATOM
                            pdb_arr.append(pdb.HETATM(line));
                             
                        except:
                            print("cif.ATOM_SITE: Error reading line:\n%s" % line);
                            err_arr.append("HETATOM");
                            
                            
            try:
                line  = "ENDMDL";
                pdb_arr.append(pdb.ENDMDL(line));
            except:
                sys.stderr.write("cif.ATOM_SITE: Error reading line:\n%s\n" % line);
                err_arr.append("ENDMDL");
            
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
    
    if(struct_conn is None or atoms is None):
        return pdb_arr, err_arr;

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

def TITLE(block):
    
    title_arr = [];
    title_err = [];
    
    struct_obj = block.getObj("struct");
    
    title_string = struct_obj.getValue("title");
    title_chunk    = int(ceil(len(title_string)/70.0));
    
    for i in range(title_chunk):
        line = "TITLE  ";
        line += " "*(2-len(str(i+1))) + str(i+1) if(i+1>1) else "  ";
        line += title_string[ (i*70) : minimum(len(title_string), (i+1)*70)];
        try:
            title_arr.append(pdb.TITLE(line));
        except:
            sys.stderr.write("cif.TITLE:    Error parsing line:\n%s" % line);
            title_err.append("TITLE");
            
    return title_arr, title_err;

def COMPND(block):
    
    compnd_arr = [];
    compnd_err = [];
    
    entity_obj = block.getObj("entity");
    
    cont = 1;
    for i in range(entity_obj.getRowCount()):
        line1 = "COMPND "
        line1 += " "*(3 - len(str(cont))) + str(cont) if(cont>1) else "   ";
        line1 += "MOL_ID: " + str(entity_obj.getValue("id", i)) + ";";
        try:
            compnd_arr.append(pdb.COMPND(line1));
        except:
            sys.stderr.write("cif.COMPND:    Error parsing line:\n%s\n" % line1);
            compnd_err.append("COMPND");
        
        cont += 1;
        
        line2 = "COMPND "
        line2 += " "*(3 - len(str(cont))) + str(cont) if(cont>1) else "   ";
        line2 += "MOLECULE: " + entity_obj.getValue("pdbx_description", i) + ";";
        try:
            compnd_arr.append(pdb.COMPND(line2));
        except:
            sys.stderr.write("cif.COMPND:    Error parsing line:\n%s\n" % line2);
            compnd_err.append("COMPND");
        
        cont += 1;
        
    return compnd_arr, compnd_err;

def SOURCE(block):
    
    src_arr = [];
    src_err = [];
    
    src_obj = block.getObj("entity_src_gen");
    
    if(src_obj is None):
        return src_arr, src_err;
    
    cont = 1;
    for i in range(src_obj.getRowCount()):
        
        if(src_obj.getValue("entity_id", 0) != "?"):
            line =  "SOURCE ";
            line += " "*(3 - len(str(cont))) + str(cont) if(cont>1) else "   ";
            line += "MOL_ID: " + str(src_obj.getValue("entity_id", i)) + ";";
            cont += 1;
            try:
                src_arr.append(pdb.SOURCE(line));
            except:
                sys.stderr.write("cif.SOURCE:    Error parsing line:\n%s\n" % line);
                src_err.append("SOURCE");
                
        if(src_obj.getValue("pdbx_gene_src_scientific_name", i) != "?"):
            line = "SOURCE ";
            line += " "*(3 - len(str(cont))) + str(cont) if(cont>1) else "   ";
            line += "ORGANISIM_SCIENTIFIC: " + src_obj.getValue("pdbx_gene_src_scientific_name", i) + ";";
            cont += 1;
            try:
                src_arr.append(pdb.SOURCE(line));
            except:
                sys.stderr.write("cif.SOURCE:    Error parsing line:\n%s\n" % line);
                src_err.append("SOURCE");
                
        if(src_obj.getValue("gene_src_common_name", i) != "?"):
            line = "SOURCE ";
            line += " "*(3 - len(str(cont))) + str(cont) if(cont>1) else "   ";
            line += "ORGANISM_COMMON: " + src_obj.getValue("gene_src_common_name", i) + ";";
            cont += 1;
            try:
                src_arr.append(pdb.SOURCE(line));
            except:
                sys.stderr.write("cif.SOURCE:    Error parsing line:\n%s\n" % line);
                src_err.append("SOURCE");
                
        if(src_obj.getValue("pdbx_gene_src_ncbi_taxonomy_id", i) != "?"):
            line = "SOURCE ";
            line += " "*(3 - len(str(cont))) + str(cont) if(cont>1) else "   ";
            line += "ORGANISM_TAXID: " + src_obj.getValue("pdbx_gene_src_ncbi_taxonomy_id", i) + ";";
            cont += 1;
            try:
                src_arr.append(pdb.SOURCE(line));
            except:
                sys.stderr.write("cif.SOURCE:    Error parsing line:\n%s\n" % line);
                src_err.append("SOURCE");
    
    return src_arr, src_err;

def KEYWDS(block):
    
    key_arr = [];
    key_err = [];
    
    key_obj = block.getObj("struct_keywords");
    
    key_string = key_obj.getValue("text");
    key_chunk    = int(ceil(len(key_string)/69.0));
    
    
    for i in range(key_chunk):
        line = "KEYWDS  ";
        line += " "*(2-len(str(i+1))) + str(i+1) if(i+1>1) else "  ";
        line += key_string[ (i*69) : minimum(len(key_string), (i+1)*69)];
        try:
            key_arr.append(pdb.KEYWDS(line));
        except:
            sys.stderr.write("cif.KEYWDS:    Error parsing line:\n%s" % line);
            key_err.append("KEYWDS");
            
    return key_arr, key_err;
 
def EXPDTA(block):
    
    ex_arr = [];
    ex_err = [];
    
    ex_obj = block.getObj("exptl");
    
    line = "EXPDTA  ";
    line += "  ";
    line += ex_obj.getValue("method", 0);
    
    try:
        ex_arr.append(pdb.EXPDTA(line));
    except:
        sys.stderr.write("cif.EXPDTA:    Error parsing line:\n%s\n" % line);
        ex_err.append("EXPDTA");
        
    return ex_arr, ex_err;

def AUTHOR(block):
    
    aut_arr = [];
    aut_err = [];
    
    aut_obj = block.getObj("audit_author");
    
    for i in range(aut_obj.getRowCount()):
        line =  "AUTHOR  ";
        line += "  "*(2 - len(str(aut_obj.getValue("pdbx_ordinal", i)))) + str(aut_obj.getValue("pdbx_ordinal", i));
        line += aut_obj.getValue("name", i);
        
        try:
            aut_arr.append(pdb.AUTHOR(line));
        except:
            sys.stderr.write("cif.AUTHOR:    Error parsing line:\n%s\n" % line);
            aut_err.append("AUTHOR");
            
    return aut_arr, aut_err;

def CRYST1(block):
    
    cry_arr = [];
    cry_err = [];
    
    cry_obj = block.getObj("cell");
    sym_obj = block.getObj("symmetry");
    
    line = "CRYST1";
    line += " "*(9 - len(str(cry_obj.getValue("length_a", 0)))) + cry_obj.getValue("length_a", 0);
    line += " "*(9 - len(str(cry_obj.getValue("length_b", 0)))) + cry_obj.getValue("length_b", 0);
    line += " "*(9 - len(str(cry_obj.getValue("length_c", 0)))) + cry_obj.getValue("length_c", 0); 
    line += " "*(7 - len(str(cry_obj.getValue("angle_alpha", 0)))) + cry_obj.getValue("angle_alpha", 0);
    line += " "*(7 - len(str(cry_obj.getValue("angle_beta", 0)))) + cry_obj.getValue("angle_beta", 0);
    line += " "*(7 - len(str(cry_obj.getValue("angle_gamma", 0)))) + cry_obj.getValue("angle_gamma", 0);
    line += " "*(11 - len(str(sym_obj.getValue("space_group_name_H-M", 0)))) + sym_obj.getValue("space_group_name_H-M", 0);
    line += " "*(4 - len(str(cry_obj.getValue("Z_PDB", 0)))) + cry_obj.getValue("Z_PDB", 0);
    
    try:
        cry_arr.append(pdb.CRYST1(line));
    except:
        sys.stderr.write("cif.CRYST1:    Error parsing line:\n%s\n" % line);
        cry_err.append(CRYST1);
        
    return cry_arr, cry_err;

def SCALEn(block):
    
    sc_arr = [];
    sc_err = [];
    
    sc_obj = block.getObj("atom_sites");
    
    scale1 =  "";
    scale1 += "SCALE1    ";
    scale1 += " "*(10 - len(str(sc_obj.getValue("fract_transf_matrix[1][1]", 0)))) + str(sc_obj.getValue("fract_transf_matrix[1][1]", 0));
    scale1 += " "*(10 - len(str(sc_obj.getValue("fract_transf_matrix[1][2]", 0)))) + str(sc_obj.getValue("fract_transf_matrix[1][2]", 0));
    scale1 += " "*(10 - len(str(sc_obj.getValue("fract_transf_matrix[1][3]", 0)))) + str(sc_obj.getValue("fract_transf_matrix[1][3]", 0));
    scale1 += "     ";
    scale1 += " "*(10 - len(str(sc_obj.getValue("fract_transf_vector[1]", 0)))) + str(sc_obj.getValue("fract_transf_vector[1]", 0));
    
    scale2 =  "";
    scale2 += "SCALE2    ";
    scale2 += " "*(10 - len(str(sc_obj.getValue("fract_transf_matrix[2][1]", 0)))) + str(sc_obj.getValue("fract_transf_matrix[2][1]", 0));
    scale2 += " "*(10 - len(str(sc_obj.getValue("fract_transf_matrix[2][2]", 0)))) + str(sc_obj.getValue("fract_transf_matrix[2][2]", 0));
    scale2 += " "*(10 - len(str(sc_obj.getValue("fract_transf_matrix[2][3]", 0)))) + str(sc_obj.getValue("fract_transf_matrix[2][3]", 0));
    scale2 += "     ";
    scale2 += " "*(10 - len(str(sc_obj.getValue("fract_transf_vector[2]", 0)))) + str(sc_obj.getValue("fract_transf_vector[2]", 0));
    
    scale3 =  "";
    scale3 += "SCALE3    ";
    scale3 += " "*(10 - len(str(sc_obj.getValue("fract_transf_matrix[3][1]", 0)))) + str(sc_obj.getValue("fract_transf_matrix[3][1]", 0));
    scale3 += " "*(10 - len(str(sc_obj.getValue("fract_transf_matrix[3][2]", 0)))) + str(sc_obj.getValue("fract_transf_matrix[3][2]", 0));
    scale3 += " "*(10 - len(str(sc_obj.getValue("fract_transf_matrix[3][3]", 0)))) + str(sc_obj.getValue("fract_transf_matrix[3][3]", 0));
    scale3 += "     ";
    scale3 += " "*(10 - len(str(sc_obj.getValue("fract_transf_vector[3]", 0)))) + str(sc_obj.getValue("fract_transf_vector[3]", 0));
    
    try:
        sc_arr.append(pdb.SCALE1(scale1));
    except:
        sys.stderr.write("cif.SCALEn:    Error parsing line:\n%s\n" % scale1);
        sc_err.append("SCALE1");
        
    try:
        sc_arr.append(pdb.SCALE2(scale2));
    except:
        sys.stderr.write("cif.SCALEn:    Error parsing line:\n%s\n" % scale2);
        sc_err.append("SCALE2");
    
    try:
        sc_arr.append(pdb.SCALE3(scale3));
    except:
        sys.stderr.write("cif.SCALEn:    Error parsing line:\n%s\n" % scale3);
        sc_err.append("SCALE3");
        
    return sc_arr, sc_err;

def ORIGXn(block):
    
    or_arr = [];
    or_err = [];
    
    or_obj = block.getObj("database_PDB_matrix");
    
    orig1 =  "ORIGX1    ";
    orig1 += " "*(10 - len(str(or_obj.getValue("origx[1][1]", 0)))) + str(or_obj.getValue("origx[1][1]", 0));
    orig1 += " "*(10 - len(str(or_obj.getValue("origx[1][2]", 0)))) + str(or_obj.getValue("origx[1][2]", 0));
    orig1 += " "*(10 - len(str(or_obj.getValue("origx[1][3]", 0)))) + str(or_obj.getValue("origx[1][3]", 0));
    orig1 += "     ";
    orig1 += " "*(10 - len(str(or_obj.getValue("origx_vector[1]", 0)))) + str(or_obj.getValue("origx_vector[1]", 0));
    
    orig2 =  "ORIGX2    ";
    orig2 += " "*(10 - len(str(or_obj.getValue("origx[2][1]", 0)))) + str(or_obj.getValue("origx[2][1]", 0));
    orig2 += " "*(10 - len(str(or_obj.getValue("origx[2][2]", 0)))) + str(or_obj.getValue("origx[2][2]", 0));
    orig2 += " "*(10 - len(str(or_obj.getValue("origx[2][3]", 0)))) + str(or_obj.getValue("origx[2][3]", 0));
    orig2 += "     ";
    orig2 += " "*(10 - len(str(or_obj.getValue("origx_vector[2]", 0)))) + str(or_obj.getValue("origx_vector[2]", 0));
    
    orig3 =  "ORIGX3    ";
    orig3 += " "*(10 - len(str(or_obj.getValue("origx[3][1]", 0)))) + str(or_obj.getValue("origx[3][1]", 0));
    orig3 += " "*(10 - len(str(or_obj.getValue("origx[3][2]", 0)))) + str(or_obj.getValue("origx[3][2]", 0));
    orig3 += " "*(10 - len(str(or_obj.getValue("origx[3][3]", 0)))) + str(or_obj.getValue("origx[3][3]", 0));
    orig3 += "     ";
    orig3 += " "*(10 - len(str(or_obj.getValue("origx_vector[3]", 0)))) + str(or_obj.getValue("origx_vector[3]", 0));
    
    try:
        or_arr.append(pdb.ORIGX1(orig1));
    except:
        sys.stderr.write("cif.ORIGXn:    Error parsing line:\n%s\n" % orig1);
        or_err.append("ORIGX1");
        
    try:
        or_arr.append(pdb.ORIGX2(orig2));
    except:
        sys.stderr.write("cif.ORIGXn:    Error parsing line:\n%s\n" % orig2);
        or_err.append("ORIGX2");
    
    try:
        or_arr.append(pdb.ORIGX3(orig3));
    except:
        sys.stderr.write("cif.ORIGXn:    Error parsing line:\n%s\n" % orig3);
        or_err.append("ORIGX3");
        
    return or_arr, or_err;

def CISPEP(block):
    
    cis_arr = [];
    cis_err = [];
    
    cis_obj = block.getObj("struct_mon_prot_cis");
    if(cis_obj is None):
        return cis_arr, cis_err;
    
    
    for i in range(cis_obj.getRowCount()):
        line =  "CISPEP ";
        line += " "*(3 - len(str(cis_obj.getValue("pdbx_id", i)))) + str(cis_obj.getValue("pdbx_id", i));
        line += " ";
        line += " "*(3 - len(cis_obj.getValue("auth_comp_id", i))) + cis_obj.getValue("auth_comp_id", i);
        line += " ";
        line += cis_obj.getValue("auth_asym_id", i);
        line += " ";
        line += " "*(4 - len(str(cis_obj.getValue("auth_seq_id", i)))) + str(cis_obj.getValue("auth_seq_id", i));
        line += cis_obj.getValue("pdbx_PDB_ins_code", i) if(cis_obj.getValue("pdbx_PDB_ins_code", i) != '?') else " ";
        line += "   ";
        line += " "*(3 - len(cis_obj.getValue("pdbx_auth_comp_id_2", i))) + cis_obj.getValue("pdbx_auth_comp_id_2", i);
        line += " ";
        line +=  cis_obj.getValue("pdbx_auth_asym_id_2", i);
        line += " ";
        line += " "*(4 - len(str(cis_obj.getValue("pdbx_auth_seq_id_2", i)))) + str(cis_obj.getValue("pdbx_auth_seq_id_2", i));
        line += cis_obj.getValue("pdbx_PDB_ins_code_2", i) if(cis_obj.getValue("pdbx_PDB_ins_code_2", i) != '?') else " ";
        line += " "*7;
        line += " "*(3 - len(str(cis_obj.getValue("pdbx_PDB_model_num", i)))) + str(cis_obj.getValue("pdbx_PDB_model_num", i));
        line += " "*7
        line += " "*(6 - len(str(cis_obj.getValue("pdbx_omega_angle", i)))) + str(cis_obj.getValue("pdbx_omega_angle", i));
        
        try:
            cis_arr.append(pdb.CISPEP(line));
        except:
            sys.stderr.write("cif.CISPEP:    Erro parsing line:\n%s\n" % line);
            cis_err.append("CISPEP");
        
        
    return cis_arr, cis_err; 

def SSBOND(block):
    
    ssb_arr = [];
    ssb_err = [];
    
    ssb_obj = block.getObj("struct_conn");
    
    if(ssb_obj is None):
        return ssb_arr, ssb_err;
    
    for i in range(ssb_obj.getRowCount()):
        line  = "SSBOND ";
        line += " "*(3 - len(str(ssb_obj.getValue("id", i)[-1]))) + str(ssb_obj.getValue("id", i)[-1]);
        line += " ";
        line += " "*(3 - len(ssb_obj.getValue("ptnr1_auth_comp_id", i))) + ssb_obj.getValue("ptnr1_auth_comp_id", i);
        line += " ";
        line += ssb_obj.getValue("ptnr1_auth_asym_id", i);
        line += " ";
        line += " "*(4 - len(str(ssb_obj.getValue("ptnr1_auth_seq_id", i)))) + str(ssb_obj.getValue("ptnr1_auth_seq_id", i));
        line += ssb_obj.getValue("pdbx_ptnr1_PDB_ins_code", i) if (ssb_obj.getValue("pdbx_ptnr1_PDB_ins_code", i) != '?') else " ";
        line += " "*3;
        line += " "*(3 - len(ssb_obj.getValue("ptnr2_auth_comp_id", i))) + ssb_obj.getValue("ptnr2_auth_comp_id", i);
        line += " ";
        line += ssb_obj.getValue("ptnr2_auth_asym_id", i);
        line += " ";
        line += " "*(4 - len(str(ssb_obj.getValue("ptnr2_auth_seq_id", i)))) + str(ssb_obj.getValue("ptnr2_auth_seq_id", i));
        line += ssb_obj.getValue("pdbx_ptnr2_PDB_ins_code", i) if(ssb_obj.getValue("pdbx_ptnr2_PDB_ins_code", i) != '?') else " ";
        line += " "*23;
        line += " "*(6 - len(ssb_obj.getValue("ptnr1_symmetry", i).replace("_", ""))) + ssb_obj.getValue("ptnr1_symmetry", i).replace("_", "");
        line += " ";
        line += " "*(6 - len(ssb_obj.getValue("ptnr2_symmetry", i).replace("_", ""))) + ssb_obj.getValue("ptnr2_symmetry", i).replace("_", "");
        line += " ";
        line += " "*(5 - len(str(ssb_obj.getValue("pdbx_dist_value", i)))) + str(ssb_obj.getValue("pdbx_dist_value", i));
        
        try:
            ssb_arr.append(pdb.SSBOND(line));
        except:
            sys.stderr.write("cif.SSBOND:    Error parsing line:\n%s\n" % line);
            ssb_err.append("SSBOND");   
        
        
    return ssb_arr, ssb_err;

def CountModels(block):
    
    atom_obj = block.getObj("atom_site");
    
    modelNum = [];
    for i in range(atom_obj.getRowCount()):
        tmp = atom_obj.getValue("pdbx_PDB_model_num", i);
        if(tmp not in modelNum):
            modelNum.append(tmp);
        else:
            pass;
        
    return modelNum;
    

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
            head_pdb, head_err   = HEADER(block);
            title_pdb, title_err = TITLE(block);
            cmpnd_pdb, cmpnd_err = COMPND(block);
            src_pdb, src_err     = SOURCE(block);
            key_pdb, key_err     = KEYWDS(block);
            ex_pdb, ex_err       = EXPDTA(block);
            aut_pdb, aut_err     = AUTHOR(block);
            ssb_pdb, ssb_err     = SSBOND(block);
            cis_pdb, cis_err     = CISPEP(block);
            cry_pdb, cry_err     = CRYST1(block);
            or_pdb, or_err       = ORIGXn(block);
            sc_pdb, sc_err       = SCALEn(block);
            ato_pdb, ato_err     = ATOM_SITE(block);
            con_pdb, con_err     = CONECT(block);
            

            pdblist = head_pdb + title_pdb + cmpnd_pdb + src_pdb + key_pdb + ex_pdb + aut_pdb + \
                    ssb_pdb + cis_pdb + cry_pdb + or_pdb + sc_pdb + ato_pdb + con_pdb;
            errlist = head_err + title_err + cmpnd_err + src_err + key_err + ex_err + aut_err + \
                    ssb_err + cis_err + cry_err + or_err + sc_err + ato_err + con_err;

        return pdblist, errlist

    else:
        print("Unknown error while reading cif file.");
        return pdblist, errlist;
