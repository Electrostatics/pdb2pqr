"""PDB parsing class

This module parses PDBs in accordance to PDB Format Description Version 2.2
(1996); it is not very forgiving.  Each class in this module corresponds
to a record in the PDB Format Description.  Much of the documentation for
the classes is taken directly from the above PDB Format Description.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
.. codeauthor::  Nathan Baker
"""

import logging

_LOGGER = logging.getLogger(__name__)


LINE_PARSERS = {}


def register_line_parser(klass):
    """Register a line parser in the global dictionary.

    :param klass:  class for line parser
    """
    LINE_PARSERS[klass.__name__] = klass
    return klass


class BaseRecord:
    """Base class for all records.

    Verifies the received record type.
    """

    def __init__(self, line):
        record = line[0:6].strip()
        if record != self.__class__.__name__:
            raise ValueError(record)
        self.original_text = line.rstrip("\r\n")

    def __str__(self):
        return self.original_text

    def record_type(self):
        """Return PDB record type as string.

        :return:  record type
        :rtype:  str
        """
        return self.original_text.split()[0]


@register_line_parser
class END(BaseRecord):
    """END class

    The END records are paired with MODEL records to group individual
    structures found in a coordinate entry.
    """

    def __init__(self, line):
        """Initialize with line.

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)


@register_line_parser
class MASTER(BaseRecord):
    """MASTER class

    The MASTER record is a control record for bookkeeping.
    It lists the number of lines in the coordinate entry or file for selected
    record types.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+------+------------+-------------------------------------+
        | COLUMNS | TYPE | FIELD      | DEFINITION                          |
        +=========+======+============+=====================================+
        | 11-15   | int  | num_remark | Number of REMARK records            |
        +---------+------+------------+-------------------------------------+
        | 21-25   | int  | num_het    | Number of HET records               |
        +---------+------+------------+-------------------------------------+
        | 26-30   | int  | numHelix   | Number of HELIX records             |
        +---------+------+------------+-------------------------------------+
        | 31-35   | int  | numSheet   | Number of SHEET records             |
        +---------+------+------------+-------------------------------------+
        | 36-40   | int  | numTurn    | Number of TURN records              |
        +---------+------+------------+-------------------------------------+
        | 41-45   | int  | numSite    | Number of SITE records              |
        +---------+------+------------+-------------------------------------+
        | 46-50   | int  | numXform   | Number of coordinate transformation |
        |         |      |            | records (ORIGX+SCALE+MTRIX)         |
        +---------+------+------------+-------------------------------------+
        | 51-55   | int  | numCoord   | Number of atomic coordinate records |
        |         |      |            | (ATOM+HETATM)                       |
        +---------+------+------------+-------------------------------------+
        | 56-60   | int  | numTer     | Number of TER records               |
        +---------+------+------------+-------------------------------------+
        | 61-65   | int  | numConect  | Number of CONECT records            |
        +---------+------+------------+-------------------------------------+
        | 66-70   | int  | numSeq     | Number of SEQRES records            |
        +---------+------+------------+-------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.num_remark = int(line[10:15].strip())
        self.num_het = int(line[20:25].strip())
        self.num_helix = int(line[25:30].strip())
        self.num_sheet = int(line[30:35].strip())
        self.num_turn = int(line[35:40].strip())
        self.num_site = int(line[40:45].strip())
        self.num_xform = int(line[45:50].strip())
        self.num_coord = int(line[50:55].strip())
        self.num_ter = int(line[55:60].strip())
        self.num_conect = int(line[60:65].strip())
        self.num_seq = int(line[65:70].strip())


@register_line_parser
class CONECT(BaseRecord):
    """CONECT class

    The CONECT records specify connectivity between atoms for which
    coordinates are supplied. The connectivity is described using the atom
    serial number as found in the entry. CONECT records are mandatory for
    HET groups (excluding water) and for other bonds not specified in the
    standard residue connectivity table which involve atoms in standard
    residues (see Appendix 4 for the list of standard residues). These
    records are generated by the PDB.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+------+----------+---------------------------------------+
        | COLUMNS | TYPE | FIELD    | DEFINITION                            |
        +=========+======+==========+=======================================+
        | 7-11    | int  | serial   | Atom serial number                    |
        +---------+------+----------+---------------------------------------+
        | 12-16   | int  | serial1  | Serial number of bonded atom          |
        +---------+------+----------+---------------------------------------+
        | 17-21   | int  | serial2  | Serial number of bonded atom          |
        +---------+------+----------+---------------------------------------+
        | 22-26   | int  | serial3  | Serial number of bonded atom          |
        +---------+------+----------+---------------------------------------+
        | 27-31   | int  | serial4  | Serial number of bonded atom          |
        +---------+------+----------+---------------------------------------+
        | 32-36   | int  | serial5  | Serial number of hydrogen bonded atom |
        +---------+------+----------+---------------------------------------+
        | 37-41   | int  | serial6  | Serial number of hydrogen bonded atom |
        +---------+------+----------+---------------------------------------+
        | 42-46   | int  | serial7  | Serial number of salt bridged atom    |
        +---------+------+----------+---------------------------------------+
        | 47-51   | int  | serial8  | Serial number of hydrogen bonded atom |
        +---------+------+----------+---------------------------------------+
        | 52-56   | int  | serial9  | Serial number of hydrogen bonded atom |
        +---------+------+----------+---------------------------------------+
        | 57-61   | int  | serial10 | Serial number of salt bridged atom    |
        +---------+------+----------+---------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.serial = int(line[6:11].strip())
        try:
            self.serial1 = int(line[11:16].strip())
        except ValueError:
            self.serial1 = None
        try:
            self.serial2 = int(line[16:21].strip())
        except ValueError:
            self.serial2 = None
        try:
            self.serial3 = int(line[21:26].strip())
        except ValueError:
            self.serial3 = None
        try:
            self.serial4 = int(line[26:31].strip())
        except ValueError:
            self.serial4 = None
        try:
            self.serial5 = int(line[31:36].strip())
        except ValueError:
            self.serial5 = None
        try:
            self.serial6 = int(line[36:41].strip())
        except ValueError:
            self.serial6 = None
        try:
            self.serial7 = int(line[41:46].strip())
        except ValueError:
            self.serial7 = None
        try:
            self.serial8 = int(line[46:51].strip())
        except ValueError:
            self.serial8 = None
        try:
            self.serial9 = int(line[51:56].strip())
        except ValueError:
            self.serial9 = None
        try:
            self.serial10 = int(line[56:61].strip())
        except ValueError:
            self.serial10 = None


@register_line_parser
class NUMMDL(BaseRecord):
    """NUMMDL class

    The NUMMDL record indicates total number of models in a PDB entry.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+------+-------------+-------------------+
        | COLUMNS | TYPE | FIELD       | DEFINITION        |
        +=========+======+=============+===================+
        | 11-14   | int  | modelNumber | Number of models. |
        +---------+------+-------------+-------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        try:
            self.model_number = int(line[10:14].strip())
        except ValueError:
            self.model_number = None


@register_line_parser
class ENDMDL(BaseRecord):
    """ENDMDL class

    The ENDMDL records are paired with MODEL records to group individual
    structures found in a coordinate entry.
    """

    def __init__(self, line):
        super().__init__(line)


@register_line_parser
class TER(BaseRecord):
    """TER class

    The TER record indicates the end of a list of ATOM/HETATM records for a
    chain.
    """

    def __init__(self, line):
        """Initialize by parsing line:

        +---------+--------+----------+--------------------------+
        | COLUMNS | TYPE   | FIELD    | DEFINITION               |
        +=========+========+==========+==========================+
        | 7-11    | int    | serial   | Serial number.           |
        +---------+--------+----------+--------------------------+
        | 18-20   | string | res_name | Residue name.            |
        +---------+--------+----------+--------------------------+
        | 22      | string | chain_id | Chain identifier.        |
        +---------+--------+----------+--------------------------+
        | 23-26   | int    | res_seq  | Residue sequence number. |
        +---------+--------+----------+--------------------------+
        | 27      | string | ins_code | Insertion code.          |
        +---------+--------+----------+--------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        try:  # Not really needed
            self.serial = int(line[6:11].strip())
            self.res_name = line[17:20].strip()
            self.chain_id = line[21].strip()
            self.res_seq = int(line[22:26].strip())
            self.ins_code = line[26].strip()
        except (IndexError, ValueError):
            self.serial = None
            self.res_name = None
            self.chain_id = None
            self.res_seq = None
            self.ins_code = None


@register_line_parser
class SIGUIJ(BaseRecord):
    """SIGUIJ class

    The SIGUIJ records present the anisotropic temperature factors.
    """

    def __init__(self, line):
        """Initialize by parsing line:

        +---------+--------+----------+-------------------------------------+
        | COLUMNS | TYPE   | FIELD    | DEFINITION                          |
        +=========+========+==========+=====================================+
        | 7-11    | int    | serial   | Atom serial number.                 |
        +---------+--------+----------+-------------------------------------+
        | 13-16   | string | name     | Atom name.                          |
        +---------+--------+----------+-------------------------------------+
        | 17      | string | alt_loc  | Alternate location indicator.       |
        +---------+--------+----------+-------------------------------------+
        | 18-20   | string | res_name | Residue name.                       |
        +---------+--------+----------+-------------------------------------+
        | 22      | string | chain_id | Chain identifier.                   |
        +---------+--------+----------+-------------------------------------+
        | 23-26   | int    | res_seq  | Residue sequence number.            |
        +---------+--------+----------+-------------------------------------+
        | 27      | string | ins_code | Insertion code.                     |
        +---------+--------+----------+-------------------------------------+
        | 29-35   | int    | sig11    | Sigma U(1,1)                        |
        +---------+--------+----------+-------------------------------------+
        | 36-42   | int    | sig22    | Sigma U(2,2)                        |
        +---------+--------+----------+-------------------------------------+
        | 43-49   | int    | sig33    | Sigma U(3,3)                        |
        +---------+--------+----------+-------------------------------------+
        | 50-56   | int    | sig12    | Sigma U(1,2)                        |
        +---------+--------+----------+-------------------------------------+
        | 57-63   | int    | sig13    | Sigma U(1,3)                        |
        +---------+--------+----------+-------------------------------------+
        | 64-70   | int    | sig23    | Sigma U(2,3)                        |
        +---------+--------+----------+-------------------------------------+
        | 73-76   | string | seg_id   | Segment identifier, left-justified. |
        +---------+--------+----------+-------------------------------------+
        | 77-78   | string | el.ment  | Element symbol, right-justified.    |
        +---------+--------+----------+-------------------------------------+
        | 79-80   | string | charge   | Charge on the atom.                 |
        +---------+--------+----------+-------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.serial = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.alt_loc = line[16].strip()
        self.res_name = line[17:20].strip()
        self.chain_id = line[21].strip()
        self.res_seq = int(line[22:26].strip())
        self.ins_code = line[26].strip()
        self.sig11 = int(line[28:35].strip())
        self.sig22 = int(line[35:42].strip())
        self.sig33 = int(line[42:49].strip())
        self.sig12 = int(line[49:56].strip())
        self.sig13 = int(line[56:63].strip())
        self.sig23 = int(line[63:70].strip())
        self.seg_id = line[72:76].strip()
        self.element = line[76:78].strip()
        self.charge = line[78:80].strip()


@register_line_parser
class ANISOU(BaseRecord):
    """ANISOU class

    The ANISOU records present the anisotropic temperature factors.
    """

    def __init__(self, line):
        """Initialize by parsing line:

        +---------+--------+----------+-------------------------------------+
        | COLUMNS | TYPE   | FIELD    | DEFINITION                          |
        +=========+========+==========+=====================================+
        | 7-11    | int    | serial   | Atom serial number.                 |
        +---------+--------+----------+-------------------------------------+
        | 13-16   | string | name     | Atom name.                          |
        +---------+--------+----------+-------------------------------------+
        | 17      | string | alt_loc  | Alternate location indicator.       |
        +---------+--------+----------+-------------------------------------+
        | 18-20   | string | res_name | Residue name.                       |
        +---------+--------+----------+-------------------------------------+
        | 22      | string | chain_id | Chain identifier.                   |
        +---------+--------+----------+-------------------------------------+
        | 23-26   | int    | res_seq  | Residue sequence number.            |
        +---------+--------+----------+-------------------------------------+
        | 27      | string | ins_code | Insertion code.                     |
        +---------+--------+----------+-------------------------------------+
        | 29-35   | int    | u00      | U(1,1)                              |
        +---------+--------+----------+-------------------------------------+
        | 36-42   | int    | u11      | U(2,2)                              |
        +---------+--------+----------+-------------------------------------+
        | 43-49   | int    | u22      | U(3,3)                              |
        +---------+--------+----------+-------------------------------------+
        | 50-56   | int    | u01      | U(1,2)                              |
        +---------+--------+----------+-------------------------------------+
        | 57-63   | int    | u02      | U(1,3)                              |
        +---------+--------+----------+-------------------------------------+
        | 64-70   | int    | u12      | U(2,3)                              |
        +---------+--------+----------+-------------------------------------+
        | 73-76   | string | seg_id   | Segment identifier, left-justified. |
        +---------+--------+----------+-------------------------------------+
        | 77-78   | string | element  | Element symbol, right-justified.    |
        +---------+--------+----------+-------------------------------------+
        | 79-80   | string | charge   | Charge on the atom.                 |
        +---------+--------+----------+-------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.serial = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.alt_loc = line[16].strip()
        self.res_name = line[17:20].strip()
        self.chain_id = line[21].strip()
        self.res_seq = int(line[22:26].strip())
        self.ins_code = line[26].strip()
        self.u00 = int(line[28:35].strip())
        self.u11 = int(line[35:42].strip())
        self.u22 = int(line[42:49].strip())
        self.u01 = int(line[49:56].strip())
        self.u02 = int(line[56:63].strip())
        self.u12 = int(line[63:70].strip())
        self.seg_id = line[72:76].strip()
        self.element = line[76:78].strip()
        self.charge = line[78:80].strip()


@register_line_parser
class SIGATM(BaseRecord):
    """SIGATM class

    The SIGATM records present the standard deviation of atomic parameters
    as they appear in ATOM and HETATM records.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+----------+-------------------------------------+
        | COLUMNS | TYPE   | FIELD    | DEFINITION                          |
        +=========+========+==========+=====================================+
        | 7-11    | int    | serial   | Atom serial number.                 |
        +---------+--------+----------+-------------------------------------+
        | 13-16   | string | name     | Atom name.                          |
        +---------+--------+----------+-------------------------------------+
        | 17      | string | alt_loc  | Alternate location indicator.       |
        +---------+--------+----------+-------------------------------------+
        | 18-20   | string | res_name | Residue name.                       |
        +---------+--------+----------+-------------------------------------+
        | 22      | string | chain_id | Chain identifier.                   |
        +---------+--------+----------+-------------------------------------+
        | 23-26   | int    | res_seq  | Residue sequence number.            |
        +---------+--------+----------+-------------------------------------+
        | 27      | string | ins_code | Code for insertion of residues.     |
        +---------+--------+----------+-------------------------------------+
        | 31-38   | float  | sig_x    | Standard deviation of orthogonal    |
        |         |        |          | coordinates for X in Angstroms.     |
        +---------+--------+----------+-------------------------------------+
        | 39-46   | float  | sig_y    | Standard deviation of orthogonal    |
        |         |        |          | coordinates for Y in Angstroms.     |
        +---------+--------+----------+-------------------------------------+
        | 47-54   | float  | sig_z    | Standard deviation of orthogonal    |
        |         |        |          | coordinates for Z in Angstroms.     |
        +---------+--------+----------+-------------------------------------+
        | 55-60   | float  | sig_occ  | Standard deviation of occupancy.    |
        +---------+--------+----------+-------------------------------------+
        | 61-66   | float  | sig_temp | Standard deviation of temperature   |
        |         |        |          | factor.                             |
        +---------+--------+----------+-------------------------------------+
        | 73-76   | string | seg_id   | Segment identifier, left-justified. |
        +---------+--------+----------+-------------------------------------+
        | 77-78   | string | element  | Element symbol, right-justified.    |
        +---------+--------+----------+-------------------------------------+
        | 79-80   | string | charge   | Charge on the atom.                 |
        +---------+--------+----------+-------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.serial = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.alt_loc = line[16].strip()
        self.res_name = line[17:20].strip()
        self.chain_id = line[21].strip()
        self.res_seq = int(line[22:26].strip())
        self.ins_code = line[26].strip()
        self.sig_x = float(line[30:38].strip())
        self.sig_y = float(line[38:46].strip())
        self.sig_z = float(line[46:54].strip())
        self.sig_occ = float(line[54:60].strip())
        self.sig_temp = float(line[60:66].strip())
        self.seg_id = line[72:76].strip()
        self.element = line[76:78].strip()
        self.charge = line[78:80].strip()


@register_line_parser
class HETATM(BaseRecord):
    """HETATM class

    The HETATM records present the atomic coordinate records for atoms
    within "non-standard" groups. These records are used for water
    molecules and atoms presented in HET groups.
    """

    def __init__(
        self, line, sybyl_type="A.aaa", l_bonds=[], l_bonded_atoms=[]
    ):
        """Initialize by parsing line

        +---------+--------+-------------+-----------------------------------+
        | COLUMNS | TYPE   | FIELD       | DEFINITION                        |
        +=========+========+=============+===================================+
        | 7-11    | int    | serial      | Atom serial number.               |
        +---------+--------+-------------+-----------------------------------+
        | 13-16   | string | name        | Atom name.                        |
        +---------+--------+-------------+-----------------------------------+
        | 17      | string | alt_loc     | Alternate location indicator.     |
        +---------+--------+-------------+-----------------------------------+
        | 18-20   | string | res_name    | Residue name.                     |
        +---------+--------+-------------+-----------------------------------+
        | 22      | string | chain_id    | Chain identifier.                 |
        +---------+--------+-------------+-----------------------------------+
        | 23-26   | int    | res_seq     | Residue sequence number.          |
        +---------+--------+-------------+-----------------------------------+
        | 27      | string | ins_code    | Code for insertion of residues.   |
        +---------+--------+-------------+-----------------------------------+
        | 31-38   | float  | x           | Orthogonal coordinates for X in   |
        |         |        |             | Angstroms.                        |
        +---------+--------+-------------+-----------------------------------+
        | 39-46   | float  | y           | Orthogonal coordinates for Y in   |
        |         |        |             | Angstroms.                        |
        +---------+--------+-------------+-----------------------------------+
        | 47-54   | float  | z           | Orthogonal coordinates for Z in   |
        |         |        |             | Angstroms.                        |
        +---------+--------+-------------+-----------------------------------+
        | 55-60   | float  | occupancy   | Occupancy.                        |
        +---------+--------+-------------+-----------------------------------+
        | 61-66   | float  | temp_factor | Temperature factor.               |
        +---------+--------+-------------+-----------------------------------+
        | 73-76   | string | seg_id      | Segment identifier, left-         |
        |         |        |             | justified.                        |
        +---------+--------+-------------+-----------------------------------+
        | 77-78   | string | element     | Element symbol, right-justified.  |
        +---------+--------+-------------+-----------------------------------+
        | 79-80   | string | charge      | Charge on the atom.               |
        +---------+--------+-------------+-----------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.serial = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.alt_loc = line[16].strip()
        try:
            self.res_name = line[17:20].strip()
            self.chain_id = line[21].strip()
            self.res_seq = int(line[22:26].strip())
            self.ins_code = line[26].strip()
        except IndexError as e:
            raise ValueError(
                "Residue name must be less than 4 characters!"
            ) from e
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())

        self.sybyl_type = sybyl_type
        self.l_bonded_atoms = l_bonded_atoms
        self.l_bonds = l_bonds
        self.radius = 1.0
        self.is_c_term = 0
        self.is_n_term = 0
        self.mol2charge = None

        try:
            self.occupancy = float(line[54:60].strip())
            self.temp_factor = float(line[60:66].strip())
            self.seg_id = line[72:76].strip()
            self.element = line[76:78].strip()
            self.charge = line[78:80].strip()
        except (ValueError, IndexError):
            self.occupancy = 0.00
            self.temp_factor = 0.00
            self.seg_id = ""
            self.element = ""
            self.charge = ""


@register_line_parser
class ATOM(BaseRecord):
    """ATOM class

    The ATOM records present the atomic coordinates for standard residues.
    They also present the occupancy and temperature factor for each atom.
    Heterogen coordinates use the HETATM record type. The element symbol is
    always present on each ATOM record; segment identifier and charge are
    optional.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+-------------+-----------------------------------+
        | COLUMNS | TYPE   | FIELD       | DEFINITION                        |
        +=========+========+=============+===================================+
        | 7-11    | int    | serial      | Atom serial number.               |
        +---------+--------+-------------+-----------------------------------+
        | 13-16   | string | name        | Atom name.                        |
        +---------+--------+-------------+-----------------------------------+
        | 17      | string | alt_loc     | Alternate location indicator.     |
        +---------+--------+-------------+-----------------------------------+
        | 18-20   | string | res_name    | Residue name.                     |
        +---------+--------+-------------+-----------------------------------+
        | 22      | string | chain_id    | Chain identifier.                 |
        +---------+--------+-------------+-----------------------------------+
        | 23-26   | int    | res_seq     | Residue sequence number.          |
        +---------+--------+-------------+-----------------------------------+
        | 27      | string | ins_code    | Code for insertion of residues.   |
        +---------+--------+-------------+-----------------------------------+
        | 31-38   | float  | x           | Orthogonal coordinates for X in   |
        |         |        |             | Angstroms.                        |
        +---------+--------+-------------+-----------------------------------+
        | 39-46   | float  | y           | Orthogonal coordinates for Y in   |
        |         |        |             | Angstroms.                        |
        +---------+--------+-------------+-----------------------------------+
        | 47-54   | float  | z           | Orthogonal coordinates for Z in   |
        |         |        |             | Angstroms.                        |
        +---------+--------+-------------+-----------------------------------+
        | 55-60   | float  | occupancy   | Occupancy.                        |
        +---------+--------+-------------+-----------------------------------+
        | 61-66   | float  | temp_factor | Temperature factor.               |
        +---------+--------+-------------+-----------------------------------+
        | 73-76   | string | seg_id      | Segment identifier,               |
        |         |        |             | left-justified.                   |
        +---------+--------+-------------+-----------------------------------+
        | 77-78   | string | element     | Element symbol, right-justified.  |
        +---------+--------+-------------+-----------------------------------+
        | 79-80   | string | charge      | Charge on the atom.               |
        +---------+--------+-------------+-----------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.serial = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.alt_loc = line[16].strip()
        self.res_name = line[17:20].strip()
        self.chain_id = line[21].strip()
        self.res_seq = int(line[22:26].strip())
        self.ins_code = line[26].strip()
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())
        try:
            self.occupancy = float(line[54:60].strip())
            self.temp_factor = float(line[60:66].strip())
            self.seg_id = line[72:76].strip()
            self.element = line[76:78].strip()
            self.charge = line[78:80].strip()
        except (ValueError, IndexError):
            self.occupancy = 0.00
            self.temp_factor = 0.00
            self.seg_id = ""
            self.element = ""
            self.charge = ""


@register_line_parser
class MODEL(BaseRecord):
    """MODEL class

    The MODEL record specifies the model serial number when multiple
    structures are presented in a single coordinate entry, as is often the
    case with structures determined by NMR.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+------+--------+----------------------+
        | COLUMNS | TYPE | FIELD  | DEFINITION           |
        +=========+======+========+======================+
        | 11-14   | int  | serial | Model serial number. |
        +---------+------+--------+----------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.serial = int(line[10:14].strip())


@register_line_parser
class TVECT(BaseRecord):
    """TVECT class

    The TVECT records present the translation vector for infinite
    covalently connected structures.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+--------+----------------------------------+
        | COLUMNS | TYPE   | FIELD  | DEFINITION                       |
        +=========+========+========+==================================+
        | 8-10    | int    | serial | Serial number                    |
        +---------+--------+--------+----------------------------------+
        | 11-20   | float  | t1     | Components of translation vector |
        +---------+--------+--------+----------------------------------+
        | 21-30   | float  | t2     | Components of translation vector |
        +---------+--------+--------+----------------------------------+
        | 31-40   | float  | t2     | Components of translation vector |
        +---------+--------+--------+----------------------------------+
        | 41-70   | string | text   | Comments                         |
        +---------+--------+--------+----------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.serial = int(line[7:10].strip())
        self.trans1 = float(line[10:20].strip())
        self.trans2 = float(line[20:30].strip())
        self.trans3 = float(line[30:40].strip())
        self.text = line[40:70].strip()


class MTRIXn(BaseRecord):
    """MTRIXn baseclass

    The MTRIXn (n = 1, 2, or 3) records present transformations expressing
    non-crystallographic symmetry.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+-------+---------+----------------------------------------+
        | COLUMNS | TYPE  | FIELD   | DEFINITION                             |
        +=========+=======+=========+========================================+
        | 8-10    | int   | serial  | Serial number                          |
        +---------+-------+---------+----------------------------------------+
        | 11-20   | float | mn1     | M31                                    |
        +---------+-------+---------+----------------------------------------+
        | 21-30   | float | mn2     | M32                                    |
        +---------+-------+---------+----------------------------------------+
        | 31-40   | float | mn3     | M33                                    |
        +---------+-------+---------+----------------------------------------+
        | 46-55   | float | vn      | V3                                     |
        +---------+-------+---------+----------------------------------------+
        | 60      | int   | i_given | 1 if coordinates for the               |
        |         |       |         | representations which are approximately|
        |         |       |         | related by the transformations of the  |
        |         |       |         | molecule are contained in the entry.   |
        |         |       |         | Otherwise, blank.                      |
        +---------+-------+---------+----------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.serial = int(line[7:10].strip())
        self.mn1 = float(line[10:20].strip())
        self.mn2 = float(line[20:30].strip())
        self.mn3 = float(line[30:40].strip())
        self.vecn = float(line[45:55].strip())
        try:
            self.i_given = int(line[59].strip())
        except (ValueError, IndexError):
            self.i_given = None


@register_line_parser
class MTRIX3(MTRIXn):
    """MATRIX3 PDB entry"""


@register_line_parser
class MTRIX2(MTRIXn):
    """MATRIX2 PDB entry"""


@register_line_parser
class MTRIX1(MTRIXn):
    """MATRIX1 PDB entry"""


class SCALEn(BaseRecord):
    """SCALEn baseclass

    The SCALEn (n = 1, 2, or 3) records present the transformation from the
    orthogonal coordinates as contained in the entry to fractional
    crystallographic coordinates. Non-standard coordinate systems should be
    explained in the remarks.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+-------+-------+------------+
        | COLUMNS | TYPE  | FIELD | DEFINITION |
        +=========+=======+=======+============+
        | 11-20   | float | sn1   | S31        |
        +---------+-------+-------+------------+
        | 21-30   | float | sn2   | S32        |
        +---------+-------+-------+------------+
        | 31-40   | float | sn3   | S33        |
        +---------+-------+-------+------------+
        | 46-55   | float | un    | U3         |
        +---------+-------+-------+------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.sn1 = float(line[10:20].strip())
        self.sn2 = float(line[20:30].strip())
        self.sn3 = float(line[30:40].strip())
        self.unif = float(line[45:55].strip())


@register_line_parser
class SCALE3(SCALEn):
    """SCALE3 PDB entry"""


@register_line_parser
class SCALE2(SCALEn):
    """SCALE2 PDB entry"""


@register_line_parser
class SCALE1(SCALEn):
    """SCALE2 PDB entry"""


class ORIGXn(BaseRecord):
    """ORIGXn class

    The ORIGXn (n = 1, 2, or 3) records present the transformation from the
    orthogonal coordinates contained in the entry to the submitted
    coordinates.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+-------+-------+------------+
        | COLUMNS | TYPE  | FIELD | DEFINITION |
        +=========+=======+=======+============+
        | 11-20   | float | on1   | O21        |
        +---------+-------+-------+------------+
        | 21-30   | float | on2   | O22        |
        +---------+-------+-------+------------+
        | 31-40   | float | on3   | O23        |
        +---------+-------+-------+------------+
        | 46-55   | float | tn    | T2         |
        +---------+-------+-------+------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.on1 = float(line[10:20].strip())
        self.on2 = float(line[20:30].strip())
        self.on3 = float(line[30:40].strip())
        self.tn = float(line[45:55].strip())


@register_line_parser
class ORIGX2(ORIGXn):
    """ORIGX2 PDB entry"""


@register_line_parser
class ORIGX3(ORIGXn):
    """ORIGX3 PDB entry"""


@register_line_parser
class ORIGX1(ORIGXn):
    """ORIGX3 PDB entry"""


@register_line_parser
class CRYST1(BaseRecord):
    """CRYST1 class

    The CRYST1 record presents the unit cell parameters, space group, and Z
    value. If the structure was not determined by crystallographic means,
    CRYST1 simply defines a unit cube.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+-------------+------------------+
        | COLUMNS | TYPE   | FIELD       | DEFINITION       |
        +=========+========+=============+==================+
        | 7-15    | float  | a           | a (Angstroms).   |
        +---------+--------+-------------+------------------+
        | 16-24   | float  | b           | b (Angstroms).   |
        +---------+--------+-------------+------------------+
        | 25-33   | float  | c           | c (Angstroms).   |
        +---------+--------+-------------+------------------+
        | 34-40   | float  | alpha       | alpha (degrees). |
        +---------+--------+-------------+------------------+
        | 41-47   | float  | beta        | beta (degrees).  |
        +---------+--------+-------------+------------------+
        | 48-54   | float  | gamma       | gamma (degrees). |
        +---------+--------+-------------+------------------+
        | 56-66   | string | space_group | Space group.     |
        +---------+--------+-------------+------------------+
        | 67-70   | int    | z           | Z value.         |
        +---------+--------+-------------+------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.a = float(line[6:15].strip())
        self.b = float(line[15:24].strip())
        self.c = float(line[24:33].strip())
        self.alpha = float(line[33:40].strip())
        self.beta = float(line[40:47].strip())
        self.gamma = float(line[47:54].strip())
        self.space_group = line[55:65].strip()
        self.z = int(line[66:70].strip())


@register_line_parser
class SITE(BaseRecord):
    """SITE class

    The SITE records supply the identification of groups comprising
    important sites in the macromolecule.
    """

    def __init__(self, line):
        """Initialize by parsing the line

        +---------+--------+-----------+-------------------------------------+
        | COLUMNS | TYPE   | FIELD     | DEFINITION                          |
        | 8-10    | int    | seq_num   | Sequence number.                    |
        +---------+--------+-----------+-------------------------------------+
        | 12-14   | string | site_id   | Site name.                          |
        +---------+--------+-----------+-------------------------------------+
        | 16-17   | int    | num_res   | Number of residues comprising site. |
        +---------+--------+-----------+-------------------------------------+
        | 19-21   | string | res_name1 | Residue name for first residue      |
        |         |        |           | comprising site.                    |
        +---------+--------+-----------+-------------------------------------+
        | 23      | string | chain_id1 | Chain identifier for first residue  |
        |         |        |           | comprising site.                    |
        +---------+--------+-----------+-------------------------------------+
        | 24-27   | int    | seq1      | Residue sequence number for first   |
        |         |        |           | residue comprising site.            |
        +---------+--------+-----------+-------------------------------------+
        | 28      | string | ins_code1 | Insertion code for first residue    |
        |         |        |           | comprising site.                    |
        +---------+--------+-----------+-------------------------------------+
        | 30-32   | string | res_name2 | Residue name for second residue     |
        |         |        |           | comprising site.                    |
        +---------+--------+-----------+-------------------------------------+
        | 34      | string | chain_id2 | Chain identifier for second residue |
        |         |        |           | comprising site.                    |
        +---------+--------+-----------+-------------------------------------+
        | 35-38   | int    | seq2      | Residue sequence number for second  |
        |         |        |           | residue comprising site.            |
        +---------+--------+-----------+-------------------------------------+
        | 39      | string | ins_code2 | Insertion code for second residue   |
        |         |        |           | comprising site.                    |
        +---------+--------+-----------+-------------------------------------+
        | 41-43   | string | res_name3 | Residue name for third residue      |
        |         |        |           | comprising site.                    |
        +---------+--------+-----------+-------------------------------------+
        | 45      | string | chain_id3 | Chain identifier for third residue  |
        |         |        |           | comprising site.                    |
        +---------+--------+-----------+-------------------------------------+
        | 46-49   | int    | seq3      | Residue sequence number for third   |
        |         |        |           | residue comprising site.            |
        +---------+--------+-----------+-------------------------------------+
        | 50      | string | ins_code3 | Insertion code for third residue    |
        |         |        |           | comprising site.                    |
        +---------+--------+-----------+-------------------------------------+
        | 52-54   | string | res_name4 | Residue name for fourth residue     |
        |         |        |           | comprising site.                    |
        +---------+--------+-----------+-------------------------------------+
        | 56      | string | chain_id4 | Chain identifier for fourth residue |
        |         |        |           | comprising site.                    |
        +---------+--------+-----------+-------------------------------------+
        | 57-60   | int    | seq4      | Residue sequence number for fourth  |
        |         |        |           | residue comprising site.            |
        +---------+--------+-----------+-------------------------------------+
        | 61      | string | ins_code4 | Insertion code for fourth residue   |
        |         |        |           | comprising site.                    |
        +---------+--------+-----------+-------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.seq_num = int(line[7:10].strip())
        self.site_id = line[11:14].strip()
        self.num_res = int(line[15:17].strip())
        self.res_name1 = line[18:21].strip()
        self.chain_id1 = line[22].strip()
        self.seq1 = int(line[23:27].strip())
        self.ins_code1 = line[27].strip()
        self.res_name2 = line[29:32].strip()
        self.chain_id2 = line[33].strip()
        self.seq2 = int(line[34:38].strip())
        self.ins_code2 = line[38].strip()
        self.res_name3 = line[40:43].strip()
        self.chain_id3 = line[44].strip()
        self.seq3 = int(line[45:49].strip())
        self.ins_code3 = line[49].strip()
        self.res_name4 = line[51:54].strip()
        self.chain_id4 = line[55].strip()
        self.seq4 = int(line[56:60].strip())
        try:
            self.ins_code4 = line[60].strip()
        except IndexError:
            self.ins_code4 = None


@register_line_parser
class CISPEP(BaseRecord):
    """CISPEP field

    CISPEP records specify the prolines and other peptides found to be in
    the cis conformation. This record replaces the use of footnote records
    to list cis peptides.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+-----------+----------------------------------+
        | COLUMNS | TYPE   | FIELD     | DEFINITION                       |
        +=========+========+===========+==================================+
        | 8-10    | int    | ser_num   | Record serial number.            |
        +---------+--------+-----------+----------------------------------+
        | 12-14   | string | pep1      | Residue name.                    |
        +---------+--------+-----------+----------------------------------+
        | 16      | string | chain_id1 | Chain identifier.                |
        +---------+--------+-----------+----------------------------------+
        | 18-21   | int    | seq_num1  | Residue sequence number.         |
        +---------+--------+-----------+----------------------------------+
        | 22      | string | icode1    | Insertion code.                  |
        +---------+--------+-----------+----------------------------------+
        | 26-28   | string | pep2      | Residue name.                    |
        +---------+--------+-----------+----------------------------------+
        | 30      | string | chain_id2 | Chain identifier.                |
        +---------+--------+-----------+----------------------------------+
        | 32-35   | int    | seq_num2  | Residue sequence number.         |
        +---------+--------+-----------+----------------------------------+
        | 36      | string | icode2    | Insertion code.                  |
        +---------+--------+-----------+----------------------------------+
        | 44-46   | int    | mod_num   | Identifies the specific model.   |
        +---------+--------+-----------+----------------------------------+
        | 54-59   | float  | measure   | Measure of the angle in degrees. |
        +---------+--------+-----------+----------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.ser_num = int(line[7:10].strip())
        self.pep1 = line[11:14].strip()
        self.chain_id1 = line[15].strip()
        self.seq_num1 = int(line[17:21].strip())
        self.icode1 = line[21].strip()
        self.pep2 = line[25:28].strip()
        self.chain_id2 = line[29].strip()
        self.seq_num2 = int(line[31:35].strip())
        self.icode2 = line[35].strip()
        self.mod_num = int(line[43:46].strip())
        self.measure = float(line[53:59].strip())


@register_line_parser
class SLTBRG(BaseRecord):
    """SLTBRG field

    The SLTBRG records specify salt bridges in the entry.
    records and is provided here for convenience in searching.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+-----------+---------------------------------+
        | COLUMNS | TYPE   | FIELD     | DEFINITION                      |
        +=========+========+===========+=================================+
        | 13-16   | string | name1     | Atom name.                      |
        +---------+--------+-----------+---------------------------------+
        | 17      | string | alt_loc1  | Alternate location indicator.   |
        +---------+--------+-----------+---------------------------------+
        | 18-20   | string | res_name1 | Residue name.                   |
        +---------+--------+-----------+---------------------------------+
        | 22      | string | chain_id1 | Chain identifier.               |
        +---------+--------+-----------+---------------------------------+
        | 23-26   | int    | res_seq1  | Residue sequence number.        |
        +---------+--------+-----------+---------------------------------+
        | 27      | string | ins_code1 | Insertion code.                 |
        +---------+--------+-----------+---------------------------------+
        | 43-46   | string | name2     | Atom name.                      |
        +---------+--------+-----------+---------------------------------+
        | 47      | string | alt_loc2  | Alternate location indicator.   |
        +---------+--------+-----------+---------------------------------+
        | 48-50   | string | res_name2 | Residue name.                   |
        +---------+--------+-----------+---------------------------------+
        | 52      | string | chain_id2 | Chain identifier.               |
        +---------+--------+-----------+---------------------------------+
        | 53-56   | int    | res_seq2  | Residue sequence number.        |
        +---------+--------+-----------+---------------------------------+
        | 57      | string | ins_code2 | Insertion code.                 |
        +---------+--------+-----------+---------------------------------+
        | 60-65   | string | sym1      | Symmetry operator for 1st atom. |
        +---------+--------+-----------+---------------------------------+
        | 67-72   | string | sym2      | Symmetry operator for 2nd atom. |
        +---------+--------+-----------+---------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.name1 = line[12:16].strip()
        self.alt_loc1 = line[16].strip()
        self.res_name1 = line[17:20].strip()
        self.chain_id1 = line[21].strip()
        self.res_seq1 = int(line[22:26].strip())
        self.ins_code1 = line[26].strip()
        self.name2 = line[42:46].strip()
        self.alt_loc2 = line[46].strip()
        self.res_name2 = line[47:50].strip()
        self.chain_id2 = line[51].strip()
        self.res_seq2 = int(line[52:56].strip())
        self.ins_code2 = line[56].strip()
        self.sym1 = line[59:65].strip()
        self.sym2 = line[66:72].strip()


@register_line_parser
class HYDBND(BaseRecord):
    """HYDBND field

    The HYDBND records specify hydrogen bonds in the entry.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+-----------+-------------------------------------+
        | COLUMNS | TYPE   | FIELD     | DEFINITION                          |
        +=========+========+===========+=====================================+
        | 13-16   | string | name1     | Atom name.                          |
        +---------+--------+-----------+-------------------------------------+
        | 17      | string | alt_loc1  | Alternate location indicator.       |
        +---------+--------+-----------+-------------------------------------+
        | 18-20   | string | res_name1 | Residue name.                       |
        +---------+--------+-----------+-------------------------------------+
        | 22      | string | chain1    | Chain identifier.                   |
        +---------+--------+-----------+-------------------------------------+
        | 23-27   | int    | res_seq1  | Residue sequence number.            |
        +---------+--------+-----------+-------------------------------------+
        | 28      | string | i_code1   | Insertion code.                     |
        +---------+--------+-----------+-------------------------------------+
        | 30-33   | string | name_h    | Hydrogen atom name.                 |
        +---------+--------+-----------+-------------------------------------+
        | 34      | string | alt_loc_h | Alternate location indicator.       |
        +---------+--------+-----------+-------------------------------------+
        | 36      | string | chain_h   | Chain identifier.                   |
        +---------+--------+-----------+-------------------------------------+
        | 37-41   | int    | res_seq_h | Residue sequence number.            |
        +---------+--------+-----------+-------------------------------------+
        | 42      | string | ins_codeH | Insertion code.                     |
        +---------+--------+-----------+-------------------------------------+
        | 44-47   | string | name2     | Atom name.                          |
        +---------+--------+-----------+-------------------------------------+
        | 48      | string | alt_loc2  | Alternate location indicator.       |
        +---------+--------+-----------+-------------------------------------+
        | 49-51   | string | res_name2 | Residue name.                       |
        +---------+--------+-----------+-------------------------------------+
        | 53      | string | chain_id2 | Chain identifier.                   |
        +---------+--------+-----------+-------------------------------------+
        | 54-58   | int    | res_seq2  | Residue sequence number.            |
        +---------+--------+-----------+-------------------------------------+
        | 59      | string | ins_code2 | Insertion code.                     |
        +---------+--------+-----------+-------------------------------------+
        | 60-65   | string | sym1      | Symmetry operator for 1st           |
        |         |        |           | non-hydrogen atom.                  |
        +---------+--------+-----------+-------------------------------------+
        | 67-72   | string | sym2      | Symmetry operator for 2nd           |
        |         |        |           | non-hydrogen atom.                  |
        +---------+--------+-----------+-------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.name1 = line[12:16].strip()
        self.alt_loc1 = line[16].strip()
        self.res_name1 = line[17:20].strip()
        self.chain1 = line[21].strip()
        self.res_seq1 = line[22:27].strip()
        self.i_code1 = line[27].strip()
        self.name_h = line[29:33].strip()
        self.alt_loc_h = line[33].strip()
        self.chain_h = line[35].strip()
        self.res_seq_h = line[36:41].strip()
        self.i_code_h = line[41].strip()
        self.name2 = line[43:47].strip()
        self.alt_loc2 = line[47].strip()
        self.res_name2 = line[48:51].strip()
        self.chain2 = line[52].strip()
        self.res_seq2 = line[53:58].strip()
        self.i_code2 = line[58].strip()
        self.sym1 = line[59:65].strip()
        self.sym2 = line[66:72].strip()


@register_line_parser
class LINK(BaseRecord):
    """LINK field

    The LINK records specify connectivity between residues that is not
    implied by the primary structure. Connectivity is expressed in terms of
    the atom names. This record supplements information given in CONECT
    records and is provided here for convenience in searching.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+-----------+---------------------------------+
        | COLUMNS | TYPE   | FIELD     | DEFINITION                      |
        +=========+========+===========+=================================+
        | 13-16   | string | name1     | Atom name.                      |
        +---------+--------+-----------+---------------------------------+
        | 17      | string | alt_loc1  | Alternate location indicator.   |
        +---------+--------+-----------+---------------------------------+
        | 18-20   | string | res_name1 | Residue name.                   |
        +---------+--------+-----------+---------------------------------+
        | 22      | string | chain_id1 | Chain identifier.               |
        +---------+--------+-----------+---------------------------------+
        | 23-26   | int    | res_seq1  | Residue sequence number.        |
        +---------+--------+-----------+---------------------------------+
        | 27      | string | ins_code1 | Insertion code.                 |
        +---------+--------+-----------+---------------------------------+
        | 43-46   | string | name2     | Atom name.                      |
        +---------+--------+-----------+---------------------------------+
        | 47      | string | alt_loc2  | Alternate location indicator.   |
        +---------+--------+-----------+---------------------------------+
        | 48-50   | string | res_name2 | Residue name.                   |
        +---------+--------+-----------+---------------------------------+
        | 52      | string | chain_id2 | Chain identifier.               |
        +---------+--------+-----------+---------------------------------+
        | 53-56   | int    | res_seq2  | Residue sequence number.        |
        +---------+--------+-----------+---------------------------------+
        | 57      | string | ins_code2 | Insertion code.                 |
        +---------+--------+-----------+---------------------------------+
        | 60-65   | string | sym1      | Symmetry operator for 1st atom. |
        +---------+--------+-----------+---------------------------------+
        | 67-72   | string | sym2      | Symmetry operator for 2nd atom. |
        +---------+--------+-----------+---------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.name1 = line[12:16].strip()
        self.alt_loc1 = line[16].strip()
        self.res_name1 = line[17:20].strip()
        self.chain_id1 = line[21].strip()
        self.res_seq1 = int(line[22:26].strip())
        self.ins_code1 = line[26].strip()
        self.name2 = line[42:46].strip()
        self.alt_loc2 = line[46].strip()
        self.res_name2 = line[47:50].strip()
        self.chain_id2 = line[51].strip()
        self.res_seq2 = int(line[52:56].strip())
        self.ins_code2 = line[56].strip()
        self.sym1 = line[59:65].strip()
        self.sym2 = line[66:72].strip()


@register_line_parser
class SSBOND(BaseRecord):
    """SSBOND field

    The SSBOND record identifies each disulfide bond in protein and
    polypeptide structures by identifying the two residues involved in the
    bond.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+-----------+------------------------------------+
        | COLUMNS | TYPE   | FIELD     | DEFINITION                         |
        +=========+========+===========+====================================+
        | 8-10    | int    | ser_num   | Serial number.                     |
        +---------+--------+-----------+------------------------------------+
        | 16      | string | chain_id1 | Chain identifier.                  |
        +---------+--------+-----------+------------------------------------+
        | 18-21   | int    | seq_num1  | Residue sequence number.           |
        +---------+--------+-----------+------------------------------------+
        | 22      | string | icode1    | Insertion code.                    |
        +---------+--------+-----------+------------------------------------+
        | 30      | string | chain_id2 | Chain identifier.                  |
        +---------+--------+-----------+------------------------------------+
        | 32-35   | int    | seq_num2  | Residue sequence number.           |
        +---------+--------+-----------+------------------------------------+
        | 36      | string | icode2    | Insertion code.                    |
        +---------+--------+-----------+------------------------------------+
        | 60-65   | string | sym1      | Symmetry operator for 1st residue. |
        +---------+--------+-----------+------------------------------------+
        | 67-72   | string | sym2      | Symmetry operator for 2nd residue. |
        +---------+--------+-----------+------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.ser_num = int(line[7:10].strip())
        self.chain_id1 = line[15].strip()
        self.seq_num1 = int(line[17:21].strip())
        self.icode1 = line[21].strip()
        self.chain_id2 = line[29].strip()
        self.seq_num2 = int(line[31:35].strip())
        self.icode2 = line[35].strip()
        self.sym1 = line[59:65].strip()
        self.sym2 = line[66:72].strip()


@register_line_parser
class TURN(BaseRecord):
    """TURN field

    The TURN records identify turns and other short loop turns which normally
    connect other secondary structure segments.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+---------------+---------------------------------+
        | COLUMNS | TYPE   | FIELD         | DEFINITION                      |
        +=========+========+===============+=================================+
        | 8-10    | int    | seq           | Turn number; starts with 1 and  |
        |         |        |               | increments by one.              |
        +---------+--------+---------------+---------------------------------+
        | 12-14   | string | turn_id       | Turn identifier.                |
        +---------+--------+---------------+---------------------------------+
        | 16-18   | string | init_res_name | Residue name of initial residue |
        |         |        |               | in turn.                        |
        +---------+--------+---------------+---------------------------------+
        | 20      | string | init_chain_id | Chain identifier for the chain  |
        |         |        |               | containing this turn.           |
        +---------+--------+---------------+---------------------------------+
        | 21-24   | int    | init_seq_num  | Sequence number of initial      |
        |         |        |               | residue in turn.                |
        +---------+--------+---------------+---------------------------------+
        | 25      | string | init_i_code   | Insertion code of initial       |
        |         |        |               | residue in turn.                |
        +---------+--------+---------------+---------------------------------+
        | 27-29   | string | end_res_name  | Residue name of terminal residue|
        |         |        |               | of turn.                        |
        +---------+--------+---------------+---------------------------------+
        | 31      | string | end_chain_id  | Chain identifier for the chain  |
        |         |        |               | containing this turn.           |
        +---------+--------+---------------+---------------------------------+
        | 32-35   | int    | end_seq_num   | Sequence number of terminal     |
        |         |        |               | residue of turn.                |
        +---------+--------+---------------+---------------------------------+
        | 36      | string | end_i_code    | Insertion code of terminal      |
        |         |        |               | residue of turn.                |
        +---------+--------+---------------+---------------------------------+
        | 41-70   | string | comment       | Associated comment.             |
        +---------+--------+---------------+---------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.seq = int(line[7:10].strip())
        self.turn_id = line[11:14].strip()
        self.init_res_name = line[15:18].strip()
        self.init_chain_id = line[19].strip()
        self.init_seq_num = int(line[20:24].strip())
        self.init_i_code = line[24].strip()
        self.end_res_name = line[26:29].strip()
        self.end_chain_id = line[30].strip()
        self.end_seq_num = int(line[31:35].strip())
        self.end_i_code = line[35].strip()
        self.comment = line[40:70].strip()


@register_line_parser
class SHEET(BaseRecord):
    """SHEET field

    SHEET records are used to identify the position of sheets in the
    molecule. Sheets are both named and numbered. The residues where the
    sheet begins and ends are noted.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+---------------+---------------------------------+
        | COLUMNS | TYPE   | FIELD         | DEFINITION                      |
        +=========+========+===============+=================================+
        | 8-10    | int    | strand        | Strand number which starts at 1 |
        |         |        |               | for each strand within a sheet  |
        |         |        |               | and increases by one.           |
        +---------+--------+---------------+---------------------------------+
        | 12-14   | string | sheet_id      | Sheet identifier.               |
        +---------+--------+---------------+---------------------------------+
        | 15-16   | int    | num_strands   | Number of strands in sheet.     |
        +---------+--------+---------------+---------------------------------+
        | 18-20   | string | init_res_name | Residue name of initial residue.|
        +---------+--------+---------------+---------------------------------+
        | 22      | string | init_chain_id | Chain identifier of initial     |
        |         |        |               | residue in strand.              |
        +---------+--------+---------------+---------------------------------+
        | 23-26   | int    | init_seq_num  | Sequence number of initial      |
        |         |        |               | residue in strand.              |
        +---------+--------+---------------+---------------------------------+
        | 27      | string | init_i_code   | Insertion code of initial       |
        |         |        |               | residue in strand.              |
        +---------+--------+---------------+---------------------------------+
        | 29-31   | string | end_res_name  | Residue name of terminal        |
        |         |        |               | residue.                        |
        +---------+--------+---------------+---------------------------------+
        | 33      | string | end_chain_id  | Chain identifier of terminal    |
        |         |        |               | residue.                        |
        +---------+--------+---------------+---------------------------------+
        | 34-37   | int    | end_seq_num   | Sequence number of terminal     |
        |         |        |               | residue.                        |
        +---------+--------+---------------+---------------------------------+
        | 38      | string | end_i_code    | Insertion code of terminal      |
        |         |        |               | residue.                        |
        +---------+--------+---------------+---------------------------------+
        | 39-40   | int    | sense         | Sense of strand with respect to |
        |         |        |               | previous strand in the sheet. 0 |
        |         |        |               | if first strand, 1 if parallel, |
        |         |        |               | -1 if anti-parallel.            |
        +---------+--------+---------------+---------------------------------+
        | 42-45   | string | cur_atom      | Registration. Atom name in      |
        |         |        |               | current strand.                 |
        +---------+--------+---------------+---------------------------------+
        | 46-48   | string | curr_res_name | Registration. Residue name in   |
        |         |        |               | current strand.                 |
        +---------+--------+---------------+---------------------------------+
        | 50      | string | curChainId    | Registration. Chain identifier  |
        |         |        |               | in current strand.              |
        +---------+--------+---------------+---------------------------------+
        | 51-54   | int    | curr_res_seq  | Registration. Residue sequence  |
        |         |        |               | number in current strand.       |
        +---------+--------+---------------+---------------------------------+
        | 55      | string | curr_ins_code | Registration. Insertion code in |
        |         |        |               | current strand.                 |
        +---------+--------+---------------+---------------------------------+
        | 57-60   | string | prev_atom     | Registration. Atom name in      |
        |         |        |               | previous strand.                |
        +---------+--------+---------------+---------------------------------+
        | 61-63   | string | prev_res_name | Registration. Residue name in   |
        |         |        |               | previous strand.                |
        +---------+--------+---------------+---------------------------------+
        | 65      | string | prevChainId   | Registration. Chain identifier  |
        |         |        |               | in previous strand.             |
        +---------+--------+---------------+---------------------------------+
        | 66-69   | int    | prev_res_seq  | Registration. Residue sequence  |
        |         |        |               | number in previous strand.      |
        +---------+--------+---------------+---------------------------------+
        | 70      | string | prev_ins_code | Registration. Insertion code in |
        |         |        |               | previous strand.                |
        +---------+--------+---------------+---------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.strand = int(line[7:10].strip())
        self.sheet_id = line[11:14].strip()
        self.num_strands = int(line[14:16].strip())
        self.init_res_name = line[17:20].strip()
        self.init_chain_id = line[21].strip()
        self.init_seq_num = int(line[22:26].strip())
        self.init_i_code = line[26].strip()
        self.end_res_name = line[28:31].strip()
        self.end_chain_id = line[32].strip()
        self.end_seq_num = int(line[33:37].strip())
        self.end_i_code = line[37].strip()
        self.sense = int(line[38:40].strip())
        try:
            self.cur_atom = line[41:45].strip()
            self.curr_res_name = line[45:48].strip()
            self.curr_chain_id = line[49].strip()
            try:
                self.curr_res_seq = int(line[50:54].strip())
            except ValueError:
                self.curr_res_seq = None
            self.curr_ins_code = line[54].strip()
            self.prev_atom = line[56:60].strip()
            self.prev_res_name = line[60:63].strip()
            self.prev_chain_id = line[64].strip()
            try:
                self.prev_res_seq = int(line[65:69].strip())
            except ValueError:
                self.prev_res_seq = None
            self.prev_ins_code = line[69].strip()
        except IndexError:
            self.cur_atom = None
            self.curr_res_name = None
            self.curr_chain_id = None
            self.curr_res_seq = None
            self.curr_ins_code = None
            self.prev_atom = None
            self.prev_res_name = None
            self.prev_chain_id = None
            self.prev_res_seq = None
            self.prev_ins_code = None


@register_line_parser
class HELIX(BaseRecord):
    """HELIX field

    HELIX records are used to identify the position of helices in the
    molecule. Helices are both named and numbered. The residues where the
    helix begins and ends are noted, as well as the total length.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+---------------+---------------------------------+
        | COLUMNS | TYPE   | FIELD         | DEFINITION                      |
        +=========+========+===============+=================================+
        | 8-10    | int    | ser_num       | Serial number of the helix. This|
        |         |        |               | starts at 1 and increases       |
        |         |        |               | incrementally.                  |
        +---------+--------+---------------+---------------------------------+
        | 12-14   | string | helix_id      | Helix identifier.  In addition  |
        |         |        |               | to a serial number, each helix  |
        |         |        |               | is given an alphanumeric        |
        |         |        |               | character helix identifier.     |
        +---------+--------+---------------+---------------------------------+
        | 16-18   | string | init_res_name | Name of the initial residue.    |
        +---------+--------+---------------+---------------------------------+
        | 20      | string | init_chain_id | Chain identifier for the chain  |
        |         |        |               | containing this helix.          |
        +---------+--------+---------------+---------------------------------+
        | 22-25   | int    | init_seq_num  | Sequence number of the initial  |
        |         |        |               | residue.                        |
        +---------+--------+---------------+---------------------------------+
        | 26      | string | init_i_code   | Insertion code of the initial   |
        |         |        |               | residue.                        |
        +---------+--------+---------------+---------------------------------+
        | 28-30   | string | end_res_name  | Name of the terminal residue of |
        |         |        |               | the helix.                      |
        +---------+--------+---------------+---------------------------------+
        | 32      | string | end_chain_id  | Chain identifier for the chain  |
        |         |        |               | containing this helix.          |
        +---------+--------+---------------+---------------------------------+
        | 34-37   | int    | end_seq_num   | Sequence number of the terminal |
        |         |        |               | residue.                        |
        +---------+--------+---------------+---------------------------------+
        | 38      | string | end_i_code    | Insertion code of the terminal  |
        |         |        |               | residue.                        |
        +---------+--------+---------------+---------------------------------+
        | 39-40   | int    | helix_class   | Helix class (see below).        |
        +---------+--------+---------------+---------------------------------+
        | 41-70   | string | comment       | Comment about this helix.       |
        +---------+--------+---------------+---------------------------------+
        | 72-76   | int    | length        | Length of this helix.           |
        +---------+--------+---------------+---------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.ser_num = int(line[7:10].strip())
        self.helix_id = line[11:14].strip()
        self.init_res_name = line[15:18].strip()
        self.init_chain_id = line[19].strip()
        self.init_seq_num = int(line[21:25].strip())
        self.init_i_code = line[25].strip()
        self.end_res_name = line[27:30].strip()
        self.end_chain_id = line[31].strip()
        self.end_seq_num = int(line[33:37].strip())
        self.end_i_code = line[37].strip()
        try:
            self.helix_class = int(line[38:40].strip())
        except ValueError:
            self.helix_class = None
        self.comment = line[40:70].strip()
        try:
            self.length = int(line[71:76].strip())
        except ValueError:
            self.length = None


@register_line_parser
class FORMUL(BaseRecord):
    """FORMUL field

    The FORMUL record presents the chemical formula and charge of a
    non-standard group.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+------------+------------------+
        | COLUMNS | TYPE   | FIELD      | DEFINITION       |
        +=========+========+============+==================+
        | 9-10    | int    | comp_num   | Component number |
        +---------+--------+------------+------------------+
        | 13-15   | string | hetatm_id  | Het identifier   |
        +---------+--------+------------+------------------+
        | 19      | string | asterisk * | for water        |
        +---------+--------+------------+------------------+
        | 20-70   | string | text       | Chemical formula |
        +---------+--------+------------+------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.comp_num = int(line[8:10].strip())
        self.hetatm_id = line[12:15].strip()
        self.asterisk = line[19].strip()
        self.text = line[19:70].strip()


@register_line_parser
class HETSYN(BaseRecord):
    """HETSYN field

    This record provides synonyms, if any, for the compound in the
    corresponding (i.e., same hetatm_id) HETNAM record. This is to allow
    greater flexibility in searching for HET groups.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+-----------------+-------------------+
        | COLUMNS | TYPE   | FIELD           | DEFINITION        |
        +=========+========+=================+===================+
        | 12-14   | string | hetatm_id       | Het identifier,   |
        |         |        |                 | right-justified.  |
        +---------+--------+-----------------+-------------------+
        | 16-70   | string | hetatm_synonyms | List of synonyms  |
        +---------+--------+-----------------+-------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.hetatm_id = line[11:14].strip()
        self.hetatm_synonyms = line[15:70].strip()


@register_line_parser
class HETNAM(BaseRecord):
    """HETNAM field

    This record gives the chemical name of the compound with the
    given hetatm_id.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+-----------+----------------------------------+
        | COLUMNS | TYPE   | FIELD     | DEFINITION                       |
        +=========+========+===========+==================================+
        | 12-14   | string | hetatm_id | Het identifier, right-justified. |
        +---------+--------+-----------+----------------------------------+
        | 16-70   | string | text      | Chemical name.                   |
        +---------+--------+-----------+----------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.hetatm_id = line[11:14].strip()
        self.text = line[15:70].strip()


@register_line_parser
class HET(BaseRecord):
    """HET field

    HET records are used to describe non-standard residues, such as
    prosthetic groups, inhibitors, solvent molecules, and ions for which
    coordinates are supplied. Groups are considered HET if they are:

    * not one of the standard amino acids, and
    * not one of the nucleic acids (C, G, A, T, U, and I), and
    * not one of the modified versions of nucleic acids (+C, +G, +A, +T, +U,
      and +I), and
    * not an unknown amino acid or nucleic acid where UNK is used to indicate
      the unknown residue name.

    Het records also describe heterogens for which the chemical identity is
    unknown, in which case the group is assigned the hetatm_id UNK.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+---------------+---------------------------------+
        | COLUMNS | TYPE   | FIELD         | DEFINITION                      |
        +=========+========+===============+=================================+
        | 8-10    | string | hetatm_id     | Het identifier, right-justified.|
        +---------+--------+---------------+---------------------------------+
        | 13      | string | ChainID       | Chain identifier.               |
        +---------+--------+---------------+---------------------------------+
        | 14-17   | int    | seq_num       | Sequence number.                |
        +---------+--------+---------------+---------------------------------+
        | 18      | string | ins_code      | Insertion code.                 |
        +---------+--------+---------------+---------------------------------+
        | 21-25   | int    | num_het_atoms | Number of HETATM records.       |
        +---------+--------+---------------+---------------------------------+
        | 31-70   | string | text          | Text describing Het group.      |
        +---------+--------+---------------+---------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.hetatm_id = line[7:10].strip()
        self.chain_id = line[12].strip()
        try:
            self.seq_num = int(line[13].strip())
        except ValueError:
            self.seq_num = None
        self.ins_code = line[17].strip()
        self.num_het_atoms = int(line[20:25].strip())
        self.text = line[30:70].strip()


@register_line_parser
class MODRES(BaseRecord):
    """MODRES field

    The MODRES record provides descriptions of modifications (e.g.,
    chemical or post-translational) to protein and nucleic acid residues.
    Included are a mapping between residue names given in a PDB entry and
    standard residues.
    """

    def __init__(self, line):
        """Initialize by parsing a line

        +---------+--------+----------+--------------------------------------+
        | COLUMNS | TYPE   | FIELD    | DEFINITION                           |
        +=========+========+==========+======================================+
        | 8-11    | string | id_code  | ID code of this entry.               |
        +---------+--------+----------+--------------------------------------+
        | 13-15   | string | res_name | Residue name used in this entry.     |
        +---------+--------+----------+--------------------------------------+
        | 17      | string | chain_id | Chain identifier.                    |
        +---------+--------+----------+--------------------------------------+
        | 19-22   | int    | seq_num  | Sequence number.                     |
        +---------+--------+----------+--------------------------------------+
        | 23      | string | ins_code | Insertion code.                      |
        +---------+--------+----------+--------------------------------------+
        | 25-27   | string | stdRes   | Standard residue name.               |
        +---------+--------+----------+--------------------------------------+
        | 30-70   | string | comment  | Description of the residue           |
        |         |        |          | modification.                        |
        +---------+--------+----------+--------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.id_code = line[7:11].strip()
        self.res_name = line[12:15].strip()
        self.chain_id = line[16].strip()
        self.seq_num = int(line[18:22].strip())
        self.ins_code = line[22].strip()
        self.stdRes = line[24:27].strip()
        self.comment = line[29:70].strip()


@register_line_parser
class SEQRES(BaseRecord):
    """SEQRES field

    SEQRES records contain the amino acid or nucleic acid sequence of
    residues in each chain of the macromolecule that was studied.
    """

    def __init__(self, line):
        """Initialize by parsing a line

        +---------+--------+----------+--------------------------------------+
        | COLUMNS | TYPE   | FIELD    | DEFINITION                           |
        +=========+========+==========+======================================+
        | 9-10    | int    | ser_num  | Serial number of the SEQRES record   |
        |         |        |          | for the current chain.  Starts at 1  |
        |         |        |          | and increments by one each line.     |
        |         |        |          | Reset to 1 for each chain.           |
        +---------+--------+----------+--------------------------------------+
        | 12      | string | chain_id | Chain identifier.  This may be any   |
        |         |        |          | single legal character, including a  |
        |         |        |          | blank which is used if there is only |
        |         |        |          | one chain.                           |
        +---------+--------+----------+--------------------------------------+
        | 14-17   | int    | num_res  | Number of residues in the chain. This|
        |         |        |          | value is repeated on every record.   |
        +---------+--------+----------+--------------------------------------+
        | 20-22   | string | res_name | Residue name.                        |
        +---------+--------+----------+--------------------------------------+
        | 24-26   | string | res_name | Residue name.                        |
        +---------+--------+----------+--------------------------------------+
        | 28-30   | string | res_name | Residue name.                        |
        +---------+--------+----------+--------------------------------------+
        | 32-34   | string | res_name | Residue name.                        |
        +---------+--------+----------+--------------------------------------+
        | 36-38   | string | res_name | Residue name.                        |
        +---------+--------+----------+--------------------------------------+
        | 40-42   | string | res_name | Residue name.                        |
        +---------+--------+----------+--------------------------------------+
        | 44-46   | string | res_name | Residue name.                        |
        +---------+--------+----------+--------------------------------------+
        | 48-50   | string | res_name | Residue name.                        |
        +---------+--------+----------+--------------------------------------+
        | 52-54   | string | res_name | Residue name.                        |
        +---------+--------+----------+--------------------------------------+
        | 56-58   | string | res_name | Residue name.                        |
        +---------+--------+----------+--------------------------------------+
        | 60-62   | string | res_name | Residue name.                        |
        +---------+--------+----------+--------------------------------------+
        | 64-66   | string | res_name | Residue name.                        |
        +---------+--------+----------+--------------------------------------+
        | 68-70   | string | res_name | Residue name.                        |
        +---------+--------+----------+--------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.ser_num = int(line[8:10].strip())
        self.chain_id = line[11].strip()
        self.num_res = int(line[13:17].strip())
        self.res_name = [line[19:22].strip()]
        self.res_name.append(line[23:26].strip())
        self.res_name.append(line[27:30].strip())
        self.res_name.append(line[31:34].strip())
        self.res_name.append(line[35:38].strip())
        self.res_name.append(line[39:42].strip())
        self.res_name.append(line[43:46].strip())
        self.res_name.append(line[47:50].strip())
        self.res_name.append(line[51:54].strip())
        self.res_name.append(line[55:58].strip())
        self.res_name.append(line[59:62].strip())
        self.res_name.append(line[63:66].strip())
        self.res_name.append(line[67:70].strip())


@register_line_parser
class SEQADV(BaseRecord):
    """SEQADV field

    The SEQADV record identifies conflicts between sequence information in
    the ATOM records of the PDB entry and the sequence database entry given
    on DBREF. Please note that these records were designed to identify
    differences and not errors. No assumption is made as to which database
    contains the correct data. PDB may include REMARK records in the entry
    that reflect the depositor's view of which database has the correct
    sequence.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+------------+------------------------------------+
        | COLUMNS | TYPE   | FIELD      | DEFINITION                         |
        +=========+========+============+====================================+
        | 8-11    | string | id_code    | ID code of this entry.             |
        +---------+--------+------------+------------------------------------+
        | 13-15   | string | res_name   | Name of the PDB residue in         |
        |         |        |            | conflict.                          |
        +---------+--------+------------+------------------------------------+
        | 17      | string | chain_id   | PDB chain identifier.              |
        +---------+--------+------------+------------------------------------+
        | 19-22   | int    | seq_num    | PDB sequence number.               |
        +---------+--------+------------+------------------------------------+
        | 23      | string | ins_code   | PDB insertion code.                |
        +---------+--------+------------+------------------------------------+
        | 25-28   | string | database   | Sequence database name.            |
        +---------+--------+------------+------------------------------------+
        | 30-38   | string | db_id_code | Sequence database accession number.|
        +---------+--------+------------+------------------------------------+
        | 40-42   | string | db_res     | Sequence database residue name.    |
        +---------+--------+------------+------------------------------------+
        | 44-48   | int    | db_seq     | Sequence database sequence number. |
        +---------+--------+------------+------------------------------------+
        | 50-70   | string | conflict   | Conflict comment.                  |
        +---------+--------+------------+------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.id_code = line[7:11].strip()
        self.res_name = line[12:15].strip()
        self.chain_id = line[16].strip()
        try:
            self.seq_num = int(line[19:22].strip())
        except ValueError:
            self.seq_num = None
        self.ins_code = line[22].strip()
        self.database = line[24:28].strip()
        self.db_id_code = line[29:38].strip()
        self.db_res = line[39:42].strip()
        try:
            self.db_seq = int(line[43:48].strip())
        except ValueError:
            self.db_seq = None
        self.conflict = line[49:70].strip()


@register_line_parser
class DBREF(BaseRecord):
    """DBREF field

    The DBREF record provides cross-reference links between PDB sequences
    and the corresponding database entry or entries. A cross reference to
    the sequence database is mandatory for each peptide chain with a length
    greater than ten (10) residues. For nucleic acid entries a DBREF record
    pointing to the Nucleic Acid Database (NDB) is mandatory when the
    corresponding entry exists in NDB.
    """

    def __init__(self, line):
        """Initialize by parsing a line.

        +---------+--------+--------------+----------------------------------+
        | COLUMNS | TYPE   | FIELD        | DEFINITION                       |
        +=========+========+==============+==================================+
        | 8-11    | string | id_code      | ID code of this entry.           |
        +---------+--------+--------------+----------------------------------+
        | 13      | string | chain_id     | Chain identifier.                |
        +---------+--------+--------------+----------------------------------+
        | 15-18   | int    | seq_begin    | Initial sequence number of the   |
        |         |        |              | PDB sequence segment.            |
        +---------+--------+--------------+----------------------------------+
        | 19      | string | insert_begin | Initial insertion code of the    |
        |         |        |              | PDB sequence segment.            |
        +---------+--------+--------------+----------------------------------+
        | 21-24   | int    | seq_end      | Ending sequence number of the    |
        |         |        |              | PDB sequence segment.            |
        +---------+--------+--------------+----------------------------------+
        | 25      | string | insert_end   | Ending insertion code of the     |
        |         |        |              | PDB sequence segment.            |
        +---------+--------+--------------+----------------------------------+
        | 27-32   | string | database     | Sequence database name.  "PDB"   |
        |         |        |              | when a corresponding sequence    |
        |         |        |              | database entry has not been      |
        |         |        |              | identified.                      |
        +---------+--------+--------------+----------------------------------+
        | 34-41   | string | db_accession | Sequence database accession code.|
        |         |        |              | For GenBank entries, this is the |
        |         |        |              | NCBI gi number.                  |
        +---------+--------+--------------+----------------------------------+
        | 43-54   | string | db_id_code   | Sequence database identification |
        |         |        |              | code. For GenBank entries, this  |
        |         |        |              | is the accession code.           |
        +---------+--------+--------------+----------------------------------+
        | 56-60   | int    | db_seq_begin | Initial sequence number of the   |
        |         |        |              | database seqment.                |
        +---------+--------+--------------+----------------------------------+
        | 61      | string | db_ins_begin | Insertion code of initial residue|
        |         |        |              | of the segment, if PDB is the    |
        |         |        |              | reference.                       |
        +---------+--------+--------------+----------------------------------+
        | 63-67   | int    | dbseq_end    | Ending sequence number of the    |
        |         |        |              | database segment.                |
        +---------+--------+--------------+----------------------------------+
        | 68      | string | db_ins_end   | Insertion code of the ending     |
        |         |        |              | residue of the                   |
        |         |        |              | segment, if PDB is the reference.|
        +---------+--------+--------------+----------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.id_code = line[7:11].strip()
        self.chain_id = line[12].strip()
        self.seq_begin = int(line[14:18].strip())
        self.insert_begin = line[18].strip()
        self.seq_end = int(line[20:24].strip())
        self.insert_end = line[24].strip()
        self.database = line[26:32].strip()
        self.db_accession = line[33:41].strip()
        self.db_id_code = line[42:54].strip()
        self.db_seq_begin = int(line[55:60].strip())
        self.db_ins_begin = line[60].strip()
        self.dbseq_end = int(line[62:67].strip())
        try:
            self.db_ins_end = line[67].strip()
        except IndexError:
            self.db_ins_end = None


@register_line_parser
class REMARK(BaseRecord):
    """REMARK field

    REMARK records present experimental details, annotations, comments, and
    information not included in other records. In a number of cases,
    REMARKs are used to expand the contents of other record types. A new
    level of structure is being used for some REMARK records. This is
    expected to facilitate searching and will assist in the conversion to a
    relational database.
    """

    def __init__(self, line):
        """Initialize by parsing line.

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.remark_num = int(line[7:10].strip())
        self.remark_dict = {}

        if self.remark_num == 1:
            subfield = line[11:20].strip()
            if subfield == "REFERENCE":
                self.remark_dict["refNum"] = int(line[21:70].strip())
            elif subfield == "AUTH":
                self.remark_dict["author_list"] = line[19:70].strip()
            elif subfield == "TITL":
                self.remark_dict["title"] = line[19:70].strip()
            elif subfield == "EDIT":
                self.remark_dict["editorList"] = line[19:70].strip()
            elif subfield == "REF":
                self.remark_dict["ref"] = line[19:66].strip()
            elif subfield == "PUBL":
                self.remark_dict["pub"] = line[19:70].strip()
            elif subfield == "REFN":
                self.remark_dict["refn"] = line[19:70].strip()
        elif self.remark_num == 2:
            restr = line[22:27].strip()
            try:
                self.remark_dict["resolution"] = float(restr)
            except ValueError:
                self.remark_dict["comment"] = line[11:70].strip()
        else:
            self.remark_dict["text"] = line[11:70].strip()


@register_line_parser
class JRNL(BaseRecord):
    """JRNL field

    The JRNL record contains the primary literature citation that describes
    the experiment which resulted in the deposited coordinate set. There is
    at most one JRNL reference per entry. If there is no primary reference,
    then there is no JRNL reference. Other references are given in REMARK 1.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+-------+---------------------+
        | COLUMNS | TYPE   | FIELD | DEFINITION          |
        +=========+========+=======+=====================+
        | 13-70   | string | text  | See details on web. |
        +---------+--------+-------+---------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.text = line[12:70].strip()


@register_line_parser
class SPRSDE(BaseRecord):
    """SPRSDE field

    The SPRSDE records contain a list of the ID codes of entries that were
    made obsolete by the given coordinate entry and withdrawn from the PDB
    release set. One entry may replace many. It is PDB policy that only the
    principal investigator of a structure has the authority to withdraw it.
    """

    def __init__(self, line):
        """Initialize by parsing line

        +---------+--------+------------+------------------------------------+
        | COLUMNS | TYPE   | FIELD      | DEFINITION                         |
        +=========+========+============+====================================+
        | 12-20   | string | super_date | Date this entry superseded the     |
        |         |        |            | listed entries.                    |
        +---------+--------+------------+------------------------------------+
        | 22-25   | string | id_code    | ID code of this entry.             |
        +---------+--------+------------+------------------------------------+
        | 32-35   | string | sid_code   | ID code of a superseded entry.     |
        +---------+--------+------------+------------------------------------+
        | 37-40   | string | sid_code   | ID code of a superseded entry.     |
        +---------+--------+------------+------------------------------------+
        | 42-45   | string | sid_code   | ID code of a superseded entry.     |
        +---------+--------+------------+------------------------------------+
        | 47-50   | string | sid_code   | ID code of a superseded entry.     |
        +---------+--------+------------+------------------------------------+
        | 52-55   | string | sid_code   | ID code of a superseded entry.     |
        +---------+--------+------------+------------------------------------+
        | 57-60   | string | sid_code   | ID code of a superseded entry.     |
        +---------+--------+------------+------------------------------------+
        | 62-65   | string | sid_code   | ID code of a superseded entry.     |
        +---------+--------+------------+------------------------------------+
        | 67-70   | string | sid_code   | ID code of a superseded entry.     |
        +---------+--------+------------+------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.super_date = line[11:20].strip()
        self.id_code = line[21:25].strip()
        self.super_id_codes = [line[31:35].strip()]
        self.super_id_codes.append(line[36:40].strip())
        self.super_id_codes.append(line[41:45].strip())
        self.super_id_codes.append(line[46:50].strip())
        self.super_id_codes.append(line[51:55].strip())
        self.super_id_codes.append(line[56:60].strip())
        self.super_id_codes.append(line[61:65].strip())
        self.super_id_codes.append(line[66:70].strip())


@register_line_parser
class REVDAT(BaseRecord):
    """REVDAT field

    REVDAT records contain a history of the modifications made to an entry
    since its release.
    """

    def __init__(self, line):
        """Initialize by parsing a line.

        .. todo::

           If multiple modifications are present, only the last one in the
           file is preserved.

        +---------+--------+----------+--------------------------------------+
        | COLUMNS | TYPE   | FIELD    | DEFINITION                           |
        +=========+========+==========+======================================+
        | 8-10    | int    | mod_num  | Modification number.                 |
        +---------+--------+----------+--------------------------------------+
        | 14-22   | string | mod_date | Date of modification (or release for |
        |         |        |          | new entries).                        |
        +---------+--------+----------+--------------------------------------+
        | 24-28   | string | mod_id   | Identifies this particular           |
        |         |        |          | modification. It links to the archive|
        |         |        |          | used internally by PDB.              |
        +---------+--------+----------+--------------------------------------+
        | 32      | int    | mod_type | An integer identifying the type of   |
        |         |        |          | modification. In case of revisions   |
        |         |        |          | with more than one possible mod_type,|
        |         |        |          | the highest value applicable will be |
        |         |        |          | assigned.                            |
        +---------+--------+----------+--------------------------------------+
        | 40-45   | string | record   | Name of the modified record.         |
        +---------+--------+----------+--------------------------------------+
        | 47-52   | string | record   | Name of the modified record.         |
        +---------+--------+----------+--------------------------------------+
        | 54-59   | string | record   | Name of the modified record.         |
        +---------+--------+----------+--------------------------------------+
        | 61-66   | string | record   | Name of the modified record.         |
        +---------+--------+----------+--------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.mod_num = int(line[7:10].strip())
        self.mod_date = line[13:22].strip()
        self.mod_id = line[23:28].strip()
        mod_type = line[31].strip()
        if mod_type:
            self.mod_type = int(mod_type)
        self.records = [line[39:45].strip()]
        self.records.append(line[46:52].strip())
        self.records.append(line[53:59].strip())
        self.records.append(line[60:66].strip())


@register_line_parser
class AUTHOR(BaseRecord):
    """AUTHOR field

    The AUTHOR record contains the names of the people responsible for the
    contents of the entry.
    """

    def __init__(self, line):
        """Initialize by parsing a line

        +---------+--------+-------------+-----------------------------------+
        | COLUMNS | TYPE   | FIELD       | DEFINITION                        |
        +=========+========+=============+===================================+
        | 11-70   | string | author_list | List of the author names,         |
        |         |        |             | separated by commas               |
        +---------+--------+-------------+-----------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.author_list = line[10:70].strip()


@register_line_parser
class EXPDTA(BaseRecord):
    """EXPDTA field

    The EXPDTA record identifies the experimental technique used. This may
    refer to the type of radiation and sample, or include the spectroscopic
    or modeling technique. Permitted values include:

    * ELECTRON DIFFRACTION
    * FIBER DIFFRACTION
    * FLUORESCENCE TRANSFER
    * NEUTRON DIFFRACTION
    * NMR
    * THEORETICAL MODEL
    * X-RAY DIFFRACTION
    """

    def __init__(self, line):
        """Initialize by parsing a line

        +---------+--------+-----------+-------------------------------------+
        | COLUMNS | TYPE   | FIELD     | DEFINITION                          |
        +=========+========+===========+=====================================+
        | 11-70   | string | technique | The experimental technique(s) with  |
        |         |        |           | optional comment describing the     |
        |         |        |           | sample or experiment                |
        +---------+--------+-----------+-------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.technique = line[10:70].strip()


@register_line_parser
class KEYWDS(BaseRecord):
    """KEYWDS field

    The KEYWDS record contains a set of terms relevant to the entry. Terms
    in the KEYWDS record provide a simple means of categorizing entries and
    may be used to generate index files. This record addresses some of the
    limitations found in the classification field of the HEADER record. It
    provides the opportunity to add further annotation to the entry in a
    concise and computer-searchable fashion.
    """

    def __init__(self, line):
        """Initialize by parsing a line

        +---------+--------+--------+----------------------------------------+
        | COLUMNS | TYPE   | FIELD  | DEFINITION                             |
        +=========+========+========+========================================+
        | 11-70   | string | keywds | Comma-separated list of keywords       |
        |         |        |        | relevant to the entry                  |
        +---------+--------+--------+----------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.keywds = line[10:70].strip()


@register_line_parser
class SOURCE(BaseRecord):
    """SOURCE field

    The SOURCE record specifies the biological and/or chemical source of
    each biological molecule in the entry. Sources are described by both
    the common name and the scientific name, e.g., genus and species.
    Strain and/or cell-line for immortalized cells are given when they help
    to uniquely identify the biological entity studied.
    """

    def __init__(self, line):
        """Initialize by parsing a line

        +---------+--------+--------+----------------------------------------+
        | COLUMNS | TYPE   | FIELD  | DEFINITION                             |
        +=========+========+========+========================================+
        | 11-70   | string | source | Identifies the source of the           |
        |         |        |        | macromolecule in a token: value format |
        +---------+--------+--------+----------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.source = line[10:70].strip()


@register_line_parser
class COMPND(BaseRecord):
    """COMPND field

    The COMPND record describes the macromolecular contents of an entry.
    Each macromolecule found in the entry is described by a set of token:
    value pairs, and is referred to as a COMPND record component. Since the
    concept of a molecule is difficult to specify exactly, PDB staff may
    exercise editorial judgment in consultation with depositors in
    assigning these names.

    For each macromolecular component, the molecule name, synonyms, number
    assigned by the Enzyme Commission (EC), and other relevant details are
    specified.
    """

    def __init__(self, line):
        """Initialize by parsing a line

        +---------+--------+----------+--------------------------------------+
        | COLUMNS | TYPE   | FIELD    | DEFINITION                           |
        +=========+========+==========+======================================+
        | 11-70   | string | compound | Description of the molecular list    |
        |         |        |          | components.                          |
        +---------+--------+----------+--------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.compound = line[10:70].strip()


@register_line_parser
class CAVEAT(BaseRecord):
    """CAVEAT field

    CAVEAT warns of severe errors in an entry. Use caution when using an entry
    containing this record.
    """

    def __init__(self, line):
        """Initialize by parsing line.

        +---------+--------+---------+---------------------------------------+
        | COLUMNS | TYPE   | FIELD   | DEFINITION                            |
        +=========+========+=========+=======================================+
        | 12-15   | string | id_code | PDB ID code of this entry.            |
        +---------+--------+---------+---------------------------------------+
        | 20-70   | string | comment | Free text giving the reason for the   |
        |         |        |         | CAVEAT.                               |
        +---------+--------+---------+---------------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.id_code = line[11:15].strip()
        self.comment = line[19:70].strip()


@register_line_parser
class TITLE(BaseRecord):
    """TITLE field

    The TITLE record contains a title for the experiment or analysis that
    is represented in the entry. It should identify an entry in the PDB in
    the same way that a title identifies a paper.
    """

    def __init__(self, line):
        """Initialize by parsing a line.

        +---------+--------+-------+--------------------------+
        | COLUMNS | TYPE   | FIELD | DEFINITION               |
        +=========+========+=======+==========================+
        | 11-70   | string | title | Title of the experiment  |
        +---------+--------+-------+--------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.title = line[10:70].strip()


@register_line_parser
class OBSLTE(BaseRecord):
    """OBSLTE field

    This record acts as a flag in an entry which has been withdrawn from the
    PDB's full release. It indicates which, if any, new entries have replaced
    the withdrawn entry.

    The format allows for the case of multiple new entries replacing one
    existing entry.
    """

    def __init__(self, line):
        """Initialize by parsing a line.

        +---------+--------+--------------+----------------------------------+
        | COLUMNS | TYPE   | FIELD        | DEFINITION                       |
        +=========+========+==============+==================================+
        | 12-20   | string | replace_date | Date that this entry was         |
        |         |        |              | replaced.                        |
        +---------+--------+--------------+----------------------------------+
        | 22-25   | string | id_code      | ID code of this entry.           |
        +---------+--------+--------------+----------------------------------+
        | 32-35   | string | rid_code     | ID code of entry that replaced   |
        |         |        |              | this one.                        |
        +---------+--------+--------------+----------------------------------+
        | 37-40   | string | rid_code     | ID code of entry that replaced   |
        |         |        |              | this one.                        |
        +---------+--------+--------------+----------------------------------+
        | 42-45   | string | rid_code     | ID code of entry that replaced   |
        |         |        |              | this one.                        |
        +---------+--------+--------------+----------------------------------+
        | 47-50   | string | rid_code     | ID code of entry that replaced   |
        |         |        |              | this one.                        |
        +---------+--------+--------------+----------------------------------+
        | 52-55   | string | rid_code     | ID code of entry that replaced   |
        |         |        |              | this one.                        |
        +---------+--------+--------------+----------------------------------+
        | 57-60   | string | rid_code     | ID code of entry that replaced   |
        |         |        |              | this one.                        |
        +---------+--------+--------------+----------------------------------+
        | 62-65   | string | rid_code     | ID code of entry that replaced   |
        |         |        |              | this one.                        |
        +---------+--------+--------------+----------------------------------+
        | 67-70   | string | rid_code     | ID code of entry that replaced   |
        |         |        |              | this one.                        |
        +---------+--------+--------------+----------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.replace_date = line[11:20].strip()
        self.id_code = line[21:25].strip()
        self.replace_id_codes = [line[31:35].strip()]
        self.replace_id_codes.append(line[36:40].strip())
        self.replace_id_codes.append(line[41:45].strip())
        self.replace_id_codes.append(line[46:50].strip())
        self.replace_id_codes.append(line[51:55].strip())
        self.replace_id_codes.append(line[56:60].strip())
        self.replace_id_codes.append(line[61:65].strip())
        self.replace_id_codes.append(line[67:70].strip())


@register_line_parser
class HEADER(BaseRecord):
    """HEADER field

    The HEADER record uniquely identifies a PDB entry through the id_code
    field. This record also provides a classification for the entry. Finally,
    it contains the date the coordinates were deposited at the PDB.
    """

    def __init__(self, line):
        """Initialize by parsing a line.

        +---------+--------+----------------+--------------------------------+
        | COLUMNS | TYPE   | FIELD          | DEFINITION                     |
        +=========+========+================+================================+
        | 11-50   | string | classification | Classifies the molecule(s)     |
        +---------+--------+----------------+--------------------------------+
        | 51-59   | string | dep_date       | Deposition date.  This is the  |
        |         |        |                | date the coordinates were      |
        |         |        |                | received by the PDB            |
        +---------+--------+----------------+--------------------------------+
        | 63-66   | string | id_code        | This identifier is unique wihin|
        |         |        |                | within PDB                     |
        +---------+--------+----------------+--------------------------------+

        :param line:  line with PDB class
        :type line:  str
        """
        super().__init__(line)
        self.classification = line[10:50].strip()
        self.dep_date = line[50:59].strip()
        self.id_code = line[62:66].strip()


def read_atom(line):
    """If the ATOM/HETATM is not column-formatted, try to get some information
    by parsing whitespace from the right.  Look for five floating point
    numbers followed by the residue number.

    :param line:  the line to parse
    :type line:  str
    """
    # Try to find 5 consecutive floats
    words = str.split(line)
    size = len(words) - 1
    consec = 0
    iword = 0
    for i in range(size):
        entry = words[size - i]
        try:
            _ = float(entry)
            consec += 1
            if consec == 5:
                iword = i
                break
        except ValueError:
            consec = 0

    record = line[0:6].strip()
    newline = line[0:22]
    newline = newline + str.rjust(words[size - iword - 1], 4)
    newline = newline + str.rjust("", 3)
    newline = newline + str.rjust(words[size - iword], 8)
    newline = newline + str.rjust(words[size - iword + 1], 8)
    newline = newline + str.rjust(words[size - iword + 2], 8)
    newline = newline + str.rjust(words[size - iword + 3], 6)
    newline = newline + str.rjust(words[size - iword + 4], 6)
    klass = LINE_PARSERS[record]
    return klass(newline)


def read_pdb(file_):
    """Parse PDB-format data into array of Atom objects.

    :param file_:  open File-like object
    :type file_:  file
    :return:  (a list of objects from this module, a list of record names that
        couldn't be parsed)
    :rtype:  (list, list)
    """
    pdblist = []  # Array of parsed lines (as objects)
    errlist = []  # List of records we can't parse

    # We can come up with nothing if can't get our file off the web.
    if file_ is None:
        return pdblist, errlist

    while True:
        line = file_.readline().strip()
        if line == "":
            break

        # We assume we have a method for each PDB record and can therefore
        # parse them automatically
        record = ""
        try:
            record = line[0:6].strip()
            if record not in errlist:
                klass = LINE_PARSERS[record]
                obj = klass(line)
                pdblist.append(obj)
        except (KeyError, ValueError) as details:
            if record not in ["HETATM", "ATOM"]:
                errlist.append(record)
                _LOGGER.error(f"Error parsing line: {details}")
                _LOGGER.error(f"<{line.strip()}>")
                _LOGGER.error(
                    f"Truncating remaining errors for record type:{record}"
                )
            else:
                raise details
        except IndexError as details:
            if record in ["ATOM", "HETATM"]:
                try:
                    obj = read_atom(line)
                    pdblist.append(obj)
                except IndexError as details:
                    _LOGGER.error(f"Error parsing line: {details},")
                    _LOGGER.error(f"<{line.strip()}>")
            elif record in ["SITE", "TURN"]:
                pass
            elif record in ["SSBOND", "LINK"]:
                _LOGGER.error("Warning -- ignoring record:")
                _LOGGER.error(f"<{line.strip()}>")
            else:
                _LOGGER.error(f"Error parsing line: {details},")
                _LOGGER.error(f"<{line.strip()}>")
    return pdblist, errlist
