"""Routines for PDB2PQR

This module contains the protein object used in PDB2PQR and associated methods

Authors:  Todd Dolinsky, Yong Huang
"""
import string
import logging
from .aa import ALA, ARG, ASN, ASP, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU
from .aa import LIG, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL, WAT
from .na import Nucleic
from .structures import Chain, Residue
from .pdb import TER, ATOM, HETATM, END, MODEL
from .forcefield import Forcefield


_LOGGER = logging.getLogger(__name__)


class Protein(object):
    """Protein class

    The protein class represents the parsed PDB, and provides a hierarchy of
    information - each Protein contains a list of Chain objects as provided in
    the PDB file.  Each Chain then contains its associated list of Residue
    objects, and each Residue contains a list of Atom objects, completing the
    hierarchy.
    """

    def __init__(self, pdblist, definition):
        """Initialize using parsed PDB file

        Args:
            pdblist: List of Classes of PDB lines as created
        """
        self.chainmap = {}
        self.chains = []
        self.residues = []
        self.referencemap = definition.map
        self.patchmap = definition.patches

        chain_dict = {}
        previous_atom = None
        residue = []
        num_models = 0
        num_chains = 1
        count = 0

        for record in pdblist: # Find number of chains
            if isinstance(record, TER):
                num_chains += 1

        for record in pdblist:
            if isinstance(record, (ATOM, HETATM)):
                if record.chain_id == "":
                    if num_chains > 1 and record.res_name not in ["WAT", "HOH"]:
                        # Assign a chain ID
                        record.chain_id = string.ascii_uppercase[count]

                chain_id = record.chain_id
                res_seq = record.res_seq
                ins_code = record.ins_code

                if previous_atom is None:
                    previous_atom = record

                if chain_id not in chain_dict:
                    my_chain = Chain(chain_id)
                    chain_dict[chain_id] = my_chain

                if res_seq != previous_atom.res_seq or \
                      ins_code != previous_atom.ins_code or \
                      chain_id != previous_atom.chain_id:
                    my_residue = self.create_residue(residue, previous_atom.res_name)
                    chain_dict[previous_atom.chain_id].add_residue(my_residue)
                    residue = []

                residue.append(record)
                previous_atom = record

            elif isinstance(record, END):
                my_residue = self.create_residue(residue, previous_atom.res_name)
                chain_dict[previous_atom.chain_id].add_residue(my_residue)
                residue = []

            elif isinstance(record, MODEL):
                num_models += 1
                if residue == []: continue
                if num_models > 1:
                    my_residue = self.create_residue(residue,
                                                     previous_atom.res_name)
                    chain_dict[previous_atom.chain_id].add_residue(my_residue)
                    break

            elif isinstance(record, TER):
                count += 1

        if residue != [] and num_models <= 1:
            my_residue = self.create_residue(residue, previous_atom.res_name)
            chain_dict[previous_atom.chain_id].add_residue(my_residue)

        # Keep a map for accessing chains via chain_id
        self.chainmap = chain_dict.copy()

        # Make a list for sequential ordering of chains
        if "" in chain_dict:
            chain_dict["ZZ"] = chain_dict[""]
            del chain_dict[""]

        keys = list(chain_dict.keys())
        keys.sort()

        for key in keys:
            self.chains.append(chain_dict[key])

        for chain in self.chains:
            for residue in chain.residues:
                self.residues.append(residue)

    def create_residue(self, residue, resname):
        """Create a residue object.

        If the resname is a known residue type, try to make that specific
        object, otherwise just make a standard residue object.

        Args:
            residue:  A list of atoms (list)
            resname:  The name of the residue (string)
        Returns:
            residue:  The residue object (Residue)
        """
        try:
            refobj = self.referencemap[resname]
            if refobj.name != resname:
                klass = globals()[refobj.name]
                residue = klass(residue, refobj)
                residue.reference = refobj
            else:
                klass = globals()[resname]
                residue = klass(residue, refobj)
        except (KeyError, NameError) as err:
            _LOGGER.debug("Parsing %s as new residue", resname)
            residue = Residue(residue)
        return residue

    def print_atoms(self, atomlist, chainflag=False, pdbfile=False):
        """Get the text for the entire protein

        Args:
            atomlist:  The list of atoms to include (list)
            chainflag:  Flag whether to print chainid or not
        Returns:
            text:  list of (stringed) atoms (list)
        """
        self.reserialize()
        text = []
        currentchain_id = None
        for atom in atomlist:
            # Print the "TER" records between chains
            if currentchain_id is None:
                currentchain_id = atom.chain_id
            elif atom.chain_id != currentchain_id:
                currentchain_id = atom.chain_id
                text.append("TER\n")

            if pdbfile is True:
                text.append("%s\n" % atom.get_pdb_string())
            else:
                text.append("%s\n" % atom.get_pqr_string(chainflag=chainflag))
        text.append("TER\nEND")
        return text

    def create_html_typemap(self, definition, outfilename):
        """Create an HTML typemap file at the desired location.

        If a type cannot be found for an atom a blank is listed.

        Args:
            definition: The definition objects.
            outfilename:  The name of the file to write (string)
        """
        # Cache the initial atom numbers
        numcache = {}
        for atom in self.atoms:
            numcache[atom] = atom.serial
        self.reserialize()

        amberff = Forcefield("amber", definition, None)
        charmmff = Forcefield("charmm", definition, None)

        with open(outfilename, "w") as file_:
            file_.write("<HTML>\n")
            file_.write("<HEAD>\n")
            file_.write("<TITLE>PQR Typemap (beta)</TITLE>\n")
            file_.write("</HEAD>\n")
            file_.write("<BODY>\n")
            file_.write(("<H3>This is a developmental page including the atom "
                         "type for the atoms in the PQR file.</H3><P>\n"))
            file_.write("<TABLE CELLSPACING=2 CELLPADDING=2 BORDER=1>\n")
            file_.write(("<tr><th>Atom Number</th><th>Atom Name</th><th>Residue "
                         "Name</th><th>Chain ID</th><th>AMBER Atom Type</th><th>"
                         "CHARMM Atom Type</th></tr>\n"))

            for atom in self.atoms:
                if isinstance(atom.residue, (Amino, WAT, Nucleic)):
                    resname = atom.residue.ffname
                else:
                    resname = atom.residue.name
                ambergroup = amberff.get_group(resname, atom.name)
                charmmgroup = charmmff.get_group(resname, atom.name)
                file_.write(("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td>"
                             "<td>%s</td><td>%s</td></tr>\n") % (atom.serial, atom.name,
                                                                 resname, atom.chain_id,
                                                                 ambergroup, charmmgroup))

            file_.write("</table>\n")
            file_.write("</BODY></HTML>\n")

        # Return the original numbers back
        for atom in self.atoms:
            atom.serial = numcache[atom]

    def reserialize(self):
        """Generate new serial numbers for atoms in the protein"""
        count = 1
        for atom in self.atoms:
            atom.serial = count
            count += 1

    @property
    def atoms(self):
        """Return all Atom objects in list format

        Returns:
            atomlist:  List of Atom objects in the protein (list)
        """
        atomlist = []
        for chain in self.chains:
            for atom in chain.atoms:
                atomlist.append(atom)
        return atomlist

    @property
    def charge(self):
        """Get the total charge on the protein

        NOTE:  Since the misslist is used to identify incorrect charge
        assignments, this routine does not list the 3 and 5 termini of nucleic
        acid chains as having non-integer charge even though they are
        (correctly) non-integer.

        Returns:
            misslist:  List of residues with non-integer charges (list)
            charge:  The total charge on the protein (float)
        """
        charge = 0.0
        misslist = []
        for chain in self.chains:
            for residue in chain.residues:
                rescharge = residue.charge
                charge += rescharge
                if isinstance(residue, Nucleic):
                    if residue.is3term or residue.is5term:
                        continue
                if float("%i" % rescharge) != rescharge:
                    misslist.append(residue)
        return misslist, charge

    def __str__(self):
        output = []
        for chain in self.chains:
            output.append(chain.get_summary())
        return ' '.join(output)
