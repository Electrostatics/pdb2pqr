PDB2PQR NAMES files
=====================

Much of the difficulty in adding a new forcefield to PDB2PQR depends on the naming scheme used in that forcefield.

===============
XML file format
===============

To start, either a flat file or XML file containing the desired forcefield's parameters should be made - see :file:`AMBER.DAT` and :file:`AMBER.xml` for examples.
If the forcefield's naming scheme matches the canonical naming scheme, that's all that is necessary.
If the naming schemes differ, however, conversions must be made. These are made in the :file:`*.names` file (see :file:`CHARMM.names`, for example).
In this file you will see sections like:</p>

.. code-block:: xml

   <residue>
     <name>WAT</name>
     <useresname>TP3M</useresname>
     <atom>
       <name>O</name>
       <useatomname>OH2</useatomname>
     </atom>
   </residue>

This section tells PDB2PQR that for the oxygen atom O in WAT, CHARMM uses the names OH2 and TP3M, respectively.
When the XML file is read in, PDB2PQR ensures that the WAT/O pair points to TP3M/OH2 such that the appropriate parameters are returned.
But for naming schemes that greatly differ from the PDB2PQR canonical naming scheme, this could get really ugly.
As a result, PDB2PQR can use regular expressions to simplify the renaming process, i.e.:

.. code-block:: xml

   <residue>
     <name>[NC]?...$</name>
     <atom>
       <name>H</name>
       <useatomname>HN</useatomname>
     </atom>
   </residue>

This section of code will ensure that the H atom of all canonical residue names that match the :regexp:`[NC]?...$` regular expression point to HN instead.
This regular expression matches all three-letter residue names, residue names with an 'N' prepended (N-Termini), and residue names with a 'C' prepended (C-Termini).
For twenty amino acids, sixty residue name changes can all be done by a single section.
The use of regular expressions is therefore a much more powerful method of handling naming scheme differences than working on a one to one basis.

There are a few other additional notes when using the :file:`.names` file.
First, the ``$group`` variable is used to denote the matching group of a regular expression, for instance:

.. code-block:: xml

   <residue>
     <name>HI([PDE])$</name>
     <useresname>HS$group</useresname>
   </residue>

This section replaces HIP/HID/HIE with HSP/HSD/HSE by first matching the HI([PDE])$ regular expression and then using the group that is enclosed by parantheses to fill in the name to use.

Second, sections are cumulative - since CHARMM, for instance, has a patch-based naming scheme, one single canonical residue name can map to multiple forcefield-scheme names. Let's look at how to map an SS-bonded Cysteine (canonical name CYX) to the CHARMM naming scheme:

.. code-block:: xml

   <residue>
     <name>CYX</name>
     <useresname>CYS</useresname>
   </residue>
   <residue>
     <name>CYX</name>
     <useresname>DISU</useresname>
     <atom>
       <name>CB</name>
       <useatomname>1CB</useatomname>
     </atom>
     <atom>
       <name>SG</name>
       <useatomname>1SG</useatomname>
     </atom>
   </residue>

The CYX residue is first mapped to CHARMM's CYS, and then to CHARMM's DISU object.
All atom names that are found in DISU overwrite those found in CYS - in effect, the DISU patch is applied to CYS, yielding the desired CYX.
This cumulative can be repeated as necessary.

=========================
Caveats about atom naming
=========================

In an ideal world each individual residue and atom would have a standard, distinct name.
Unfortunately `several naming schemes for atoms exist <http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl>`_, particularly for hydrogens.
As such, in order to detect the presence/absence of atoms in a biomolecule, an internal canonical naming scheme is used.
The naming scheme used in PDB2PQR is the one recommended by the PDB itself, and derives from the IUPAC naming recommendations [#naming]_

.. [#naming] J. L. Markley, et al., "Recommendations for the Presentation of NMR Structures of Proteins and Nucleic Acids," Pure & Appl. Chem., 70 (1998): 117-142.  DOI:`10.1046/j.1432-1327.1998.2560001.x <http://dx.doi.org/10.1046/j.1432-1327.1998.2560001.x>`_

This canonical naming scheme is used as the default PDB2PQR output.
All conversions in PDB2PQR use the internal canonical naming scheme to determine distinct atom names.
In previous versions of PDB2PQR, these conversions were stored in long lists of if statements, but for transparency and editing this is a bad thing.
Instead, all conversions can now be found in XML as described above.

There are a few additions to the canonical naming scheme, mirrored after the AMBER naming scheme (chosen since for the most part it follows the IUPAC recommendations).
These changes are made in :file:`PATCHES.xml`, and allow any of the following to be patched as necessary as well as detected on input:

:regexp:`N*`
  N-Terminal Residue (i.e. NALA, NLEU)
:regexp:`NEUTRAL-N*`
  Neutral N-Terminal Residue
:regexp:`C*`
  C-Terminal Residue (i.e. CLYS, CTYR)
:regexp:`NEUTRAL-C*`
  Neutral C-Terminal Residue
:regexp:`*5`
  5-Terminus for Nucleic Acids (i.e. DA5)
:regexp:`*3`
  3-Terminus for Nucleic Acids (i.e. DA3)
``ASH``
  Neutral ASP
``CYX``
  SS-bonded CYS
``CYM``
  Negative CYS
``GLH``
  Neutral GLU
``HIP``
  Positive HIS
``HID``
  Neutral HIS, proton HD1 present
``HIE``
  Neutral HIS, proton HE2 present
``LYN``
  Neutral LYS
``TYM``
  Negative TYR
