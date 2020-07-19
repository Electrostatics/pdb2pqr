.. _pqr:

PQR molecular structure format
==============================

This format is a modification of the PDB format which allows users to add charge and radius parameters to existing PDB data while keeping it in a format amenable to visualization with standard molecular graphics programs.
The origins of the PQR format are somewhat uncertain, but has been used by several computational biology software programs, including MEAD and AutoDock.
UHBD uses a very similar format called QCD.

APBS reads very loosely-formatted PQR files: all fields are whitespace-delimited rather than the strict column formatting mandated by the PDB format.
This more liberal formatting allows coordinates which are larger/smaller than ± 999 Å.
APBS reads data on a per-line basis from PQR files using the following format:::

  Field_name Atom_number Atom_name Residue_name Chain_ID Residue_number X Y Z Charge Radius

where the whitespace is the most important feature of this format.
The fields are:

``Field_name``
  A string which specifies the type of PQR entry and should either be ATOM or HETATM in order to be parsed by APBS.

``Atom_number``
  An integer which provides the atom index.

``Atom_name``
  A string which provides the atom name.

``Residue_name``
  A string which provides the residue name.

``Chain_ID``
  An optional string which provides the chain ID of the atom.
  Note that chain ID support is a new feature of APBS 0.5.0 and later versions.

``Residue_number``
  An integer which provides the residue index.

``X Y Z``
  3 floats which provide the atomic coordinates (in Å)

``Charge``
  A float which provides the atomic charge (in electrons).

``Radius``
  A float which provides the atomic radius (in Å).

Clearly, this format can deviate wildly from PDB due to the use of whitespaces rather than specific column widths and alignments.
This deviation can be particularly significant when large coordinate values are used.
However, in order to maintain compatibility with most molecular graphics programs, the PDB2PQR program and the utilities provided with APBS attempt to preserve the PDB format as much as possible.
