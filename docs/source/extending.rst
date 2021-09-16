=================
Extending PDB2PQR
=================

--------------------------------
Adding new forcefield parameters
--------------------------------

If you are just adding the parameters of a few residues and atoms to an existing forcefield (e.g., AMBER), you can open the forcefield data file distributed with PDB2PQR (:file:`dat/AMBER.DAT`) directly and add your parameters.
After the parameter addition, save the force field data file with your changes.
You should also update the corresponding .names file (:file:`dat/AMBER.names`) if your added residue or atom naming scheme is different from the PDB2PQR canonical naming scheme.

---------------------------------
Adding an entirely new forcefield
---------------------------------

The following steps outline how to add a new force field to PDB2PQR.

You will need to generate a forcefield data file (e.g., :file:`myff.DAT`) and, if your atom naming scheme of the forcefield is different from the PDB2PQR canonical naming scheme, you will also need to provide a names files (:file:`myFF.names`).
The format of the names file is described in :doc:`/formats/xml-names`.
It is recommended to build your own forcefield data and names files based on existing PDB2PQR :file:`.DAT` and :file:`.names` examples provided with PDB2PQR in the :file:`dat` directory.
After finishing your forcefield data file and names file, these can be used with either the command line or the web server versions of PDB2PQR.

------------------------
Helping with development
------------------------

^^^^^^^^^^^^^^^^^^^^^^^^
Adding new functionality
^^^^^^^^^^^^^^^^^^^^^^^^

PDB2PQR welcomes new contributions; the software API is documented in :ref:`api-label`.
To contribute code, submit a *pull request* against the master branch in the `PDB2PQR repository <https://github.com/Electrostatics/pdb2pqr>`_.
Please be sure to run PDB2PQR tests, as described in :ref:`testing-label`, before submitting new code.

^^^^^^^^^^^^^^^^^^^^^^^^
Helping with to-do items
^^^^^^^^^^^^^^^^^^^^^^^^

A list of "to-do" items for the code is available in `GitHub Issues <https://github.com/Electrostatics/pdb2pqr/issues>`_.
A loosely maintained list auto-generated from the documentation is also presented below.

.. todolist::
