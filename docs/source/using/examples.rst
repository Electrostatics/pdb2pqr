========
Examples
========

In order to perform electrostatics calculations on your biomolecular structure of interest, you need to provide atomic charge and radius information to APBS.
Charges are used to form the biomolecular charge distribution for the Poisson-Boltzmann (PB) equation while the radii are used to construct the dielectric and ionic accessibility functions.

The PDB2PQR web service and software will convert most PDB files into PQR format with some caveats.
Although PDB2PQR can fix some missing heavy atoms in sidechains, it does not currently have the (nontrivial) capability to model in large regions of missing backbone and sidechain coordinates.
Be patient and make certain that the job you submitted to the PDB2PQR website has finished and you have downloaded the resulting PQR file correctly.
It usually takes less than 10 minutes for the job to finish.

These examples assume that you have `registered <http://eepurl.com/by4eQr>`_ and have access to `the PDB2PQR web server <http://server.poissonboltzmann.org>`_.

---------------------------------------------------------------------
Adding hydrogens and assigning parameters with the PDB2PQR web server
---------------------------------------------------------------------

This example uses http://server.poissonboltzmann.org.

^^^^^^^^^^^^^^^^
Pick a structure
^^^^^^^^^^^^^^^^

Start by choosing a PDB file to process.
Either enter the 4-character PDB ID into PDB2PQR or accession number (`1FAS <http://www.rcsb.org/pdb/explore.do?structureId=1FAS>`_ is a good starting choice) or upload your own PDB file.
Note that, if you choose to enter a 4-character PDB ID, PDB2PQR will process all recognizable chains of PDB file as it was deposited in the PDB (e.g., not the biological unit, any related transformations, etc.).

^^^^^^^^^^^^^^^^^
Pick a forcefield
^^^^^^^^^^^^^^^^^

For most applications, the choice is easy: PARSE.
This forcefield has been optimized for implicit solvent calculation and is probably the best choice for visualization of biomolecular electrostatics and many common types of energetic calculations for biomolecules.
However, AMBER and CHARMM may be more appropriate if you are attempting to compare directly to simulations performed with those force fields, require nucleic acid support, are simulating ligands parameterized with those force fields, etc.

It is also possible to upload a user-defined forcefield (e.g., to define radii and charges for ligands or unusual residues).
Please see :doc:`/extending` for more information.

^^^^^^^^^^^^^^^^^^^^
Pick a naming scheme
^^^^^^^^^^^^^^^^^^^^

This choice is largely irrelevant to electrostatics calculations but may be important for some visualization programs.
When in doubt, choose the "Internal naming scheme" which attempts to conform to IUPAC standards.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Reconstruct missing atoms (hydrogens)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under options, be sure the "Ensure that new atoms are not rebuilt too close to existing atoms" and "Optimize the hydrogen bonding network" options are selected.
You can select other options as well, if interested.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Download and view the results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download the resulting PQR file and visualize in a molecular graphics program to examine how the hydrogens were added and how hydrogen bonds were optimized.

--------------------------------------------------
Parameterizing ligands with the PDB2PQR web server
--------------------------------------------------

This section outlines the parameterization of ligands using the PEOE_PB methods (see `Czodrowski et al, 2006 <http://dx.doi.org/10.1002/prot.21110>`_ for more information).

.. warning::

   PDB2PQR is unable to process more than one ligand per structure.
   This limitation includes structures with multiple copies of the same ligand--only one copy of the ligand will be processed.

The PDB structure `1HPX <http://www.rcsb.org/pdb/explore.do?structureId=1hpx>`_ includes HIV-1 protease complexed with an inhibitor at 2.0 Å resolution.
HIV-1 protease has two chains; residue D25 is anionic on one chain and neutral on the other -- these titration states are important in the role of D25 as an acid in the catalytic mechanism.

^^^^^^^^^^^^^^^^^^^
Ignoring the ligand
^^^^^^^^^^^^^^^^^^^

If we don't want to include the ligand, then the process is straightforward:

#. From the `PDB2PQR server web page <http://server.poissonboltzmann.org>`_, enter ``1HPX`` into the PDB ID field.

#. Choose whichever forcefield and naming schemes you prefer.

#. Under options, be sure the "Ensure that new atoms are not rebuilt too close to existing atoms", "Optimize the hydrogen bonding network", and "Use PROPKA to assign protonation states at pH" options are selected. Choose pH 7 for your initial calculations. You can select other options as well, if interested.

#. Hit the "Submit" button.

#. Once the calculations are complete, you should see a web page with a link to the PROPKA output, a new PQR file, and warnings about the ligand KNI (since we didn't choose to parameterize it in this calculation). For comparison, you might download the the `original PDB file <http://www.pdb.org/pdb/explore.do?structureId=1HPX>`_ and compare the PDB2PQR-generated structure with the original to see where hydrogens were placed.

^^^^^^^^^^^^^^^^^^^^^^^^^
Parameterizing the ligand
^^^^^^^^^^^^^^^^^^^^^^^^^

This section outlines the parameterization of ligands using the PEOE_PB methods (see DOI:`10.1002/prot.21110 <http://dx.doi.org/10.1002/prot.21110>`_).

Ligand parameterization currently requires a :doc:`MOL2-format </formats/mol2>` representation of the ligand to provide the necessary bonding information.
MOL2-format files can be obtained through the `PRODRG web server <http://davapc1.bioch.dundee.ac.uk/cgi-bin/prodrg>`_ or some molecular modeling software packages.
PRODRG provides documentation as well as several examples on ligand preparation on its web page.

We're now ready to look at the 1HPV crystal structure from above and parameterize its ligand, KNI-272.

#. From the PDB2PQR server web page, enter ``1HPX`` into the PDB ID field.

#. Choose whichever forcefield and naming schemes you prefer.

#. Under options, be sure the "Ensure that new atoms are not rebuilt too close to existing atoms", "Optimize the hydrogen bonding network", and "Assign charges to the ligand specified in a MOL2 file" options are selected. You can select other options as well, if interested.

#. Hit the "Submit" button.

#. Once the calculations are complete, you should see a web page with a link to the new PQR file with a warning about debumping P81 (but no warnings about ligand parameterization!).

As a second example, we use the PDB structure `1ABF <http://www.rcsb.org/pdb/explore.do?structureId=1abf>`_ of L-arabinose binding protein in complex with a sugar ligand at 1.90 Å resolution.
To parameterize both this protein and its ligand:

#. From the PDB2PQR server web page, enter `1ABF` into the PDB ID field.

#. Choose whichever forcefield and naming schemes you prefer.

#. Under options, be sure the "Ensure that new atoms are not rebuilt too close to existing atoms", "Optimize the hydrogen bonding network", and "Assign charges to the ligand specified in a MOL2 file" options are selected. You can select other options as well, if interested.

#. Hit the "Submit" button.

#. Once the calculations are complete, you should see a web page with a link to the new PQR file with a warning about debumping P66, K295, and K306 (but no warnings about ligand parameterization!).
