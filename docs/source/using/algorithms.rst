--------------------------
Algorithms used by PDB2PQR
--------------------------

^^^^^^^^^
Debumping
^^^^^^^^^

Unless otherwise instructed with ``--nodebump``, PDB2PQR will attempt to remove steric clashes (debump) between residues.

To determine if a residue needs to be debumped, PDB2PQR compares its atoms to all nearby atoms.
With the exception of donor/acceptor pairs and CYS residue SS bonded pairs, a residue needs to be debumped if any of its atoms are within cutoff distance of any other atoms.
The cutoff is 1.0 angstrom for hydrogen/hydrogen collisions, 1.5 angstrom for hydrogen/heavy collisions, and 2.0 angstrom otherwise.

Considering the atoms that are conflicted, PDB2PQR changes selected dihedral angle configurations in increments of 5.0 degrees, looking for positions where the residue does not conflict with other atoms.
If modifying a dihedral angle does not result in a debumped configuration then the dihedral angle is reset and the next one is tried.
If 10 angles are tried without success the algorithm reports failure.

.. warning::

   It should be noted that this is not an optimal solution.
   This method is not guaranteed to find a solution if it exists and will accept the first completely debumped state found, not the optimal state.

   Additionally, PDB2PQR does not consider water atoms when looking for conflicts.

^^^^^^^^^^^^^^^^^^^^^^^^^^
Hydrogen bond optimization
^^^^^^^^^^^^^^^^^^^^^^^^^^

Unless otherwise indicated with ``--noopts``, PDB2PQR will attempt to add hydrogens in a way that optimizes hydrogen bonding.

The hydrogen bonding network optimization seeks, as the name suggests, to optimize the hydrogen bonding network of the biomolecule.
Currently this entails manipulating the following residues:

* Flipping the side chains of HIS (including user defined HIS states), ASN, and GLN residues;
* Rotating the sidechain hydrogen on SER, THR, TYR, and CYS (if available);
* Determining the best placement for the sidechain hydrogen on neutral HIS, protonated GLU, and protonated ASP;
* Optimizing all water hydrogens.

^^^^^^^^^^^^^^^^^^^^^^^^^^
Titration state assignment
^^^^^^^^^^^^^^^^^^^^^^^^^^

Versions 2.1 and earlier of PDB2PQR offered the following methods to assign titration states to molecules at a specified pH.

* PROPKA. This method uses the `PROPKA software <https://github.com/jensengroup/propka>`_ to assign titration states.  More information about PROPKA can be found `on its website <https://github.com/jensengroup/propka>`_.

* PDB2PKA. Uses a Poisson-Boltzmann method to assign titration states. This approach is loosely related to the method described by Nielsen and Vriend (2001) doi:`10.1002/prot.10 <https://doi.org/10.1002/prot.1053>`_.

The current version
(|release|)
of PDB2PQR currently only supports PROPKA while we address portability issues in PDB2PKA.

PDB2PQR has the ability to recognize certain protonation states and keep them fixed during optimization.
To use this feature manually rename the residue name in the PDB file as follows:

* Neutral ASP: ASH
* Negative CYS: CYM
* Neutral GLU: GLH
* Neutral HIS:   HIE or HSE (epsilon-protonated); HID or HSD (delta-protonated)
* Positive HIS: HIP or HSP
* Neutral LYS: LYN
* Negative TYR: TYM

PDB2PQR is unable to assign charges and radii when they are not available in the forcefield - thus this warning message will occur for most ligands unless a ``MOL2`` file is provided for the ligand with the ``--ligand`` option.
Occasionally this message will occur in error for a standard amino acid residue where an atom or residue may be misnamed.
However, some of the protonation states derived from the PROPKA results are not supported in the requested forcefield and thus PDB2PQR is unable to get charges and radii for that state.
PDB2PQR currently supports the following states as derived from PROPKA:

================== ============= ============== =============
Protonation State  AMBER Support CHARMM Support PARSE Support
================== ============= ============== =============
Neutral N-Terminus No            No             Yes
Neutral C-Terminus No            No             Yes
Neutral ARG        No            No             No
Neutral ASP        Yes [#but]_   Yes            Yes
Negative CYS       Yes [#but]_   No             Yes
Neutral GLU        Yes [#but]_   Yes            Yes
Neutral HIS        Yes           Yes            Yes
Neutral LYS        Yes [#but]_   No             Yes
Negative TYR       No            No             Yes
================== ============= ============== =============

.. [#but] Only if residue is not a terminal residue; if the residue is terminal it will not be set to this state.