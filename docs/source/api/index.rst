.. _api-label:

=============
API Reference
=============

.. currentmodule:: pdb2pqr
.. module:: pdb2pqr


The :program:`pdb2pqr` command provides a command-line interface to PDB2PQR's functionality.
It is built on classes and functions in the :mod:`pdb2pqr` module.
The API of :mod:`pdb2pqr` is documented here for developers who might want to directly use the PDB2PQR code.

.. Note::

   The API is still changing and there is currently no guarantee that
   it will remain stable between minor releases.

--------------------------
Forcefield support modules
--------------------------

.. toctree::
   :maxdepth: 2

   definitions
   forcefield

---------------------------
Molecular structure modules
---------------------------

.. toctree::
   :maxdepth: 2

   aa
   hydrogens
   ligand
   na
   biomolecule
   residue
   structures
   topology

-----------
I/O modules
-----------

.. toctree::
   :maxdepth: 2

   cif
   io
   inputgen
   pdb

-------------
Other modules
-------------

.. toctree::
   :maxdepth: 2

   cells
   config
   debump
   main
   psize
   quatfit
   run
   utilities
