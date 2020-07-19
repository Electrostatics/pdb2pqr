.. _api-label:

=============
API Reference
=============

The :program:`pdb2pqr30` command provides a command-line interface to PDB2PQR's functionality.
It is built on classes and functions in the :mod:`pdb2pqr` module.
The API of :mod:`pdb2pqr` is documented here for developers who might want to directly use the PDB2PQR code.

.. Note::

   The API is still changing and there is currently no guarantee that
   it will remain stable between minor releases.

.. currentmodule:: pdb2pqr

.. module:: pdb2pqr

------------------
Forcefield support
------------------

.. autosummary::
   :toctree: api

   definitions
   forcefield

--------------------
Molecular structures
--------------------

.. autosummary::
   :toctree: api

   aa
   hydrogens
   hydrogens.optimize
   hydrogens.structures
   ligand
   ligand.mol2
   ligand.peoe
   ligand.topology
   na
   protein
   residue
   structures
   topology

-----------
I/O support
-----------

.. autosummary::
   :toctree: api

   cif
   io
   inputgen
   pdb

----------------------------------------------
Other modules that need to be better organized
----------------------------------------------

.. autosummary::
   :toctree: api

   cells
   debump
   main
   psize
   quatfit
   run
   utilities
