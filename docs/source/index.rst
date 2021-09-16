.. pdb2pqr documentation master file, created by
   sphinx-quickstart on Thu Jul 16 13:30:31 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:Release: |release|
:Date: |today|

PDB2PQR
=======

========
Overview
========

The use of continuum solvation methods such as `APBS <http://www.poissonboltzmann.org>`_ requires accurate and complete structural data as well as force field parameters such as atomic charges and radii.
Unfortunately, the limiting step in continuum electrostatics calculations is often the addition of missing atomic coordinates to molecular structures from the Protein Data Bank and the assignment of parameters to these structures.
To address this problem, we have developed PDB2PQR.
This software automates many of the common tasks of preparing structures for continuum solvation calculations as well as many other types of biomolecular structure modeling, analysis, and simulation. These tasks include:

* Adding a limited number of missing heavy (non-hydrogen) atoms to biomolecular structures.
* Estimating titration states and protonating biomolecules in a manner consistent with favorable hydrogen bonding.
* Assigning charge and radius parameters from a variety of force fields.
* Generating “PQR” output compatible with several popular computational modeling and analysis packages.
* This service is intended to facilitate the setup and execution of electrostatics calculations for both experts and non-experts and thereby broaden the accessibility of biomolecular solvation and electrostatics analyses to the biomedical community.

========
Contents
========

.. toctree::
   :maxdepth: 2

   getting
   using/index
   help
   supporting
   extending
   formats/index
   api/index
   releases

------------------
Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
