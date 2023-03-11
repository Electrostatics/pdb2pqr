===============
Getting PDB2PQR
===============

.. note::

   *Before you begin!* PDB2PQR funding is dependent on your help for continued development and support. Please `register <http://eepurl.com/by4eQr>`_ before using the software so we can accurately report the number of users to our funding agencies.

-----------
Web servers
-----------

Most functionality is available through our online web servers.

The PDB2PQR web server offers a simple way to use both APBS and PDB2PQR without the need to download and install additional programs.

After `registering <http://eepurl.com/by4eQr>`_, please visit http://server.poissonboltzmann.org/ to access the web server.

------------------------------
Python Package Installer (PIP)
------------------------------

Most users who want to use the software offline should install it via :program:`pip` with the following command

.. code-block:: bash

   pip install pdb2pqr

from within your favorite `virtual environment <https://docs.python.org/3/tutorial/venv.html>`_.

The PIP package provides the :mod:`pdb2pqr` Python module as well as the :program:`pdb2pqr` program which can be used from the command line.

-----------------------------
Installation from source code
-----------------------------

You can also download the source code from `GitHub <https://github.com/Electrostatics/pdb2pqr>`_ (we recommend using `a tagged release <https://github.com/Electrostatics/pdb2pqr/releases>`_) and building the code yourself with:

.. code-block:: bash

   pip install .

from the top-level of the source code directory.
Note that developers may want to install the code in "editable" mode with

.. code-block:: bash

   pip install -e .

.. _testing-label:

^^^^^^^
Testing
^^^^^^^

The software can be tested for correct functioning via `Coverage.py <https://coverage.readthedocs.io/en/coverage-5.2/>`_ with

.. code-block::

    coverage run -m pytest

or `pytest <https://docs.pytest.org/en/stable/>`_

.. code-block::

    python -m pytest

from the top level of the source directory. 
These commands run basic tests; more extensive testing can be performed by adding the ``--run-long`` option to these commands.
