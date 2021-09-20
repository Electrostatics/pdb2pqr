"""PDB2PQR

This package takes a PDB file as input and performs optimizations before
yielding a new PDB-style file as output.

For more information, see http://www.poissonboltzmann.org/
"""
import logging
from sys import version_info
from pdb2pqr.main import main

assert version_info >= (3, 5)


_LOGGER = logging.getLogger(__name__)
logging.captureWarnings(True)


if __name__ == "__main__":
    main()
