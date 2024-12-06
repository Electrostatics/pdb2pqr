"""PDB2PQR

This package takes a PDB file as input and performs optimizations before
yielding a new PDB-style file as output.

For more information, see http://www.poissonboltzmann.org/
"""

import logging

from pdb2pqr.main import main

_LOGGER = logging.getLogger(__name__)
logging.captureWarnings(capture=True)


if __name__ == "__main__":
    main()
