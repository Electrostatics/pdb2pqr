"""PDB2PQR

This package takes a PDB file as input and performs optimizations before
yielding a new PDB-style file as output.

For more information, see http://www.poissonboltzmann.org/
"""

import logging

from ._version import __version__
from .main import main as pdb2pqr_main
from .main import run_pdb2pqr

_LOGGER = logging.getLogger(__name__)
logging.captureWarnings(capture=True)


if __name__ == "__main__":
    pdb2pqr_main()
