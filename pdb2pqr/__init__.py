"""PDB2PQR

This package takes a PDB file as input and performs optimizations before
yielding a new PDB-style file as output.

For more information, see http://www.poissonboltzmann.org/
"""
import logging
from sys import version_info
from .main import main_driver, build_main_parser
from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

assert version_info >= (3, 5)

_LOGGER = logging.getLogger(__name__)
logging.captureWarnings(True)

if __name__ == "__main__":
    PARSER = build_main_parser()
    ARGS = PARSER.parse_args()
    main_driver(ARGS)
