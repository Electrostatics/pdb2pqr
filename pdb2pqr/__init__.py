"""PDB2PQR

This package takes a PDB file as input and performs optimizations before
yielding a new PDB-style file as output.

For more information, see http://www.poissonboltzmann.org/
"""
import logging
from sys import version_info
assert version_info >= (3, 5)
from .main import main, build_parser


_LOGGER = logging.getLogger(__name__)
logging.captureWarnings(True)


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    logging.captureWarnings(True)
    parser = build_parser()
    args = parser.parse_args()
    main(args)
