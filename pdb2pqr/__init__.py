"""PDB2PQR

This package takes a PDB file as input and performs optimizations before
yielding a new PDB-style file as output.

For more information, see http://www.poissonboltzmann.org/
"""
import logging
from sys import version_info
from pathlib import Path
assert version_info >= (3, 5)
from .main import main, build_parser


_LOGGER = logging.getLogger(__name__)
logging.captureWarnings(True)


if __name__ == "__main__":
    PARSER = build_parser()
    ARGS = PARSER.parse_args()
    LOG_CACHE = Path(Path(ARGS.output_pqr).parent, Path(ARGS.output_pqr).stem + '.log')
    logging.basicConfig(
        filename=LOG_CACHE,
        level=logging.DEBUG
    )
    logging.captureWarnings(True)
    main(ARGS)
