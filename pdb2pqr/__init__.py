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
    parser = build_parser()
    args = parser.parse_args()

    # Get the output logging location
    output_path = Path(args.output_pqr)
    log_file = Path(output_path.parent, output_path.stem + '.log')
    logging.basicConfig(
        filename=log_file,
        level=logging.DEBUG)
    logging.captureWarnings(True)
    _LOGGER.info(args)

    main(args)
