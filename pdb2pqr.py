"""PDB2PQR

This package takes a PDB file as input and performs optimizations before
yielding a new PDB-style file as output.

For more information, see http://www.poissonboltzmann.org/
"""
from pdb2pqr import main
from pdb2pqr import cli


if __name__ == "__main__":
    parser = cli.build_parser()
    args = parser.parse_args()
    main(args)
