"""PDB2PQR

This package takes a PDB file as input and performs optimizations before
yielding a new PDB-style file as output.

For more information, see http://www.poissonboltzmann.org/
"""
from main import mainCommand


if __name__ == "__main__":
    """ Determine if called from command line or CGI """
    mainCommand()
