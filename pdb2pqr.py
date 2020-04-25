"""Driver for PDB2PQR

This module takes a PDB file as input and performs optimizations before yielding
a new PDB-style file as output.
"""
from main import mainCommand


__version__ = "3.0"


if __name__ == "__main__":
    """ Determine if called from command line or CGI """
    mainCommand()
