#!/usr/bin/python
"""Get dimensions and other information from a PQR file.

.. todo:: This code could be combined with :mod:`inputgen`.

.. todo:: This code should be moved to the APBS code base.

.. codeauthor:: Dave Sept
.. codeauthor:: Nathan Baker
.. codeauthor:: Todd Dolinksy
.. codeauthor:: Yong Huang
"""

import argparse
import logging
from math import log

from .config import TITLE_STR

#: The number of Angstroms added to the molecular dimensions to determine the
#: find grid dimensions
FADD = 20.0
#: The fine grid dimensions are multiplied by this constant to calculate the
#: coarse grid dimensions
CFAC = 1.7
#: Desired fine grid spacing (in Angstroms)
SPACE = 0.50
#: Approximate memory usage (in bytes) can be estimated by multiplying the
#: number of grid points by this constant
GMEMFAC = 200
#: Maxmimum memory (in MB) to be used for a calculation
GMEMCEIL = 400
#: The fractional overlap between grid partitions in a parallel focusing
#: calculation
OFRAC = 0.1
#: The maximum factor by which a domain can be "shrunk" during a focusing
#: calculation
REDFAC = 0.25
_LOGGER = logging.getLogger(__name__)


class Psize:
    """Class for parsing input files and suggesting settings."""

    def __init__(
        self,
        cfac=CFAC,
        fadd=FADD,
        space=SPACE,
        gmemfac=GMEMFAC,
        gmemceil=GMEMCEIL,
        ofrac=OFRAC,
        redfac=REDFAC,
    ):
        """Initialize.

        :param cfac:  factor by which to expand molecular dimensions to get
            coarse grid dimensions
        :type cfac:  float
        :param fadd:  amount (in Angstroms) to add to molecular dimensions to
            get the fine grid dimensions
        :type fadd:  float
        :param space:  desired fine mesh resolution (in Angstroms)
        :type space:  float
        :param gmemfac:  number of bytes per grid point required for
            a sequential multigrid calculation
        :type gmemfac:  float
        :param gmemceil:  maximum memory (in MB) allowed for sequential
            multigrid calculation. Adjust this value to force the script to
            perform faster calculations (which require more parallelism).
        :type gmemceil:  float
        :param ofrac:  overlap factor between mesh partitions in parallel
            focusing calculation
        :type ofrac:  float
        :param redfac:  the maximum factor by which a domain dimension can be
            reduced during focusing
        :type redfac:  float
        """
        self.minlen = [None, None, None]
        self.maxlen = [None, None, None]
        self.cfac = cfac
        self.fadd = fadd
        self.space = space
        self.gmemfac = gmemfac
        self.gmemceil = gmemceil
        self.ofrac = ofrac
        self.redfac = redfac
        self.charge = 0.0
        self.gotatom = 0
        self.gothet = 0
        self.mol_length = [0.0, 0.0, 0.0]
        self.center = [0.0, 0.0, 0.0]
        self.coarse_length = [0.0, 0.0, 0.0]
        self.fine_length = [0.0, 0.0, 0.0]
        self.ngrid = [0, 0, 0]
        self.proc_grid = [0.0, 0.0, 0.0]
        self.nsmall = [0, 0, 0]
        self.nfocus = 0

    def parse_string(self, structure):
        """Parse the input structure as a string in PDB or PQR format.

        :param structure:  input structure as string in PDB or PQR format.
        :type structure:  str
        """
        lines = structure.split("\n")
        self.parse_lines(lines)

    def parse_input(self, filename):
        """Parse input structure file in PDB or PQR format.

        :param filename:  string with path to PDB- or PQR-format file.
        :type filename:  str
        """
        with open(filename, encoding="utf-8") as file_:
            self.parse_lines(file_.readlines())

    def parse_lines(self, lines):
        """Parse the PQR/PDB lines.

        .. todo::
           This is messed up. Why are we parsing the PQR manually here when
           we already have other routines to do that?  This function should
           be replaced by a call to existing routines.

        :param lines:  PDB/PQR lines to parse
        :type lines:  [str]
        """
        for line in lines:
            if line.find("ATOM") == 0:
                self.gotatom += 1
            elif line.find("HETATM") == 0:
                self.gothet = self.gothet + 1
            subline = line[30:].replace("-", " -")
            words = subline.split()
            if len(words) < 5:
                continue
            self.charge = self.charge + float(words[3])
            rad = float(words[4])
            center = [float(word) for word in words[0:3]]
            for i in range(3):
                if self.minlen[i] is None or center[i] - rad < self.minlen[i]:
                    self.minlen[i] = center[i] - rad
                if self.maxlen[i] is None or center[i] + rad > self.maxlen[i]:
                    self.maxlen[i] = center[i] + rad

    def set_length(self, maxlen, minlen):
        """Compute molecular dimensions, adjusting for zero-length values.

        .. todo:: Replace hard-coded values in this function.

        :param maxlen:  maximum dimensions from molecule
        :type maxlen:  [float, float, float]
        :param minlen:  minimum dimensions from molecule
        :type minlen:  [float, float, float]
        :return:  molecular dimensions
        :rtype:  [float, float, float]
        """
        for i in range(3):
            self.mol_length[i] = maxlen[i] - minlen[i]
            self.mol_length[i] = max(self.mol_length[i], 0.1)
        return self.mol_length

    def set_coarse_grid_dims(self, mol_length):
        """Compute coarse mesh lengths.

        :param mol_length:  input molecule lengths
        :type mol_length:  [float, float, float]
        :return:  coarse grid dimensions
        :rtype:  [float, float, float]
        """
        for i in range(3):
            self.coarse_length[i] = self.cfac * mol_length[i]
        return self.coarse_length

    def set_fine_grid_dims(self, mol_length, coarse_length):
        """Compute fine mesh lengths.

        :param mol_length:  input molecule lengths
        :type mol_length:  [float, float, float]
        :param coarse_length:  coarse grid lengths
        :type coarse_length:  [float, float, float]
        :return:  fine grid lengths
        :rtype:  [float, float, float]
        """
        for i in range(3):
            self.fine_length[i] = mol_length[i] + self.fadd
            self.fine_length[i] = min(self.fine_length[i], coarse_length[i])
        return self.fine_length

    def set_center(self, maxlen, minlen):
        """Compute molecular center.

        :param maxlen:  maximum molecule lengths
        :type maxlen:  [float, float, float]
        :param minlen:  minimum molecule lengths
        :type minlen:  [float, float, float]
        :return:  center of molecule
        :rtype:  [float, float, float]
        """
        for i in range(3):
            self.center[i] = (maxlen[i] + minlen[i]) / 2
        return self.center

    def set_fine_grid_points(self, fine_length):
        """Compute mesh grid points, assuming 4 levels in multigrid hierarchy.

        .. todo:: remove hard-coded values from this function.

        :param fine_length:  lengths of the fine grid
        :type fine_length:  [float, float, float]
        :return:  number of grid points in each direction
        :rtype:  [int, int, int]
        """
        temp_num = [0, 0, 0]
        for i in range(3):
            temp_num[i] = int(fine_length[i] / self.space + 0.5)
            self.ngrid[i] = 32 * (int((temp_num[i] - 1) / 32.0 + 0.5)) + 1
            self.ngrid[i] = max(self.ngrid[i], 33)
        return self.ngrid

    def set_smallest(self, ngrid):
        """Set smallest dimensions.

        Compute parallel division of domain in case the memory requirements
        for the calculation are above the memory ceiling. Find the smallest
        dimension and see if the number of grid points in that dimension will
        fit below the memory ceiling Reduce nsmall until an nsmall^3 domain
        will fit into memory.

        .. todo:: Remove hard-coded values from this function.

        :param ngrid:  number of grid points
        :type ngrid:  [int, int, int]
        :return:  smallest number of grid points in each direction to fit in
            memory
        :rtype:  [int, int, int]
        """
        nsmall = [ngrid[i] for i in range(3)]
        while 1:
            nsmem = 200.0 * nsmall[0] * nsmall[1] * nsmall[2] / 1024 / 1024
            if nsmem < self.gmemceil:
                break
            else:
                i = nsmall.index(max(nsmall))
                nsmall[i] = 32 * ((nsmall[i] - 1) / 32 - 1) + 1
                if nsmall[i] <= 0:
                    _LOGGER.error(
                        "You picked a memory ceiling that is too small"
                    )
                    raise ValueError(nsmall[i])
        self.nsmall = nsmall
        return nsmall

    def set_proc_grid(self, ngrid, nsmall):
        """Calculate the number of processors required in a parallel focusing
        calculation to span each dimension of the grid given the grid size
        suitable for memory constraints.

        :param ngrid:  number of needed grid points
        :type ngrid:  [int, int, int]
        :param nsmall:  number of grid points that will fit in memory
        :type nsmall:  [int, int, int]
        :return:  number of processors needed in each direction
        :rtype:  [int, int, int]
        """
        zofac = 1 + 2 * self.ofrac
        for i in range(3):
            self.proc_grid[i] = ngrid[i] / float(nsmall[i])
            if self.proc_grid[i] > 1:
                self.proc_grid[i] = int(zofac * ngrid[i] / nsmall[i] + 1.0)
        return self.proc_grid

    def set_focus(self, fine_length, nproc, coarse_length):
        """Calculate the number of levels of focusing required for each processor subdomain.

        :param fine_length:  fine grid length
        :type fine_length:  [float, float, float]
        :param nproc:  number of processors in each dimension
        :type nproc:  [int, int, int]
        :param coarse_length:  coarse grid length
        :type coarse_length:  [float, float, float]
        """
        nfoc = [0, 0, 0]
        for i in range(3):
            nfoc[i] = int(
                log((fine_length[i] / nproc[i]) / coarse_length[i])
                / log(self.redfac)
                + 1.0
            )
        nfocus = nfoc[0]
        nfocus = max(nfoc[1], nfocus)
        nfocus = max(nfoc[2], nfocus)
        if nfocus > 0:
            nfocus = nfocus + 1
        self.nfocus = nfocus

    def set_all(self):
        """Set up all of the things calculated individually above."""
        maxlen = self.maxlen
        minlen = self.minlen
        self.set_length(maxlen, minlen)
        mol_length = self.mol_length
        self.set_coarse_grid_dims(mol_length)
        coarse_length = self.coarse_length
        self.set_fine_grid_dims(mol_length, coarse_length)
        fine_length = self.fine_length
        self.set_center(maxlen, minlen)
        self.set_fine_grid_points(fine_length)
        ngrid = self.ngrid
        self.set_smallest(ngrid)
        nsmall = self.nsmall
        self.set_proc_grid(ngrid, nsmall)
        nproc = self.proc_grid
        self.set_focus(fine_length, nproc, coarse_length)

    def run_psize(self, filename):
        """Parse input PQR file and set parameters.

        :param filename:  path of PQR file
        :type filename:  str
        """
        self.parse_input(filename)
        self.set_all()

    def __str__(self):
        """Return a string with the formatted results.

        :return:  string with formatted results
        :rtype:  str
        """
        str_ = "\n"
        if self.gotatom > 0:
            maxlen = self.maxlen
            minlen = self.minlen
            charge = self.charge
            mol_length = self.mol_length
            coarse_length = self.coarse_length
            fine_length = self.fine_length
            center = self.center
            ngrid = self.ngrid
            nsmall = self.nsmall
            nproc = self.proc_grid
            nfocus = self.nfocus
            # Compute memory requirements
            nsmem = 200.0 * nsmall[0] * nsmall[1] * nsmall[2] / 1024 / 1024
            gmem = 200.0 * ngrid[0] * ngrid[1] * ngrid[2] / 1024 / 1024
            # Print the calculated entries
            str_ += "######## MOLECULE INFO ########\n"
            str_ += f"Number of ATOM entries = {self.gotatom}\n"
            str_ += f"Number of HETATM entries (ignored) = {self.gothet}\n"
            str_ += f"Total charge = {charge:.3f} e\n"
            str_ += f"Dimensions = {mol_length[0]:.3f} Å x "
            str_ += f"{mol_length[1]:.3f} Å x {mol_length[2]:.3f} Å\n"
            str_ += f"Center = {center[0]:.3f} Å x {center[1]:.3f} Å x "
            str_ += f"{center[2]:.3f} Å\n"
            str_ += f"Lower corner = {minlen[0]:.3f} Å x {minlen[1]:.3f} Å x "
            str_ += f"{minlen[2]:.3f} Å\n"
            str_ += f"Upper corner = {maxlen[0]:.3f} Å x {maxlen[1]:.3f} Å x "
            str_ += f"{maxlen[2]:.3f} Å\n"
            str_ += "\n"
            str_ += "######## GENERAL CALCULATION INFO ########\n"
            str_ += f"Course grid dims = {coarse_length[0]:.3f} Å x "
            str_ += f"{coarse_length[1]:.3f} Å x "
            str_ += f"{coarse_length[2]:.3f} Å\n"
            str_ += f"Fine grid dims = {fine_length[0]:.3f} Å x "
            str_ += f"{fine_length[1]:.3f} Å x "
            str_ += f"{fine_length[2]:.3f} Å\n"
            str_ += f"Num. fine grid pts. = {ngrid[0]:d} Å x "
            str_ += f"{ngrid[1]:d} Å x "
            str_ += f"{ngrid[2]:d} Å\n"
            str_ += "\n"
            if gmem > self.gmemceil:
                str_ += f"Parallel solve required ({gmem:.3f} MB > "
                str_ += f"{self.gmemceil:.3f} MB)\n"
                str_ += "Total processors required = "
                str_ += f"{nproc[0] * nproc[1] * nproc[2]}\n"
                str_ += f"Proc. grid = {nproc[0]:d} x {nproc[1]:d} x "
                str_ += f"{nproc[2]:d}\n"
                str_ += f"Grid pts. on each proc. = {nsmall[0]:d} x "
                str_ += f"{nsmall[1]:d} x {nsmall[2]:d}\n"
                xglob = nproc[0] * round(
                    nsmall[0] / (1 + 2 * self.ofrac - 0.001)
                )
                yglob = nproc[1] * round(
                    nsmall[1] / (1 + 2 * self.ofrac - 0.001)
                )
                zglob = nproc[2] * round(
                    nsmall[2] / (1 + 2 * self.ofrac - 0.001)
                )
                if nproc[0] == 1:
                    xglob = nsmall[0]
                if nproc[1] == 1:
                    yglob = nsmall[1]
                if nproc[2] == 1:
                    zglob = nsmall[2]
                str_ += "Fine mesh spacing = "
                str_ += f"{fine_length[0] / (xglob - 1):g} x "
                str_ += f"{fine_length[1] / (yglob - 1):g} x "
                str_ += f"{fine_length[2] / (zglob - 1):g} A\n"
                str_ += "Estimated mem. required for parallel solve = "
                str_ += f"{nsmem:.3f} MB/proc.\n"
                ntot = nsmall[0] * nsmall[1] * nsmall[2]
            else:
                str_ += "Fine mesh spacing = "
                str_ += f"{fine_length[0] / (ngrid[0] - 1):g} x "
                str_ += f"{fine_length[1] / (ngrid[1] - 1):g} x "
                str_ += f"{fine_length[2] / (ngrid[2] - 1):g} A\n"
                str_ += "Estimated mem. required for sequential solve = "
                str_ += f"{gmem:.3f} MB\n"
                ntot = ngrid[0] * ngrid[1] * ngrid[2]
            str_ += f"Number of focusing operations = {nfocus}\n"
            str_ += "\n"
            str_ += "######## ESTIMATED REQUIREMENTS ########\n"
            str_ += "Memory per processor = "
            str_ += f"{200.0 * ntot / 1024 / 1024:.3f} MB\n"
            str_ += "Grid storage requirements (ASCII) = "
            grid_storage_req = (
                8.0 * 12 * nproc[0] * nproc[1] * nproc[2] * ntot / 1024 / 1024
            )
            str_ += f"{grid_storage_req:.3f} MB\n"
            str_ += "\n"
        else:
            str_ = "No ATOM entries in file!\n\n"
        return str_


def build_parser():
    """Build argument parser.

    :return:  argument parser
    :rtype:  argparse.ArgumentParser
    """
    desc = f"{TITLE_STR}\npsize: figuring out the size of electrostatics "
    desc += "calculations since (at least) 2002."

    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--cfac",
        default=CFAC,
        type=float,
        help=(
            "Factor by which to expand molecular dimensions to "
            "get coarse grid dimensions"
        ),
    )
    parser.add_argument(
        "--fadd",
        default=FADD,
        type=float,
        help="Amount to add to mol dims to get fine grid dims",
    )
    parser.add_argument(
        "--space",
        default=SPACE,
        type=float,
        help="Desired fine mesh resolution",
    )
    parser.add_argument(
        "--gmemfac",
        default=GMEMFAC,
        type=int,
        help=(
            "Number of bytes per grid point required for "
            "sequential MG calculation"
        ),
    )
    parser.add_argument(
        "--gmemceil",
        default=GMEMCEIL,
        type=int,
        help=(
            "Max MB allowed for sequential MG calculation. "
            "Adjust this to force the script to perform faster "
            "calculations (which require more parallelism)."
        ),
    )
    parser.add_argument(
        "--ofrac",
        default=OFRAC,
        type=float,
        help="Overlap factor between mesh partitions",
    )
    parser.add_argument(
        "--redfac",
        default=REDFAC,
        type=float,
        help=(
            "The maximum factor by which a domain dimension "
            "can be reduced during focusing"
        ),
    )
    parser.add_argument("mol_path", help="Path to PQR file.")
    return parser
