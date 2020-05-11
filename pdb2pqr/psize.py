#!/usr/bin/python
"""psize

Get dimensions and other information from a PQR file.

Authors:  Dave Sept, Nathan Baker, Todd Dolinksy, Yong Huang
"""
from math import log
import logging
import argparse


CFAC = 1.7
FADD = 20.0
SPACE = 0.50
GMEMFAC = 200
GMEMCEIL = 400
OFRAC = 0.1
REDFAC = 0.25
_LOGGER = logging.getLogger(__name__)


class Psize(object):
    """Master class for parsing input files and suggesting settings"""
    def __init__(self, cfac=CFAC, fadd=FADD, space=SPACE, gmemfac=GMEMFAC,
                 gmemceil=GMEMCEIL, ofrac=OFRAC, redfac=REDFAC):
        """Initialize Psize.

        Args:
            cfac:  Factor by which to expand molecular dimensions to get coarse
                    grid dimensions
            fadd:  Amount to add to mol dims to get fine grid dims
            space:  Desired fine mesh resolution
            gmemfac:  Number of bytes per grid point required for sequential MG
                        calculation
            gmemceil: Max MB allowed for sequential MG calculation. Adjust this
                        to force the script to perform faster calculations (which
                        require more parallelism).
            ofrac:  Overlap factor between mesh partitions
            redfac:  The maximum factor by which a domain dimension can be reduced
                     during focusing
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
        """ Parse the input structure as a string in PDB or PQR format """
        # TODO - modernize str functions in this module
        lines = str.split(structure, "\n")
        self.parse_lines(lines)

    def parse_input(self, filename):
        """ Parse input structure file in PDB or PQR format """
        with open(filename, 'rt', encoding="utf-8") as file_:
            self.parse_lines(file_.readlines())

    def parse_lines(self, lines):
        """ Parse the lines """
        for line in lines:
            # TODO -- This is messed up.
            # Why are we parsing the PQR manually here when we have routines to do that?
            if str.find(line, "ATOM") == 0:
                subline = str.replace(line[30:], "-", " -")
                words = str.split(subline)
                if len(words) < 5:
                    continue
                self.gotatom += 1
                self.charge = self.charge + float(words[3])
                rad = float(words[4])
                center = []
                for word in words[0:3]:
                    center.append(float(word))
                for i in range(3):
                    if self.minlen[i] is None or center[i]-rad < self.minlen[i]:
                        self.minlen[i] = center[i]-rad
                    if self.maxlen[i] is None or center[i]+rad > self.maxlen[i]:
                        self.maxlen[i] = center[i]+rad
            elif str.find(line, "HETATM") == 0:
                self.gothet = self.gothet + 1
                # Special handling for no ATOM entries in the pqr file, only HETATM entries
                if self.gotatom == 0:
                    subline = str.replace(line[30:], "-", " -")
                    words = str.split(subline)
                    if len(words) < 5:
                        continue
                    self.charge = self.charge + float(words[3])
                    rad = float(words[4])
                    center = []
                    for word in words[0:3]:
                        center.append(float(word))
                    for i in range(3):
                        if self.minlen[i] is None or center[i]-rad < self.minlen[i]:
                            self.minlen[i] = center[i]-rad
                        if self.maxlen[i] is None or center[i]+rad > self.maxlen[i]:
                            self.maxlen[i] = center[i]+rad

    def set_length(self, maxlen, minlen):
        """ Compute molecule dimensions """
        for i in range(3):
            self.mol_length[i] = maxlen[i] - minlen[i]
            if self.mol_length[i] < 0.1:
                self.mol_length[i] = 0.1
        return self.mol_length

    def set_coarse_grid_dims(self, mol_length):
        """ Compute coarse mesh dimensions """
        for i in range(3):
            self.coarse_length[i] = self.cfac * mol_length[i]
        return self.coarse_length

    def set_fine_grid_dims(self, mol_length, coarse_length):
        """ Compute fine mesh dimensions """
        for i in range(3):
            self.fine_length[i] = mol_length[i] + self.fadd
            if self.fine_length[i] > coarse_length[i]:
                self.fine_length[i] = coarse_length[i]
        return self.fine_length

    def set_center(self, maxlen, minlen):
        """ Compute molecule center """
        for i in range(3):
            self.center[i] = (maxlen[i] + minlen[i]) / 2
        return self.center

    def set_fine_grid_points(self, fine_length):
        """ Compute mesh grid points, assuming 4 levels in MG hierarchy """
        temp_num = [0, 0, 0]
        for i in range(3):
            temp_num[i] = int(fine_length[i]/self.space + 0.5)
            self.ngrid[i] = 32*(int((temp_num[i] - 1) / 32.0 + 0.5)) + 1
            if self.ngrid[i] < 33:
                self.ngrid[i] = 33
        return self.ngrid

    def set_smallest(self, ngrid):
        """ Compute parallel division in case memory requirement above ceiling
        Find the smallest dimension and see if the number of grid points in
        that dimension will fit below the memory ceiling
        Reduce nsmall until an nsmall^3 domain will fit into memory """
        nsmall = []
        for i in range(3):
            nsmall.append(ngrid[i])
        while 1:
            nsmem = 200.0 * nsmall[0] * nsmall[1] * nsmall[2] / 1024 / 1024
            if nsmem < self.gmemceil:
                break
            else:
                i = nsmall.index(max(nsmall))
                nsmall[i] = 32 * ((nsmall[i] - 1)/32 - 1) + 1
                if nsmall[i] <= 0:
                    _LOGGER.error("You picked a memory ceiling that is too small")
                    raise ValueError(nsmall[i])
        self.nsmall = nsmall
        return nsmall

    def set_proc_grid(self, ngrid, nsmall):
        """ Calculate the number of processors required to span each
        dimension """
        zofac = 1 + 2 * self.ofrac
        for i in range(3):
            self.proc_grid[i] = ngrid[i]/float(nsmall[i])
            if self.proc_grid[i] > 1:
                self.proc_grid[i] = int(zofac*ngrid[i]/nsmall[i] + 1.0)
        return self.proc_grid

    def set_focus(self, fine_length, nproc, coarse_length):
        """ Calculate the number of levels of focusing required for each
        processor subdomain """
        nfoc = [0, 0, 0]
        for i in range(3):
            nfoc[i] = int(log((fine_length[i]/nproc[i])/coarse_length[i])/log(self.redfac) + 1.0)
        nfocus = nfoc[0]
        if nfoc[1] > nfocus:
            nfocus = nfoc[1]
        if nfoc[2] > nfocus:
            nfocus = nfoc[2]
        if nfocus > 0:
            nfocus = nfocus + 1
        self.nfocus = nfocus

    def set_all(self):
        """ Set up all of the things calculated individually above """
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
        """ Parse input PQR file and set parameters """
        self.parse_input(filename)
        self.set_all()

    def __str__(self):
        """ Return a string with the formatted results """
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
            str_ = str_ + "################# MOLECULE INFO ####################\n"
            str_ = str_ + "Number of ATOM entries = %i\n" % self.gotatom
            str_ = str_ + "Number of HETATM entries (ignored) = %i\n" % self.gothet
            str_ = str_ + "Total charge = %.3f e\n" % charge
            str_ = str_ + "Dimensions = %.3f x %.3f x %.3f A\n" % (mol_length[0],
                                                                   mol_length[1],
                                                                   mol_length[2])
            str_ = str_ + "Center = %.3f x %.3f x %.3f A\n" % (center[0], center[1], center[2])
            str_ = str_ + "Lower corner = %.3f x %.3f x %.3f A\n" % (minlen[0],
                                                                     minlen[1],
                                                                     minlen[2])
            str_ = str_ + "Upper corner = %.3f x %.3f x %.3f A\n" % (maxlen[0],
                                                                     maxlen[1],
                                                                     maxlen[2])

            str_ = str_ + "\n"
            str_ = str_ + "############## GENERAL CALCULATION INFO #############\n"
            str_ = str_ + "Coarse grid dims = %.3f x %.3f x %.3f A\n" % (coarse_length[0],
                                                                         coarse_length[1],
                                                                         coarse_length[2])
            str_ = str_ + "Fine grid dims = %.3f x %.3f x %.3f A\n" % (fine_length[0],
                                                                       fine_length[1],
                                                                       fine_length[2])
            str_ = str_ + "Num. fine grid pts. = %i x %i x %i\n" % (ngrid[0],
                                                                    ngrid[1],
                                                                    ngrid[2])

            if gmem > self.gmemceil:
                str_ = str_ + ("Parallel solve required "
                               "(%.3f MB > %.3f MB)\n") % (gmem, self.gmemceil)
                str_ = str_ + "Total processors required = %i\n" % (nproc[0]*nproc[1]*nproc[2])
                str_ = str_ + "Proc. grid = %i x %i x %i\n" % (nproc[0], nproc[1], nproc[2])
                str_ = str_ + ("Grid pts. on each proc. = "
                               "%i x %i x %i\n") % (nsmall[0], nsmall[1], nsmall[2])
                xglob = nproc[0]*round(nsmall[0]/(1 + 2*self.ofrac - 0.001))
                yglob = nproc[1]*round(nsmall[1]/(1 + 2*self.ofrac - 0.001))
                zglob = nproc[2]*round(nsmall[2]/(1 + 2*self.ofrac - 0.001))
                if nproc[0] == 1:
                    xglob = nsmall[0]
                if nproc[1] == 1:
                    yglob = nsmall[1]
                if nproc[2] == 1:
                    zglob = nsmall[2]
                str_ = str_ + "Fine mesh spacing = %g x %g x %g A\n" % (fine_length[0]/(xglob-1),
                                                                        fine_length[1]/(yglob-1),
                                                                        fine_length[2]/(zglob-1))
                str_ = str_ + "Estimated mem. required for parallel solve = %.3f MB/proc.\n" % nsmem
                ntot = nsmall[0]*nsmall[1]*nsmall[2]

            else:
                str_ = str_ + "Fine mesh spacing = %g x %g x %g A\n" % (fine_length[0]/(ngrid[0]-1),
                                                                        fine_length[1]/(ngrid[1]-1),
                                                                        fine_length[2]/(ngrid[2]-1))
                str_ = str_ + "Estimated mem. required for sequential solve = %.3f MB\n" % gmem
                ntot = ngrid[0]*ngrid[1]*ngrid[2]

            str_ = str_ + "Number of focusing operations = %i\n" % nfocus
            str_ = str_ + "\n"
            str_ = str_ + "################# ESTIMATED REQUIREMENTS ####################\n"
            str_ = str_ + "Memory per processor = %.3f MB\n" % (200.0*ntot/1024/1024)
            str_ = str_ + ("Grid storage requirements (ASCII) "
                           "= %.3f MB\n") % (8.0*12*nproc[0]*nproc[1]*nproc[2]*ntot/1024/1024)
            str_ = str_ + "\n"
        else:
            str_ = str_ + "No ATOM entries in file!\n\n"
        return str_


def build_parser():
    """Build argument parser.

    Returns:
        ArgumentParser
    """
    parser = argparse.ArgumentParser(description="Set size parameters for APBS",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--cfac", default=CFAC, type=float,
                        help=("Factor by which to expand molecular dimensions to "
                              "get coarse grid dimensions"))
    parser.add_argument("--fadd", default=FADD, type=float,
                        help="Amount to add to mol dims to get fine grid dims")
    parser.add_argument("--space", default=SPACE, type=float,
                        help="Desired fine mesh resolution")
    parser.add_argument("--gmemfac", default=GMEMFAC, type=int,
                        help=("Number of bytes per grid point required for "
                              "sequential MG calculation"))
    parser.add_argument("--gmemceil", default=GMEMCEIL, type=int,
                        help=("Max MB allowed for sequential MG calculation. "
                              "Adjust this to force the script to perform faster "
                              "calculations (which require more parallelism)."))
    parser.add_argument("--ofrac", default=OFRAC, type=float,
                        help="Overlap factor between mesh partitions")
    parser.add_argument("--redfac", default=REDFAC, type=float,
                        help=("The maximum factor by which a domain dimension "
                              "can be reduced during focusing"))
    parser.add_argument("mol_path", help="Path to PQR file.")
    return parser


def main():
    """Main driver for module."""
    parser = build_parser()
    args = parser.parse_args()

    psize = Psize(cfac=args.cfac, fadd=args.fadd, space=args.space,
                  gmemfac=args.gmemfac, gmemceil=args.gmemceil,
                  ofrac=args.ofrac, redfac=args.redfac)
    print(psize)


if __name__ == "__main__":
    main()
