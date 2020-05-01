#!/usr/bin/python
"""psize

Get dimensions and other information from a PQR file.

Authors:  Dave Sept, Nathan Baker, Todd Dolinksy, Yong Huang
"""
import string
from math import log
import logging
import argparse


CFAC = 1.7
FADD = 20.0
SPACE = 0.50
GMEMFAC = 200
GMEMCEIL = 400
OFRAC = 0.1
REDFRAC = 0.25
_LOGGER = logging.getLogger(__name__)


class Psize(object):
    """Master class for parsing input files and suggesting settings"""
    def __init__(self):
        self.constants = {"cfac": CFAC, "fadd": FADD, "space": SPACE,
                          "gmemfac": GMEMFAC, "gmemceil": GMEMCEIL,
                          "ofrac": OFRAC, "redfac": REDFRAC}
        self.minlen = [None, None, None]
        self.maxlen = [None, None, None]
        self.charge = 0.0
        self.gotatom = 0
        self.gothet = 0
        self.olen = [0.0, 0.0, 0.0]
        self.cen = [0.0, 0.0, 0.0]
        self.clen = [0.0, 0.0, 0.0]
        self.flen = [0.0, 0.0, 0.0]
        self.ngrid = [0, 0, 0]
        self.proc_grid = [0.0, 0.0, 0.0]
        self.nsmall = [0, 0, 0]
        self.nfocus = 0

    def parse_string(self, structure):
        """ Parse the input structure as a string in PDB or PQR format """
        # TODO - modernize str functions in this module
        lines = str.split(structure, "\n")
        self.parseLines(lines)

    def parse_input(self, filename):
        """ Parse input structure file in PDB or PQR format """
        with open(filename, 'rt', encoding="utf-8") as file_:
            self.parseLines(file_.readlines())

    def parseLines(self, lines):
        """ Parse the lines """
        for line in lines:
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
                    if self.minlen[i] == None or center[i]-rad < self.minlen[i]:
                        self.minlen[i] = center[i]-rad
                    if self.maxlen[i] == None or center[i]+rad > self.maxlen[i]:
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
                        if self.minlen[i] == None or center[i]-rad < self.minlen[i]:
                            self.minlen[i] = center[i]-rad
                        if self.maxlen[i] == None or center[i]+rad > self.maxlen[i]:
                            self.maxlen[i] = center[i]+rad

    def setLength(self, maxlen, minlen):
        """ Compute molecule dimensions """
        for i in range(3):
            self.olen[i] = maxlen[i] - minlen[i]
            if self.olen[i] < 0.1:
                self.olen[i] = 0.1
        return self.olen

    def setCoarseGridDims(self, olen):
        """ Compute coarse mesh dimensions """
        for i in range(3):
            self.clen[i] = self.constants["cfac"] * olen[i]
        return self.clen

    def setFineGridDims(self, olen, clen):
        """ Compute fine mesh dimensions """
        for i in range(3):
            self.flen[i] = olen[i] + self.constants["fadd"]
            if self.flen[i] > clen[i]:
                self.flen[i] = clen[i]
        return self.flen

    def setCenter(self, maxlen, minlen):
        """ Compute molecule center """
        for i in range(3):
            self.cen[i] = (maxlen[i] + minlen[i]) / 2
        return self.cen

    def setFineGridPoints(self, flen):
        """ Compute mesh grid points, assuming 4 levels in MG hierarchy """
        tn = [0, 0, 0]
        for i in range(3):
            tn[i] = int(flen[i]/self.constants["space"] + 0.5)
            self.ngrid[i] = 32*(int((tn[i] - 1) / 32.0 + 0.5)) + 1
            if self.ngrid[i] < 33:
                self.ngrid[i] = 33
        return self.ngrid

    def setSmallest(self, n):
        """ Compute parallel division in case memory requirement above ceiling
        Find the smallest dimension and see if the number of grid points in
        that dimension will fit below the memory ceiling
        Reduce nsmall until an nsmall^3 domain will fit into memory """
        nsmall = []
        for i in range(3):
            nsmall.append(n[i])
        while 1:
            nsmem = 200.0 * nsmall[0] * nsmall[1] * nsmall[2] / 1024 / 1024
            if nsmem < self.constants["gmemceil"]: break
            else:
                i = nsmall.index(max(nsmall))
                nsmall[i] = 32 * ((nsmall[i] - 1)/32 - 1) + 1
                if nsmall[i] <= 0:
                    _LOGGING.error("You picked a memory ceiling that is too small")
                    raise ValueError(nsmall[i])

        self.nsmall = nsmall
        return nsmall

    def setProcGrid(self, n, nsmall):
        """ Calculate the number of processors required to span each
        dimension """

        zofac = 1 + 2 * self.constants["ofrac"]
        for i in range(3):
            self.proc_grid[i] = n[i]/float(nsmall[i])
            if self.proc_grid[i] > 1: self.proc_grid[i] = int(zofac*n[i]/nsmall[i] + 1.0)
        return self.proc_grid

    def setFocus(self, flen, np, clen):
        """ Calculate the number of levels of focusing required for each
        processor subdomain """
        nfoc = [0, 0, 0]
        for i in range(3):
            nfoc[i] = int(log((flen[i]/np[i])/clen[i])/log(self.constants["redfac"]) + 1.0)
        nfocus = nfoc[0]
        if nfoc[1] > nfocus: nfocus = nfoc[1]
        if nfoc[2] > nfocus: nfocus = nfoc[2]
        if nfocus > 0: nfocus = nfocus + 1
        self.nfocus = nfocus

    def setAll(self):
        """ Set up all of the things calculated individually above """
        maxlen = self.maxlen
        minlen = self.getMin()
        self.setLength(maxlen, minlen)
        olen = self.getLength()

        self.setCoarseGridDims(olen)
        clen = self.getCoarseGridDims()

        self.setFineGridDims(olen, clen)
        flen = self.getFineGridDims()

        self.setCenter(maxlen, minlen)
        cen = self.getCenter()

        self.setFineGridPoints(flen)
        n = self.getFineGridPoints()

        self.setSmallest(n)
        nsmall = self.getSmallest()

        self.setProcGrid(n, nsmall)
        np = self.proc_grid

        self.setFocus(flen, np, clen)
        nfocus = self.getFocus()

    def getMax(self): return self.maxlen
    def getMin(self): return self.minlen
    def getLength(self): return self.olen
    def getCoarseGridDims(self): return self.clen
    def getFineGridDims(self): return self.flen
    def getCenter(self): return self.cen
    def getFineGridPoints(self): return self.ngrid
    def getSmallest(self): return self.nsmall
    def getFocus(self): return self.nfocus

    def run_pize(self, filename):
        """ Parse input PQR file and set parameters """
        self.parse_input(filename)
        self.setAll()

    def printResults(self):
        """ Return a string with the formatted results """

        str_ = "\n"

        if self.gotatom > 0:

            maxlen = self.getMax()
            minlen = self.getMin()
            q = self.charge
            olen = self.getLength()
            clen = self.getCoarseGridDims()
            flen = self.getFineGridDims()
            cen = self.getCenter()
            n = self.getFineGridPoints()
            nsmall = self.getSmallest()
            np = self.proc_grid
            nfocus = self.getFocus()

            # Compute memory requirements

            nsmem = 200.0 * nsmall[0] * nsmall[1] * nsmall[2] / 1024 / 1024
            gmem = 200.0 * n[0] * n[1] * n[2] / 1024 / 1024

            # Print the calculated entries
            str_ = str_ + "################# MOLECULE INFO ####################\n"
            str_ = str_ + "Number of ATOM entries = %i\n" % self.gotatom
            str_ = str_ + "Number of HETATM entries (ignored) = %i\n" % self.gothet
            str_ = str_ + "Total charge = %.3f e\n" % q
            str_ = str_ + "Dimensions = %.3f x %.3f x %.3f A\n" % (olen[0], olen[1],
                                                                   olen[2])
            str_ = str_ + "Center = %.3f x %.3f x %.3f A\n" % (cen[0], cen[1], cen[2])
            str_ = str_ + "Lower corner = %.3f x %.3f x %.3f A\n" % (minlen[0], minlen[1],
                                                                     minlen[2])
            str_ = str_ + "Upper corner = %.3f x %.3f x %.3f A\n" % (maxlen[0], maxlen[1],
                                                                     maxlen[2])

            str_ = str_ + "\n"
            str_ = str_ + "############## GENERAL CALCULATION INFO #############\n"
            str_ = str_ + "Coarse grid dims = %.3f x %.3f x %.3f A\n" % (clen[0], clen[1], clen[2])
            str_ = str_ + "Fine grid dims = %.3f x %.3f x %.3f A\n" % (flen[0], flen[1], flen[2])
            str_ = str_ + "Num. fine grid pts. = %i x %i x %i\n" % (n[0], n[1], n[2])

            if gmem > self.constants["gmemceil"]:
                str_ = str_ + ("Parallel solve required "
                               "(%.3f MB > %.3f MB)\n") % (gmem, self.constants["gmemceil"])
                str_ = str_ + "Total processors required = %i\n" % (np[0]*np[1]*np[2])
                str_ = str_ + "Proc. grid = %i x %i x %i\n" % (np[0], np[1], np[2])
                str_ = str_ + ("Grid pts. on each proc. = "
                               "%i x %i x %i\n") % (nsmall[0], nsmall[1], nsmall[2])
                xglob = np[0]*round(nsmall[0]/(1 + 2*self.constants["ofrac"]) - .001)
                yglob = np[1]*round(nsmall[1]/(1 + 2*self.constants["ofrac"]) - .001)
                zglob = np[2]*round(nsmall[2]/(1 + 2*self.constants["ofrac"]) - .001)
                if np[0] == 1:
                    xglob = nsmall[0]
                if np[1] == 1:
                    yglob = nsmall[1]
                if np[2] == 1:
                    zglob = nsmall[2]
                str_ = str_ + "Fine mesh spacing = %g x %g x %g A\n" % (flen[0]/(xglob-1),
                                                                        flen[1]/(yglob-1),
                                                                        flen[2]/(zglob-1))
                str_ = str_ + "Estimated mem. required for parallel solve = %.3f MB/proc.\n" % nsmem
                ntot = nsmall[0]*nsmall[1]*nsmall[2]

            else:
                str_ = str_ + "Fine mesh spacing = %g x %g x %g A\n" % (flen[0]/(n[0]-1),
                                                                        flen[1]/(n[1]-1),
                                                                        flen[2]/(n[2]-1))
                str_ = str_ + "Estimated mem. required for sequential solve = %.3f MB\n" % gmem
                ntot = n[0]*n[1]*n[2]

            str_ = str_ + "Number of focusing operations = %i\n" % nfocus
            str_ = str_ + "\n"
            str_ = str_ + "################# ESTIMATED REQUIREMENTS ####################\n"
            str_ = str_ + "Memory per processor = %.3f MB\n" % (200.0*ntot/1024/1024)
            str_ = str_ + ("Grid storage requirements (ASCII) "
                           "= %.3f MB\n") % (8.0*12*np[0]*np[1]*np[2]*ntot/1024/1024)
            str_ = str_ + "\n"
        else:
            str_ = str_ + "No ATOM entires in file!\n\n"
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
    parser.add_argument("--redfac", default=REDFRAC, type=float,
                        help=("The maximum factor by which a domain dimension "
                              "can be reduced during focusing"))
    return parser


def main():
    """Main driver for module."""
    parser = build_parser()
    args = parser.parse_args()

    psize = Psize()
    psize.setConstant()

    psize.cfac = args.cfac
    psize.fadd = args.fadd
    psize.space = args.space
    psize.gmemfac = args.gmemfac
    psize.gmemceil = args.gmemceil
    psize.ofrac = args.ofrac
    psize.redfac = args.redfac

    psize.run_pize(filename)

    _LOGGING.info("# Constants used: \n")
    for key in psize.constants.keys():
        _LOGGING.info("# \t%s: %s\n" % (key, psize.constants[key]))
    _LOGGING.info("# Run:\n")
    _LOGGING.info("#    `%s --help`\n" % argv[0])
    _LOGGING.info("# for more information on these default values\n")
    _LOGGING.info(psize.printResults())


if __name__ == "__main__":
    main()
