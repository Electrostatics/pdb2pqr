"""Create an APBS input file using psize data

Authors: Todd Dolinsky based on original sed script by Nathan Baker
"""
import pickle
import os.path
import logging
import argparse
from . import psize
from . import utilities


_LOGGER = logging.getLogger(__name__)


class Elec(object):
    """An object for the ELEC section of an APBS input file"""

    def __init__(self, pqrpath, size, method, asyncflag, istrng=0, potdx=False):
        # If this is an async or parallel calc, we want to use
        # the per-grid dime rather than the global dime.
        self.dime = size.ngrid
        gmem = 200.0 * self.dime[0] * self.dime[1] * self.dime[2] / 1024.0 / 1024.0
        if method == "": # method not named - use ceiling
            if gmem > size.gmemceil:
                method = "mg-para"
            else:
                method = "mg-auto"

        if method == "mg-para":
            self.dime = size.getSmallest()

        self.method = method
        self.istrng = istrng
        self.glen = size.coarse_length
        self.cglen = size.coarse_length
        self.fglen = size.fine_length
        self.pdime = size.proc_grid

        self.label = ""
        self.nlev = 4
        self.ofrac = 0.1
        self.async_ = 0
        self.asyncflag = asyncflag
        self.cgcent = "mol 1"
        self.fgcent = "mol 1"
        self.gcent = "mol 1"
        self.mol = 1
        self.lpbe = 1
        self.npbe = 0
        self.bcfl = "sdh"
        # TODO - where did these very arbitrary numbers come from?
        self.ion = [[-1, 1.815], [1, 1.875]] # Multiple ions possible
        self.pdie = 2.0
        self.sdie = 78.54
        self.srfm = "smol"
        self.chgm = "spl2"
        self.sdens = 10.0
        self.srad = 1.4
        self.swin = 0.3
        self.temp = 298.15
        self.gamma = 0.105
        self.calcenergy = "total"
        self.calcforce = "no"
        if potdx:
            self.write = [["pot", "dx", pqrpath]]
        else:
            self.write = [["pot", "dx", "pot"]] # Multiple write statements possible

    def __str__(self):
        text = "elec %s\n" % self.label
        text += "    %s\n" % self.method
        text += "    dime %i %i %i\n" % (self.dime[0], self.dime[1], self.dime[2])
        if self.method == "mg-manual":
            text += "    glen %.3f %.3f %.3f\n" % (self.glen[0], self.glen[1], self.glen[2])
            text += "    gcent %s\n" % self.gcent
        elif self.method == "mg-auto":
            text += "    cglen %.4f %.4f %.4f\n" % (self.cglen[0], self.cglen[1], self.cglen[2])
            text += "    fglen %.4f %.4f %.4f\n" % (self.fglen[0], self.fglen[1], self.fglen[2])
            text += "    cgcent %s\n" % self.cgcent
            text += "    fgcent %s\n" % self.fgcent
        elif self.method == "mg-para":
            text += "    pdime %i %i %i\n" % (self.pdime[0], self.pdime[1], self.pdime[2])
            text += "    ofrac %.1f\n" % self.ofrac
            text += "    cglen %.4f %.4f %.4f\n" % (self.cglen[0], self.cglen[1], self.cglen[2])
            text += "    fglen %.4f %.4f %.4f\n" % (self.fglen[0], self.fglen[1], self.fglen[2])
            text += "    cgcent %s\n" % self.cgcent
            text += "    fgcent %s\n" % self.fgcent
            if self.asyncflag:
                text += "    async %i\n" % self.async_
        text += "    mol %i\n" % self.mol
        if self.lpbe:
            text += "    lpbe\n"
        else:
            text += "    npbe\n"
        text += "    bcfl %s\n" % self.bcfl
        if self.istrng > 0:
            for ion in self.ion:
                text += "    ion charge %.2f conc %.3f radius %.4f\n" % (ion[0],
                                                                         self.istrng,
                                                                         ion[1])
        text += "    pdie %.4f\n" % self.pdie
        text += "    sdie %.4f\n" % self.sdie
        text += "    srfm %s\n" % self.srfm
        text += "    chgm %s\n" % self.chgm
        text += "    sdens %.2f\n" % self.sdens
        text += "    srad %.2f\n" % self.srad
        text += "    swin %.2f\n" % self.swin
        text += "    temp %.2f\n" % self.temp
        text += "    calcenergy %s\n" % self.calcenergy
        text += "    calcforce %s\n" % self.calcforce
        for write in self.write:
            text += "    write %s %s %s\n" % (write[0], write[1], write[2])
        text += "end\n"
        return text

class Input(object):
    """The input class.  Each input object is one APBS input file."""

    def __init__(self, pqrpath, size, method, asyncflag, istrng=0, potdx=False):
        """Initialize the input file class.

        Each input file contains a PQR name, a list of elec objects, and a
        list of strings containing print statements.  For starters assume two
        ELEC statements are needed, one for the inhomgenous and the other for
        the homogenous dielectric calculations.

        Users can edit the elec statements and the print statements.

        This assumes you have already run psize, either by
            size.run_pize(/path/to/pqr) or

            size.parse_string(string)
            size.set_all()

        Args:
            pqrpath:   The path to the PQR file (string)
            size:      The Psize object (psize)
            method:    The method (para, auto, manual, async) to use
            asyncflag: 1 if async is desired, 0 otherwise
        """
        self.pqrpath = pqrpath
        self.asyncflag = asyncflag

        # Initialize variables to default elec values
        elec1 = Elec(pqrpath, size, method, asyncflag, istrng, potdx)
        if not potdx:
            elec2 = Elec(pqrpath, size, method, asyncflag, istrng, potdx)
            setattr(elec2, "sdie", 2.0)
            setattr(elec2, "write", [])
        else:
            elec2 = ""
        self.elecs = [elec1, elec2]

        self.pqrname = os.path.basename(pqrpath)

        if not potdx:
            self.prints = ["print elecEnergy 2 - 1 end"]
        else:
            self.prints = ["print elecEnergy 1 end"]

    def __str__(self):
        text = "read\n"
        text += "    mol pqr %s\n" % self.pqrname
        text += "end\n"
        for elec in self.elecs:
            text += str(elec)
        for prints in self.prints:
            text += prints
        text += "\nquit\n"
        return text

    def print_input_files(self):
        """Make the input file(s) associated with this object"""
        base_pqr_name = utilities.getPQRBaseFileName(self.pqrpath)
        if self.asyncflag == 1:
            outname = base_pqr_name + "-para.in"

            # Temporarily disable async flag
            for elec in self.elecs:
                elec.asyncflag = 0
            file_ = open(outname, "w")
            file_.write(str(self))
            file_.close()

            # Now make the async files
            elec = self.elecs[0]
            nproc = elec.pdime[0] * elec.pdime[1] * elec.pdime[2]
            for i in range(int(nproc)):
                outname = base_pqr_name + "-PE%i.in" % i
                for elec in self.elecs:
                    elec.asyncflag = 1
                    elec.async_ = i
                file_ = open(outname, "w")
                file_.write(str(self))
                file_.close()

        else:
            outname = base_pqr_name + ".in"
            file_ = open(outname, "w")
            file_.write(str(self))
            file_.close()

    def dump_pickle(self):
        """Make a Python pickle associated with the APBS input parameters"""
        # TODO - is this function still useful?
        base_pqr_name = utilities.getPQRBaseFileName(self.pqrpath)
        outname = base_pqr_name + "-input.p"
        pfile = open(outname, "wb")
        pickle.dump(self, pfile)
        pfile.close()


def split_input(filename):
    """Split the parallel input file into multiple async file names

    Args:
        filename:  The path to the original parallel input file (string)
    """
    nproc = 0
    with open(filename, "rt") as file_:
        text = ""
        while True:
            line = file_.readline()
            if not line:
                break
            text += line
            line = line.strip()
            if line.startswith("pdime"): # Get # Procs
                words = line.split()
                nproc = int(words[1]) * int(words[2]) * int(words[3])

    if nproc == 0:
        errstr = "%s is not a valid APBS parallel input file!\n" % filename
        errstr = errstr + "The inputgen script was unable to asynchronize this file!"
        raise RuntimeError(errstr)

    base_pqr_name = utilities.getPQRBaseFileName(filename)
    for iproc in range(nproc):
        outname = base_pqr_name + "-PE%i.in" % iproc
        outtext = text.replace("mg-para\n", "mg-para\n    async %i\n" % iproc)
        outfile = open(outname, "w")
        outfile.write(outtext)
        outfile.close()


def build_parser():
    """Build argument parser."""
    parse = argparse.ArgumentParser(description=("Use this script to generate new APBS input "
                                                 "files or split an existing parallel input "
                                                 "file into multiple async files"),
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parse.add_argument("--asynch", action="store_true", 
                       help="Perform an asynchronous parallel calculation.")
    parse.add_argument("--split", action="store_true",
                       help=("Split an existing parallel input file to multiple "
                             "async input files."))
    parse.add_argument("--potdx", action="store_true",
                       help=("Create an input to compute an electrostatic potential map."))
    parse.add_argument("--method",
                       help=("Force output file to write a specific APBS ELEC method."),
                       choices=["para", "auto", "manual", "async"])
    parse.add_argument("--cfac", type=float, default=psize.CFAC,
                       help=("Factor by which to expand molecular dimensions to "
                             "get coarse grid dimensions."))
    parse.add_argument("--fadd", type=float, default=psize.FADD,
                       help=("Amount to add to molecular dimensions to get fine "
                             "grid dimensions."))
    parse.add_argument("--space", type=float, default=psize.SPACE,
                       help="Desired fine mesh resolution")
    parse.add_argument("--gmemfac", type=int, default=psize.GMEMFAC,
                       help=("Number of bytes per grid point required for sequential "
                             "MG calculation"))
    parse.add_argument("--gmemceil", type=int, default=psize.GMEMCEIL,
                       help=("Max MB allowed for sequential MG calculation. Adjust "
                             "this to force the script to perform faster calculations "
                             "(which require more parallelism)"))
    parse.add_argument("--ofrac", type=float, default=psize.OFRAC,
                       help="Overlap factor between mesh partitions (parallel)")
    parse.add_argument("--redfac", type=float, default=psize.REDFAC,
                       help=("The maximum factor by which a domain dimension can "
                             "be reduced during focusing"))
    parse.add_argument("--istrng", help="Ionic strength (M). Na+ anc Cl- ions will be used")
    parse.add_argument("filename")
    return parse


def main():
    """Main driver"""
    parser = build_parser()
    args = parser.parse_args()
    size = psize.Psize()
    filename = args.filename

    if args.split:
        split_input(filename)
    else:
        size.run_pize(filename)
        input_ = Input(filename, size, args.method, args.asynch, args.istrng,
                       args.potdx)
        input_.print_input_files()


if __name__ == "__main__":
    main()
