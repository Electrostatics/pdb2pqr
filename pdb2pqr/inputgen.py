"""Create an APBS input file using :mod:`psize` data.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Nathan Baker
"""

import argparse
import logging
import pickle
from pathlib import Path

from . import psize
from .config import TITLE_STR

_LOGGER = logging.getLogger(__name__)


class Elec:
    """An object for the ELEC section of an APBS input file."""

    def __init__(
        self, pqrpath, size, method, asyncflag, istrng=0, *, potdx=False
    ):
        """Initialize object.

        .. todo::  Remove hard-coded parameters.

        :param pqrpath:  path to PQR file
        :type pqrpath:  str
        :param size:  parameter sizing object
        :type size:  Psize
        :param method:  solution method (e.g., mg-para, mg-auto, etc.)
        :type method:  str
        :param asyncflag:  perform an asynchronous parallel focusing
            calculation
        :type asyncflag:  bool
        :param istrng:  ionic strength/concentration (M)
        :type istring:  float
        :param potdx:  whether to write out potential information in DX format
        :type potdx:  bool
        """
        # If this is an async or parallel calc, we want to use
        # the per-grid dime rather than the global dime.
        self.dime = size.ngrid
        gmem = (
            200.0
            * self.dime[0]
            * self.dime[1]
            * self.dime[2]
            / 1024.0
            / 1024.0
        )
        if method == "":  # method not named - use ceiling
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
        self.ion = [[-1, 1.815], [1, 1.875]]  # Multiple ions possible
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
            # Multiple write statements possible
            self.write = [["pot", "dx", "pot"]]

    def __str__(self):
        text = f"elec {self.label}\n"
        text += f"    {self.method}\n"
        text += f"    dime {int(self.dime[0])} {int(self.dime[1])} "
        text += f"{int(self.dime[2])}\n"
        if self.method == "mg-auto":
            text += f"    cglen {self.cglen[0]:.4f} {self.cglen[1]:.4f} "
            text += f"{self.cglen[2]:.4f}\n"
            text += f"    fglen {self.fglen[0]:.4f} {self.fglen[1]:.4f} "
            text += f"{self.fglen[2]:.4f}\n"
            text += f"    cgcent {self.cgcent}\n"
            text += f"    fgcent {self.fgcent}\n"
        elif self.method == "mg-manual":
            text += f"    glen {self.glen[0]:.3f} {self.glen[1]:.3f} "
            text += f"{self.glen[2]:.3f}\n"
            text += f"    gcent {self.gcent}\n"
        elif self.method == "mg-para":
            text += f"    pdime {int(self.pdime[0])} {int(self.pdime[1])} "
            text += f"{int(self.pdime[2])}\n"
            text += f"    ofrac {self.ofrac:.1f}\n"
            text += f"    cglen {self.cglen[0]:.4f} {self.cglen[1]:.4f} "
            text += f"{self.cglen[2]:.4f}\n"
            text += f"    fglen {self.fglen[0]:.4f} {self.fglen[1]:.4f} "
            text += f"{self.fglen[2]:.4f}\n"
            text += f"    cgcent {self.cgcent}\n"
            text += f"    fgcent {self.fgcent}\n"
            if self.asyncflag:
                text += f"    async {self.async_}\n"
        text += f"    mol {int(self.mol)}\n"
        text += "    lpbe\n" if self.lpbe else "    npbe\n"
        text += f"    bcfl {self.bcfl}\n"
        if self.istrng > 0:
            for ion in self.ion:
                text += f"    ion charge {ion[0]:.2f} conc {self.istrng:.3f} "
                text += f"radius {ion[1]:.4f}\n"
        text += f"    pdie {self.pdie:.4f}\n"
        text += f"    sdie {self.sdie:.4f}\n"
        text += f"    srfm {self.srfm}\n"
        text += f"    chgm {self.chgm}\n"
        text += f"    sdens {self.sdens:.2f}\n"
        text += f"    srad {self.srad:.2f}\n"
        text += f"    swin {self.swin:.2f}\n"
        text += f"    temp {self.temp:.2f}\n"
        text += f"    calcenergy {self.calcenergy}\n"
        text += f"    calcforce {self.calcforce}\n"
        for write in self.write:
            text += f"    write {write[0]} {write[1]} {write[2]}\n"
        text += "end\n"
        return text


class Input:
    """Each object of this class is one APBS input file."""

    def __init__(
        self, pqrpath, size, method, asyncflag, istrng=0, *, potdx=False
    ):
        """Initialize the input file class.

        Each input file contains a PQR name, a list of elec objects, and a
        list of strings containing print statements.
        For starters, assume two ELEC statements are needed, one for the
        inhomgenous and the other for the homogenous dielectric calculations.

        .. note::
            This assumes you have already run psize, either by
            :func:`size.run_psize(...)` or :func:`size.parse_string(...)`
            followed by :func:`size.set_all()`.

        :param pqrpath:  path to PQR file
        :type pqrpath:  str
        :param size:  parameter sizing object
        :type size:  Psize
        :param method:  solution method (e.g., mg-para, mg-auto, etc.)
        :type method:  str
        :param asyncflag:  perform an asynchronous parallel focusing
            calculation
        :type asyncflag:  bool
        :param istrng:  ionic strength/concentration (M)
        :type istring:  float
        :param potdx:  whether to write out potential information in DX format
        :type potdx:  bool
        """
        self.pqrpath = Path(pqrpath)
        self.pqrname = self.pqrpath.name
        self.asyncflag = asyncflag
        # Initialize variables to default elec values
        elec1 = Elec(pqrpath, size, method, asyncflag, istrng, potdx=potdx)
        if not potdx:
            elec2 = Elec(pqrpath, size, method, asyncflag, istrng, potdx=potdx)
            elec2.sdie = 2.0
            elec2.write = []
        else:
            elec2 = ""
        self.elecs = [elec1, elec2]
        if not potdx:
            self.prints = ["print elecEnergy 2 - 1 end"]
        else:
            self.prints = ["print elecEnergy 1 end"]

    def __str__(self):
        text = "read\n"
        text += f"    mol pqr {self.pqrname}\n"
        text += "end\n"
        for elec in self.elecs:
            text += str(elec)
        for prints in self.prints:
            text += prints
        text += "\nquit\n"
        return text

    def print_input_files(self, output_path):
        """Generate the input file(s) associated with this object.

        :param output_path:  location for generated files
        :type output_path:  str
        """
        path = Path(output_path)
        base_name = path.stem
        if self.asyncflag:
            outname = base_name + "-para.in"
            # Temporarily disable async flag
            for elec in self.elecs:
                elec.asyncflag = False
            with open(outname, "w") as out_file:
                out_file.write(str(self))
            # Now make the async files
            elec = self.elecs[0]
            nproc = elec.pdime[0] * elec.pdime[1] * elec.pdime[2]
            for i in range(int(nproc)):
                outname = base_name + f"-PE{i}.in"
                for elec in self.elecs:
                    elec.asyncflag = True
                    elec.async_ = i
                with open(outname, "w") as out_file:
                    out_file.write(str(self))
        else:
            with open(path, "w") as out_file:
                out_file.write(str(self))

    def dump_pickle(self):
        """Make a Python pickle associated with the APBS input parameters.

        .. todo::  is this function still useful?
        """
        base_pqr_name = self.pqrpath.stem
        outname = base_pqr_name + "-input.p"
        with open(outname, "wb") as pfile:
            pickle.dump(self, pfile)


def split_input(filename):
    """Split the parallel input file into multiple async file names.

    :param filename:  the path to the original parallel input file
    :type filename:  str
    """
    nproc = 0
    with open(filename) as file_:
        text = ""
        while True:
            line = file_.readline()
            if not line:
                break
            text += line
            line = line.strip()
            if line.startswith("pdime"):  # Get # Procs
                words = line.split()
                nproc = int(words[1]) * int(words[2]) * int(words[3])
    if nproc == 0:
        errstr = f"{filename} is not a valid APBS parallel input file!\n"
        errstr = errstr + (
            "The inputgen script was unable to asynchronize this file!"
        )
        raise RuntimeError(errstr)
    base_pqr_name = Path(filename).stem
    for iproc in int(range(nproc)):
        outname = base_pqr_name + f"-PE{iproc}.in"
        outtext = text.replace("mg-para\n", f"mg-para\n    async {iproc}\n")
        with open(outname, "w") as outfile:
            outfile.write(outtext)


def build_parser():
    """Build argument parser."""
    desc = f"{TITLE_STR}\ninputgen: generating APBS input files since "
    desc += "(at least) 2004"
    parse = argparse.ArgumentParser(
        description=desc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parse.add_argument(
        "--asynch",
        action="store_true",
        help="perform an asynchronous parallel calculation.",
    )
    parse.add_argument(
        "--split",
        action="store_true",
        help=(
            "split an existing parallel input file to multiple "
            "async input files."
        ),
    )
    parse.add_argument(
        "--potdx",
        action="store_true",
        help=("create an input to compute an electrostatic potential map."),
    )
    parse.add_argument(
        "--method",
        help=("force output file to write a specific APBS ELEC method."),
        choices=["para", "auto", "manual", "async"],
    )
    parse.add_argument(
        "--cfac",
        type=float,
        default=psize.CFAC,
        help=(
            "factor by which to expand molecular dimensions to "
            "get coarse grid dimensions."
        ),
    )
    parse.add_argument(
        "--fadd",
        type=float,
        default=psize.FADD,
        help=(
            "amount to add to molecular dimensions to get fine "
            "grid dimensions."
        ),
    )
    parse.add_argument(
        "--space",
        type=float,
        default=psize.SPACE,
        help="desired fine mesh resolution",
    )
    parse.add_argument(
        "--gmemfac",
        type=int,
        default=psize.GMEMFAC,
        help=(
            "number of bytes per grid point required for sequential "
            "MG calculation"
        ),
    )
    parse.add_argument(
        "--gmemceil",
        type=int,
        default=psize.GMEMCEIL,
        help=(
            "max MB allowed for sequential MG calculation; adjust "
            "this to force the script to perform faster calculations "
            "(which require more parallelism)"
        ),
    )
    parse.add_argument(
        "--ofrac",
        type=float,
        default=psize.OFRAC,
        help="overlap factor between mesh partitions (parallel)",
    )
    parse.add_argument(
        "--redfac",
        type=float,
        default=psize.REDFAC,
        help=(
            "the maximum factor by which a domain dimension can "
            "be reduced during focusing"
        ),
    )
    parse.add_argument(
        "--istrng", help="Ionic strength (M); Na+ and Cl- ions will be used"
    )
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
        size.run_psize(filename)
        input_ = Input(
            filename,
            size,
            args.method,
            args.asynch,
            args.istrng,
            potdx=args.potdx,
        )
        path = Path(filename)
        output_path = f"{path.parent}/{path.stem}.in"
        input_.print_input_files(output_path)


if __name__ == "__main__":
    main()
