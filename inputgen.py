#!/usr/bin/env python2
# You may need to edit the above to point to your version of Python 2.0



#
# inputgen.py 
# Create an input file for APBS using psize data
#
# Written by Todd Dolinsky based on original sed script by Nathan Baker
#
# Version:  $Id$
#
# APBS -- Adaptive Poisson-Boltzmann Solver
#
# Nathan A. Baker (baker@biochem.wustl.edu)
# Dept. of Biochemistry and Molecular Biophysics
# Center for Computational Biology
# Washington University in St. Louis
#
# Additional contributing authors listed in the code documentation.
#
# Copyright (c) 2002-2004.  Washington University in St. Louis.
# All Rights Reserved.
# Portions Copyright (c) 1999-2002.  The Regents of the University of
# California.
# Portions Copyright (c) 1995.  Michael Holst.
#
# This file is part of APBS.
#
# APBS is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# APBS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with APBS; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
#


# User - Definable Variables: Default values

# CFAC = 1.7                  # Factor by which to expand mol dims to
                              # get coarse grid dims
# FADD = 20                   # Amount to add to mol dims to get fine
                              # grid dims
# SPACE = 0.50                # Desired fine mesh resolution
# GMEMFAC = 160               # Number of bytes per grid point required 
                              # for sequential MG calculation 
# GMEMCEIL = 400              # Max MB allowed for sequential MG
                              # calculation.  Adjust this to force the
                              # script to perform faster calculations (which
                              # require more parallelism).
# OFAC = 0.1                  # Overlap factor between mesh partitions
# REDFAC = 0.25               # The maximum factor by which a domain
                              # dimension can be reduced during focusing

import string, sys
import psize

class inputgen:
    def __init__(self):
        self.method = ""
        self.async = 0

    def printInput(self, pqrname, coarsedim, finedim, procgrid, finegridpoints, outname, id):

        file = open(outname, "w")

        file.write("read\n")
        file.write("    mol pqr %s\n" % pqrname)
        file.write("end\n")
        file.write("elec\n")
        file.write("    %s\n" % self.method)
        if self.method == "mg-manual":
            file.write("    dime %i %i %i\n" % (finegridpoints[0], finegridpoints[1], finegridpoints[2]))
            file.write("    nlev 4\n")
            file.write("    glen %.3f %.3f %.3f\n" % (coarsedim[0], coarsedim[1], coarsedim[2]))
            file.write("    gcent mol 1\n")
        elif self.method == "mg-auto":
            file.write("    dime %i %i %i\n" % (finegridpoints[0], finegridpoints[1], finegridpoints[2]))
            file.write("    cglen %.4f %.4f %.4f\n" % (coarsedim[0], coarsedim[1], coarsedim[2]))
            file.write("    fglen %.4f %.4f %.4f\n" % (finedim[0], finedim[1], finedim[2]))
            file.write("    cgcent mol 1\n")
            file.write("    fgcent mol 1\n")  
        elif self.method == "mg-para":
            file.write("    pdime %i %i %i\n" % (procgrid[0], procgrid[1], procgrid[2]))
            file.write("    ofrac 0.1\n")
            file.write("    dime %i %i %i\n" % (finegridpoints[0], finegridpoints[1], finegridpoints[2]))
            file.write("    cglen %.3f %.3f %.3f\n" % (coarsedim[0], coarsedim[1], coarsedim[2]))
            file.write("    fglen %.3f %.3f %.3f\n" % (finedim[0], finedim[1], finedim[2]))
            file.write("    cgcent mol 1\n")
            file.write("    fgcent mol 1\n")
            if self.async == 1:
                file.write("    async %i\n" % id)
        file.write("    mol 1\n")                             
        file.write("    lpbe\n")                              
        file.write("    bcfl sdh\n")                           
        file.write("    ion 1 0.150 2.0\n")             
        file.write("    ion -1 0.150 2.0\n")           
        file.write("    pdie 2.0\n")                  
        file.write("    sdie 78.54\n")               
        file.write("    srfm smol\n")                  
        file.write("    chgm spl2\n")
        file.write("    srad 1.4\n")         
        file.write("    swin 0.3\n")        
        file.write("    temp 298.15\n")    
        file.write("    gamma 0.105\n")   
        file.write("    calcenergy total\n")
        file.write("    calcforce no\n")
        #file.write("    write pot dx pot%s\n" % pqrname[0:4])
        #file.write("    write smol dx acc%s\n" % pqrname[0:4])
        file.write("end\n")
        
        file.write("elec\n")
        file.write("    %s\n" % self.method)
        if self.method == "mg-manual":
            file.write("    dime %i %i %i\n" % (finegridpoints[0], finegridpoints[1], finegridpoints[2]))
            file.write("    nlev 4\n")
            file.write("    glen %.3f %.3f %.3f\n" % (coarsedim[0], coarsedim[1], coarsedim[2]))
            file.write("    gcent mol 1\n")
        elif self.method == "mg-auto":
            file.write("    dime %i %i %i\n" % (finegridpoints[0], finegridpoints[1], finegridpoints[2]))
            file.write("    cglen %.4f %.4f %.4f\n" % (coarsedim[0], coarsedim[1], coarsedim[2]))
            file.write("    fglen %.4f %.4f %.4f\n" % (finedim[0], finedim[1], finedim[2]))
            file.write("    cgcent mol 1\n")
            file.write("    fgcent mol 1\n")  
        elif self.method == "mg-para":
            file.write("    pdime %i %i %i\n" % (procgrid[0], procgrid[1], procgrid[2]))
            file.write("    ofrac 0.1\n")
            file.write("    dime %i %i %i\n" % (finegridpoints[0], finegridpoints[1], finegridpoints[2]))
            file.write("    cglen %.3f %.3f %.3f\n" % (coarsedim[0], coarsedim[1], coarsedim[2]))
            file.write("    fglen %.3f %.3f %.3f\n" % (finedim[0], finedim[1], finedim[2]))
            file.write("    cgcent mol 1\n")
            file.write("    fgcent mol 1\n")
            if self.async == 1:
                file.write("    async %i\n" % id)
        file.write("    mol 1\n")                             
        file.write("    lpbe\n")                              
        file.write("    bcfl sdh\n")                           
        file.write("    ion 1 0.150 2.0\n")             
        file.write("    ion -1 0.150 2.0\n")           
        file.write("    pdie 2.0\n")                  
        file.write("    sdie 2.00\n")               
        file.write("    srfm smol\n")                  
        file.write("    chgm spl2\n")
        file.write("    srad 1.4\n")         
        file.write("    swin 0.3\n")        
        file.write("    temp 298.15\n")    
        file.write("    gamma 0.105\n")   
        file.write("    calcenergy total\n")
        file.write("    calcforce no\n")
        file.write("end\n")
        file.write("\nprint energy 2 - 1 end\n")
        file.write("\nquit\n")
        
        file.close()

    def runinputgen(self, filename, size):
        size.runPsize(filename)
        coarsedim = size.getCoarseGridDims()
        finedim = size.getFineGridDims()
        procgrid = size.getProcGrid()
        n = size.getFineGridPoints()
        gmem = 160.0 * n[0] * n[1] * n[2] / 1024 / 1024

        i = string.rfind(filename, "/") + 1
        pqrname = filename[i:]

        if self.async == 1:
            returns = []
            self.method = "mg-para"
            self.async = 0
            j = string.find(filename,".")
            outname = filename[0:j] + "-para.in"
            returns.append(outname)
            self.printInput(pqrname, coarsedim, finedim, procgrid, n, outname, 0)
            self.async = 1
            nproc = procgrid[0] * procgrid[1] * procgrid[2]
    
            for i in range(nproc):
                j = string.find(filename,".")
                outname = filename[0:j] + "-PE%s.in" % str(i)
                returns.append(outname)
                self.printInput(pqrname, coarsedim, finedim, procgrid, n, outname, i)
            return returns
        
        else:
            if self.method == "": #method not named
                if gmem > size.getConstant("GMEMCEIL"):
                    #sys.stderr.write("***Estimated mem. required larger than allowed - Parallel solve required\n")
                    self.method = "mg-para"
                else:
                    #sys.stderr.write("***Sequential solve allowed\n")
                    self.method = "mg-auto"
            j = string.find(filename,".")
            outname = filename[0:j] + ".in"
            self.printInput(pqrname, coarsedim, finedim, procgrid, n, outname, 0)
            return outname

def usage():
    size = psize.Psize()
    usage = "\n"
    usage = usage + "Inputgen script\n"
    usage = usage + "Usage: inputgen.py [opts] <filename>\n"
    usage = usage + "Optional Arguments:\n"
    usage = usage + "  --help               : Display this text\n"
    usage = usage + "  --METHOD=<value>     : Force output file to write for parallel\n"
    usage = usage + "                         (para), automatic (auto) or asynchronous (async)\n"
    usage = usage + "                         solve.  async will result in multiple input files.\n"
    usage = usage + "  --CFAC=<value>       : Factor by which to expand mol dimsto\n"
    usage = usage + "                         get coarse grid dims\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("CFAC")
    usage = usage + "  --FADD=<value>       : Amount to add to mol dims to get fine\n"
    usage = usage + "                         grid dims\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("FADD")
    usage = usage + "  --SPACE=<value>      : Desired fine mesh resolution\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("SPACE")
    usage = usage + "  --GMEMFAC=<value>    : Number of bytes per grid point required\n"
    usage = usage + "                         for sequential MG calculation\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("GMEMFAC")
    usage = usage + "  --GMEMCEIL=<value>   : Max MB allowed for sequential MG\n"
    usage = usage + "                         calculation.  Adjust this to force the\n"
    usage = usage + "                         script to perform faster calculations (which\n"
    usage = usage + "                         require more parallelism).\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("GMEMCEIL")
    usage = usage + "  --OFAC=<value>       : Overlap factor between mesh partitions\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("OFAC")
    usage = usage + "  --REDFAC=<value>     : The maximum factor by which a domain\n"
    usage = usage + "                         dimension can be reduced during focusing\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("REDFAC")
    sys.stderr.write(usage)
    sys.exit(2)

def main():

    import getopt
    filename = ""
    shortOptList = ""
    longOptList = ["help","METHOD=","CFAC=","SPACE=","GMEMCEIL=","GMEMFAC=","OFAC=","REDFAC="]

    try:
        opts, args = getopt.getopt(sys.argv[1:], shortOptList, longOptList)
    except getopt.GetoptError, details:
        sys.stderr.write("Option error (%s)!\n" % details)
        usage()
        
    if len(args) != 1:
        sys.stderr.write("Invalid argument list!\n")
        usage()
    else:
        filename = args[0]

    igen = inputgen()
    size = psize.Psize()
    async = 0
    
    for o, a in opts:
        if o == "--help":
            usage()
        if o == "--METHOD":
            if a == "para":
                sys.stdout.write("Forcing a parallel calculation\n")
                igen.method = "mg-para"
            elif a == "auto":
                sys.stdout.write("Forcing a sequential calculation\n")
                igen.method = "mg-auto"
            elif a == "async":
                sys.stdout.write("Forcing an asynchronous calculation\n")
                igen.method = "mg-para"
                igen.async = 1
            else:
                sys.stdout.write("Incorrect method argument: %s\n" % a)
                sys.stdout.write("Defaulting to memory dependent result\n")
        if o == "--CFAC":
            size.setConstant("CFAC", float(a))
        if o == "--SPACE":
            size.setConstant("SPACE", float(a))
        if o == "--GMEMFAC":
            size.setConstant("GMEMFAC", int(a))
        if o == "--GMEMCEIL":
            size.setConstant("GMEMCEIL",  int(a))
        if o == "--OFAC":
            size.setConstant("OFAC", float(a))
        if o == "--REDFAC":
            size.setConstant("REDFAC", float(a))
            
    igen.runinputgen(filename, size)

if __name__ == "__main__": main()
