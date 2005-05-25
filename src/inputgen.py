#!/usr/bin/env python
# You may need to edit the above to point to your version of Python

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
# GMEMFAC = 200               # Number of bytes per grid point required 
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

class inputGen:
    """
        The inputGen class is used for creating input files for APBS
    """
    
    def __init__(self, pqrpath, size, method, async):
        """
            Initialize the inputGen class
            
            Parameters
                pqrpath: The full path to the PQR file (string)
                size:    The psize object with appropriate constant values (Psize)
                method:  The APBS ELEC method (string)
                async:   A flag to denote whether to generate
                         asynchronous input files or not (int)     
        """
        self.method = method
        self.async = async
        self.size = size
        self.coarsedim = []
        self.finedim = []
        self.procgrid = []
        self.finegridpoints = []
        self.id = None
        self.center = []
        self.fullpath = pqrpath
       
        i = string.rfind(pqrpath, "/") + 1
        self.pqrname = pqrpath[i:]
   
        self.setup()

    def getCenter(self):
        """
            Get the desired center for the APBS input file.
            Uses an arbitrary center if the input file is intended for a pKa
            calculation, otherwise centers on the center of the molecule.

            Returns
                center:  The center for the input file (string)
        """
        center = ""
        if self.center != []:
            if self.method == "mg-auto" or self.method == "mg-para":
                center =  "    cgcent %.3f %.3f %.3f\n" % \
                         (self.center[0], self.center[1], self.center[2])
                center += "    fgcent %.3f %.3f %.3f\n" % \
                         (self.center[0], self.center[1], self.center[2])
            elif self.method == "mg-manual":
                center =  "    gcent %.3f %.3f %.3f\n" % \
                         (self.center[0], self.center[1], self.center[2])
        else:
            if self.method == "mg-auto" or self.method == "mg-para":
                center  = "    cgcent mol 1\n"
                center += "    fgcent mol 1\n"
            elif self.method == "mg-manual":
                center  = "    gcent mol 1\n"
        return center
    
    def setCenter(self, center):
        """
            Set the center of the inputGen to a specific point

            Parameters
                center:  The desired center for the input file (string)
        """
        self.center = center
        
    def getText(self):
        """
            Get the text associated with the inputgen object

            Returns
                text:  The input file (string)
        """
        
        text  = "read\n"
        text += "    mol pqr %s\n" % self.pqrname
        text += "end\n"
        text += "elec\n"
        text += "    %s\n" % self.method
        if self.method == "mg-manual":
            text += "    dime %i %i %i\n" % (self.finegridpoints[0], self.finegridpoints[1], self.finegridpoints[2])
            text += "    nlev 4\n"
            text += "    glen %.3f %.3f %.3f\n" % (self.coarsedim[0], self.coarsedim[1], self.coarsedim[2])
        elif self.method == "mg-auto":
            text += "    dime %i %i %i\n" % (self.finegridpoints[0], self.finegridpoints[1], self.finegridpoints[2])
            text += "    cglen %.4f %.4f %.4f\n" % (self.coarsedim[0], self.coarsedim[1], self.coarsedim[2])
            text += "    fglen %.4f %.4f %.4f\n" % (self.finedim[0], self.finedim[1], self.finedim[2])
        elif self.method == "mg-para":
            text += "    pdime %i %i %i\n" % (self.procgrid[0], self.procgrid[1], self.procgrid[2])
            text += "    ofrac 0.1\n"
            text += "    dime %i %i %i\n" % (self.finegridpoints[0], self.finegridpoints[1], self.finegridpoints[2])
            text += "    cglen %.3f %.3f %.3f\n" % (self.coarsedim[0], self.coarsedim[1], self.coarsedim[2])
            text += "    fglen %.3f %.3f %.3f\n" % (self.finedim[0], self.finedim[1], self.finedim[2])
            if self.async == 1:
                text += "    async %i\n" % self.id
        text += self.getCenter()
        text += "    mol 1\n"                            
        text += "    lpbe\n"                             
        text += "    bcfl sdh\n"                           
        text += "    ion 1 0.150 2.0\n"            
        text += "    ion -1 0.150 2.0\n"           
        text += "    pdie 2.0\n"                
        text += "    sdie 78.54\n"                
        text += "    srfm smol\n"                   
        text += "    chgm spl2\n"
        text += "    srad 1.4\n"          
        text += "    swin 0.3\n"         
        text += "    temp 298.15\n"     
        text += "    gamma 0.105\n"    
        text += "    calcenergy total\n"
        text += "    calcforce no\n"
        text += "    write pot dx pot\n"
        text += "    write smol dx acc\n"
        text += "end\n"
        
        text += "elec\n"
        text += "    %s\n" % self.method
        if self.method == "mg-manual":
            text += "    dime %i %i %i\n" % (self.finegridpoints[0], self.finegridpoints[1], self.finegridpoints[2])
            text += "    nlev 4\n"
            text += "    glen %.3f %.3f %.3f\n" % (self.coarsedim[0], self.coarsedim[1], self.coarsedim[2])
        elif self.method == "mg-auto":
            text += "    dime %i %i %i\n" % (self.finegridpoints[0], self.finegridpoints[1], self.finegridpoints[2])
            text += "    cglen %.4f %.4f %.4f\n" % (self.coarsedim[0], self.coarsedim[1], self.coarsedim[2])
            text += "    fglen %.4f %.4f %.4f\n" % (self.finedim[0], self.finedim[1], self.finedim[2])
        elif self.method == "mg-para":
            text += "    pdime %i %i %i\n" % (self.procgrid[0], self.procgrid[1], self.procgrid[2])
            text += "    ofrac 0.1\n"
            text += "    dime %i %i %i\n" % (self.finegridpoints[0], self.finegridpoints[1], self.finegridpoints[2])
            text += "    cglen %.3f %.3f %.3f\n" % (self.coarsedim[0], self.coarsedim[1], self.coarsedim[2])
            text += "    fglen %.3f %.3f %.3f\n" % (self.finedim[0], self.finedim[1], self.finedim[2])
            if self.async == 1:
                text += "    async %i\n" % self.id
        text += self.getCenter()
        text += "    mol 1\n"                              
        text += "    lpbe\n"                               
        text += "    bcfl sdh\n"                            
        text += "    ion 1 0.150 2.0\n"              
        text += "    ion -1 0.150 2.0\n"            
        text += "    pdie 2.0\n"                  
        text += "    sdie 2.00\n"               
        text += "    srfm smol\n"                  
        text += "    chgm spl2\n"
        text += "    srad 1.4\n"        
        text += "    swin 0.3\n"        
        text += "    temp 298.15\n"    
        text += "    gamma 0.105\n"   
        text += "    calcenergy total\n"
        text += "    calcforce no\n"
        text += "end\n"
        text += "\nprint energy 2 - 1 end\n"
        text += "\nquit\n"
        return text
    
    def setup(self):
        """
            Do some setting up for input generation
        """
        size = self.size
        size.runPsize(self.fullpath)
        self.coarsedim = size.getCoarseGridDims()
        self.finedim = size.getFineGridDims()
        self.procgrid = size.getProcGrid()
        self.finegridpoints = size.getFineGridPoints()

        if self.async == 1:
            self.finegridpoints = size.getSmallest()
        else:
            n = self.finegridpoints
            gmem = 200.0 * n[0] * n[1] * n[2] / 1024 / 1024
            if self.method == "": #method not named
                if gmem > size.getConstant("GMEMCEIL"): self.method = "mg-para"
                else: self.method = "mg-auto"
            if self.method == "mg-para": self.finegridpoints = size.getSmallest()

    def printInput(self):
        """
            Print the input file(s) to the correct locations
        """
        if self.async == 1:
            self.method = "mg-para"
            self.async = 0
            period = string.find(self.fullpath,".")
            outname = self.fullpath[0:period] + "-para.in"

            file = open(outname, "w")
            file.write(self.getText())
            file.close()
            
            self.async = 1
            nproc = self.procgrid[0] * self.procgrid[1] * self.procgrid[2]
            for i in range(nproc):
                period = string.find(self.fullpath,".")
                outname = self.fullpath[0:period] + "-PE%i.in" % i
                self.id = i
                file = open(outname, "w")
                file.write(self.getText())
                file.close()
        
        else:
            period = string.find(self.fullpath,".")
            if period > 0:
                outname = self.fullpath[0:period] + ".in"
            else:
                outname = self.fullpath + ".in"
            file = open(outname, "w")
            file.write(self.getText())
            file.close()

def splitInput(filename):
    """
        Split the parallel input file into multiple async file names
    """
    nproc = 0
    file = open(filename)
    text = ""
    while 1:
        line = file.readline()
        if line == "": break
        text += line
        line = string.strip(line)
        if line.startswith("pdime"): # Get # Procs
            words = string.split(line)
            nproc = int(words[1]) * int(words[2]) * int(words[3])

    if nproc == 0:
        sys.stderr.write("%s is not a valid APBS parallel input file!\n" % filename)
        sys.stderr.write("The inputgen script was unable to asynchronize this file!\n")
        sys.exit(2)

    period = string.find(filename,".")
    for i in range(nproc):
        outname = filename[0:period] + "-PE%i.in" % i
        outtext = string.replace(text, "mg-para\n","mg-para\n    async %i\n" % i)
        outfile = open(outname, "w")
        outfile.write(outtext)
        outfile.close()
          
def usage():
    """
        Display the usage information for this script
    """
    size = psize.Psize()
    usage = "\n"
    usage = usage + "Use this script to generate new APBS input files or split an existing\n"
    usage = usage + "parallel input file into multiple async files.\n\n"
    usage = usage + "Usage: inputgen.py [opts] <filename>\n"
    usage = usage + "Optional Arguments:\n"
    usage = usage + "  --help               : Display this text\n"
    usage = usage + "  --split              : Split an existing parallel input file to multiple\n"
    usage = usage + "                         async input files.\n"
    usage = usage + "  --METHOD=<value>     : Force output file to write a specific APBS ELEC\n"
    usage = usage + "                         method.  Options are para (parallel), auto\n"
    usage = usage + "                         (automatic), manual (manual), or async (asynchronous).\n"
    usage = usage + "                         solve.  async will result in multiple input files.\n"
    usage = usage + "  --CFAC=<value>       : Factor by which to expand molecular dimensions to\n"
    usage = usage + "                         get coarse grid dimensions.\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("CFAC")
    usage = usage + "  --FADD=<value>       : Amount to add to molecular dimensions to get fine\n"
    usage = usage + "                         grid dimensions.\n"
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
    longOptList = ["help","split","METHOD=","CFAC=","SPACE=","GMEMCEIL=","GMEMFAC=","OFAC=","REDFAC="]

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

    method = ""
    size = psize.Psize()
    async = 0
    split = 0
    
    for o, a in opts:
        if o == "--help":
            usage()
        if o == "--split": split = 1
        if o == "--METHOD":
            if a == "para":
                sys.stdout.write("Forcing a parallel calculation\n")
                method = "mg-para"
            elif a == "auto":
                sys.stdout.write("Forcing a sequential calculation\n")
                method = "mg-auto"
            elif a == "async":
                sys.stdout.write("Forcing an asynchronous calculation\n")
                method = "mg-para"
                async = 1
            elif a == "manual":
                sys.stdout.write("Forcing a manual calculation\n")
                method = "mg-manual"
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

    if split == 1:
        splitInput(filename)
    else:
        igen = inputGen(filename, size, method, async)
        igen.printInput()

if __name__ == "__main__": main()
