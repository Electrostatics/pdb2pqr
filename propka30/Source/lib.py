#!/usr/bin/python
#
# * This library is free software; you can redistribute it and/or
# * modify it under the terms of the GNU Lesser General Public
# * License as published by the Free Software Foundation; either
# * version 2.1 of the License, or (at your option) any later version.
# *
# * This library is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# * Lesser General Public License for more details.
#

#propka3.0, revision 182                                                                      2011-08-09
#-------------------------------------------------------------------------------------------------------
#--                                                                                                   --
#--                                   PROPKA: A PROTEIN PKA PREDICTOR                                 --
#--                                                                                                   --
#--                              VERSION 3.0,  01/01/2011, COPENHAGEN                                 --
#--                              BY MATS H.M. OLSSON AND CHRESTEN R. SONDERGARD                       --
#--                                                                                                   --
#-------------------------------------------------------------------------------------------------------
#
#
#-------------------------------------------------------------------------------------------------------
# References:
#
#   Very Fast Empirical Prediction and Rationalization of Protein pKa Values
#   Hui Li, Andrew D. Robertson and Jan H. Jensen
#   PROTEINS: Structure, Function, and Bioinformatics 61:704-721 (2005)
#
#   Very Fast Prediction and Rationalization of pKa Values for Protein-Ligand Complexes
#   Delphine C. Bas, David M. Rogers and Jan H. Jensen
#   PROTEINS: Structure, Function, and Bioinformatics 73:765-783 (2008)
#
#   PROPKA3: Consistent Treatment of Internal and Surface Residues in Empirical pKa predictions
#   Mats H.M. Olsson, Chresten R. Sondergard, Michal Rostkowski, and Jan H. Jensen
#   Journal of Chemical Theory and Computation, 7, 525-537 (2011)
#-------------------------------------------------------------------------------------------------------


import string, sys, copy, math
import output


def loadOptions():
    """
    load the arguments parser with options
    """
    from optparse import OptionParser

    # printing out header before parsing input
    output.printHeader()

    # printing out header before parsing input
    checkPythonVersion()

    # defining a 'usage' message
    usage = "usage: %prog [options] filename"

    # creating a parser
    parser = OptionParser(usage)

    # loading the parser
    parser.add_option("-f", "--file", action="append", dest="filenames", 
           help="read data from <filename>, i.e. <filename> is added to arguments")
    parser.add_option("-r", "--reference", dest="reference", default="neutral", 
           help="setting which reference to use for stability calculations [neutral/low-pH]")
    parser.add_option("-c", "--chain", action="append", dest="chains", 
           help="creating the protein with only a specified chain, note, chains without ID are labeled 'A' [all]")
    parser.add_option("-t", "--thermophile", action="append", dest="thermophiles", 
           help="defining a thermophile filename; usually used in 'alignment-mutations'")
    parser.add_option("-a", "--alignment", action="append", dest="alignment", 
           help="alignment file connecting <filename> and <thermophile> [<thermophile>.pir]")
    parser.add_option("-m", "--mutation", action="append", dest="mutations", 
           help="specifying mutation labels which is used to modify <filename> according to, e.g. N25R/N181D")
    parser.add_option("-v", "--version", dest="version_label", default="Nov30", 
           help="specifying the sub-version of propka [Jan15/Dec19]")
    parser.add_option("-z", "--verbose", dest="verbose", action="store_true", default=True, 
           help="sleep during calculations")
    parser.add_option("-q", "--quiet", dest="verbose", action="store_false",
           help="sleep during calculations")
    parser.add_option(      "--mute", dest="verbose", action="store_false",
           help="sleep during calculations")
    parser.add_option("-s", "--silent",  dest="verbose", action="store_false", 
           help="not activated yet")
    parser.add_option("--verbosity",  dest="verbosity", action="store_const", 
           help="level of printout - not activated yet")
    parser.add_option("--protonation", dest="protonation", default="old-school", 
           help="setting protonation scheme")
    parser.add_option("-p", "--pH", dest="pH", type="float", default=7.0, 
           help="setting pH-value used in e.g. stability calculations [7.0]")
    parser.add_option("--window", dest="window", nargs=3, type="float", default=(0.0, 14.0, 1.0),
           help="setting the pH-window to show e.g. stability profiles [0.0, 14.0, 1.0]")
    parser.add_option("--grid",   dest="grid",   nargs=3, type="float", default=(0.0, 14.0, 0.1),
           help="setting the pH-grid to calculate e.g. stability related properties [0.0, 14.0, 0.1]")
    parser.add_option("--mutator", dest="mutator", 
           help="setting approach for mutating <filename> [alignment/scwrl/jackal]")
    parser.add_option("--mutator-option", dest="mutator_options", action="append",
           help="setting property for mutator [e.g. type=\"side-chain\"]")

    parser.add_option("-d","--display-coupled-residues", dest="display_coupled_residues", action="store_true",
           help="Displays alternative pKa values due to coupling of titratable groups")
    parser.add_option("--print-iterations", dest="print_iterations", action="store_true",
           help="Displays the pKa iterations in the Tanford-Roxby scheme")


    # parsing and returning options and arguments
    options, args = parser.parse_args()

    # adding specified filenames to arguments
    if options.filenames:
      for filename in options.filenames:
        args.append(filename)

    # checking at early stage that there is at least one pdbfile to work with
    if len(args) == 0:
      pka_print("Warning: no pdbfile provided")
      #sys.exit(9)


    # --- post-processing; interpreting some of the arguments ---

    # interpreting 'options.mutator' and switching to object
    interpretMutator(options)

    # interpreting 'options.mutations' and switching to dictionary
    if False:
      if options.mutations != None:
        interpretDictionaryMutations(options)

    # setting/checking default alignment files
    setDefaultAlignmentFiles(options)

    # done!
    return  options, args

proPKA_verbose = False
def setVerbose(value):
    global proPKA_verbose
    proPKA_verbose = value
    
def pka_print(txt):
    if proPKA_verbose:
        print(txt)

def checkPythonVersion():
    """
    checking that this is python 2.6.0 or later
    """

    if sys.hexversion < 0x02060000:
      pka_print("propka does not run under python %s, please use version 2.6 or later" % (sys.version.split()[0]))
      sys.exit(8)

    return


def interpretMutator(options):
    """
    setting the mutator defined by options.mutator
    """
    # creating a mutator object
    myMutator = Mutator(label=options.mutator)

    # setting mutator-options
    if options.mutator_options == None:
      """ do nothing """
    else:
      for item in options.mutator_options:
        if item[:4] == "type":
          property = "type"; value = "\"%s\"" % item[5:]
        elif item[:3] == "prm":
          property = "prm"; value = item[4:]
        elif item[:3] == "min":
          property = "min"; value = item[4:]
        elif item[:3] == "ini":
          property = "ini"; value = item[4:]
        elif item[:3] == "rtm":
          property = "rtm"; value = item[4:]
        cmd = "myMutator.setProperty(%s=%s)" % (property, value)
        pka_print(cmd)
        exec(cmd)

    # resetting mutator
    options.mutator = myMutator


def setDefaultAlignmentFiles(options):
    """
    setting the default alignmentfiles to [<thermophile.pir>, ...]
    """
    if options.mutator.label in ["alignment", "overlap"] and options.alignment == None:
      options.alignment = []
      for mutation in options.mutations:
        if isinstance(mutation, str):
          pdbcode = mutation[:4]
          filename = "%s.pir" % (pdbcode)
          if filename not in options.alignment:
            options.alignment.append( filename )
        else:
          for code in mutation.keys():
            filename = "%s.pir" % ( extractName(code) )
            if filename not in options.alignment:
              options.alignment.append( filename )


def interpretMutationsDictionary_old(options):
    """
    interprets the mutations in options; i.e. separating pdb-key and mutation label 
    in e.g. '2vuj:N25R/N181D' - trying to use dictionary
    """
    mutations = []
    for mutation_line in options.mutations:
      separator = None
      for i in range(len(mutation)):
        if mutation[i] == ":":
          separator = i
      if separator:
        code  = mutation[:separator]
        label = mutation[separator+1:]
      else:
        # pdbcode None, trying to resolve ambigous situation
        if options.mutator.label == "alignment":
          if len(options.thermophiles) == 1:
            code  = extractName(options.thermophiles[0])
            label = mutation
          else:
            pka_print("cannot assign pdbcode to mutation %s; specify pdbcode or use a different mutator" % (mutation))
            sys.exit(9)
        else:
          # not using alignment anyway
          code  = None
          label = mutation

      if code in mutations:
        mutations[code].append(label)
      else:
        mutations[code] = [label]
    
    # resetting the content of 'options.mutations to dictionary
    options.mutations = mutations


def interpretDictionaryMutations(options):
    """
    interprets the mutations in options; i.e. trying to understand "pdb:chainID:mutation" or any combinations
    like that; e.g. '2vuj:A:N25R/N181D', '2vuj:N25R/N181D' or 'N25R/N181D'
    """
    mutations = []
    
    for mutation_line in options.mutations:
      separator = []
      for i in range(len(mutation_line)):
        if mutation_line[i] == ":":
          separator.append(i)
      if   len(separator) == 0:
        code    = None
        chainID = None
        label   = mutation_line
      elif len(separator) == 1:
        # trying to resolve ambigous situation
        if len( mutation_line[:separator[0]] ) == 1:
          code    = None
          chainID = mutation_line[:separator[0]]
        else:
          code    = mutation_line[:separator[0]]
          chainID = None
        label   = mutation_line[separator[-1]+1:]
      elif len(separator) == 2:
        code    = mutation_line[:separator[0]]
        chainID = mutation_line[separator[0]+1:separator[1]]
        label   = mutation_line[separator[-1]+1:]
      else:
        # too difficult to figure out what user want
        pka_print("cannot interpret mutation '%s'; specify as '2vuj:A:N25R/N181D'" % (mutation))
        sys.exit(9)

      mutation = {code: {chainID: label}}
      mutations.append( mutation )
      #if code != None and code not in options.thermophiles:
      #  options.thermophiles.append(code)

    pka_print("interpreted mutations as:")
    for mutation in mutations:
      pka_print(mutation)

    # resetting the content of 'options.mutations to new list
    options.mutations = mutations


def interpretMutationsList(options):
    """
    interprets the mutations in options; i.e. trying to understand "pdb:chainID:mutation" or any combinations
    like that; e.g. '2vuj:A:N25R/N181D', '2vuj:N25R/N181D' or 'N25R/N181D'
    """
    mutations = []
    for mutation in options.mutations:
      separator = []
      for i in range(len(mutation)):
        if mutation[i] == ":":
          separator.append(i)
      if   len(separator) == 0:
        code    = None
        chainID = None
        label   = mutation
      elif len(separator) == 1:
        # trying to resolve ambigous situation
        if len( mutation[:separator[0]] ) == 1:
          code    = None
          chainID = mutation[:separator[0]]
        else:
          code    = mutation[:separator[0]]
          chainID = None
        label   = mutation[separator[-1]+1:]
      elif len(separator) == 2:
        code    = mutation[:separator[0]]
        chainID = mutation[separator[0]+1:separator[1]]
        label   = mutation[separator[-1]+1:]
      else:
        # too difficult to figure out what user want
        pka_print("cannot interpret mutation '%s'; specify as '2vuj:A:N25R/N181D'" % (mutation))
        sys.exit(9)
      mutations.append( [code, chainID, label] )
      if code != None and code not in options.thermophiles:
        options.thermophiles.append(code)
    
    # resetting the content of 'options.mutations to new list
    options.mutations = mutations


def residueList(name):
    """
    Creates a list of residue labels for use n various parts of the program
    """
    if name == "all":
        residue_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', \
                        'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'C- ', 'N+ ']
    elif name == "standard":
        residue_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', \
                        'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    elif name == "coulomb":
        residue_list = ['ASP', 'GLU', 'HIS', 'CYS', 'TYR', 'LYS', 'ARG', 'C- ', 'N+ ']
    elif name == "propka1":
        residue_list = ['ASP', 'GLU', 'C- ', 'HIS', 'CYS', 'TYR', 'LYS', 'ARG', 'N+ ']
    elif name == "acids":
        residue_list = ['ASP', 'GLU', 'CYS', 'TYR', 'C- ']
    elif name == "bases":
        residue_list = ['ARG', 'LYS', 'HIS', 'N+ ']
    elif name == "protonated":
        residue_list = ['ASN', 'GLN', 'TRP', 'ARG', 'HIS']
    elif name == "excluded":
        residue_list = ["H2O", "HOH", "SO4", "PO4", "HEM", "cu", "zn", "GTT", "PMS", "MYG"]
    else:
        pka_print('cannot understand \"%s\" residueList' % (name))
        sys.exit(0)

    return residue_list



def residueInteractionList(grpName):
    """
    returns a list of residues to determine the pKa determinants.
    """
    residue_list = []; str = ""
    if   grpName == "COO":
        str = "SER THR CYS TYR LYS ASN GLN ARG TRP HIS"
    elif grpName == "CYS":
        str = "SER THR COO TYR LYS ASN GLN ARG TRP HIS"
    elif grpName == "TYR":
        str = "SER THR CYS TYR LYS ASN GLN ARG TRP HIS"
    elif grpName == "ALL":
        #str = "ASP GLU SER THR ASN GLN TRP HIS CYS TYR LYS ARG 'C- ' "
        residue_list = ['ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'TRP', 'HIS', 'CYS', 'TYR', 'LYS', 'ARG', 'C- ', 'N+ ']
    else:
        pka_print("Don't understand %s" % (grpName))
        sys.exit(0)

    if len(residue_list) == 0:
        return  string.split(str)
    else:
        return  residue_list


def residueCenterAtomList(resName):
    """
    Creates a list of atom names that constitutes the 'residue center'
    """
    names = {}
    names["GLU"] = ["OE1", "OE2"]
    names["ASP"] = ["OD1", "OD2"]
    names["HIS"] = ["CG" , "ND1", "CD2", "NE2", "CE1"]
    names["CYS"] = ["SG"]
    names["TYR"] = ["OH"]
    names["SER"] = ["OG"]
    names["THR"] = ["OG1"]
    names["LYS"] = ["NZ"]
    names["ARG"] = ["CZ"]
    names["GLN"] = ["OE1", "NE2"]
    names["ASN"] = ["OD1", "ND2"]
    names["TRP"] = ["NE1"]
    names["N+ "] = ["N"]
    names["C- "] = ["O", "OXT"]

    if resName in names:
      return names[resName]
    else:
      return []


def atomList(resName):
    """
    Creates a list of heavy atoms for each residue - used, e.g., to write the atoms in a 'correct order'
    in the new pdbfile.
    """
    atom_list = []
    str = ""

    if resName == "GLY":
        str = "N CA C O"
    elif resName == "ALA":
        str = "N CA C O CB"
    elif resName == "ASP":
        str = "N CA C O CB CG OD1 OD2"
    elif resName == "ASN":
        str = "N CA C O CB CG OD1 ND2"
    elif resName == "ARG":
        str = "N CA C O CB CG CD NE CZ NH1 NH2"
    elif resName == "CYS":
        str = "N CA C O CB SG"
    elif resName == "GLN":
        str = "N CA C O CB CG CD OE1 NE2"
    elif resName == "GLU":
        str = "N CA C O CB CG CD OE1 OE2"
    elif resName == "HIS":
        str = "N CA C O CB CG ND1 CD2 CE1 NE2"
    elif resName == "ILE":
        str = "N CA C O CB CG1 CG2 CD1"
    elif resName == "LEU":
        str = "N CA C O CB CG CD1 CD2"
    elif resName == "LYS":
        str = "N CA C O CB CG CD CE NZ"
    elif resName == "MET":
        str = "N CA C O CB CG SD CE"
    elif resName == "PHE":
        str = "N CA C O CB CG CD1 CD2 CE1 CE2 CZ"
    elif resName == "PRO":
        str = "N CA C O CB CG CD"
    elif resName == "SER":
        str = "N CA C O CB OG"
    elif resName == "THR":
        str = "N CA C O CB OG1 CG2"
    elif resName == "TRP":
        str = "N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2"
    elif resName == "TYR":
        str = "N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH"
    elif resName == "VAL":
        str = "N CA C O CB CG1 CG2"
    elif resName == "N+ ":
        str = ""
    elif resName == "C- ":
        str = ""
    else:
        pka_print("Don't understand %s in atomList(resName)" % (resName))
        sys.exit(0)

    atom_list = str.split()
    return atom_list



def extractName(filename):
    """
    Creates name without initial directory and extension, i.e. '/home/molsson/1xnb.pdb' gives '1xnb'
    """
    root, extension = splitFileName(filename)

    return root


def splitFileName(filename):
    """
    splits a filename into root and extension, e.g. '1xnb.pdb' gives ['1xnb', 'pdb']
    """
    name = getFileName(filename)
    dot = None
    for i in range(len(name)):
      if name[i] == '.':
        dot = i

    if dot == None:
      root = name; extension = None
    else:
      root = name[:dot]; extension = name[dot+1:]

    return root, extension


def getFileName(name):
    """
    reads the actual file name, e.g. '/home/molsson/1xnb.pdb' gives '1xnb.pdb'
    """
    start = None
    for i in range(len(name)):
      if name[i] == '/':
        start = i+1

    return name[start:]



def checkBuried(Nmass1, Nmass2):
    """
    returns True if an interaction is buried
    """

    if (Nmass1 + Nmass2 <= 900) and (Nmass1 <= 400 or Nmass2 <= 400):
        return False
    else:
        return True


def makeResidueLabel(resName, resNumb, chainID):
    """
    just returning a residue label
    """
    label = "%s" % (resName)
    for i in range (len(label), 3):
      label += " "
    label += "%4d%2s" % (resNumb, chainID)

    return label


def sortConfigurationKeys(keys):
    """
    find the 'default' configuration key
    """
    sorted_keys = []
    for key in keys:
      insert = False
      for i in range( len(sorted_keys) ):
        if   int(key[1:-2]) < int(sorted_keys[i][1:-2]):
          insert = True
        elif int(key[1:-2]) == int(sorted_keys[i][1:-2]) and ord(key[-1]) < ord(sorted_keys[i][-1]):
          insert = True
        if insert == True:
          sorted_keys.insert(i, key)
          break
      if insert == False:
        sorted_keys.append(key)

    return  sorted_keys


def get_sorted_configurations(configuration_keys):
    """
    extract and sort configurations
    """
    configurations = list(configuration_keys)
    configurations.sort(key=configuration_compare)
    return configurations

def configuration_compare(conf):
    return 100*int(conf[1:-2]) + ord(conf[-1])



def groupName(resName):
    """
    returns a group interaction label
    """

    if   resName == "ASP" or resName == "GLU" or resName == "C- ":
        name = "COO"
    elif resName == "ASN" or resName == "GLN":
        name = "AMD"
    elif resName == "SER" or resName == "THR":
        name = "ROH"
    else:
        name = resName

    return name


def groupList(grpName):
    """
    returns a list of resNames corresponding to grpName
    """
    grpList = {
               'ASP': ["ASP"],
               'GLU': ["GLU"],
               'C- ': ["C- "],
               'CYS': ["CYS"],
               'TYR': ["TYR"],
               'HIS': ["HIS"],
               'LYS': ["LYS"],
               'ARG': ["ARG"],
               'COO': ["ASP", "GLU", "C- "],
              }

    return grpList[grpName]


def convertResidueCode(code=None, resName=None):
    """
    returns the code and resName label
    """
    data = [["A", "ALA"],
            ["R", "ARG"],
            ["N", "ASN"],
            ["D", "ASP"],
            ["C", "CYS"],
            ["Q", "GLN"],
            ["E", "GLU"],
            ["G", "GLY"],
            ["H", "HIS"],
            ["I", "ILE"],
            ["L", "LEU"],
            ["K", "LYS"],
            ["M", "MET"],
            ["F", "PHE"],
            ["P", "PRO"],
            ["S", "SER"],
            ["T", "THR"],
            ["W", "TRP"],
            ["Y", "TYR"],
            ["V", "VAL"]]

    for test1, test2 in data:
      if code == test1 or resName == test2:
        return test1, test2

    pka_print("could not figure out code=%s resName=%s" % (code, resName) )
    sys.exit(9)



def extractResidueNumbers(description_line):
    """
    extracts the first and last residue numbers
    """
    stops = []
    for i in range(0, len(description_line)):
        if description_line[i] == ':':
            stops.append(i)
    first = int(description_line[ (stops[1]+1) : (stops[2]) ])
    last  = int(description_line[ (stops[3]+1) : (stops[4]) ])
    #print 'residues: %4d to %4d' % (first, last)

    return first, last



def readAlignments(filename, name):
    """
    reads alignment from alignment file
    """
    FOUND = False
    alignment  = ""
    alignments = []
    text       = []
    file = open(filename)
    # searching for correct protein entry
    while True:
      line = file.readline()
      if line == "":
        break
      if line[:4] == ">P1;":
        if name == line[4:].strip():
          pka_print( "%s=%s" % (name, line[4:].strip()) )
          FOUND = True
          text.append(line)
          line = file.readline()
          text.append(line)
          first, last = extractResidueNumbers(line)
          break

    # getting the alignment data
    if FOUND == False:
      pka_print( "could not find alignment for %s" % (name) )
    else:
      for line in file.readlines():
        if line == "\n":
            break
        positions = len(line)-1
        if line[0] == '/':
          alignments.append(alignment)
          alignment = ""
        else:
          alignment += line[:-1]
      alignments.append(alignment)

    file.close()
    pka_print(alignments)

    return first, last, alignments


def writeFile(filename, lines):
    """
    Writes a new file
    """
    file = open(filename, 'w')

    for line in lines:
        file.write( "%s\n" % (line) )
    file.close()


def extractResidueType(labels, restype=None, sort=False):
    """
    extracts all labels, e.g. 'HIS', from labels list (used in protein.compareWithExperiment)
    """
    
    # extracting residue types
    if restype == None or restype == "ALL":
      extracted = labels
    else:
      extracted = []
      wantedLabels = groupList(restype)
      for label in labels:
        if label[:3] in wantedLabels:
          extracted.append(label)

    # sorting according to residue number
    """ not there yet """

    return extracted


def examineNeighbours(file, protein):
    """
    extracts all labels, e.g. 'HIS', from labels list (used in protein.compareWithExperiment)
    """
    import data
    experimental = data.getExperiment(name=protein.name)
    residues = []
    str = "# #%s" % (protein.name)
    file.write( "%s\n" % (str) )
    for resName in ["ASP", "GLU"]:
      for residue in protein.residue_dictionary[resName]:
        if residue.label in experimental:
          residues.append(residue)

    dpKa = None
    for residue in residues:
      dpKa = 0.00
      for atom2, atom3 in protein.COlist:
        center = [residue.x, residue.y, residue.z]
        distance, f_angle, nada = calculateAngleFactor(None, atom2, atom3, center)
        if distance <  6.0 and f_angle > 0.001:
          value = 1.0-(distance-3.0)/(6.0-3.0)
          dpKa += 1.2*min(1.0, value)

      str = "%6.2lf %6.2lf" % (experimental[residue.label]-residue.pKa_mod, dpKa)
      file.write( "%s\n" % (str) )


def int2roman(number):
    numerals = { 1 : "I", 4 : "IV", 5 : "V", 9 : "IX", 10 : "X", 40 : "XL",
        50 : "L", 90 : "XC", 100 : "C", 400 : "CD", 500 : "D", 900 : "CM", 1000 : "M" }
    result = ""
    for value, numeral in sorted(numerals.items(), reverse=True):
        while number >= value:
            result += numeral
            number -= value

    return  result


class Mutator:
    """
      mutator object, contains information for mutator
    """

    def __init__(self, label=None):
        """
        Contructer of determinant object - simple, but helps in creating structure!
        """
        self.label = label

        # alignment mutation
        self.type  = "side-chain" # options: [all/side-chain/back-track]

        # overlap mutation
        # no options so far
        self.iterations = 100

        # Jackal mutation
        self.prm   = 2            # force-field:   [CHARM-AA/AMBER-AA/CHARM-UA/AMBER-UA]
        self.min   = 1            # minimization:  [0/1]
        self.ini   = 3            # initial tries: []
        self.rtm   = 1            # rotamer:       [large/medium/mixed/small]

        # Scwrl mutation
        # no options so far


    def setProperty(self, type=None, prm=None, min=None, ini=None, rtm=None, iterations=None):
        """
        setting the secondary mutator information
        """
        if type != None:
          self.type  = type
        if prm  != None:
          self.prm   = prm
        if min  != None:
          self.min   = min
        if ini  != None:
          self.ini   = ini
        if rtm  != None:
          self.rtm   = rtm
        if iterations  != None:
          self.iterations   = iterations


