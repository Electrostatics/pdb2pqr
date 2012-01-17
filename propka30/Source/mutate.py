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
import math, os, sys, re

import lib
import output
import pdb
#import debug


def optimizeMutationDeterminants(protein, mutation=None, atoms=None, alignment=None, version=None, options=None):
    """
    Finding the most stable set of determinants around a stabilizing central mutation.
    """
    print(" ----- Finding most stable determinant mutation -----")

    # preliminaries with "WT"
    if protein.status['done pka'] == False:
      protein.calculatePKA(version=version, options=options)
    dG_ref = protein.calculateFoldingEnergy(options=options)
    # dictionary to store proteins and energies
    proteins = {"WT": [dG_ref, protein, None]}

    # reading alignment if needed
    if options.mutator.label in ["alignment", "overlap"] and alignment == None:
      print("reading alignment from optimizeMutationDeterminants()")
      alignment = readAlignmentFiles(filenames=options.alignment, mesophile=protein.name, options=options)

    # creating a list of all mutation combinations
    if isinstance(mutation, str):
      chains                  = generateChainList(mutation, protein, atoms)
      dict_mutation           = recastIntoDictionary(mutation, alignment, protein, chains)
      combinatorial_mutations = generateCombinatorialDictionaryMutations(dict_mutation)
    else:
      combinatorial_mutations = generateCombinatorialDictionaryMutations(mutation)

    # finding the optimal mutation set.
    for combinatorial_mutation in combinatorial_mutations:
      newProtein = makeMutatedProtein(protein, mutation=combinatorial_mutation, atoms=atoms, alignment=alignment, options=options)
      newProtein.calculatePKA(version=version, options=options)
      dG_mut = newProtein.calculateFoldingEnergy(options=options)
      mutation_label  = extractMutationLabel(combinatorial_mutation)
      proteins[mutation_label] = [dG_mut, newProtein, combinatorial_mutation]

    # printing out result of determinant combinatorial search
    best_mutation = postProcessCombinatorialMutation(proteins, options=options)

    # returned in string-format
    return  proteins[best_mutation][2]


def optimizeMultipleMutations(protein, mutations=None, atoms=None, alignment=None, version=None, options=None):
    """
    making permutations of base-mutations and getting the optimal combination, e.g. ["2vuj::A:N25R/A:N181D", "2vuj::A:N151D/A:F48R"]
    """
    print(" ----- Finding most stable meta mutation -----")

    # preliminaries with "WT"
    if protein.status['done pka'] == False:
      protein.calculatePKA(version=version, options=options)
    dG_ref = protein.calculateFoldingEnergy(options=options)
    proteins = {"WT": [dG_ref, protein]}

    combinatorial_mutations = generateCombinatorialDictionaryMutations(mutations)

    for combinatorial_mutation in combinatorial_mutations:
      newProtein = makeMutatedProtein(protein, mutation=combinatorial_mutation, atoms=atoms, options=options)
      newProtein.calculatePKA(version=version, options=options)
      dG_mut = newProtein.calculateFoldingEnergy(options=options)
      label  = extractMutationLabel(combinatorial_mutation)
      proteins[ label ] = [dG_mut, newProtein, combinatorial_mutation]

    best_mutation = postProcessCombinatorialMutation(proteins, options=options)

    return  proteins[best_mutation][2]


def postProcessCombinatorialMutation(proteins, options=None):
    """
    post-processing of optimizeMultipleMutation() and optimizeMutationDeterminants()
    """
    best_mutation = "WT"
    dG_ref = proteins["WT"][0]
    print("%8s %6s\n%8.2lf %6.2lf  %s" % ("dG_fold", "ddG_fold", dG_ref, dG_ref-dG_ref, "WT"))
    for label in proteins.keys():
      if label != "WT":
        dG_mut, protein, full_mutation = proteins[label]
        print("%8.2lf %6.2lf  %s" % (dG_mut, dG_mut-dG_ref, label))
        if dG_mut < proteins[best_mutation][0]:
          best_mutation = label

    print(" \nbest mutation:   %s\n" % (best_mutation))
    if best_mutation != "WT":
      # creating best mutated protein (not WT)
      newProtein = proteins[best_mutation][1]
      newProtein.name = "%s_best" % (proteins["WT"][1].name)
      newProtein.writePDB(options=options)
      newProtein.writePKA(options=options)

    return best_mutation


def generateChainList(mutation, protein, atoms=None):
    """
    Finding the most stable set of determinants around a stabilizing central mutation.
    """
    chains = []
    mesophile         = protein.name
    thermophile, rest = splitStringMutationInTwo(mutation)
    mesophile_chains   = sorted( protein.atoms.keys() )
    if thermophile == None or atoms == None:
      thermophile_chains = []
    else:
      thermophile_chains = sorted( atoms[thermophile].keys() )

    for i in range( len(mesophile_chains) ):
      if len(thermophile_chains) > i:
        chains.append( [mesophile_chains[i], thermophile_chains[i]] )
      else:
        chains.append( [mesophile_chains[i], None] )

    return  chains


def makeMutatedProtein(protein, mutation=None, atoms=None, alignment=None, options=None):
    """
    returning a new protein, mutation approach determined from options.mutator
    """
    from protein import Protein
    print("mutator: %s" % (options.mutator.label))
    print("type: %s\n" % (options.mutator.type))

    # checking mutation format
    # making sure string-mutation has chainID in it, and then recasting into dictionary-mutation.
    if   isinstance(mutation, str):
      if alignment == None:
        alignment = readAlignmentFiles(filenames=options.alignment, mesophile=protein.name, options=options)
      pdbcode, rest = splitStringMutationInTwo(mutation)
      if pdbcode == None:
        mutation = remakeStringMutation(mutation, protein.atoms.keys())
      else:
        mutation = remakeStringMutation(mutation, protein.atoms.keys(), atoms[pdbcode].keys())
      dict_mutation = recastIntoDictionary(mutation, alignment=alignment, protein=protein)
    elif isinstance(mutation, dict):
      # this is already in dictionary-format
      dict_mutation = mutation


    pdbfile = "%s.pdb" % (protein.name)
    newfile = "%s_mut.pdb" % (protein.name)
    if   options.mutator.label == "alignment":
        newAtoms = mutateAtomsDictionaryUsingAlignment(protein, mutation=dict_mutation, atoms=atoms, options=options)
        printMutation(dict_mutation)
        name = "%s_mut" % (protein.name)
        newProtein = Protein(atoms=newAtoms, name=name, options=options)
        newProtein.writePDB(filename=newfile, options=options)
        newProtein = Protein(pdbfile=newfile, name=name, options=options)
        #fixProtonationIssues(protein, newProtein)
    elif options.mutator.label in ["overlap"]:
        newAtoms = mutateAtomsDictionaryUsingOverlap(protein, mutation=dict_mutation, atoms=atoms, options=options)
        printMutation(dict_mutation)
        name = "%s_mut" % (protein.name)
        newProtein = Protein(atoms=newAtoms, name=name, options=options)
        newProtein.writePDB(filename=newfile, options=options)
        newProtein = Protein(pdbfile=newfile, name=name, options=options)
        #fixProtonationIssues(protein, newProtein)
    elif options.mutator.label in ["scwrl3", "scwrl4", "scwrl"]:
        seqfile = "x-ray.seq"
        sequence = protein.makeSequence(mutation=dict_mutation, program="scwrl", options=options)
        output.writeScwrlSequenceFile(sequence, filename=seqfile, options=options)
        executeScwrl(pdbfile, newfile, seqfile, options=options)
        newProtein = Protein(pdbfile=newfile, options=options)
        cmd = "rm %s" % (newfile); os.system(cmd)
    elif options.mutator.label in ["jackal", "scap"]:
        scapfile = "%s_scap.list" % protein.name
        mutation_data = splitIntoMutationData(mutation)
        output.writeJackalScapFile(mutationData=mutation_data, filename=scapfile, options=options)
        executeJackal(pdbfile, scapfile, newfile, options=options)
        newProtein = Protein(pdbfile=newfile, options=options)
        cmd = "rm %s %s" % (scapfile, newfile); os.system(cmd)
    elif options.mutator.label in ["whatif", "WhatIf"]:
        mutation_data = splitIntoMutationData(mutation)
        output.writeWhatIfFile(mutationData=mutation_data, pdbfile=pdbfile, newfile=newfile, options=options)
        executeWhatIf(newfile=newfile, options=options)
        newProtein = Protein(pdbfile=newfile, options=options)
        cmd = "rm %s" % (newfile); os.system(cmd)
    else:
        print("don't know how to use mutator \"%s\"" % (options.mutator.label))
        sys.exit(9)

    return  newProtein


def printMutation(mutation):
    """
    print out the mutation in formatted form
    """
    print(" ----- mutation ----- \n", "  residue    rmsd")
    for key in mutation.keys():
      print("  %s %6.2lf" % (mutation[key]['label'], mutation[key]['rmsd']))


def executeJackal(pdbfile, scapfile, newfile=None, options=None):
    """
    excute Jackal
    """
    prm = options.mutator.prm
    min = options.mutator.min
    ini = options.mutator.ini
    rtm = options.mutator.rtm

    cmd  = "scap -prm %d -min %d -ini %d -rtm %d  %s %s" % (prm, min, ini, rtm, pdbfile, scapfile)

    if newfile != None:
      # chaning the name of resulting pdbfile
      for i in range(len(pdbfile)):
        if pdbfile[i] == ".":
          stop = i
          break
      scap_pdbfile = "%s_scap.pdb" % (pdbfile[:stop])
      cmd += "; mv %s %s" % (scap_pdbfile, newfile)

    os.system(cmd)


def executeScwrl(pdbfile, filename, seqfile, options=None):
    """
    excute Scwrl
    """
    if options.mutator.label in ["scwrl4", "scwrl4.0"]:
      cmd = "scwrl4  -i %s -s %s -o %s" % (pdbfile, seqfile, filename)
    else:
      cmd = "scwrl3  -i %s -s %s -o %s" % (pdbfile, seqfile, filename)
    os.system(cmd)


def executeWhatIf(exefile="whatif.sh", newfile=None, options=None):
    """
    excute WhatIf
    """
    cmd  = ""
    cmd += "chmod u+x %s" % (exefile)
    cmd += "; %s" % (exefile)
    os.system(cmd)


def readAlignmentFiles(filenames=None, mesophile=None, options=None):
    """
    reading several alignment files and returning the resulting dictionary
    alignment[thermophile][thermophile, mesophile]
    """
    if filenames == None:
      return  None
    else:
      alignments = {}
      for filename in filenames:
    
        print("reading alignmentfile %s" % (filename))
        file = open(filename, 'r')
        alignment = {}
        sequence  = ""
      
        while True:
          line = file.readline()
          if   line == "":
            break
          elif line == "\n":
            alignment[name] = {"name": name,
                               "chainID": chainID,
                               "resNumb": resNumb,
                               "sequence": sequence}
            if re.search(mesophile, name):
              mesophile_key = name
            sequence = ""
          elif line[:3] == ">P1":
            name = line[5:-1]
          elif line[:10] == "structureX":
            separators = []
            for i in range(len(line)):
              if line[i] == ":":
                separators.append(i)
            resNumb = int(line[separators[1]+1:separators[2]])
            chainID =     line[separators[2]+1:separators[3]]
          else:
            sequence += line[:-1]
      
        file.close()
    
    
        # store it as alignment-pairs in 'alignments'
        for key in alignment.keys():
          if re.search(mesophile, key):
            """ do nothing, this matches the mesophile """
            # print("match,  don't do \"%s\" \"%s\":" % (mesophile, key))
          else:
            # print("no match, please do \"%s\" \"%s\":" % (mesophile, key))
            alignments[key] = {mesophile: alignment[mesophile_key]}
            alignments[key][key] = alignment[key]
    
      return alignments


def setupBackBoneDictionary(target=None, template=None, options=None):
    """
    returning a dictionary for target and template back-bone atoms; used for easy search and comparison
    """
    backbone = {'N':  [],
                'CA': [],
                'C':  [],
                'O':  []}

    for atoms in [target, template]:
      for atom in atoms:
        if atom.name == "CA" or atom.name == "N" or atom.name == "C" or atom.name == "O":
          backbone[atom.name].append(atom)
          
    return  backbone


def getBackTrackTranslation(backbone, options=None):
    """
    get the translation for the 'back-tracking' in the 'alignment-mutation'
    """
    rmsd = 0.00; number_of_atoms = 0; translation = [0.00, 0.00, 0.00]
    for key in backbone.keys():
      number_of_atoms += 1
      atom1, atom2 = backbone[key]
      dX = atom2.x - atom1.x;
      dY = atom2.y - atom1.y;
      dZ = atom2.z - atom1.z;
      translation[0] += dX
      translation[1] += dY
      translation[2] += dZ
      rmsd += dX*dX + dY*dY + dZ*dZ
    rmsd = math.sqrt(rmsd/number_of_atoms)
    for i in range(3):
      translation[i] = translation[i]/number_of_atoms

    return translation, rmsd


def copyAtomsDictionary(atoms=None, options=None):
    """
    creating a copy of the atoms dictionary
    """
    newAtoms = {}
    for chainID in atoms.keys():
      # creating new chain dictionary
      newAtoms[chainID] = {'keys': []}
      for key in atoms[chainID]['keys']:
        newAtoms[chainID]['keys'].append(key)
        newAtoms[chainID][key] = []
        for atom in atoms[chainID][key]:
          newAtoms[chainID][key].append(atom)

    return  newAtoms


def mutateAtomsDictionaryUsingAlignment(protein, atoms=None, mutation=None, options=None):
    """
    creating a new atoms dictionary with swapped residues, from other 'atoms' dictionaries 
    note, mutation should be in dictionary format
    """
    # create a copy of target 'atoms-dictionary'
    newAtoms = copyAtomsDictionary(atoms=protein.atoms, options=options)
    configs  = [ protein.configurations[0] ]

    for key in mutation.keys():

      singleMutateAtomsDictionary(atoms=newAtoms, 
                                  template=atoms, 
                                  mutation=mutation[key], 
                                  type=options.mutator.type, 
                                  configs=configs,
                                  options=options)

    return  newAtoms


def singleMutateAtomsDictionary(atoms=None, template=None, mutation=None, type=None, configs=None, options=None):
    """
    replaces a single 'atoms'-list with that of a template
    """
    # extracting preliminary information from 'mutation'
    code             = mutation['pdb']
    target_label     = mutation['target']
    target_chainID   = target_label[-1]
    template_label   = mutation['template']
    template_chainID = template_label[-1]
    label            = mutation['label']
    resNumb          = int( label[3:7] )
    resName          = label[:3]

    new_atoms = []
    backbone_names = ["N", "CA", "C", "O"]
    # getting sorted references to back-bone atoms
    backbone = setupBackBoneDictionary(target=atoms[target_chainID][target_label], template=template[code][template_chainID][template_label])
    # getting 'back-track' translation based on back-bone atoms
    translation, rmsd = getBackTrackTranslation(backbone, options=options)
    mutation['rmsd'] = rmsd
    # going through back-bone atoms
    for name in backbone_names:
      target_atom, template_atom = backbone[name]
      if   type == "side-chain":
        new_atom = target_atom.makeCopy(resName=resName, configs=configs)
      elif type == "back-track":
        new_atom = target_atom.makeCopy(resName=resName, configs=configs)
        new_atom.translate(translation)
      else:
        new_atom = template_atom.makeCopy(resNumb=resNumb, chainID=target_chainID, configs=configs)
      new_atoms.append(new_atom)
    # going through template
    for atom in template[code][template_chainID][template_label]:
      if atom.name not in backbone_names:
        new_atom = atom.makeCopy(resNumb=resNumb, chainID=target_chainID, configs=configs)
        new_atoms.append(new_atom)

    # updating atoms dictionary
    del atoms[target_chainID][target_label]; atoms[target_chainID][label] = new_atoms
    for i in range( len(atoms[target_chainID]['keys']) ):
      if atoms[target_chainID]['keys'][i] == target_label:
        atoms[target_chainID]['keys'][i] = label
        break

    return


def mutateAtomsDictionary_old(protein, atoms=None, mutation=None, options=None):
    """
    creating an atoms dictionary with swapped residues, note, mutation should be in dictionary format
    """
    print(" ----- mutation -----\n", mutation)
    mesophile = protein.name

    newAtoms = {}
    # creating a new 'atoms dictionary', mostly from the mesophile, to make the mutant protein
    for chainID in sorted(protein.atoms.keys()):
      # creating new chain dictionary
      newAtoms[chainID] = {"keys": []}
      for key in protein.atoms[chainID]["keys"]:
        atomList = []
        if key in mutation.keys():
          # adding residue from thermophile data
          code = mutation[key]['pdb']
          key1 = mutation[key]['template']
          key2 = mutation[key]['label']
          newAtoms[chainID]["keys"].append(key2)
          if options.mutator.type == "back-track":
            dX, dY, dZ = getBackTrackTranslation(protein.atoms[chainID][key], atoms[code][chainID][key1])
          resNumb = int( key[3:7] )
          resName = key2[:3]
          template_chainID = key1[-1]
          for atom in atoms[code][template_chainID][key1]:
            if atom.name in ["N", "CA", "C", "O"]:
              # hmm, this is a back-bone atom
              if   options.mutator.type == "all":
                newAtom = atom.makeCopy(resNumb=resNumb, chainID=chainID)
                atomList.append(newAtom)
              elif options.mutator.type == "back-track":
                x=atom.x+dX; y=atom.y+dY; z=atom.z+dZ
                newAtom = atom.makeCopy(resNumb=resNumb, x=x, y=y, z=z, chainID=chainID)
                atomList.append(newAtom)
            else:
              # copy the thermophile side-chain atom
              newAtom = atom.makeCopy(resNumb=resNumb, chainID=chainID)
              newAtom.checkConfigurations(configurations=protein.configurations)
              atomList.append(newAtom)
          # adding mesophile back-bone atoms if "side-chain"
          if options.mutator.type == "side-chain":
            for atom in protein.atoms[chainID][key]:
              if atom.name in ["N", "CA", "C", "O"]:
                newAtom = atom.makeCopy(resName=resName)
                atomList.append(newAtom)
          newAtoms[chainID][key2] = atomList
        else:
          # non-mutated residue, using mesophile data
          newAtoms[chainID]["keys"].append(key)
          for atom in protein.atoms[chainID][key]:
            atomList.append(atom)
          newAtoms[chainID][key] = atomList

    return newAtoms


def splitIntoMutationData(generic_mutation):
    """
    split 'mutation' into multiple mutation data, e.g. 'A:N25R/A:N181D' -> [[A, N, 25, R], [A, N, 181, D]]
    """
    # making sure mutation is used properly
    mutation_data = []
    if isinstance(generic_mutation, str):
      mutation_list = [generic_mutation]
    elif isinstance(generic_mutation, list):
      mutation_list = generic_mutation
    elif isinstance(generic_mutation, dict):
      for key in generic_mutation.keys():
        chainID1 = generic_mutation[key]['target'][-1]
        code1, resName = lib.convertResidueCode(resName=generic_mutation[key]['target'][:3])
        resNumb = int(generic_mutation[key]['target'][3:7])
        code2, resName = lib.convertResidueCode(resName=generic_mutation[key]['label'][:3])
        mutation_data.append([chainID1, code1, resNumb, code2])

    if len(mutation_data) > 0:
      return  mutation_data
    else:
      for mutation in mutation_list:
        thermophile, mutation_string =  splitStringMutationInTwo(mutation)
        for single_mutation in splitIntoSingleMutations(mutation_string):
          chainID1, chainID2, mutation_label = splitSingleMutation(single_mutation)
          code1, resNumb, code2 = extractSingleMutationInformation(mutation_label)
          mutation_data.append([chainID1, code1, resNumb, code2])
      return  mutation_data


def splitIntoSingleMutations(mutation):
    """
    split multiple mutation into single mutation list, e.g. 'A:N25R/A:N181D' -> ['A:N25R', 'A:N181D']
    """
    single_mutations = []
    start = 0

    for i in range( len(mutation) ):
      if mutation[i] == '/':
        single_mutations.append(mutation[start:i])
        start = i+1
    single_mutations.append(mutation[start:])

    return single_mutations


def splitStringMutationInTwo(mutation):
    """
    split a mutation in 'string-format', i.e. '2vuj::AB:N25R/AB:N181D' into '2vuj' & 'AB:N25R/AB:N181D'
    """
    separator = None
    for i in range( len(mutation) ):
      if mutation[i:i+2] == "::":
        separator = i

    if separator == None:
      return  None, mutation
    else:
      return  mutation[:separator], mutation[separator+2:]


def extractSingleMutationInformation(mutation_label):
    """
    split single mutation label into code1, resNumb, & code2, e.g. 'N25R' -> 'N', 25, 'R'
    """
    if mutation_label[1] == ':':
      print("this should be called through splitSingleMutation() instead")
      sys.exit(8)
    else:
      code1   = mutation_label[:1]
      resNumb = int(mutation_label[1:-1])
      code2   = mutation_label[-1]

    return [code1, resNumb, code2]


def generateCombinatorialMutations(mutations):
    """
    generating combinatorial mutationlist with all possible permutations/combinations from 'mutations'
    """
    base = []
    composite_mutations = []
    maxCombinations = 1
    for mutation in mutations:
      base.append(maxCombinations)
      maxCombinations *= 2
    base.reverse()

    for combination in range(1, maxCombinations):
      rest = combination
      str = ""
      short = ""
      for i in range(len(base)):
        if rest >= base[i]:
          rest -= base[i]
          short = "Y%s" % (short)
          if str == "":
            str += "%s" % (mutations[i])
          else:
            str += "/%s" % (mutations[i])
        else:
          short = "N%s" % (short)

      if False:
        print( "%3d %s:  \"%s\"" % (combination, short, str) )
      composite_mutations.append(str)

    return composite_mutations


def generateCombinatorialStringMutations(string_mutation):
    """
    generating all combinations with one string-mutation
    """
    pdbcode, mutation = splitStringMutationInTwo(string_mutation)
    single_mutations  = splitIntoSingleMutations(mutation)
    combinations      = generateCombinatorialMutations(single_mutations)

    mutation_combinations = []
    for combination in combinations:
      if pdbcode == None:
        mutation_combinations.append( combination )
      else:
        mutation_combinations.append( "%s::%s" % (pdbcode, combination) )

    return  mutation_combinations


def generateCombinatorialListMutations(mutations):
    """
    generating combinatorial mutationlist with all possible permutations/combinations from 'mutations'
    """
    base = []
    composite_mutations = []
    maxCombinations = 1
    for mutation in mutations:
      base.append(maxCombinations)
      maxCombinations *= 2
    base.reverse()

    print(" combinatorial mutation list ")
    for combination in range(1, maxCombinations):
      composite_mutation = []
      rest  = combination
      str   = ""
      short = ""
      for i in range(len(base)):
        if rest >= base[i]:
          rest -= base[i]
          short = "Y%s" % (short)
          composite_mutation.append( mutations[i] )
        else:
          short = "N%s" % (short)

      if False:
        print( "%3d %s:  \"%s\"" % (combination, short, composite_mutation) )
      composite_mutations.append(composite_mutation)

    return composite_mutations


def generateCombinatorialInformation(number_of_mutations):
    """
    generating the number of combinations and the base vector to generate combinatorial mutations.
    """
    base = []
    max_combinations = 1
    for counter in range(number_of_mutations):
      base.append(max_combinations)
      max_combinations *= 2
    base.reverse()

    return  max_combinations, base


def generateCombinatorialDictionaryMutations(mutations=None):
    """
    generating combinatorial mutationlist with all possible permutations/combinations from 'mutations'
    """

    max_combinations, base = generateCombinatorialInformation( len(mutations) )

    if isinstance(mutations, dict):
      keys = sorted(mutations.keys())

    combinatorial_mutations = []

    for combination in range(1, max_combinations):
      mutation = {}
      rest = combination
      str = ""
      short = ""
      for i in range(len(base)):
        if rest >= base[i]:
          rest -= base[i]
          short += "Y"
          if   isinstance(mutations, dict):
            mutation[ keys[i] ] = mutations[ keys[i] ]
          elif isinstance(mutations, list):
            for key in mutations[i].keys():
              mutation[ key ] = mutations[i][ key ]
        else:
          short += "N"

      str = "%3d %s  %s" % (combination, short, extractMutationLabel(mutation))
      print(str)
      combinatorial_mutations.append(mutation)

    return combinatorial_mutations


def recastIntoDictionary(input_mutation=None, alignment=None, protein=None, chains=None):
    """
    recasting the mutation from string or list format into dictionary format; needed for more complicated mutations
    """

    if protein == None:
      mesophile = None
    else:
      mesophile = protein.name

    # making sure mutation is used properly
    if isinstance(input_mutation, str):
      mutation_list = [input_mutation]
    elif isinstance(input_mutation, list):
      mutation_list = input_mutation
    elif isinstance(input_mutation, dict):
      return  input_mutation

    dictionary_mutation = {}
    for mutation in mutation_list:
      thermophile, mutation_string =  splitStringMutationInTwo(mutation)
      for single_mutation in splitIntoSingleMutations(mutation_string):
        chainID1, chainID2, mutation_label = splitSingleMutation(single_mutation)
        code1, resNumb, code2     = extractSingleMutationInformation(mutation_label)
        code, resName_mesophile   = lib.convertResidueCode(code=code1)
        code, resName_thermophile = lib.convertResidueCode(code=code2)
        if alignment == None or mesophile == None:
          # No alignment mutation as far as I can judge, just reformatting
          thermophile = None
          template    = None
          # repeat for all chains
          if chainID1 == None:
            chainIDs = chains
          else:
            chainIDs = [[chainID1, None]]
          for chainID, chainID2 in chainIDs:
            key                      = "%s%4d%2s" % (resName_mesophile, resNumb, chainID)
            new_label                = "%s%4d%2s" % (resName_thermophile, resNumb, chainID)
            dictionary_mutation[key] = {'pdb': thermophile, 'target': key, 'template': template, 'label': new_label}
        else:
          # preparing for alignment or overlap mutation
          # thermophile
          resNumb_mesophile   = alignment[thermophile][mesophile]["resNumb"]
          resNumb_thermophile = alignment[thermophile][thermophile]["resNumb"]
          for i in range(len(alignment[thermophile][mesophile]["sequence"])):
            if alignment[thermophile][mesophile]["sequence"][i] in   "ARNDCQEGHILKMFPSTWYV":
              resNumb_mesophile += 1
            if alignment[thermophile][thermophile]["sequence"][i] in "ARNDCQEGHILKMFPSTWYV":
              resNumb_thermophile += 1
            if resNumb_mesophile == resNumb:
              if chainID1 == None and chainID2 == None:
                chainIDs = chains
              else:
                chainIDs = [[chainID1, chainID2]]
              for chainID1, chainID2 in chainIDs:
                key       = "%s%4d%2s" % (resName_mesophile, resNumb, chainID1)
                template  = "%s%4d%2s" % (resName_thermophile, resNumb_thermophile, chainID2)
                new_label = "%s%4d%2s" % (resName_thermophile, resNumb_mesophile, chainID1)
                dictionary_mutation[key] = {'pdb': thermophile, 'target': key, 'template': template, 'label': new_label}

    return  dictionary_mutation


def splitSingleMutation(single_mutation):
    """
    splits a mutation into chain and mutation components, e.g. 'AB:N25R' -> 'A', 'B', 'N25R'
    """
    separator = None
    for i in range( len(single_mutation) ):
      if single_mutation[i] == ":":
        separator = i

    if   separator == None:
      return  None, None, single_mutation
    elif separator == 1:
      return  single_mutation[0], single_mutation[0], single_mutation[2:]
    elif separator == 2:
      return  single_mutation[0], single_mutation[1], single_mutation[3:]
    else:
      print("could not split mutation \"%s\" correctly in splitSingleMutation()" % (single_mutation))
      sys.exit(8)


def extractMutationLabel(generic_mutation):
    """
    extracts a mutation-label string from a mutation, e.g. ["2vuj::A:N25R/A:S27T/A:N181D"] -> "N25R/S27T/N181D"
    """
    if   isinstance(generic_mutation, str):
      mutations = [generic_mutation]
    elif isinstance(generic_mutation, list):
      mutations =  generic_mutation 
    elif isinstance(generic_mutation, dict):
      """ do nothing """

    label = ""
    if isinstance(generic_mutation, dict):
      for key in generic_mutation.keys():
        code1, resName = lib.convertResidueCode(resName=key[:3])
        code2, resName = lib.convertResidueCode(resName=generic_mutation[key]['label'][:3])
        label += "/%s%d%s" % (code1, int(key[3:7]), code2)
    else:
      for mutation in mutations:
        pdbcode, mutation_string = splitStringMutationInTwo(mutation)
        single_mutations         = splitIntoSingleMutations(mutation_string)
        for single_mutation in single_mutations:
          chainID1, chainID2, gotcha = splitSingleMutation(single_mutation)
          label += "/%s" % (gotcha)

    return  label[1:]


def remakeStringMutation(generic_mutation, keys1=None, keys2=None):
    """
    remaking mutation according to '2vuj::A:N25R/A:N181D', needed since Michal's GUI ommits chainID if
    there is only one chain - don't like it :o(
    """
    if   isinstance(generic_mutation, str):
      mutation  =  generic_mutation
    elif isinstance(generic_mutation, list):
      print("don't know how to deal with list-mutations in checkStringMutation()")
      sys.exit(8)
      mutations =  generic_mutation 
    elif isinstance(generic_mutation, dict):
      print("don't know how to deal with dictionary-mutations in checkStringMutation()")
      sys.exit(8)

    new_string = ""
    pdbcode, mutation_string = splitStringMutationInTwo(mutation)
    single_mutations         = splitIntoSingleMutations(mutation_string)
    for single_mutation in single_mutations:
      chainID1, chainID2, gotcha = splitSingleMutation(single_mutation)
      if chainID1 == None and chainID2 == None:
        # try to interpret chainIDs
        if   keys1 == None and keys2 == None:
          chain_tag = None
        elif keys1 == None and keys2 != None:
          chain_tag = "%s:" % (sorted(keys2)[0])
        elif keys1 != None and keys2 == None:
          chain_tag = "%s:" % (sorted(keys1)[0])
        else:
          chain_tag = "%s%s:" % (sorted(keys1)[0], sorted(keys2)[0])
      else:
        chain_tag = "%s%s:" % (chainID1, chainID2)
      new_string += "/%s%s" % (chain_tag, gotcha)

    if pdbcode == None:
      remade_mutation = "%s" % (new_string[1:])
    else:
      remade_mutation = "%s::%s" % (pdbcode, new_string[1:])

    return  remade_mutation


def makeDictionaryMutation(mutation=None, label=None, pdb=None, template=None):
    """
    updating the dictionary mutation; making sure all elements are there
    """
    # new[label] = {'pdb': pdb, 'template': label2, 'label': new_label}
    if   isinstance(mutation, str):
      print("not coded this conversion in makeDictionaryMutation()")
      sys.exit(8)
    elif isinstance(mutation, list):
      new = {}
      for target_label, template_label in mutation:
        new[target_label]             = {'pdb': pdb}
        new[target_label]['label']    = "%s%s" % (template_label[:3], target_label[3:])
        new[target_label]['template'] = template_label
        new[target_label]['target']   = target_label
      mutation = new
    elif isinstance(mutation, dict):
      for key in mutation.keys():
        mutation[key]['target'] = key
        if "label" not in mutation[key] and "template" in mutation[key]:
          # generating new label
          mutation[key]['label'] = "%s%s" % (mutation[key]['template'][:3], label[3:])

    return  mutation



def calculateRMSDs(position):
    """
    calculating rmsd for new and old position
    """
    rmsd_new = 0.00
    rmsd_old = 0.00
    number_of_points = len(position['target'])
    for key in position['target'].keys():
      for i in range(3):
        rmsd_new += pow( (position['new'][key][i] - position['target'][key][i]), 2)
        rmsd_old += pow( (position['old'][key][i] - position['target'][key][i]), 2)
    rmsd_new = rmsd_new/number_of_points
    rmsd_old = rmsd_old/number_of_points
    rmsd_new = math.sqrt(rmsd_new)
    rmsd_old = math.sqrt(rmsd_old)

    return  rmsd_new, rmsd_old


def makePositionDictionaryFromResidues(original_residue=None, template_residue=None):
    """
    creating a position dictionary for test-rotations and rmsd calculations
    """
    x=0; y=1; z=2
    labels    = []
    translate = {}
    position  = {'target':   {},
                 'new':      {},
                 'old':      {},
                 'template': {}}

    # find 'CA' and set up translation
    ref_to_atom = original_residue.getAtom(name="CA")
    translate['target']   = [-ref_to_atom.x, -ref_to_atom.y, -ref_to_atom.z]
    ref_to_atom = template_residue.getAtom(name="CA")
    translate['template'] = [-ref_to_atom.x, -ref_to_atom.y, -ref_to_atom.z]
    # print(translate)

    # creating position dictionary with translated 'target' and 'template' positions
    for atom in template_residue.atoms:
      if atom.name[0] != "H":
        ref_to_atom = original_residue.getAtom(name=atom.name)
        if ref_to_atom == None:
          """ do nothing """
        else:
          labels.append(atom.name)
          position['target'][atom.name]   = [ref_to_atom.x + translate['target'][x],  
                                             ref_to_atom.y + translate['target'][y],  
                                             ref_to_atom.z + translate['target'][z]]
          position['template'][atom.name] = [atom.x + translate['template'][x],  
                                             atom.y + translate['template'][y],  
                                             atom.z + translate['template'][z]]

    return  labels, translate, position


def getAtomCoordinates(atoms=None, name=None):
    """
    finding and extracting coordinates from 'atoms'-list
    """
    for atom in atoms:
      if name == atom.name:
        return  [atom.x, atom.y, atom.z]


def makePositionDictionaryFromAtoms(target=None, template=None, corresponding_names=None):
    """
    creating a position dictionary for test-rotations and rmsd calculations
    """
    labels    = []
    translate = {'target':   [0.00, 0.00, 0.00],
                 'template': [0.00, 0.00, 0.00]}
    position  = {'target':   {},
                 'new':      {},
                 'old':      {},
                 'template': {}}

    # extracting current atom positions
    for name1, name2 in corresponding_names:
      position['template'][name1] = getAtomCoordinates(atoms=template, name=name1)
      position['target'][name1]   = getAtomCoordinates(atoms=target,   name=name2)

    # setting translation vector
    for key in ['target', 'template']:
      for i in range(3):
        translate[key][i] = -position[key]['CA'][i]

    # translating according to translation vector
    for key in ['target', 'template']:
      for name1, name2 in corresponding_names:
        for i in range(3):
          position[key][name1][i] += translate[key][i]
      
    return  labels, translate, position


def mutateAtomsDictionaryUsingOverlap(protein=None, atoms=None, mutation=None, options=None):
    """
    This routine overlaps two residues based on an array of atom labels, 'center'
    """
    from corresponding_atoms import makeCorrespondingAtomNames
    from rotate import rotatePosition, translatePosition, makeCrossProduct, calculateVectorLength, \
                              makeScalarProduct, generateRandomDisplacement, generateRandomRotation, rotateAtoms, translateAtoms

    # create copy of target
    newAtoms = copyAtomsDictionary(atoms=protein.atoms, options=options)
    # create dictionary with 'corresponding atom names'.
    corresponding_names = makeCorrespondingAtomNames()
    configs  = [ protein.configurations[0] ]
    dR       = [None, None, None]
    center   = ['CA']

    for key in mutation.keys():

      # mutate the new atoms dictionary (position unchanged)
      singleMutateAtomsDictionary(atoms=newAtoms,
                                  template=atoms,
                                  mutation=mutation[key],
                                  type="all",
                                  configs=configs,
                                  options=options)

      # generate position-dictionary
      code             = mutation[key]['pdb']
      target_label     = mutation[key]['target']
      target_chainID   = target_label[-1]
      template_label   = mutation[key]['template']
      template_chainID = template_label[-1]
      label            = mutation[key]['label']
      center, translate, position = makePositionDictionaryFromAtoms(target=protein.atoms[target_chainID][target_label], 
                                                                    template=atoms[code][template_chainID][template_label], 
                                                                    corresponding_names=corresponding_names[ template_label[:3] ][ target_label[:3] ])
      translateAtoms(newAtoms[target_chainID][label], translate['template'])

      # debug printout
      if False:
        for atom in newAtoms[target_chainID][label]:
          print(atom.makePDBLine())

      # -- rotation 1 ---
      cross_product = makeCrossProduct(position['target']['N'], position['template']['N'])
      template_length      = calculateVectorLength(position['template']['N'])
      target_length        = calculateVectorLength(position['target']['N'])
      cross_product_length = calculateVectorLength(cross_product)
      theta1 = math.asin(cross_product_length/(target_length*template_length))
      # print("REMARK: theta1 = %6.3lf" % (theta1*180/math.pi))
      rotatePosition(position['template'], cross_product, -theta1, center=['CA'])
      rotateAtoms(newAtoms[target_chainID][label], cross_product, -theta1, center=['CA'])

      #  ---  rotation 2  ---
      target_cross_product   = makeCrossProduct(position['target']['N'], position['target']['C'])
      template_cross_product = makeCrossProduct(position['template']['N'], position['template']['C'])
      scalar_product  = makeScalarProduct(target_cross_product, template_cross_product)
      target_length   = calculateVectorLength(target_cross_product)
      template_length = calculateVectorLength(template_cross_product)
      if makeScalarProduct(position['template']['C'], makeCrossProduct(position['target']['C'], position['target']['N'])) < 0.00:
        theta2 =  math.acos(scalar_product/(target_length*template_length))
      else:
        theta2 = -math.acos(scalar_product/(target_length*template_length))
      # print("REMARK: theta2 = %6.3lf" % (theta2*180/math.pi))
      rotatePosition(position['template'], position['target']['N'], -theta2, center=['CA'])
      rotateAtoms(newAtoms[target_chainID][label], position['target']['N'], -theta2, center=['CA'])

      # initializing iterations: copy over 'template' to 'new' & 'old'
      chainNumb = ord("A")
      center    = []
      for name in position['template'].keys():
        center.append(name)
        position['new'][name] = [position['template'][name][0], position['template'][name][1], position['template'][name][2]]
        position['old'][name] = [position['template'][name][0], position['template'][name][1], position['template'][name][2]]

      # debug printout
      if False:
        for atom in newAtoms[target_chainID][label]:
          print(atom.makePDBLine())


      #  ---  iterations  ---
      for iter in range(1, options.mutator.iterations+1):

        theta, axis = generateRandomRotation(0.10)
        dR          = generateRandomDisplacement(0.10)

        rotatePosition(position['new'], axis, theta, center=center)
        translatePosition(position['new'], dR)

        rmsd_new, rmsd_old = calculateRMSDs(position)

        if rmsd_new < rmsd_old:
          # print("REMARK   update structure (iter %3d (%8.3lf%8.3lf))" % (iter, rmsd_new, rmsd_old))
          chainNumb += 1
          for name in position['old'].keys():
            for i in range(3):
              position['old'][name][i] = position['new'][name][i]
            if False:
              Rx, Ry, Rz = position['new'][name]
              chainNumb = ord("B")
              print("ATOM    669  %-3s GLU %s  89    %8.3lf%8.3lf%8.3lf" % ("H", chr(chainNumb), Rx, Ry, Rz))
          rotateAtoms(newAtoms[target_chainID][label], axis, theta, center=center)
          translateAtoms(newAtoms[target_chainID][label], dR)
        else:
          for name in position['new'].keys():
            for i in range(3):
              position['new'][name][i] = position['old'][name][i]


      # translate back to 'target', i.e. 'mutation-site'
      dR[0] = -translate['target'][0]; dR[1] = -translate['target'][1]; dR[2] = -translate['target'][2]
      translateAtoms(newAtoms[target_chainID][label], dR)
      mutation[key]['rmsd'] = rmsd_old

      # debug printout
      if False:
        print("REMARK   Final structure")
        for atom in newAtoms[target_chainID][label]:
          print(atom.makePDBLine())
      #sys.exit(8)

    return  newAtoms


def fixProtonationIssues(original_protein=None, mutant_protein=None):
    """
    This routine fixes the problems with protonation with 'alignment' and 'overlap' mutations, i.e. the 
    'new' protonation scheme does not protonate references to atoms.
    """
    for chain in mutant_protein.chains:
      for residue in chain.residues:
        original_residue = original_protein.getResidue(label=residue.label)
        if original_residue == None:
          """ did not find residue, doing nothing """
        else:
          str = "protons for %s:" % (residue.label)
          for original_atom in original_residue.atoms:
            if original_atom.element == 'H':
              mutant_atom = residue.getAtom(name=original_atom.name)
              if mutant_atom == None:
                str += " %s" % (original_atom.name)
                residue.atoms.append(original_atom)

    return


def makeCompositeAtomsDictionary(protein=None, pdbfiles=None, options=None):
    """
    This routine creates a composite 'atoms' dictionary.
    """
    if   options.mutations == None:
      # not making mutations, we don't need 'atoms-dictionary'
      return  None
    elif options.mutator.label in ["alignment", "overlap"]:
      atoms = {protein.name: protein.atoms}
      if options.thermophiles == None:
        # read the pdbcode from mutations and get the pdb atoms
        for mutation in options.mutations:
          if isinstance(mutation, str):
            # keeping string notation of mutation
            pdbcode = mutation[:4]
            if pdbcode not in atoms:
              atoms[pdbcode] = pdb.readPDB(filename=pdbcode)
          else:
            # extracting pdbcode if mutation is in dictionary format
            for key in mutation.keys():
              pdbcode = mutation[key]['pdb']
              if pdbcode not in atoms:
                atoms[pdbcode] = pdb.readPDB(filename=pdbcode)
      else:
        # Match thermophile pdb files to the mutation code
        if len(options.thermophiles) != 0:
          import re
          for pdbname in options.thermophiles:
            for mutation in options.mutations:
              code, rest = splitStringMutationInTwo(mutation)
              print(mutation, code, rest)
              if re.search(code, pdbname):
                atoms[code] = pdb.readPDB(filename=pdbname)
      return  atoms
    else:
      # not using alignment or overlap mutator, we don't need 'atoms dictionary'
      return  None



def makeMutationLight(original_protein=None, template_protein=None, mutation=None, iterations=500):
    """
    This routine overlaps two residues based on an array of atom labels, 'center'
    Note, each rotation & translation has to be perfomed twice: position & residue. This is a choice to make things simpler!
    position  = dictionary with atoms used for the overlap (rmsd)
    translate = dictionary for translating target and template residues to the origin
    center    = array with atom names in 'position-dictionary': used for 'rotation-center'
    """
    from rotate import rotatePosition, translatePosition, makeCrossProduct, calculateVectorLength, makeScalarProduct, generateRandomDisplacement, generateRandomRotation

    # start mutations
    for label in mutation.keys():

      print("%s ==> %s" % (label, mutation[label]['template']))

      # find corresponding residues
      original_residue = original_protein.getResidue(label=label)
      template_residue = template_protein.getResidue(label=mutation[label]['template'])
      copied_residue   = template_residue.makeCopy(chainID=original_residue.chainID, resNumb=original_residue.resNumb)

      # set up position-dictionary for corresponding atoms used for evaluating rmsd
      center, translate, position = makePositionDictionary(original_residue, template_residue)
      copied_residue.translate(translate['template'])

      # debug printout
      if False:
        for key in position['target'].keys():
          Rx, Ry, Rz = position['target'][key]
          print("ATOM    669  %-3s %s A  89    %8.3lf%8.3lf%8.3lf" % (key, original_residue.resName, Rx, Ry, Rz))
        for key in position['template'].keys():
          Rx, Ry, Rz = position['template'][key]
          print("ATOM    669  %-3s GLU B  89    %8.3lf%8.3lf%8.3lf" % ("H", Rx, Ry, Rz))

      #  ---  rotation 1  ---
      cross_product = makeCrossProduct(position['target']['N'], position['template']['N'])
      template_length      = calculateVectorLength(position['template']['N'])
      target_length        = calculateVectorLength(position['target']['N'])
      cross_product_length = calculateVectorLength(cross_product)
      theta1 = math.asin(cross_product_length/(target_length*template_length))
      print("REMARK: theta1 = %6.3lf" % (theta1*180/math.pi))
      rotatePosition(position['template'], cross_product, -theta1, center=['CA'])
      copied_residue.rotate(cross_product, -theta1, center=['CA'])

      # debug printout
      if False:
        for key in position['template'].keys():
          Rx, Ry, Rz = position['template'][key]
          print("ATOM    669  %-3s GLU C  89    %8.3lf%8.3lf%8.3lf" % ("H", Rx, Ry, Rz))

      #  ---  rotation 2  ---
      target_cross_product   = makeCrossProduct(position['target']['N'], position['target']['C'])
      template_cross_product = makeCrossProduct(position['template']['N'], position['template']['C'])
      scalar_product  = makeScalarProduct(target_cross_product, template_cross_product)
      target_length   = calculateVectorLength(target_cross_product)
      template_length = calculateVectorLength(template_cross_product)
      if makeScalarProduct(position['template']['C'], makeCrossProduct(position['target']['C'], position['target']['N'])) > 0.00:
        theta2 =  math.acos(scalar_product/(target_length*template_length))
      else:
        theta2 = -math.acos(scalar_product/(target_length*template_length))
      print("REMARK: theta2 = %6.3lf" % (theta2*180/math.pi))
      rotatePosition(position['template'], position['target']['N'], -theta2, center=['CA'])
      copied_residue.rotate(position['target']['N'], -theta2, center=['CA'])

      # debug printout
      if False:
        for key in position['template'].keys():
          Rx, Ry, Rz = position['template'][key]
          print("ATOM    669  %-3s GLU B  89    %8.3lf%8.3lf%8.3lf" % ("H", Rx, Ry, Rz))

      # initializing iterations: copy over 'template' to 'new' & 'old'
      for key in position['template'].keys():
        position['new'][key] = [position['template'][key][0], position['template'][key][1], position['template'][key][2]]
        position['old'][key] = [position['template'][key][0], position['template'][key][1], position['template'][key][2]]

      chainNumb = ord("A")
      for iter in range(1, iterations+1):

        theta, axis = generateRandomRotation(0.10)
        dR          = generateRandomDisplacement(0.10)

        rotatePosition(position['new'], axis, theta, center=center)
        translatePosition(position['new'], dR)

        rmsd_new, rmsd_old = calculateRMSDs(position)

        if rmsd_new < rmsd_old:
          # print("REMARK   update structure (iter %3d (%8.3lf%8.3lf))" % (iter, rmsd_new, rmsd_old))
          chainNumb += 1
          for key in position['old'].keys():
            for i in range(3):
              position['old'][key][i] = position['new'][key][i]
            if False:
              Rx, Ry, Rz = position['new'][key]
              chainNumb = ord("B")
              print("ATOM    669  %-3s GLU %s  89    %8.3lf%8.3lf%8.3lf" % ("H", chr(chainNumb), Rx, Ry, Rz))
          copied_residue.rotate(axis, theta, center=center)
          copied_residue.translate(dR)
        else:
          # print("no update        (iter %3d (%8.3lf%8.3lf))" % (iter, rmsd_new, rmsd_old))
          for key in position['new'].keys():
            for i in range(3):
              position['new'][key][i] = position['old'][key][i]

      # final translation and printout
      for atom in original_residue.atoms:
        str = atom.makePDBLine()
        print(str)
      dR[0] = -translate['target'][0]; dR[1] = -translate['target'][1]; dR[2] = -translate['target'][2]
      copied_residue.translate(dR)
      for atom in copied_residue.atoms:
        str = atom.makePDBLine()
        print(str)



