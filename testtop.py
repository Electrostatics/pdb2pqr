import string
import sys
import getopt
import os
import time
from src import pdb
from src import utilities
from src import structures
from src import routines
from src import protein
from src.pdb import *
from src.utilities import *
from src.structures import *
from src.definitions import *
from src.forcefield import *
from src.routines import *
from src.protein import *
from src.server import *
from src.hydrogens import *
from src.aconf import *
from src.topology import *
from StringIO import *

import src.topology

TOPOLOGYPATH = "dat/TOPOLOGY.xml"

def testTop(argv):
    """
        Main driver for running program from the command line.
    """

    # Append Numeric/Numpy path to sys.path if the user specified a non-standard location during configuration
    sys.argv=argv
    package_path = PACKAGE_PATH
    if package_path != "":
        sys.path.extend(package_path.split(":"))

    shortOptlist = "n,c"
    longOptlist = ["nterm","cterm"]
    pdblist = []
    atomlist = []
    keyatomlist = []
    keyrecords = {}
    newrecords = {}
    newconformerrecords = {}
    addrecords = {}
    missingatomlist = []
    addedhlist = []
    missingcount = 0
    keyloc = 0
    conformerloc = 0
    conformercount = 0
    recordindex = 1

    try: opts, args = getopt.getopt(sys.argv[1:], shortOptlist, longOptlist)
    except getopt.GetoptError, details:
        sys.stderr.write("GetoptError:  %s\n" % details)

#    if len(args) != 2:
#        sys.stderr.write("Incorrect number (%d) of arguments!\n" % len(args))
        #usage(2)

    options = {"nterm":0, "cterm":0}
 
    outpath = None
    ff = None
    for o,a in opts:
        undashed = o[2:]
        if o in ("-n", "--nterm", "--Nterm", "--NTERM"):
            options["nterm"] = 1                 
        elif o in ("-c", "--cterm", "--Cterm", "--CTERM"):
            options["cterm"] = 1         

    text =  "\n--------------------------\n"
    text += "TESTTOP - A Python file for testing topology\n"
    text += "--------------------------\n"
    sys.stderr.write(text)
            
    path = args[0]
    pathname = string.strip(path, string.split(path,'/')[-1])   # get the directory name
    file = getPDBFile(path)
    oldpdblist, errlist = readPDB(file)

    for record in oldpdblist:
        if not isinstance(record, REMARK):
            pdblist.append(record)
        if isinstance(record, ATOM):
            if record.resSeq == 2: 
                keyrecord = record
                keyrecords[keyrecord.name] = keyrecord
                keyres = record.resName
                keyatomlist.append(record.name)
            atomlist.append(record)
        elif isinstance(record, TER): 
            terrecord = record
        elif isinstance(record, END): 
            endrecord = record

    conformerlistreference = copy.deepcopy(pdblist)
    keyrecordscopy = copy.deepcopy(keyrecords)
    
    if len(pdblist) == 0 and len(errlist) == 0:
        try: os.remove(path)
        except OSError: pass
        raise ValueError, "Unable to find file %s!\n" % path

    if options["nterm"] == 0 and options["cterm"] == 0: 
        refmap = {}
        titrationstatemap = {}
        tautomermap = {}
        conformermap = {}
        atommap = {}
        keyconformermap = {}

        toppath = getDatFile(TOPOLOGYPATH)
        if toppath == "":
            raise ValueError, "Could not find %s!" % TOPOLOGYPATH 
     
        topfile = open(toppath)
        top = Topology(topfile)
        topfile.close()

        # reference map from TOPOLOGY.xml
        for res in top.residues:
            refmap[res.name] = res.reference
            for atom in refmap[res.name].atoms:
                atommap[res.name, atom.name] = atom
            for titrationstate in res.titrationStates:
                titrationstatemap[titrationstate.name] = titrationstate
                for tautomer in titrationstate.tautomers:
                    tautomermap[tautomer.name] = tautomer
                    for conformer in tautomer.conformers:
                        conformermap[conformer.name] = conformer
                        if res.name == keyres: #and conformer.name != res.name:
                            conformercount += 1
                            keyconformermap[conformer.name] = conformer

        # find missing atoms in key residue
        for refatom in refmap[keyres].atoms:    
            if refatom.name not in keyatomlist:
                missingatomlist.append(refatom.name)     
                missingcount += 1 

        # add missing atoms to the key residue
        newrecord = copy.deepcopy(keyrecord)
        keyloc = newrecord.serial
        conformerlocreference = newrecord.serial

        seenmap = {}
        nummissing = len(missingatomlist)

        for missatom in missingatomlist:
            newrecord = copy.deepcopy(keyrecord)
            newrecord.name = missatom
            newrecord.x = atommap[keyres, missatom].x
            newrecord.y = atommap[keyres, missatom].y
            newrecord.z = atommap[keyres, missatom].z
            newrecords[newrecord.name] = newrecord

        myDefinition = Definition()
        myProtein = Protein(pdblist, myDefinition)
        myRoutines = Routines(myProtein, 0)

        myRoutines.setTermini()
        myRoutines.updateBonds()

        myRoutines.findMissingHeavy()
        myRoutines.updateSSbridges()

        myRoutines.debumpProtein()

        myRoutines.addHydrogens()

        pdblist = []
        for residue in myProtein.getResidues():
            if residue.resSeq == 2:
                keyresidue = residue
            for atom in residue.getAtoms():
                pdblist.append(atom)
                if atom.resSeq == 2:
                    keyrecords[atom.name] = atom

        keyrecordscopyforconformer = copy.deepcopy(keyrecords)

        for keyconformer in keyconformermap:
            seenmap = {}
            conformerloc = conformerlocreference
            conformerlist = copy.deepcopy(conformerlistreference)
            conformerkeylist = copy.deepcopy(keyrecordscopyforconformer)

            for remove in keyconformermap[keyconformer].conformerRemoves:
                for removeatom in remove.atoms: 
                    del conformerkeylist[removeatom.name]

            for add in keyconformermap[keyconformer].conformerAdds:

                addedhlist = []
                addrecords = {}
                for addatom in add.atoms: 
                    addrecord = copy.deepcopy(keyrecord)
                    addrecord.name = addatom.name
                    addrecord.x = addatom.x
                    addrecord.y = addatom.y
                    addrecord.z = addatom.z
                    addrecords[addrecord.name] = addrecord
                addatomlist = addrecords.keys()
                addatomlistcopy = copy.deepcopy(addatomlist)         # making a copy of addatomlist
                numadd = len(addatomlist)

                # Coordinates translation
                while len(addatomlist) > 0:
                    coords = []
                    refcoords = []

                    atomname = addatomlist.pop(0)
                    refatomcoords = [addrecords[atomname].x, addrecords[atomname].y, addrecords[atomname].z]

                    bondlist = getNearestBonds(keyres, atomname)

                    for bond in bondlist:
                        if bond in conformerkeylist:
                            atom = conformerkeylist[bond]
                        else: 
                            atom = None

                        if atom == None: continue

                        # Get coordinates, reference coordinates
                    
                        if atom.name in addatomlistcopy:
                            coords.append([atom.x, atom.y, atom.z])
                            refcoords.append([addrecords[atom.name].x, addrecords[atom.name].y, addrecords[atom.name].z])
                        else:
                            coords.append([atom.x, atom.y, atom.z])
                            refcoords.append([atommap[keyres, atom.name].x, atommap[keyres, atom.name].y, atommap[keyres, atom.name].z])

                        if len(coords) == 3: break

                    if len(coords) < 3:
                        try:
                            seenmap[atomname] += 1
                        except KeyError: seenmap[atomname] = 1
                        if seenmap[atomname] > nummissing:
                            text = "Too few atoms present to reconstruct or cap residue %s in structure!\n" % (keyres)
                            text += "This error is generally caused by missing backbone atoms in this protein;\n"
                            text += "you must use an external program to complete gaps in the protein backbone."
                            raise ValueError, text
                        else: addatomlist.append(atomname)

                    else: # Rebuild the atom
                        newcoords = findCoordinates(3, coords, refcoords, refatomcoords)
                        if not newrecords.has_key(atomname):
                            newrecords[atomname] = addrecord
                        newrecords[atomname].x = newcoords[0]
                        newrecords[atomname].y = newcoords[1]
                        newrecords[atomname].z = newcoords[2]
                        newrecords[atomname].chainID = "" 
                        keyresidue.createAtom(atomname, newcoords)
                        conformerkeylist[atomname] = keyresidue.getAtom(atomname)
                        if "H" in atomname:
                            addedhlist.append(atomname)

            print "conformer: %s" % (keyconformer)
        
            conformerkeys = conformerkeylist.keys()
            while len(conformerkeys) > 0:
                atomname = conformerkeys.pop(0)
                if atomname in keyatomlist: pass
                else:
                    conformerlist.insert(conformerloc, conformerkeylist[atomname])
                    conformerloc += 1		

            for record in conformerlist:
                if isinstance(record, ATOM):
                    record.serial = recordindex
                    recordindex += 1


            for record in conformerlist:
                if isinstance(record, Atom):
                    record.resName = keyconformer[:3]   #recordresName
                elif isinstance(record, ATOM):
                    if record.resSeq == 2:
                        record.resName = keyconformer[:3]

            myDefinition = Definition()
            myProtein = Protein(conformerlist, myDefinition)
            myRoutines = Routines(myProtein, 0)

            for residue in myProtein.getResidues():
                if residue.resSeq == 2:
                    for atom in residue.getAtoms():
                        if atom.name in addedhlist: 
                            atom.added = 1
                        elif residue.name == "THR" and atom.name == "HG1":
                            atom.added = 1
                        elif residue.name == "CYS" and atom.name == "HG":
                            atom.added = 1

            if keyresidue.name in ["ASP", "ASN", "THR", "CYS"]:
                myRoutines.setTermini()
                myRoutines.updateBonds()

                myRoutines.updateSSbridges()
                myRoutines.debumpProteinTopology()

            for residue in myProtein.getResidues():
                for atom in residue.getAtoms():
                    for record in conformerlist:
                        if isinstance(record, Atom):
                            if record.resName == atom.resName and record.name == atom.name:
                                if record.x != atom.x:
                                    record.x = atom.x
                                if record.y != atom.y:
                                    record.y = atom.y
                                if record.z != atom.z:
                                    record.z = atom.z

            for record in conformerlist:
                if isinstance(record, Atom):
                    recordname = record.name
                    recordresName = record.resName 
                    recordx = record.x
                    recordy = record.y
                    recordz = record.z
                    recordresSeq = record.resSeq
                    recordserial = record.serial
                    conformerlist.remove(record)
                    newrecord = copy.deepcopy(conformerlist[0])
                    newrecord.name = recordname
                    if keyres == "ARG":
                        newrecord.resName = keyres
                    else:
                        newrecord.resName = keyconformer[:3]   #recordresName
                    newrecord.resSeq = recordresSeq
                    newrecord.x = recordx
                    newrecord.y = recordy
                    newrecord.z = recordz
                    conformerlist.insert(recordserial - 1, newrecord)
                elif isinstance(record, ATOM):
                    if record.resSeq == 2:
                        if keyres == "ARG":
                            record.resName = keyres
                        else:
                            record.resName = keyconformer[:3]

            recordindex = 1
            for record in conformerlist:
                if isinstance(record, ATOM):
                    record.serial = recordindex
                    recordindex += 1

            conformerfile = open(pathname + keyconformer+".pdb", 'w')
            for record in conformerlist:
                if isinstance(record, ATOM):
                    conformerfile.write(str(record)+"\n")
                elif isinstance(record, TER):
                    conformerfile.write("TER\n")
                elif isinstance(record, END):
                    conformerfile.write("END\n")
                else: 
                    pass                 
            conformerfile.close()

            recordindex = 1

    elif options["nterm"] == 1:
        refmap = {}
        titrationstatemap = {}
        tautomermap = {}
        conformermap = {}
        atommap = {}
        keyconformermap = {}
        ntermlocreference = 0
        nummissing = 0

        toppath = getDatFile(TOPOLOGYPATH)
        if toppath == "":
            raise ValueError, "Could not find %s!" % TOPOLOGYPATH 
     
        topfile = open(toppath)
        top = Topology(topfile)
        topfile.close()

        # reference map from TOPOLOGY.xml
        for res in top.residues:
            refmap[res.name] = res.reference
            for atom in refmap[res.name].atoms:
                atommap[res.name, atom.name] = atom
            for titrationstate in res.titrationStates:
                titrationstatemap[titrationstate.name] = titrationstate
                for tautomer in titrationstate.tautomers:
                    tautomermap[tautomer.name] = tautomer
                    for conformer in tautomer.conformers:
                        conformermap[conformer.name] = conformer
                        if res.name == "NTER": #and conformer.name != res.name:
                            conformercount += 1
                            keyconformermap[conformer.name] = conformer

        myDefinition = Definition()
        myProtein = Protein(pdblist, myDefinition)
        myRoutines = Routines(myProtein, 0)

        myRoutines.setTermini()
        myRoutines.updateBonds()

        myRoutines.findMissingHeavy()
        myRoutines.updateSSbridges()

        myRoutines.debumpProtein()

        myRoutines.addHydrogens()

        for residue in myProtein.getResidues():
            for atom in residue.getAtoms():
                if atom.resSeq == 1:
                    atomname = atom.name
                    if atomname in ["H", "H2", "H3"]:
                        residue.removeAtom(atomname)

        # make sure H, H2, H3 removed from N-term
        for residue in myProtein.getResidues():
            for atom in residue.getAtoms():
                if atom.resSeq == 1:
                    atomname = atom.name
                    if atomname in ["H", "H2", "H3"]:
                        residue.removeAtom(atomname)

        for residue in myProtein.getResidues():
            for atom in residue.getAtoms():
                if atom.resSeq == 1: 
                    ntermlocreference += 1

        keyrecords = {}
        pdblist = []
        for residue in myProtein.getResidues():
            if residue.resSeq == 1:
                keyresidue = residue
            for atom in residue.getAtoms():
                pdblist.append(atom)
                if atom.resSeq == 1:
                    keyrecords[atom.name] = atom

        keyrecordscopyforconformer = copy.deepcopy(keyrecords)

        for keyconformer in keyconformermap:
            seenmap = {}
            keyatomlist = [] 
            conformerloc = ntermlocreference
            conformerlist = copy.deepcopy(pdblist)
            conformerkeylist = copy.deepcopy(keyrecordscopyforconformer)
            for k in conformerkeylist:
                keyatomlist.append(k)

            for add in keyconformermap[keyconformer].conformerAdds:

                addedhlist = []
                addrecords = {}
                for addatom in add.atoms: 
                    addrecord = copy.deepcopy(keyrecord)
                    addrecord.name = addatom.name
                    addrecord.x = addatom.x
                    addrecord.y = addatom.y
                    addrecord.z = addatom.z
                    addrecords[addrecord.name] = addrecord
                addatomlist = addrecords.keys()
                if "N+1" in addatomlist:       # no need to add N+1 in the next residue
                    addatomlist.remove("N+1")
                addatomlistcopy = copy.deepcopy(addatomlist)         # making a copy of addatomlist
                numadd = len(addatomlist)

                # Coordinates translation
                while len(addatomlist) > 0:
                    coords = []
                    refcoords = []

                    atomname = addatomlist.pop(0)
                    refatomcoords = [addrecords[atomname].x, addrecords[atomname].y, addrecords[atomname].z]

                    bondlist = ["N", "CA", "HA"]

                    for bond in bondlist:
                        if bond in conformerkeylist:
                            atom = conformerkeylist[bond]
                        else: 
                            atom = None

                        if atom == None: continue

                        # Get coordinates, reference coordinates
                    
                        if atom.name in addatomlistcopy:
                            coords.append([atom.x, atom.y, atom.z])
                            refcoords.append([addrecords[atom.name].x, addrecords[atom.name].y, addrecords[atom.name].z])
                        else:
                            coords.append([atom.x, atom.y, atom.z])
                            refcoords.append([myDefinition.map[keyres].getAtom(atom.name).x, myDefinition.map[keyres].getAtom(atom.name).y, myDefinition.map[keyres].getAtom(atom.name).z])

                        if len(coords) == 3: break

                    if len(coords) < 3:
                        try:
                            seenmap[atomname] += 1
                        except KeyError: seenmap[atomname] = 1
                        if seenmap[atomname] > numadd:
                            text = "Too few atoms present to reconstruct or cap residue %s in structure!\n" % (keyres)
                            text += "This error is generally caused by missing backbone atoms in this protein;\n"
                            text += "you must use an external program to complete gaps in the protein backbone."
                            raise ValueError, text
                        else: addatomlist.append(atomname)

                    else: # Rebuild the atom
                        newcoords = findCoordinates(3, coords, refcoords, refatomcoords)

                        if not newrecords.has_key(atomname):
                            newrecords[atomname] = addrecord
                        newrecords[atomname].x = newcoords[0]
                        newrecords[atomname].y = newcoords[1]
                        newrecords[atomname].z = newcoords[2]
                        newrecords[atomname].chainID = "" 
                        keyresidue.createAtom(atomname, newcoords)
                        conformerkeylist[atomname] = keyresidue.getAtom(atomname)
                        if "H" in atomname:
                            addedhlist.append(atomname)

            print "conformer: %s" % (keyconformer)
        
            conformerkeys = conformerkeylist.keys()
            while len(conformerkeys) > 0:
                atomname = conformerkeys.pop(0)
                if atomname in keyatomlist: pass
                else:
                    conformerlist.insert(conformerloc, conformerkeylist[atomname])
                    conformerloc += 1		

            for record in conformerlist:
                if isinstance(record, ATOM):
                    record.serial = recordindex
                    recordindex += 1

            for record in conformerlist:
                if isinstance(record, Atom):
                    recordname = record.name
                    recordresName = record.resName 
                    recordx = record.x
                    recordy = record.y
                    recordz = record.z
                    recordresSeq = record.resSeq
                    recordserial = record.serial
                    conformerlist.remove(record)
                    newrecord = copy.deepcopy(keyrecord)
                    newrecord.name = recordname
                    newrecord.resSeq = recordresSeq
                    newrecord.x = recordx
                    newrecord.y = recordy
                    newrecord.z = recordz
                    conformerlist.insert(recordserial - 1, newrecord)
            
            conformerlist.append(terrecord)
            conformerlist.append(endrecord)
     
            # Re-serialize atoms
            recordindex = 1
            for record in conformerlist:
                if isinstance(record, ATOM):
                    record.serial = recordindex
                    recordindex += 1


            conformerfile = open(pathname + keyconformer+".pdb", 'w')
            for record in conformerlist:
                if isinstance(record, ATOM):
                    conformerfile.write(str(record)+"\n")
                elif isinstance(record, TER):
                    conformerfile.write("TER\n")
                elif isinstance(record, END):
                    conformerfile.write("END\n")
                else: 
                    pass                 
            conformerfile.close()

            recordindex = 1

    elif options["cterm"] == 1:  
        refmap = {}
        titrationstatemap = {}
        tautomermap = {}
        conformermap = {}
        atommap = {}
        keyconformermap = {}
        ctermlocreference = 0
        nummissing = 0

        toppath = getDatFile(TOPOLOGYPATH)
        if toppath == "":
            raise ValueError, "Could not find %s!" % TOPOLOGYPATH 
     
        topfile = open(toppath)
        top = Topology(topfile)
        topfile.close()

        # reference map from TOPOLOGY.xml
        for res in top.residues:
            refmap[res.name] = res.reference
            for atom in refmap[res.name].atoms:
                atommap[res.name, atom.name] = atom
            for titrationstate in res.titrationStates:
                titrationstatemap[titrationstate.name] = titrationstate
                for tautomer in titrationstate.tautomers:
                    tautomermap[tautomer.name] = tautomer
                    for conformer in tautomer.conformers:
                        conformermap[conformer.name] = conformer
                        if res.name == "CTER": #and conformer.name != res.name:
                            conformercount += 1
                            keyconformermap[conformer.name] = conformer

        myDefinition = Definition()
        myProtein = Protein(pdblist, myDefinition)
        myRoutines = Routines(myProtein, 0)

        myRoutines.setTermini()
        myRoutines.updateBonds()

        myRoutines.findMissingHeavy()
        myRoutines.updateSSbridges()

        myRoutines.debumpProtein()

        myRoutines.addHydrogens()

        for residue in myProtein.getResidues():
            for atom in residue.getAtoms():
                if atom.resSeq == 3:
                    atomname = atom.name
                    if atomname in ["O", "OXT"]:
                        residue.removeAtom(atomname)

        # make sure O, OXT removed from C-term
        for residue in myProtein.getResidues():
            for atom in residue.getAtoms():
                if atom.resSeq == 3:
                    atomname = atom.name
                    if atomname in ["O", "OXT"]:
                        residue.removeAtom(atomname)

        for residue in myProtein.getResidues():
            for atom in residue.getAtoms():
                ctermlocreference += 1

        keyrecords = {}
        pdblist = []
        for residue in myProtein.getResidues():
            if residue.resSeq == 3:
                keyresidue = residue
            for atom in residue.getAtoms():
                pdblist.append(atom)
                if atom.resSeq == 3:
                    keyrecords[atom.name] = atom

        keyrecordscopyforconformer = copy.deepcopy(keyrecords)

        for keyconformer in keyconformermap:
            seenmap = {}
            keyatomlist = [] 
            conformerloc = ctermlocreference
            conformerlist = copy.deepcopy(pdblist)
            conformerkeylist = copy.deepcopy(keyrecordscopyforconformer)
            for k in conformerkeylist:
                keyatomlist.append(k)

            for add in keyconformermap[keyconformer].conformerAdds:
                addedhlist = []
                addrecords = {}
                for addatom in add.atoms: 
                    addrecord = copy.deepcopy(keyrecord)
                    addrecord.name = addatom.name
                    addrecord.x = addatom.x
                    addrecord.y = addatom.y
                    addrecord.z = addatom.z
                    addrecords[addrecord.name] = addrecord
                addatomlist = addrecords.keys()
                if "C-1" in addatomlist:       # no need to add C-1 in the previous residue
                    addatomlist.remove("C-1")
                addatomlistcopy = copy.deepcopy(addatomlist)         # making a copy of addatomlist
                numadd = len(addatomlist)

                # Coordinates translation
                while len(addatomlist) > 0:
                    coords = []
                    refcoords = []

                    atomname = addatomlist.pop(0)
                    refatomcoords = [addrecords[atomname].x, addrecords[atomname].y, addrecords[atomname].z]
                    bondlist = ["C", "CA", "CB"]

                    for bond in bondlist:
                        if bond in conformerkeylist:
                            atom = conformerkeylist[bond]
                        else: 
                            atom = None

                        if atom == None: continue

                        # Get coordinates, reference coordinates
                    
                        if atom.name in addatomlistcopy:
                            coords.append([atom.x, atom.y, atom.z])
                            refcoords.append([addrecords[atom.name].x, addrecords[atom.name].y, addrecords[atom.name].z])
                        else:
                            coords.append([atom.x, atom.y, atom.z])
                            refcoords.append([myDefinition.map[keyres].getAtom(atom.name).x, myDefinition.map[keyres].getAtom(atom.name).y, myDefinition.map[keyres].getAtom(atom.name).z])

                        if len(coords) == 3: break

                    if len(coords) < 3:
                        try:
                            seenmap[atomname] += 1
                        except KeyError: seenmap[atomname] = 1
                        if seenmap[atomname] > numadd:
                            text = "Too few atoms present to reconstruct or cap residue %s in structure!\n" % (keyres)
                            text += "This error is generally caused by missing backbone atoms in this protein;\n"
                            text += "you must use an external program to complete gaps in the protein backbone."
                            raise ValueError, text
                        else: addatomlist.append(atomname)

                    else: # Rebuild the atom
                        newcoords = findCoordinates(3, coords, refcoords, refatomcoords)
                        if not newrecords.has_key(atomname):
                            newrecords[atomname] = addrecord
                        newrecords[atomname].x = newcoords[0]
                        newrecords[atomname].y = newcoords[1]
                        newrecords[atomname].z = newcoords[2]
                        newrecords[atomname].chainID = "" 
                        keyresidue.createAtom(atomname, newcoords)
                        conformerkeylist[atomname] = keyresidue.getAtom(atomname)
                        if "H" in atomname:
                            addedhlist.append(atomname)

            print "conformer: %s" % (keyconformer)
        
            conformerkeys = conformerkeylist.keys()
            while len(conformerkeys) > 0:
                atomname = conformerkeys.pop(0)
                if atomname in keyatomlist: pass
                else:
                    conformerlist.insert(conformerloc, conformerkeylist[atomname])
                    conformerloc += 1		

            for record in conformerlist:
                if isinstance(record, ATOM):
                    record.serial = recordindex
                    recordindex += 1

            for record in conformerlist:
                if isinstance(record, Atom):
                    recordname = record.name
                    recordresName = record.resName 
                    recordx = record.x
                    recordy = record.y
                    recordz = record.z
                    recordresSeq = record.resSeq
                    recordserial = record.serial
                    conformerlist.remove(record)
                    newrecord = copy.deepcopy(keyrecord)
                    newrecord.name = recordname
                    newrecord.resSeq = recordresSeq
                    newrecord.x = recordx
                    newrecord.y = recordy
                    newrecord.z = recordz
                    conformerlist.insert(recordserial - 1, newrecord)
            
            conformerlist.append(terrecord)
            conformerlist.append(endrecord)
     
            # Re-serialize atoms
            recordindex = 1
            for record in conformerlist:
                if isinstance(record, ATOM):
                    record.serial = recordindex
                    recordindex += 1

            conformerfile = open(pathname + keyconformer+".pdb", 'w')
            for record in conformerlist:
                if isinstance(record, ATOM):
                    conformerfile.write(str(record)+"\n")
                elif isinstance(record, TER):
                    conformerfile.write("TER\n")
                elif isinstance(record, END):
                    conformerfile.write("END\n")
                else: 
                    pass                 
            conformerfile.close()

            recordindex = 1

def getNearestBonds(self, atomname):
    """
        Parameters
            number:   The number of bonds to get
        Returns
            bonds:    A list of atomnames that are within three bonds of
                      the atom and present in residue (list)
    """
    bonds = []
    lev2bonds = []

    refmap = {}
    titrationstatemap = {}
    tautomermap = {}
    conformermap = {}
    conformers = {} 
    atommap = {}

    toppath = getDatFile(TOPOLOGYPATH)
    if toppath == "":
        raise ValueError, "Could not find %s!" % TOPOLOGYPATH 
     
    topfile = open(toppath)
    top = Topology(topfile)
    topfile.close()

    # reference map from TOPOLOGY.xml
    for res in top.residues:
        conformers[res.name] = []
        refmap[res.name] = res.reference
        for atom in refmap[res.name].atoms:
            atommap[res.name, atom.name] = atom
        for titrationstate in res.titrationStates:
            titrationstatemap[titrationstate.name] = titrationstate
            for tautomer in titrationstate.tautomers:
                tautomermap[tautomer.name] = tautomer
                for conformer in tautomer.conformers:
                    conformermap[conformer.name] = conformer
                    if conformer.name != res.name:
                        conformers[res.name].append(conformer.name)

    if (self, atomname) in atommap.keys():
        atom = atommap[self, atomname]
    else:    # find nearest bonds in conformers
        for conformer in conformers[self]:  
            for add in conformermap[conformer].conformerAdds: 
                for addatom in add.atoms: 
                    if atomname == addatom.name:
                        atom = addatom
    
    for bondedatom in atom.bonds:
        if bondedatom not in bonds:
            bonds.append(bondedatom)

    # Get bonded atoms 2 bond lengths away

    for bondedatom in atom.bonds:
        for bond2 in atommap[self, bondedatom].bonds:
            if bond2 not in bonds and bond2 != atomname:
                bonds.append(bond2)
                lev2bonds.append(bond2)

    # Get bonded atoms 3 bond lengths away
    for lev2atom in lev2bonds:
        for bond3 in atommap[self, lev2atom].bonds:
            if bond3 not in bonds:
                bonds.append(bond3)
     
    return bonds

if __name__ == "__main__":
    testTop(sys.argv)
