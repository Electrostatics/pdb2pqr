from itertools import ifilter, izip_longest, izip

def isAtomLine(line):
    try:
        return line[:4] == 'ATOM' or line[:6] == 'HETATM'
    except IndexError:
        return False
    
def parsePQRAtomLine(line, has_chain):
    #Parses ATOM line into a more comparable tuple
    #First peel off the element type
    #This will keep us from running into problems with tests
    # that have enough elements so that the serial runs into
    # the record type.
    recordType = line[:6].strip()
    sLine = line[6:].split()
    
    if has_chain:
        strings = (recordType, sLine[1], sLine[2],sLine[3])
        ints = (int(sLine[0]), int(sLine[4]))
        floats = tuple(float(x) for x in sLine[5:])
    else:
        strings = (recordType, sLine[1], sLine[2])
        ints = (int(sLine[0]), int(sLine[3]))
        floats = tuple(float(x) for x in sLine[4:])
    
    
    return strings,ints,floats

def compareParsedAtoms(atom1, atom2):
    return atom1[0:1] == atom2[0:1] and all(abs(x-y)<0.1 for x,y in izip(atom1[2],atom2[2]))

def ComparePQRAction(outputFileName, testFileName, correctFileName, has_chain=False):
    failure = False
    results = []
    with open(testFileName) as testFile:
        with open(correctFileName) as correctFile:
            testAtoms = ifilter(isAtomLine, testFile) 
            correctAtoms = ifilter(isAtomLine, correctFile)
            check_error = False
            correct_total = 0.0
            test_total = 0.0
            for testAtom, correctAtom in izip_longest(testAtoms, correctAtoms, fillvalue=None):
                if testAtom is None or correctAtom is None:
                    results.append('TEST ERROR: Result file is the wrong length!\n')
                    failure = True
                    break
                parsedTest = parsePQRAtomLine(testAtom, has_chain)
                parsedCorrect = parsePQRAtomLine(correctAtom, has_chain)
                
                if not compareParsedAtoms(parsedTest,parsedCorrect):
                    results.append('WARNING: Mismatch ->\n%s%s\n' % (testAtom, correctAtom))
                    check_error = True
                    
                test_total += sum(parsedTest[1]) + sum(parsedTest[2])
                correct_total += sum(parsedCorrect[1]) + sum(parsedCorrect[2])
                
            if check_error and abs(test_total - correct_total) > 20.0:
                results.append('TEST ERROR: Result file does not match target close enough!\n')
                failure = True
                
    with open(outputFileName, 'w') as outputFile:
        for line in results:
            outputFile.write(line)
            
        outputFile.write('FAILURE!' if failure else 'SUCCESS!')
            
    return failure

def getSummaryLines(sourceFile):
    while 'SUMMARY' not in sourceFile.next():
        pass
    
    #Skip header
    sourceFile.next()    
    line = sourceFile.next()
    
    results = []
    
    while '-----------' not in line:
        sLine = line.split()
        
        
        strings = tuple(sLine[:3])
        floats = tuple(float(x) for x in sLine[3:])
        
        results.append((strings,floats))
        
        line = sourceFile.next()
        
    return results
    
def ComparePROPKAAction(outputFileName, testFileName, correctFileName):
    failure = False
    results = []
    with open(testFileName) as testFile:
        with open(correctFileName) as correctFile:
            testPKAs = getSummaryLines(testFile) 
            correctPKAs = getSummaryLines(correctFile)
            correct_total = 0.0
            test_total = 0.0
            for testPKA, correctPKA in izip_longest(testPKAs, correctPKAs, fillvalue=None):
                if testPKA is None or correctPKA is None:
                    results.append('TEST ERROR: Result file is the wrong length!')
                    failure = True
                    break
                
                test_total += sum(testPKA[1])
                correct_total += sum(correctPKA[1])
                
            if abs(test_total - correct_total) > 20.0:
                results.append('TEST ERROR: Result file does not match target close enough!')
                failure = True
            
    with open(outputFileName, 'w') as outputFile:
        for line in results:
            outputFile.write(line + '\n')
            
        outputFile.write('FAILURE!\n' if failure else 'SUCCESS!\n')
            
    return failure

def CompareStringFunc(outputFileName, targetfile, sourcefile, has_chain=None):
    return 'Comparing files ("%s", "%s") -> %s' % (targetfile, sourcefile, outputFileName)






          