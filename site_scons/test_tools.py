#from itertools import ifilter, izip_longest, izip
try:
    from itertools import ifilter as filter
except ImportError:
    pass

try:
    from itertools import izip_longest as zip_longest
except ImportError:
    from itertools import zip_longest

try:
    from itertools import izip as zip
except ImportError:
    pass


from glob import glob
from os import path
import csv

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
    #return atom1[0:1] == atom2[0:1] and all(abs(x-y)<0.1 for x,y in izip(atom1[2],atom2[2]))
    return atom1[0:1] == atom2[0:1] and all(abs(x-y)<0.1 for x,y in zip(atom1[2],atom2[2]))

def ComparePQRAction(outputFileName, testFileName, correctFileName, has_chain=False):
    failure = False
    results = []
    with open(testFileName) as testFile:
        with open(correctFileName) as correctFile:
            #testAtoms = ifilter(isAtomLine, testFile)
            #correctAtoms = ifilter(isAtomLine, correctFile)
            testAtoms = filter(isAtomLine, testFile)
            correctAtoms = filter(isAtomLine, correctFile)
            check_error = False
            correct_total = 0.0
            test_total = 0.0
            #for testAtom, correctAtom in izip_longest(testAtoms, correctAtoms, fillvalue=None):
            for testAtom, correctAtom in zip_longest(testAtoms, correctAtoms, fillvalue=None):
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


def get_csv_data(in_file):
    reader = csv.reader(in_file)
    results = [(float(pH),float(titr)) for pH, titr in reader]
    return results

def get_curve_data(input_path):
    scatter_data = {}
    for file_name in glob(input_path+'/titration_curves/*.csv'):
        base_name = path.basename(file_name)

        name = base_name.rsplit('.', 1)[0]

        with open(file_name, 'rb') as in_file:
            file_data = get_csv_data(in_file)

        scatter_data[name] = file_data

    return scatter_data

def merge_curves(curve1, curve2):
    combined = dict((ph, [value, None]) for (ph, value) in curve1)

    for ph, value in curve2:
        if ph not in combined:
            combined[ph]=[None, value]
        else:
            combined[ph][1] = value

    keys = combined.keys()
    keys.sort()

    combined_list = [(combined[ph][0], combined[ph][1]) for ph in keys]
    return combined_list

def CompareTitCurvesAction(outputFileName, testDirName, correctDirName):
    results = []
    EPSILON = 0.01

    test_data = get_curve_data(testDirName)
    correct_data = get_curve_data(correctDirName)

    for name in correct_data:
        correct_curve = correct_data[name]
        test_curve = test_data.pop(name, None)

        if test_curve is None:
            results.append("ERROR: test results missing curve for residue " + name + "\n")
            continue

        combined = merge_curves(correct_curve, test_curve)

        report_extra_data = False
        report_missing_data = False
        total_error = 0.0
        bad_point_count = 0

        for correct_value, test_value in combined:
            if correct_value is None:
                report_extra_data = True
            elif test_value is None:
                report_missing_data = True
            else:
                diff = abs(correct_value - test_value)
                if diff > EPSILON:
                    bad_point_count += 1
                    total_error += diff


        if report_extra_data:
            results.append("ERROR: test curve " + name +" has extra data\n")

        if report_missing_data:
            results.append("ERROR: test curve " + name +" has missing data\n")

        if bad_point_count > 0:
            results.append("ERROR: test curve {name} has {count}"
                           " bad points with {total} cumulative error.\n".format(name=name, count=bad_point_count, total=total_error))


    for name in test_data:
        results.append("ERROR: extra curve for residue" + name + " in test results\n")

    failure = bool(results)

    with open(outputFileName, 'w') as outputFile:
        for line in results:
            outputFile.write(line)

        outputFile.write('FAILURE!' if failure else 'SUCCESS!')

    return failure



def getSummaryLines(sourceFile):
    #while 'SUMMARY' not in sourceFile.next():
    while 'SUMMARY' not in next(sourceFile):
        pass

    #Skip header
    _ = next(sourceFile)
    line = next(sourceFile)

    results = []

    while '-----------' not in line:
        sLine = line.split()


        strings = tuple(sLine[:3])
        floats = tuple(float(x) for x in sLine[3:])

        results.append((strings,floats))

        line = next(sourceFile)

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
            #for testPKA, correctPKA in izip_longest(testPKAs, correctPKAs, fillvalue=None):
            for testPKA, correctPKA in zip_longest(testPKAs, correctPKAs, fillvalue=None):
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

def CompareDirectoryFunc(outputFileName, targetdir, sourcedir):
    return 'Comparing directories ("%s", "%s") -> %s' % (targetdir, sourcedir, outputFileName)
