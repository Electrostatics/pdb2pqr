"""Utilities for PDB2PQR Suite

This module provides various utilities for the PDB2PQR suite to be imported
into other Python scripts.

Authors:  Todd Dolinsky, Yong Huang
"""
import math
import os
import io
import requests
import logging
from os.path import splitext
from sys import path as sys_path


SMALL = 1.0e-7
DIHEDRAL = 57.2958


_LOGGER = logging.getLogger(__name__)


class ExtraOptions(object):
    pass

def createPropkaOptions(pH, verbose=False, reference='neutral'):
    """
    Create a propka options object for running propka.
    """
    #build propka options
    propkaOpts = ExtraOptions()
    propkaOpts.pH = pH
    propkaOpts.reference = reference
    propkaOpts.chains = None
    propkaOpts.thermophiles = None
    propkaOpts.alignment = None
    propkaOpts.mutations = None
    propkaOpts.verbose = verbose
    propkaOpts.protonation = "old-school"
    propkaOpts.window = (0.0, 14.0, 1.0)
    propkaOpts.grid = (0.0, 14.0, 0.1)
    propkaOpts.mutator = None
    propkaOpts.mutator_options = None
    propkaOpts.display_coupled_residues = None
    propkaOpts.print_iterations = None
    propkaOpts.version_label = "Nov30"

    return propkaOpts


def getPQRBaseFileName(filename):
    root, ext = splitext(filename)
    if ext.lower() == '.pqr':
        return root
    return filename


def sortDictByValue(inputdict):
    """
        Sort a dictionary by its values

        Parameters
            inputdict:  The dictionary to sort (inputdict)
        Returns
            items: The dictionary sorted by value (list)
    """
    items = [(v, k) for k, v in inputdict.items()]
    items.sort()
    items.reverse()
    items = [ k for v, k in items]
    return items


def shortestPath(graph, start, end, path=[]):
    """
        Uses recursion to find the shortest path from one node to
        another in an unweighted graph.  Adapted from
        http://www.python.org/doc/essays/graphs.html .

        Parameters:
            graph: A mapping of the graph to analyze, of the form
                   {0: [1,2], 1:[3,4], ...} . Each key has a list
                   of edges.
            start: The ID of the key to start the analysis from
            end:   The ID of the key to end the analysis
            path:  Optional argument used during the recursive step
                   to keep the current path up to that point

        Returns:
            (variable): Returns a list of the shortest path (list)
                        Returns None if start and end are not
                        connected
    """

    path = path + [start]
    if start == end:
        return path
    if start not in graph:
        return None
    shortest = None
    for node in graph[start]:
        if node not in path:
            newpath = shortestPath(graph, node, end, path)
            if newpath:
                if not shortest or len(newpath) < len(shortest):
                    shortest = newpath
    return shortest


def analyzeConnectivity(map, key):
    """
        Analyze the connectivity of a given map using the key value.

        Parameters
            map:  The map to analyze (dict)
            key:  The key value (variable)
        Returns
            list: A list of connected values to the key (list)
    """
    list = []
    keys = [key]
    while len(keys) > 0:
        key = keys[0]
        if key not in list:
            list.append(key)
            # The following 4 lines are modified by Greg Cipriano as a bug fix
            if key in map:
                for value in map[key]:
                    if value not in list:
                        keys.append(value)

        keys.pop(keys.index(key))

    return list


def getAngle(coords1, coords2, coords3):
        """
            Get the angle between three coordinates

            Parameters
                coords1:  The first coordinate set (atom)
                coords2:  The second (vertex) coordinate set (atom)
                coords3:  The third coordinate set (atom)
            Returns
                angle:  The angle between the atoms (float)
        """
        angle = 0.0
        c1 = subtract(coords3, coords2)
        c2 = subtract(coords1, coords2)
        norm1 = normalize(c1)
        norm2 = normalize(c2)
        dotted = dot(norm1, norm2)
        if dotted > 1.0: # If normalized, this is due to rounding error
            dotted = 1.0
        rad = abs(math.acos(dotted))
        angle = rad*180.0/math.pi
        if angle > 180.0:
            angle = 360.0 - angle
        return angle


#TODO: with changes to --userff and --usernames getFFfile and get_names_file do not need to go on wild
#goose chases to find the files in question.
def getFFfile(name):
    """
        Grab the forcefield file.  May or may not residue in the dat/
        directory.
    """

    if name is None:
        return ''

    path = ""
    dirs = sys_path + ["dat"]
    if name in ["amber", "charmm", "parse", "tyl06", "peoepb", "swanson"]: name = name.upper()

    names = ["dat/%s.DAT" % name]

    names.append("%s.DAT" % name)
    names.append("%s.dat" % name)
    names.append("dat/%s" % name)
    names.append(name)

    for guess in names:
        if os.path.isfile(guess):
            return guess

        for p in dirs:
            testpath = "%s/%s" % (p, guess)
            if os.path.isfile(testpath):
                return testpath

    # If we get here return empty string

    return ""


def get_names_file(name):
    """
        Grab the *.names file that contains the XML mapping.

        Parameters
            name:  The name of the forcefield (string)
        Returns
            path:  The path to the file (string)
    """

    if name is None:
        return ''

    path = ""
    dirs = sys_path + ["dat"]
    if name in ["amber", "charmm", "parse", "tyl06", "peoepb", "swanson"]: name = name.upper()

    names = ["dat/%s.names" % name]
    names.append("%s.names" % name)

    for guess in names:
        if os.path.isfile(guess):
            return guess

        for p in dirs:
            testpath = "%s/%s" % (p, guess)
            if os.path.isfile(testpath):
                return testpath

    # If we get here return empty string

    return ""


# TODO - This is a terrible function that should be replaced.
def getDatFile(name):
    """
        Grab a data file. If the file cannot be found in the
        given directory, try the current system path.

        Parameters
            name:  The name of the file to get (string)
        Returns
            path:  The path to the file (string)
    """
    path = ""

    if os.path.isfile(name):
        path = name

    for p in sys_path:
        testpath = "%s/%s" % (p, name)
        if os.path.isfile(testpath):
            path = testpath

    return path


def getPDBFile(path):
    """
        Obtain a PDB file.  First check the path given on the command
        line - if that file is not available, obtain the file from the
        PDB webserver at http://www.rcsb.org/pdb/ .

        Parameters
            path:  Name of PDB file to obtain (string)

        Returns
            file:  File object containing PDB file (file object)
    """

    import os

    if not os.path.isfile(path):
        URLpath = "https://files.rcsb.org/download/" + path + ".pdb"
        _LOGGER.debug("Fetching PDB from %s", URLpath)
        resp = requests.get(URLpath)
        if resp.status_code != 200:
            errstr = "Got code %d while retrieving %s" % (resp.status_code, URLpath)
            raise IOError(errstr)
        return io.StringIO(resp.text)
    else:
        return open(path, 'rt', encoding="utf-8")


# TODO - this should be replaced by numpy
def distance(coords1, coords2):
    """
        Calculate the distance between two coordinates, as denoted by

            dist = sqrt((x2- x1)^2 + (y2 - y1)^2 + (z2 - z1)^2))

        Parameters
            coords1: Coordinates of form [x,y,z]
            coords2: Coordinates of form [x,y,z]
        Returns
            dist:  Distance between the two coordinates (float)
    """
    dist = 0.0
    list = []

    p = coords2[0] - coords1[0]
    q = coords2[1] - coords1[1]
    r = coords2[2] - coords1[2]
    dist = math.sqrt(p*p + q*q + r*r)

    return dist


def add(coords1, coords2):
    """
        Add one 3-dimensional point to another

        Parameters
            coords1: coordinates of form [x,y,z]
            coords2: coordinates of form [x,y,z]
        Returns
            list:  List of coordinates equal to coords2 + coords1 (list)
    """
    x = coords1[0] + coords2[0]
    y = coords1[1] + coords2[1]
    z = coords1[2] + coords2[2]
    return [x,y,z]

# TODO - replace this with numpy
def subtract(coords1, coords2):
    """
        Subtract one 3-dimensional point from another

        Parameters
            coords1: coordinates of form [x,y,z]
            coords2: coordinates of form [x,y,z]
        Returns
            list:  List of coordinates equal to coords1 - coords2 (list)
    """
    x = coords1[0] - coords2[0]
    y = coords1[1] - coords2[1]
    z = coords1[2] - coords2[2]
    return [x,y,z]

# TODO - replace with numpy
def cross(coords1, coords2):
    """
        Find the cross product of two 3-dimensional points

        Parameters
            coords1: coordinates of form [x,y,z]
            coords2: coordinates of form [x,y,z]
        Returns
            list:  Cross product coords2 and coords1 (list)
    """
    list = []
    x = coords1[1]*coords2[2] -  coords1[2]*coords2[1]
    y = coords1[2]*coords2[0] -  coords1[0]*coords2[2]
    z = coords1[0]*coords2[1] -  coords1[1]*coords2[0]
    list = [x,y,z]
    return list

# TODO - replace with numpy
def dot(coords1, coords2):
    """
        Find the dot product of two 3-dimensional points

        Parameters
            coords1: coordinates of form [x,y,z]
            coords2: coordinates of form [x,y,z]
        Returns
            value:  Dot product coords2 and coords1 (float)
    """
    value = 0.0
    for i in range(3):
        value += coords1[i]*coords2[i]
    return value


# TODO - replace with numpy
def normalize(coords):
    """
        Normalize a set of coordinates

        Parameters
            coords: coordinates of form [x,y,z]
        Returns
            list: normalized coordinates (list)
    """
    list = []
    dist = math.sqrt(pow(coords[0],2) + pow(coords[1],2) + pow(coords[2],2))
    if dist > SMALL:
        a = coords[0]/dist
        b = coords[1]/dist
        c = coords[2]/dist
        list = [a,b,c]
    else:
        list = coords
    return list

def factorial(n):
    """
        Returns the factorial of the given number n
    """
    if n <= 1 : return 1
    return n*factorial(n-1)

def getDihedral(coords1, coords2, coords3, coords4):
    """
        Calculate the angle using the four atoms

        Parameters
            coords1: First of four coordinates of form [x,y,z]
            coords2: Second of four
            coords3: Third of four
            coords4: Fourth of four
        Returns
            value: Size of the angle (float)
    """
    value = 0.0

    list43 = subtract(coords4, coords3)
    list32 = subtract(coords3, coords2)
    list12 = subtract(coords1, coords2)

    A = cross(list12, list32)
    Anorm = normalize(A)
    B = cross(list43, list32)
    Bnorm = normalize(B)

    scal = dot(Anorm, Bnorm)
    if abs(scal + 1.0) < SMALL:
        value = 180.0
    elif abs(scal - 1.0) < SMALL:
        value = 0.0
    else:
        value = DIHEDRAL * math.acos(scal)

    chiral = dot(cross(Anorm, Bnorm),list32)
    if chiral < 0:
        value = value * -1.0
    return value