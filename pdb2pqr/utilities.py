"""Utilities for the PDB2PQR software suite.

.. todo:
    The functions in this module are great examples of why PDB2PQR needs
    :mod:`numpy`.
    More efforts should be made to subsitute with :mod:`numpy` data types and
    functions wherever possible throughout the code base.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
.. codeauthor::  Nathan Baker
"""

import logging
import math

import numpy as np

from .config import CHARGE_ERROR, RADIANS_TO_DEGREES, SMALL_NUMBER

_LOGGER = logging.getLogger(__name__)


def noninteger_charge(charge, error_tol=CHARGE_ERROR) -> str:
    """Test whether a charge is an integer.

    :param float charge:  value to test
    :param float error_tol:  absolute error tolerance
    :returns:  string with descripton of problem or empty string if no problem
    """
    abs_error = abs(charge - round(charge))
    if abs_error > abs(error_tol):
        return (
            f"{charge} deviates by {abs_error} from integral, exceeding error "
            f"tolerance {error_tol}"
        )
    return ""


def sort_dict_by_value(inputdict):
    """Sort a dictionary by its values.

    :param inputdict:  the dictionary to sort
    :type inputdict:  dict
    :return:  list of keys sorted by value
    :rtype:  list
    """
    items = sorted(inputdict.items(), key=lambda x: x[1], reverse=True)
    items = [v for v, k in items]
    return items


def shortest_path(graph, start, end, path=[]):
    """Find the shortest path between two nodes.

    Uses recursion to find the shortest path from one node to another in an
    unweighted graph.
    Adapted from http://www.python.org/doc/essays/graphs.html

    :param graph:  a mapping of the graph to analyze, of the form {0: [1,2],
        1:[3,4], ...} . Each key has a list of edges.
    :type graph:  dict
    :param start:  the ID of the key to start the analysis from
    :type start:  str
    :param end:  the ID of the key to end the analysis
    :type end:  str
    :param path:  optional argument used during the recursive step to keep
        the current path up to that point
    :type path:  list

    :return:  list of the shortest path or ``None`` if start and end are not
        connected
    :rtype:  list
    """
    path = [*path, start]
    if start == end:
        return path
    if start not in graph:
        return None
    shortest = None
    for node in graph[start]:
        if node not in path:
            newpath = shortest_path(graph, node, end, path)
            if newpath and (not shortest or len(newpath) < len(shortest)):
                shortest = newpath
    return shortest


def analyze_connectivity(map_, key):
    """Analyze the connectivity of a given map using the key value.

    :param map:  map to analyze
    :type map:  dict
    :param key:  key value
    :type key:  str
    :return:  list of connected values to the key
    :rtype:  list
    """
    clist = []
    keys = [key]
    while keys:
        key = keys[0]
        if key not in clist:
            clist.append(key)
            if key in map_:
                for value in map_[key]:
                    if value not in clist:
                        keys.append(value)
        keys.pop(keys.index(key))
    return clist


def angle(coords1, coords2, coords3):
    """Get the angle between three coordinates.

    :param coords1:  first coordinate set
    :type coords1:  [float, float, float]
    :param coords2:  second (vertex) coordinate set
    :type coords2:  [float, float, float]
    :param coords3:  third coordinate set
    :type coords3:  [float, float, float]
    :return:  angle between the atoms (in degrees)
    :rtype:  float
    """
    diff32 = np.array(coords3) - np.array(coords2)
    diff12 = np.array(coords1) - np.array(coords2)
    norm1 = normalize(diff32)
    norm2 = normalize(diff12)
    dotted = np.inner(norm1, norm2)
    if dotted > 1.0:  # If normalized, this is due to rounding error
        dotted = 1.0
    elif dotted < -1.0:
        dotted = -1.0
    rad = np.absolute(np.arccos(dotted))
    value = rad * 180.0 / np.pi
    if value > 180.0:
        value = 360.0 - value
    return value


def distance(coords1, coords2):
    """Calculate the distance between two coordinates.

    :param coords1:  coordinates of form [x,y,z]
    :type coords1:  [float, float, float]
    :param coords2:  coordinates of form [x,y,z]
    :type coords2:  [float, float, float]

    :return:  distance between the two coordinates
    :rtype:  float
    """
    coords1 = np.array(coords1)
    coords2 = np.array(coords2)
    return np.linalg.norm(coords1 - coords2)


def add(coords1, coords2):
    """Add one 3-dimensional point to another.

    :param coords1:  coordinates of form [x,y,z]
    :type coords1:  [float, float, float]
    :param coords2:  coordinates of form [x,y,z]
    :type coords2:  [float, float, float]
    :return:  list of coordinates equal to coords2 + coords1
    :rtype:  numpy.ndarray
    """
    return np.array(coords1) + np.array(coords2)


def subtract(coords1, coords2):
    """Suntract one 3-dimensional point from another.

    :param coords1:  coordinates of form [x,y,z]
    :type coords1:  [float, float, float]
    :param coords2:  coordinates of form [x,y,z]
    :type coords2:  [float, float, float]
    :return:  list of coordinates equal to coords2 - coords1
    :rtype:  numpy.ndarray
    """
    return np.array(coords1) - np.array(coords2)


def cross(coords1, coords2):
    """Find the cross-product of one 3-dimensional point with another.

    :param coords1:  coordinates of form [x,y,z]
    :type coords1:  [float, float, float]
    :param coords2:  coordinates of form [x,y,z]
    :type coords2:  [float, float, float]
    :return:  list of coordinates equal to coords2 cross coords1
    :rtype:  numpy.ndarray
    """
    return np.cross(np.array(coords1), np.array(coords2))


def dot(coords1, coords2):
    """Find the dot-product of one 3-dimensional point with another.

    :param coords1:  coordinates of form [x,y,z]
    :type coords1:  [float, float, float]
    :param coords2:  coordinates of form [x,y,z]
    :type coords2:  [float, float, float]
    :return:  list of coordinates equal to the inner product of coords2 with
        coords1
    :rtype:  numpy.ndarray
    """
    return np.inner(np.array(coords1), np.array(coords2))


def normalize(coords):
    """Normalize a set of coordinates to unit vector.

    :param coords:  coordinates of form [x,y,z]
    :type coords:  [float, float, float]
    :return: normalized coordinates
    :rtype:  numpy.ndarray
    """
    return coords / np.linalg.norm(coords)


def factorial(num):
    """Returns the factorial of the given number.

    :param num:  number for which to compute factorial
    :type num:  int
    :return:  factorial of number
    :rtype:  int
    """
    if num <= 1:
        return 1
    return num * factorial(num - 1)


def dihedral(coords1, coords2, coords3, coords4):
    """Calculate the dihedral angle from four atoms' coordinates.

    :param coords1:  one of four coordinates of form [x,y,z]
    :type coords1:  [float, float, float]
    :param coords2:  one of four coordinates of form [x,y,z]
    :type coords2:  [float, float, float]
    :param coords3:  one of four coordinates of form [x,y,z]
    :type coords3:  [float, float, float]
    :param coords4:  one of four coordinates of form [x,y,z]
    :type coords4:  [float, float, float]
    :return:  the angle (in degrees)
    :rtype:  float
    """
    diff43 = np.array(coords4) - np.array(coords3)
    diff32 = np.array(coords3) - np.array(coords2)
    diff12 = np.array(coords1) - np.array(coords2)
    c12_32 = np.cross(diff12, diff32)
    c12_32 = normalize(c12_32)
    c43_32 = np.cross(diff43, diff32)
    c43_32 = normalize(c43_32)
    scal = np.inner(c12_32, c43_32)
    if np.absolute(scal + 1.0) < SMALL_NUMBER:
        value = 180.0
    elif np.absolute(scal - 1.0) < SMALL_NUMBER:
        value = 0.0
    else:
        value = RADIANS_TO_DEGREES * math.acos(scal)
    chiral = np.inner(np.cross(c12_32, c43_32), diff32)
    if chiral < 0:
        value *= -1.0
    return value
