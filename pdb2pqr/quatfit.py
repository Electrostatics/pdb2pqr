"""Quatfit routines for PDB2PQR

This module is used to find the coordinates of a new atom based on a reference
set of coordinates and a definition set of coordinates.

Original Code by David J. Heisterberg, The Ohio Supercomputer Center,
1224 Kinnear Rd., Columbus, OH  43212-1163, (614)292-6036, djh@osc.edu,
djh@ohstpy.bitnet, ohstpy::djh

Translated to C from fitest.f program and interfaced with Xmol program by
Jan Labanowski, jkl@osc.edu, jkl@ohstpy.bitnet, ohstpy::jkl

Authors:  David Heisterberg, Jan Labanowski, Jens Erik Nielsen, Todd Dolinsky
"""
import math
from .utilities import normalize


def find_coordinates(numpoints, refcoords, defcoords, defatomcoords):
    """Driver for the quaternion file.

    Provide the coordinates as inputs and obtain the coordinates for the new atom as output.

    Parameters
        numpoints:     The number of points in each list (int)
        refcoords:     The reference coordinates, a list of lists of form [x,y,z] (list)
        defcoords:     The definition coordinates, a list of lists of form [x,y,z] (list)
        defatomcoords:  The definition coordinates for the atom to be placed in the
                        reference frame (list)
    Returns
        newcoords:  The coordinates of the new atom in the reference frame (list)
    """
    refcenter, fitcenter, rotation = qfit(numpoints, refcoords, defcoords)
    newcoords = qtransform(1, defatomcoords, refcenter, fitcenter, rotation)

    # Only return the first coordinates
    return newcoords[0]


def qtransform(numpoints, defcoords, refcenter, fitcenter, rotation):
    """Transform the set of defcoords using the reference center, the fit center,
    and a rotation matrix.

    Parameters
        numpoints: The number of points in each list (int)
        defcoords: Definition coordinates (list)
        refcenter: The reference center (list)
        defcenter: The definition center (list)
        rotation:  The rotation matrix (list)
    Returns
        newcoords: The coordinates of the new point (list)
    """
    if numpoints == 1:
        defcoords = [defcoords]

    fitcoords = translate(numpoints, defcoords, fitcenter, 1)
    rotated = rotmol(numpoints, fitcoords, rotation)
    newcoords = translate(numpoints, rotated, refcenter, 2)

    return newcoords


def qfit(numpoints, refcoords, defcoords):
    """Method for getting new atom coordinates from sets of reference and definition
    coordinates.

    Parameters
        numpoints: The number of points in each list (int)
        refcoords: List of reference coordinates, with each set a list of form [x,y,z] (list)
        defcoords: List of definition coordinates, with each set a list of form [x,y,z] (list)
    """
    nrot = 30
    refcenter, refcoords = center(numpoints, refcoords)
    defcenter, defcoords = center(numpoints, defcoords)

    _, lrot = qtrfit(numpoints, defcoords, refcoords, nrot)
    rotated = rotmol(numpoints, defcoords, lrot)
    _ = translate(numpoints, rotated, refcenter, 2)

    return refcenter, defcenter, lrot


def qchichange(initcoords, refcoords, angle):
    """Change the chiangle of the reference coordinate using the initcoords and the given angle

    Parameters
        initcoords: Coordinates based on the point and basis atoms (one dimensional list)
        difchi    : The angle to use (float)
        refcoords : The atoms to analyze (list of many coordinates)
    Returns
        newcoords : The new coordinates of the atoms (list of many coords)
    """
    # Initialize
    left = []
    right = []
    for _ in range(3):
        left.append(0.0)
        right.append([0.0, 0.0, 0.0])

    # Convert to radians and normalize
    radangle = math.pi * angle/180.0
    normalized = normalize(initcoords)
    left[0] = normalized[0]
    left[1] = normalized[1]
    left[2] = normalized[2]

    # Construct the rotation matrix
    right[0][0] = math.cos(radangle) + left[0]*left[0] * (1.0 - math.cos(radangle))
    right[1][1] = math.cos(radangle) + left[1]*left[1] * (1.0 - math.cos(radangle))
    right[2][2] = math.cos(radangle) + left[2]*left[2] * (1.0 - math.cos(radangle))
    right[1][0] = left[0]*left[1]*(1.0 - math.cos(radangle)) - left[2] * math.sin(radangle)
    right[2][0] = left[0]*left[2]*(1.0 - math.cos(radangle)) + left[1] * math.sin(radangle)
    right[0][1] = left[1]*left[0]*(1.0 - math.cos(radangle)) + left[2] * math.sin(radangle)
    right[2][1] = left[1]*left[2]*(1.0 - math.cos(radangle)) - left[0] * math.sin(radangle)
    right[0][2] = left[2]*left[0]*(1.0 - math.cos(radangle)) - left[1] * math.sin(radangle)
    right[1][2] = left[2]*left[1]*(1.0 - math.cos(radangle)) + left[0] * math.sin(radangle)

    numpoints = len(refcoords)
    newcoords = rotmol(numpoints, refcoords, right)

    return newcoords


def rotmol(numpoints, coor, lrot):
    """Rotate a molecule

    Parameters
        numpoints:  The number of points in the list (int)
        coor:  The input coordinates (list)
        lrot:  The left rotation matrix (list)
    Returns
        out:  The rotated coordinates out=u * x (list)
    """
    out = []
    for i in range(numpoints):
        out.append([])
        out[i].append(lrot[0][0] * coor[i][0] + lrot[1][0] * coor[i][1] + lrot[2][0] * coor[i][2])
        out[i].append(lrot[0][1] * coor[i][0] + lrot[1][1] * coor[i][1] + lrot[2][1] * coor[i][2])
        out[i].append(lrot[0][2] * coor[i][0] + lrot[1][2] * coor[i][1] + lrot[2][2] * coor[i][2])

    return out


def qtrfit(numpoints, defcoords, refcoords, nrot):
    """Find the quaternion, q, [and left rotation matrix, u] that minimizes
        | qTXq - Y | ^ 2 [|uX - Y| ^ 2]
    This is equivalent to maximizing Re (qTXTqY)
    The left rotation matrix, u, is obtained from q by
        u = qT1q

    Parameters
        numpoints:  The number of points in each list (int)
        defcoords:  List of definition coordinates, with each set a list of form [x,y,z] (list)
        refcoords:  List of fitted coordinates, with each set a list of form [x,y,z] (list)
        nrot:  The maximum number of jacobi sweeps
    Returns
        quat:  The best-fit quaternion
        lrot:  The best-fit left rotation matrix
    """
    xxyx = 0.0
    xxyy = 0.0
    xxyz = 0.0
    xyyx = 0.0
    xyyy = 0.0
    xyyz = 0.0
    xzyx = 0.0
    xzyy = 0.0
    xzyz = 0.0

    quat = []
    cmat = []

    for i in range(numpoints):
        xxyx = xxyx + defcoords[i][0] * refcoords[i][0]
        xxyy = xxyy + defcoords[i][0] * refcoords[i][1]
        xxyz = xxyz + defcoords[i][0] * refcoords[i][2]
        xyyx = xyyx + defcoords[i][1] * refcoords[i][0]
        xyyy = xyyy + defcoords[i][1] * refcoords[i][1]
        xyyz = xyyz + defcoords[i][1] * refcoords[i][2]
        xzyx = xzyx + defcoords[i][2] * refcoords[i][0]
        xzyy = xzyy + defcoords[i][2] * refcoords[i][1]
        xzyz = xzyz + defcoords[i][2] * refcoords[i][2]

    for i in range(4):
        cmat.append([])
        for _ in range(4):
            cmat[i].append(0.0)

    cmat[0][0] = xxyx + xyyy + xzyz
    cmat[0][1] = xzyy - xyyz
    cmat[0][2] = xxyz - xzyx
    cmat[0][3] = xyyx - xxyy

    cmat[1][1] = xxyx - xyyy - xzyz
    cmat[1][2] = xxyy + xyyx
    cmat[1][3] = xzyx + xxyz

    cmat[2][2] = xyyy - xzyz - xxyx
    cmat[2][3] = xyyz + xzyy
    cmat[3][3] = xzyz - xxyx - xyyy

    _, vmat = jacobi(cmat, nrot) # diagonalize c

    for i in range(4):
        quat.append(vmat[i][3])
    lrot = q2mat(quat)
    return quat, lrot


def jacobi(amat, nrot):
    """Jacobi diagonalizer with sorted output, only good for 4x4 matrices

    Parameters
        a:    Matrix to diagonalize (4x4 list)
        nrot: Maximum number of sweeps
    Returns
        dvec:    Eigenvalues
        v:    Eigenvectors
    """
    vmat = []
    dvec = []
    for j in range(4):
        dvec.append(0)
        vmat.append([])
        for i in range(4):
            vmat[j].append(0.0)
        vmat[j][j] = 1.0
        dvec[j] = amat[j][j]

    for l in range(nrot):
        dnorm = 0.0
        onorm = 0.0
        for j in range(4):
            dnorm = dnorm + abs(dvec[j])
            for i in range(j):
                onorm = onorm + abs(amat[i][j])

        if dnorm != 0:
            if onorm/dnorm <= 1e-12:
                the_l = l
                break

        for j in range(1, 4):
            for i in range(j):
                bscl = amat[i][j]
                if abs(bscl) > 0.0:
                    dma = dvec[j] - dvec[i]
                    if abs(dma) + abs(bscl) <= abs(dma):
                        tscl = bscl / dma
                    else:
                        qscl = 0.5 * dma/bscl
                        tscl = 1.0/(abs(qscl) + math.sqrt(1 + qscl*qscl))
                        if qscl < 0:
                            tscl = tscl * -1
                    cscl = 1.0/math.sqrt(tscl*tscl + 1)
                    sscl = tscl*cscl
                    amat[i][j] = 0.0
                    for k in range(i):
                        atemp = cscl * amat[k][i] - sscl * amat[k][j]
                        amat[k][j] = sscl * amat[k][i] + cscl * amat[k][j]
                        amat[k][i] = atemp
                    for k in range(i+1, j):
                        atemp = cscl * amat[i][k] - sscl * amat[k][j]
                        amat[k][j] = sscl * amat[i][k] + cscl * amat[k][j]
                        amat[i][k] = atemp
                    for k in range(j+1, 4):
                        atemp = cscl * amat[i][k] - sscl * amat[j][k]
                        amat[j][k] = sscl * amat[i][k] + cscl * amat[j][k]
                        amat[i][k] = atemp
                    for k in range(4):
                        vtemp = cscl * vmat[k][i] - sscl * vmat[k][j]
                        vmat[k][j] = sscl * vmat[k][i] + cscl * vmat[k][j]
                        vmat[k][i] = vtemp
                    dtemp = cscl * cscl * dvec[i] + sscl * sscl * dvec[j] \
                        - 2.0 * cscl * sscl * bscl
                    dvec[j] = sscl * sscl * dvec[i] + cscl * cscl * dvec[j] \
                        +  2.0 * cscl * sscl * bscl
                    dvec[i] = dtemp

    nrot = the_l
    for j in range(3):
        k = j
        dtemp = dvec[k]
        for i in range(j+1, 4):
            if dvec[i] < dtemp:
                k = i
                dtemp = dvec[k]
        if k > j:
            dvec[k] = dvec[j]
            dvec[j] = dtemp
            for i in range(4):
                dtemp = vmat[i][k]
                vmat[i][k] = vmat[i][j]
                vmat[i][j] = dtemp
    return dvec, vmat


def q2mat(quat):
    """Generate a left rotation matrix from a normalized quaternion

    Parameters
        quat:  The normalized quaternion (list)
    Returns
        urot:  The rotation matrix (2-dimensional list)
    """
    urot = []
    for i in range(3):
        urot.append([])
        for _ in range(3):
            urot[i].append(0.0)

    urot[0][0] = quat[0] * quat[0] + quat[1] * quat[1] - quat[2] * quat[2] - quat[3] * quat[3]
    urot[0][1] = 2.0 * (quat[1] * quat[2] - quat[0] * quat[3])
    urot[0][2] = 2.0 * (quat[1] * quat[3] + quat[0] * quat[2])

    urot[1][0] = 2.0 * (quat[2] * quat[1] + quat[0] * quat[3])
    urot[1][1] = quat[0] * quat[0] - quat[1] * quat[1] + quat[2] * quat[2] - quat[3] * quat[3]
    urot[1][2] = 2.0 * (quat[2] * quat[3] - quat[0] * quat[1])

    urot[2][0] = 2.0 *(quat[3] * quat[1] - quat[0] * quat[2])
    urot[2][1] = 2.0 * (quat[3] * quat[2] + quat[0] * quat[1])
    urot[2][2] = quat[0] * quat[0] - quat[1] * quat[1] - quat[2] * quat[2] + quat[3] * quat[3]

    return urot


def center(numpoints, refcoords):
    """Center a molecule using equally weighted points

    Parameters
        numpoints:  Number of points
        refcoords:  List of reference coordinates, with each set
                    a list of form [x,y,z] (list)
    Returns
        refcenter: Center of the set of points (list)
        relcoords: Moved refcoords relative to refcenter (list)
    """
    refcenter = []
    relcoords = []

    for i in range(3):
        refcenter.append(0.0)

    for i in range(numpoints):
        refcenter[0] += refcoords[i][0]
        refcenter[1] += refcoords[i][1]
        refcenter[2] += refcoords[i][2]

    for i in range(3):
        refcenter[i] = refcenter[i] / numpoints

    for i in range(numpoints):
        relcoords.append([])
        relcoords[i].append(refcoords[i][0] - refcenter[0])
        relcoords[i].append(refcoords[i][1] - refcenter[1])
        relcoords[i].append(refcoords[i][2] - refcenter[2])

    return refcenter, relcoords


def translate(numpoints, refcoords, center_, mode):
    """Translate a molecule using equally weighted points

    Parameters
        numpoints:  Number of points
        refcoords:  List of reference coordinates, with each set a list of form [x,y,z] (list)
        center:  Center of the system(list)
        mode:  If 1, center will be subtracted from refcoords.
               If 2, center will be added to refcoords
    Returns
        relcoords: Moved refcoords relative to refcenter (list)
    """
    relcoords = []

    if mode == 1:
        modif = -1
    elif mode == 2:
        modif = 1

    for i in range(numpoints):
        relcoords.append([])
        relcoords[i].append(refcoords[i][0] + modif * center_[0])
        relcoords[i].append(refcoords[i][1] + modif * center_[1])
        relcoords[i].append(refcoords[i][2] + modif * center_[2])

    return relcoords
