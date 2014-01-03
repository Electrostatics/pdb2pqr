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
import math, random, string

from lib import pka_print


def InterAtomDistance(atom1, atom2):
    """
    calculates distance between atom1 and atom2
    """
    dX = atom1.x - atom2.x
    dY = atom1.y - atom2.y
    dZ = atom1.z - atom2.z

    return  math.sqrt( dX*dX + dY*dY + dZ*dZ )


def InterResidueDistance(residue1, residue2):
    """
    calculates distance between atom1 and atom2
    """
    dX = residue1.x - residue2.x
    dY = residue1.y - residue2.y
    dZ = residue1.z - residue2.z

    return  math.sqrt( dX*dX + dY*dY + dZ*dZ )


def AngleFactorX(atom1=None, atom2=None, atom3=None, center=None):
    """
    Calculates the distance and angle-factor from three atoms for back-bone interactions,
    IMPORTANT: you need to use atom1 to be the e.g. ASP atom if distance is reset at return: [O1 -- H2-N3]
    Also generalized to be able to be used for residue 'centers' for C=O COO interactions.
    """
    dX_32 = atom2.x - atom3.x
    dY_32 = atom2.y - atom3.y
    dZ_32 = atom2.z - atom3.z

    distance_23 = math.sqrt( dX_32*dX_32 + dY_32*dY_32 + dZ_32*dZ_32 )

    dX_32 = dX_32/distance_23
    dY_32 = dY_32/distance_23
    dZ_32 = dZ_32/distance_23

    if atom1 == None:
      dX_21 = center[0] - atom2.x
      dY_21 = center[1] - atom2.y
      dZ_21 = center[2] - atom2.z
    else:
      dX_21 = atom1.x - atom2.x
      dY_21 = atom1.y - atom2.y
      dZ_21 = atom1.z - atom2.z

    distance_12 = math.sqrt( dX_21*dX_21 + dY_21*dY_21 + dZ_21*dZ_21 )

    dX_21 = dX_21/distance_12
    dY_21 = dY_21/distance_12
    dZ_21 = dZ_21/distance_12

    f_angle = dX_21*dX_32 + dY_21*dY_32 + dZ_21*dZ_32

    return distance_12, f_angle, distance_23



def linearCoulombEnergy(distance, weight, version, options=None):
    """
    calculates the Coulomb interaction pKa shift
    """
    DIS1  = version.coulomb_cutoff[0]
    DIS2  = version.coulomb_cutoff[1]
    value = 1.0-(distance-DIS1)/(DIS2-DIS1)
    value = min(1.0, value)
    value = max(0.0, value)
    dpka  = version.coulomb_maxpka * value * weight

    return abs(dpka)


def CoulombEnergy(distance, weight, version, options=None):
    """
    calculates the Coulomb interaction pKa shift based on Coulombs law
    eps = 60.0 for the moment; to be scaled with 'weight'
    """
    # setting the dielectric constant
    if   isinstance(version.coulomb_diel, list):
      diel = version.coulomb_diel[1] - (version.coulomb_diel[1]-version.coulomb_diel[0])*weight
    elif isinstance(version.coulomb_diel, float):
      diel = version.coulomb_diel
    elif isinstance(version.coulomb_diel, int):
      diel = version.coulomb_diel

    R = max(distance, version.coulomb_cutoff[0])
    dpka  =244.12/(diel*R) - 244.12/(diel*version.coulomb_cutoff[1])
    if version.coulomb_scaled == True:
      dpka = dpka*weight

    return abs(dpka)


def distanceScaledCoulombEnergy(distance, weight, version, options=None):
    """
    calculates the Coulomb interaction pKa shift based on Coulombs law
    eps = 60.0 for the moment; to be scaled with 'weight'
    """
    # setting the dielectric constant
    if   isinstance(version.coulomb_diel, list):
      diel = version.coulomb_diel[1] - (version.coulomb_diel[1]-version.coulomb_diel[0])*weight
    elif isinstance(version.coulomb_diel, float):
      diel = version.coulomb_diel
    elif isinstance(version.coulomb_diel, int):
      diel = version.coulomb_diel

    # making sure short contacts doesn't blow up
    R = max(distance, version.coulomb_cutoff[0])
    # making sure that the Coulomb dies off at cutoff[1] in a nice way.
    scale = ( R - version.coulomb_cutoff[1] ) / ( version.coulomb_cutoff[0] - version.coulomb_cutoff[1] )
    scale = max(0.0, scale)
    scale = min(1.0, scale)
    dpka  = 244.12/(diel*R) * scale

    return abs(dpka)


def MixedCoulombEnergy(distance, weight, version, options=None):
    """
    calculates the Coulomb interaction pKa shift based on Coulombs law
    eps = 60.0 for the moment; to be scaled with 'weight'
    """
    R = max(distance, version.coulomb_cutoff[0])
    eps       = (1 + 39*(1 - math.exp(-0.18*R)))
    dpka_sur  = 244.12/(80.*R) - 244.12/(80.*version.coulomb_cutoff[1])
    dpka_bur  = 244.12/(eps*R) - 244.12/(eps*version.coulomb_cutoff[1])
    dpka      = weight*dpka_bur + (1.0-weight)*dpka_sur

    return abs(dpka)


def HydrogenBondEnergy(distance, dpka_max, cutoff, f_angle=1.0):
    """
    returns a hydrogen-bond interaction pKa shift
    """
    if   distance < cutoff[0]:
      value = 1.00
    elif distance > cutoff[1]:
      value = 0.00
    else:
      value = 1.0-(distance-cutoff[0])/(cutoff[1]-cutoff[0])

    dpKa  = dpka_max*value*f_angle

    return abs(dpKa)


def buriedRatio(Nmass):
    """
    returns the buried ratio given Nmass
    """
    Nmin =  300.0
    Nmax =  600.0
    buried_ratio = (float(Nmass) - Nmin)/(Nmax - Nmin)
    buried_ratio = max(0.00, buried_ratio)
    buried_ratio = min(1.00, buried_ratio)

    return buried_ratio



def radialVolumeDesolvation(residue, atoms, version, options=None):
    """
    calculates the desolvation according to the ScaledRadialVolumeModel
    """
    if residue.label == "BKB  50 A":
      pka_print("found %s [%6.3lf%6.3lf%6.3lf]!" % (residue.label, residue.x, residue.y, residue.z))
      pka_print("buried_cutoff_sqr = %s!" % (version.buried_cutoff_sqr))
      pka_print("desolv_cutoff_sqr = %s!" % (version.desolv_cutoff_sqr))
    scale_factor = 0.8527*1.36  # temporary weight for printing out contributions
    residue.Nlocl = 0
    residue.Nmass = 0
    residue.Elocl = 0.00
    dV            = 0.00
    volume        = 0.00
    min_distance_4th = pow(2.75, 4)
    for chainID in atoms.keys():
      for key in atoms[chainID]["keys"]:
        for atom in atoms[chainID][key]:
          if atom.element != "H":
            if atom.resNumb != residue.resNumb or atom.chainID != residue.chainID:
                # selecting atom type
                if   atom.name in ["C", "CA"]:
                  atomtype = "C"
                elif atom.name in ["N", "NE1", "NE2", "ND1", "ND2", "NZ", "NE", "NH1", "NH2"]:
                  atomtype = "N"
                elif atom.name in ["O", "OD1", "OD2", "OE1", "OE2", "OH", "OG", "OG1", "OXT"]:
                  atomtype = "O"
                elif atom.name in ["S", "SD", "SG"]:
                  atomtype = "S"
                else:
                  atomtype = "C4"
                dV = version.desolvationVolume[atomtype]
                # calculating distance (atom - residue)
                dX = atom.x - residue.x
                dY = atom.y - residue.y
                dZ = atom.z - residue.z
                distance_sqr = dX*dX + dY*dY + dZ*dZ
                if  distance_sqr < version.desolv_cutoff_sqr:
                  dV_inc  = dV/max(min_distance_4th, distance_sqr*distance_sqr)
                  volume += dV_inc
                  if residue.label in ["ASP   8 a", "ASP  10 a", "GLU 172 a", "ASP  92 a", "GLU  66 a"]:
                    # test printout
                    distance = max(2.75, math.sqrt(distance_sqr))
                    if distance < 20.0:
                      str  = "%6.2lf %8.4lf" % (distance, residue.Q * version.desolvationPrefactor * max(0.00, dV_inc)*scale_factor)
                      str += " %s" % (atomtype)
                      #str += " %s" % (residue.label)
                      pka_print(str)
                if distance_sqr < version.buried_cutoff_sqr:
                  residue.Nmass += 1
                  residue.Vmass += dV
    weight = version.calculateWeight(residue.Nmass)
    scale_factor = 1.0 - (1.0 - version.desolvationSurfaceScalingFactor)*(1.0 - weight)
    residue.buried = weight
    residue.Emass = residue.Q * version.desolvationPrefactor * max(0.00, volume-version.desolvationAllowance) * scale_factor

    return 0.00, 0.00, 0.00, 0.00



def contactDesolvation(residue, atoms, version, options=None):
    """
    calculates the desolvation according to the Contact Model, the old default
    """
    if residue.resName in version.desolvationRadii:
      local_cutoff = version.desolvationRadii[residue.resName]
    else:
      local_cutoff = 0.00
    residue.Nmass = 0
    residue.Nlocl = 0
    for chainID in atoms.keys():
      for key in atoms[chainID]["keys"]:
        for atom in atoms[chainID][key]:
          if atom.element != "H":
            if atom.resNumb != residue.resNumb or atom.chainID != residue.chainID:
                dX = atom.x - residue.x
                dY = atom.y - residue.y
                dZ = atom.z - residue.z
                distance = math.sqrt(dX*dX + dY*dY + dZ*dZ)
                if distance < local_cutoff:
                    residue.Nlocl += 1
                if distance < version.buried_cutoff:
                    residue.Nmass += 1
    if residue.Nmass > 400:
        residue.location = "BURIED "
    else:
        residue.location = "SURFACE"
    residue.Emass = residue.Q * version.desolvationPrefactor * max(0.00, residue.Nmass-version.desolvationAllowance)
    residue.Elocl = residue.Q * version.desolvationLocal * residue.Nlocl
    # Buried ratio - new feature in propka3.0
    # Note, there will be an unforseen problem: e.g. if one residue has Nmass > Nmax and
    # the other Nmass < Nmax, the Npair will not be Nmass1 + Nmass2!
    residue.buried = version.calculateWeight(residue.Nmass)

    return 0.00, 0.00, 0.00, 0.00


def originalDesolvation(residue=None, atoms=None, version=None, options=None):
    """
    calculates the desolvation according to the Contact Model, the old default
    """
    if residue.resName in version.desolvationRadii:
      local_cutoff = version.desolvationRadii[residue.resName]
    else:
      local_cutoff = 0.00
    Nlocl_his4   = 0
    Nlocl_his6   = 0
    residue.Nmass = 0
    residue.Nlocl = 0
    for chainID in atoms.keys():
      for key in atoms[chainID]["keys"]:
        for atom in atoms[chainID][key]:
          HYDROGEN_ATOM = ( (atom.name[0] == 'H') or (atom.name[0] in string.digits and atom.name[1] == 'H') )
          if HYDROGEN_ATOM == False:
            if atom.resNumb != residue.resNumb or atom.chainID != residue.chainID:
                dX = atom.x - residue.x
                dY = atom.y - residue.y
                dZ = atom.z - residue.z
                distance = math.sqrt(dX*dX + dY*dY + dZ*dZ)
                if residue.resName == "HIS":
                    # special case for HIS
                    if distance <  4.0:
                        Nlocl_his4 += 1
                    if distance <  6.0:
                        Nlocl_his6 += 1
                else:
                    # everything else
                    if distance < local_cutoff:
                        residue.Nlocl += 1
                if distance < version.buried_cutoff:
                    residue.Nmass += 1
    if residue.Nmass > 400:
        residue.location = "BURIED "
    else:
        residue.location = "SURFACE"
    if residue.resName == "HIS":
        if residue.location == "SURFACE":
            residue.Nlocl = Nlocl_his4
        else:
            residue.Nlocl = Nlocl_his6
    residue.Emass = residue.Q  * version.desolvationPrefactor * max(0.00, residue.Nmass-version.desolvationAllowance)
    residue.Elocl = residue.Q * version.desolvationLocal * residue.Nlocl
    # Buried ratio - new feature in propka3.0
    # Note, there will be an unforseen problem: e.g. if one residue has Nmass > Nmax and
    # the other Nmass < Nmax, the Npair will not be Nmass1 + Nmass2!
    residue.buried = version.calculateWeight(residue.Nmass)

    return 0.00, 0.00, 0.00, 0.00


def BackBoneReorganization(protein):
    """
    adding test stuff
    """
    residues = []
    for resName in ["ASP", "GLU"]:
      for residue in protein.residue_dictionary[resName]:
        residues.append(residue)

    for residue in residues:
      weight = residue.buried
      dpKa = 0.00
      for atom3, atom2 in protein.COlist:
        center = [residue.x, residue.y, residue.z]
        distance, f_angle, nada = AngleFactorX(atom2=atom2, atom3=atom3, center=center)
        if distance <  6.0 and f_angle > 0.001:
          value = 1.0-(distance-3.0)/(6.0-3.0)
          dpKa += 0.80*min(1.0, value)

      residue.Elocl = dpKa*weight


def TmProfile(protein, reference="neutral", grid=[0., 14., 0.1], Tm=None, Tms=None, ref=None, options=None):
    """
    Calculates the folding profile
    """
    Nres = 0
    for chain in protein.chains:
      Nres += len(chain.residues)
    dS = 0.0173*Nres
    pH_ref = 5.0; dG_ref = protein.calculateFoldingEnergy(pH_ref, reference=reference)
    if ref == None:
      Tm_ref = 0.00
      if Tms == None:
        Tm_list = [Tm]
      else:
        Tm_list = Tms
      number_of_Tms = float(len(Tm_list))
      ave_diff = 1.0
      while abs(ave_diff) > 0.005:
        ave_diff = 0.00
        for pH, Tm in Tm_list:
          dG = protein.calculateFoldingEnergy(pH, reference=reference)
          dTm = -4.187*(dG - dG_ref)/dS
          Tm_calc = Tm_ref+dTm
          ave_diff += (Tm_calc - Tm)/number_of_Tms
          #Tm_ref -= (Tm_old+dTm - Tm)/(2*number_of_Tms)
        Tm_ref -= ave_diff
        #pka_print("%6.2lf %6.2lf %6.2lf" % (Tm_ref, ave_diff, Tm_ref - Tm_old))
    else:
      dTm_ref = -4.187*(dG_ref - ref[2])/dS
      Tm_ref = ref[1] + dTm_ref

    pka_print("ref = %6.2lf%6.2lf%6.2lf" % (pH_ref, Tm_ref, dG_ref))
    profile = []
    pH, end, increment = grid
    while pH <= end:
      dG = protein.calculateFoldingEnergy(pH, reference=reference)
      dTm = -4.187*(dG - dG_ref)/dS
      profile.append([pH, Tm_ref+dTm])
      pH += increment

    return profile


def ChargeProfile(protein, options=None):
    """
    Calculates the folding profile
    """
    profile = []
    for i_pH in range(0, 15):
      pH = float(i_pH)
      Q_pro, Q_mod = protein.calculateCharge(pH)
      profile.append([pH, Q_pro, Q_mod])

    return profile


def pI(protein, pI=7.0, options=None):
    """
    Calculates the iso electric point
    """
    pI_pro = pI - 0.50
    pI_mod = pI + 0.50
    Q1_pro, Q1_mod = protein.calculateCharge(pI_pro)
    Q2_pro, Q2_mod = protein.calculateCharge(pI_mod)
    iter = 0

    while abs(Q1_pro) > 0.005 and abs(Q2_mod) > 0.005:
      if iter == 50:
        pka_print("pI iterations did not converge after %d iterations %s, switching to bracketing" % (iter, protein.name))
        pI_pro, pI_mod = bracketingPI(protein)
        break
      else:
        iter += 1
      if abs(pI_pro - pI_mod) < 0.010:
        shift_pro = random.random()*0.02 - 0.01
        shift_mod = random.random()*0.02 - 0.01
        shift = (shift_pro-shift_mod)
        pI_pro += shift_pro
        pI_mod += shift_mod
      Q1_pro, Q1_mod = protein.calculateCharge(pI_pro)
      Q2_pro, Q2_mod = protein.calculateCharge(pI_mod)
      k1 = (Q1_pro - Q2_pro)/(pI_pro - pI_mod)
      k2 = (Q2_mod - Q1_mod)/(pI_mod - pI_pro)
      shift = -Q1_pro/k1
      if abs(shift) > 4.0:
        shift = shift/abs(shift)
      pI_pro += shift
      shift = -Q2_mod/k2
      if abs(shift) > 4.0:
        shift = shift/abs(shift)
      pI_mod += shift
      #pka_print("%4d%8.3lf%8.3lf" % (iter, pI_pro, pI_mod))
    #if options.verbose == True:
    #  pka_print("%10d pI iterations" % (iter))

    return pI_pro, pI_mod


def bracketingPI(protein, bracket=[0.0, 14.0]):
    """
    Calculates the pI using 'bracketing'
    """
    iter = 0
    pI     = [0., 0.]
    Q_min  = [0., 0.]; Q_max  = [0., 0.]
    pI_min = [bracket[0], bracket[0]]; pI_max = [bracket[1], bracket[1]]
    Q_min[0], Q_min[1] = protein.calculateCharge( 0.00)
    Q_max[0], Q_max[1] = protein.calculateCharge(14.00)
    while True:
      pI[0] = random.uniform(pI_min[0], pI_max[0])
      pI[1] = random.uniform(pI_min[1], pI_max[1])
      Q = []
      Q.append(protein.calculateCharge(pI[0]))
      Q.append(protein.calculateCharge(pI[1]))
      # folded structure
      if Q[0][0] > 0.00:
        pI_min[0] = pI[0]; Q_min[0] = Q[0][0]
      else:
        pI_max[0] = pI[0]; Q_max[0] = Q[0][0]
      if True:
        if Q[1][0] > 0.00 and Q[1][0] < Q_min[0]:
          pI_min[0] = pI[1]; Q_min[0] = Q[1][0]
        elif Q[1][0] < 0.00 and Q[1][0] > Q_max[0]:
          pI_max[0] = pI[1]; Q_max[0] = Q[1][0]
      # unfolded structure
      if Q[1][1] > 0.00:
        pI_min[1] = pI[1]; Q_min[1] = Q[1][1]
      else:
        pI_max[1] = pI[1]; Q_max[1] = Q[1][1]
      if True:
        if Q[0][1] > 0.00 and Q[0][1] < Q_min[1]:
          pI_min[1] = pI[0]; Q_min[1] = Q[0][1]
        elif Q[0][1] < 0.00 and Q[0][1] > Q_max[1]:
          pI_max[1] = pI[0]; Q_max[1] = Q[0][1]
      iter += 1
      pka_print("%4d protein = %6.2lf [%6.2lf%6.2lf] [%6.2lf%6.2lf]" % (iter, pI[0], Q_min[0], Q_max[0], pI_min[0], pI_max[0]))
      if Q_min[0] <  0.005 and Q_min[1] <  0.005 and \
         Q_max[0] > -0.005 and Q_max[1] > -0.005:
        break

    return pI[0], pI[1]

