#
# $Id: pKa_utility_functions.py 361 2007-05-06 11:44:10Z nielsen $
#
# pKa utility functions
#
# Copyright (C) Jens Erik Nielsen, EMBL 2000, UCSD/HHMI 2002-2003
# University College Dublin 2003 -
# All rights reserved
#

#
# Routines for getting WHAT IF residue identifiers
#

def getWI_resid(line):
    #
    # Take a WHAT IF line of the form:
    # Residue:     1 LYS  (   1  )      (Prp= 1.00)
    #
    # and return a unique resid
    #
    split = line.split()
    residuename = split[2]
    last = line.split('(')[1]
    number_ = last.split(')')[0]
    if len(last.split(')')[1].strip()) > 0:
        chain_id = last.split(')')[1][0]
    else:
        chain_id = ''
    number_ = number_.strip()
    resid = '%s:%4s:%3s' %(chain_id, number_.zfill(4), residuename)
    return resid

def getWI_resid2(line, format_=None):
    """
     Convert a WHAT IF line of the form:
     43 ASP  (  43  )         xxxxxx
     or
     182 SER  ( 182  ) A
     and return a unique resid


    If format is pdb2pka, then just return the line"""
    if format_ == 'pdb2pka':
        return line.split(' ')[0].strip()
    #
    # WHAT IF format
    #

    fields = line.split()[:-1]
    line = fields.join()
    #
    # Line reformatted
    #
    residuename = fields[1]
    last = line.split('(')[1]
    number = last.split(')')[0]
    last_part = last.split(')')[1].strip()
    last_split = last_part.split()
    if len(last_split) == 1:
        chain_id = last_split[0]
    else:
        chain_id = ''
    number = number.strip()
    resid = '%s:%4s:%3s' % (chain_id, number.zfill(4), residuename)
    return resid

def getWI_resid3(line):
    #
    # Take a WHAT IF line of the form:
    # Residue:     1 THR      1    A     Prp= 0.00
    # and return a unique resid
    #
    split = line.split()
    residuename = split[2]
    number = split[3]
    cid = split[4]
    chain_id = ''
    if len(cid) == 1:
        chain_id = cid
    else:
        chain_id = ''
    number = number.strip()
    resid = '%s:%4s:%3s' % (chain_id, number.zfill(4), residuename)
    return resid

def getWI_resid4(line):
    #
    # Take a WHAT IF line of the form:
    #   1 THR  (   1  )       N   xxxxxxx
    # and return a unique resid
    #
    # First trim the string
    #
    line = line.strip()[:-1].join()
    #
    # Now the line looks like: 1 THR ( 1 ) N
    #
    split = line.split()
    atomname = split[-1]
    residuename = split[1]
    last = line.split('(')[1]
    number = last.split(')')[0].strip()
    # Chain ID
    # TODO: 2020/07/06 intendo - I think line 109 should is supposed to be:
    # last = last.strip()
    #   or line 102 should be:
    # last = line.strip().split('(')[1]
    # UNUSED: laststrip = last.strip()
    lastsplit = last.split(')')[1]
    split_lastsplit = lastsplit.split()
    if len(split_lastsplit) > 1:
        chain_id = split_lastsplit[0]
    else:
        chain_id = ''
    resid = '%s:%4s:%3s:%s' % (chain_id, number.zfill(4), residuename, atomname)
    return resid

#
# -----
#

def get_resid(unique_id):
    return unique_id.split(':')[:2].join(':')

def get_resnum(unique_id):
    # Given a unique_id this function returns the residue number
    return unique_id.split(':')[1]

def get_resname(unique_id):
    return unique_id.split(':')[2]

# TODO: 2020/07/06 intendo - this does not seem to be used anywhere. Can it be removed?
#def get_chainid(self, unique_id):
#    return unique_id.split(':')[0]

def is_terminal(unique_id):
    #
    # Is this residue a terminal titratable group?
    #
    if unique_id.split(':')[-1] == 'TERM':
        return 1
    return None

ACIDBASE = {'ARG': 1, 'HIS': 1, 'LYS': 1, 'TYR': -1, 'ASP': -1, 'GLU': -1,
            'CYS': -1, 'CTERM': -1, 'NTERM': 1, 'SER': -1, 'THR': -1}

CHARGED = ['ARG', 'LYS', 'HIS', 'ASP', 'GLU']

def is_normally_titratable(unique_id):
    """Does this group have a pKa value in the range 2-12"""
    type_ = unique_id.split(':')[-1]
    if istitratable(unique_id) and type_ != 'SER' and type_ != 'THR':
        return 1
    return None

def is_titratable(unique_id):
    return istitratable(unique_id)

def istitratable(unique_id):
    type_ = unique_id.split(':')[-1]
    if type_ in ACIDBASE:
        return 1
    else:
        return None

def is_charged(unique_id):
    type_ = unique_id.split(':')[-1]
    if type_ in CHARGED:
        return 1
    return None

def isacid(unique_id):
    #
    # Is the residue a base or an acid?
    #
    type_ = unique_id.split(':')[-1]
    val = ACIDBASE[type_]
    if val == -1:
        return 1
    return None

def acibas(unique_id):
    #
    # Return -1 (acid) or 1 (base)
    #
    if isacid(unique_id):
        return -1
    return 1

def reformat_name(oldname, nterm=None, format_='WHAT IF'):
    #
    # Reformat the residue names
    #
    if format_ == 'WHAT IF':
        import copy
        residue = copy.copy(oldname)
        newname = ''
        tflag = None
        if residue[0] == 'T':
            tflag = 1
            residue = residue[1:]
        number = residue[:4]
        name = residue[4:7]
        chain_id = ''
        if len(residue) == 8:
            chain_id = residue[7]
        newname = '%s:%4s:%3s' % (chain_id, number, name)
        if tflag:
            #
            # Nterm is 1 for an Nterm and 0 for a Cterm
            #
            if nterm == 1:
                newname += ':NTERM'
            else:
                newname += ':CTERM'
    else:
        return oldname.strip()
    return newname
