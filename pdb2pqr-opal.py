#!@WHICHPYTHON@
"""
    Driver for interfacing pdb2pqr web form with OPAL SOAP job submission.
    
"""

__date__  = "27 August 2007"
__author__ = "Wes Goodman"
__version__ = "1.3.0"

import string
import sys
import getopt
import os
import time
import re
from src import pdb
from src import utilities
from src import structures
from src import routines
from src import protein
from src import server
from src.pdb import *
from src.utilities import *
from src.structures import *
from src.definitions import *
from src.forcefield import *
from src.routines import *
from src.protein import *
from src.server import *
from src.hydrogens import *
from StringIO import *

import httplib
from AppService_services import AppServiceLocator, getAppMetadataRequest, launchJobRequestWrapper, \
    launchJobBlockingRequestWrapper, getOutputAsBase64ByNameRequestWrapper
from AppService_services_types import ns1
from ZSI.TC import String
import cgi
import cgitb

def printheader(pagetitle,jobid=None):
  """
        Function to print html headers
  """
  if jobid:
    print "Location: querystatus.cgi?jobid=%s\n" % jobid
  print "Content-type: text/html\n"
  print "<HTML>"
  print "<HEAD>"
  print "\t<TITLE>%s</TITLE>" % pagetitle
  print "\t<link rel=\"stylesheet\" href=\"%s\" type=\"text/css\">\n" % STYLESHEET
  print "</HEAD>"
  return
   
def mainCGI():
    """
        Opal driver for running PDB2PQR from a web page
    """
    serviceURL = "@OPALURL@"
    
    cgitb.enable()
    form = cgi.FieldStorage()
    
    options = {}
 
    ff = form["FF"].value 
    options["ff"] = ff
    fffile = None
    input = 0
  
    if form.has_key("DEBUMP"): options["debump"] = 1
    else: options["debump"] = 0
    if form.has_key("OPT"): options["opt"] = 1
    else: options["opt"] = 0
    if form.has_key("PROPKA"):
        try:
            ph = float(form["PH"].value)
            if ph < 0.0 or ph > 14.0: raise ValueError
            options["ph"] = ph
        except ValueError:
             text = "The entered pH of %.2f is invalid!  " % form["PH"].value
             text += "Please choose a pH between 0.0 and 14.0."
#             print "Content-type: text/html\n"
             print text
             sys.exit(2)
    if form.has_key("PDBID"):
        filename = form["PDBID"].value
        infile = getPDBFile(form["PDBID"].value)
    elif form.has_key("PDB"):
        filename = form["PDB"].filename
        filename=re.split(r'[/\\]',filename)[-1]
        infile = StringIO(form["PDB"].value)
    if form.has_key("INPUT"):
        input = 1
        options["apbs"] = 1
    if form.has_key("USERFF"):
#        userff = StringIO(form["USERFF"].value)
#        ff = "user-defined"
#        options["userff"] = userff
        ffname = form["USERFF"].filename
        ffname = re.split(r'[/\\]',ffname)[-1]
        fffile = StringIO(form["USERFF"].value)
        options["ff"] = ffname
    if form.has_key("FFOUT"):
        if form["FFOUT"].value != "internal":
            options["ffout"] = form["FFOUT"].value
    if form.has_key("CHAIN"):
        options["chain"] = 1
    if form.has_key("LIGAND"):
        ligandfilename=str(form["LIGAND"].filename)
        ligandfilename=re.split(r'[/\\]',ligandfilename)[-1]
        options["ligand"] = StringIO(form["LIGAND"].value)
        
    try:
#        starttime = time.time()
#        name = setID(starttime)
        name = filename
        ligandFile=None
        ffFile=None
        # begin SOAP changes
        # need to switch options from a dictionary to something resembling a command line query
        # such as --chain
        myopts=""
        for key in options:
            if key=="opt":
                if options[key]==0:
                    # user does not want optimization
                    key="noopt"
                else:
                    # pdb2pqr optimizes by default, don't bother with flag
                    continue
            elif key=="debump":
                if options[key]==0:
                    # user does not want debumping
                    key="nodebump"
                else:
                    # pdb2pqr debumps by default, so change this flag to --nodebump
                    continue
            elif key=="ph":
                val=options[key]
                key="with-ph=%s" % val
            elif key=="ffout":
                val=options[key]
                key="ffout=%s" % val
            elif key=="ligand":
                val=ligandfilename
                key="ligand=%s" % val
                ligandFile = ns1.InputFileType_Def()
                ligandFile.Set_name(val)
                ligandFileString = options["ligand"].read()
                options["ligand"].close()
                ligandFile.Set_contents(ligandFileString)
            elif key=="apbs":
                key="apbs-input"
            elif key=="chain":
                key="chain"
            elif key=="ff":
                val=options[key]
                key="ff=%s" % val
                if fffile:
                  ffFile = ns1.InputFileType_Def()
                  ffFile.Set_name(val)
                  ffFileString = fffile.read()
                  fffile.close()
                  ffFile.Set_contents(ffFileString)
            myopts+="--"+str(key)+" "
        myopts+=str(filename)+" "
        myopts+="%s.pqr" % str(name)
        appLocator = AppServiceLocator()
        appServicePort = appLocator.getAppServicePortType(serviceURL)
        # launch job
        req = launchJobRequestWrapper()
        req._argList = myopts
        inputFiles = []
        pdbFile = ns1.InputFileType_Def()
        pdbFile.Set_name(filename)
        pdbFileString = infile.read()
        infile.close()
        pdbFile.Set_contents(pdbFileString)
        inputFiles.append(pdbFile)
        if ligandFile:
          inputFiles.append(ligandFile)
        if ffFile:
          inputFiles.append(ffFile)
        req._inputFile=inputFiles
        try:
          resp=appServicePort.launchJob(req)
        except Exception, e:
          printheader("PDB2PQR Job Submission - Error")
          print "<BODY>\n<P>"
          print "There was an error with your job submission<br>"
          print "</P>\n</BODY>"
          print "</HTML>"
          sys.exit(2)
        printheader("PDB2PQR Job Submission",resp._jobID)
        
    except StandardError, details:
        print details
        createError(name, details)

# File should only be called as CGI
mainCGI()
