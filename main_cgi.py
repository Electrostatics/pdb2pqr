"""
    Driver for PDB2PQR

    This module takes a PDB file as input and performs optimizations
    before yielding a new PDB-style file as output.

    Ported to Python by Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis

    Parsing utilities provided by Nathan A. Baker (baker@biochem.wustl.edu)
    Washington University in St. Louis

	Copyright (c) 2002-2009, Jens Erik Nielsen, University College Dublin; 
	Nathan A. Baker, Washington University in St. Louis; Paul Czodrowski & 
	Gerhard Klebe, University of Marburg

	All rights reserved.

	Redistribution and use in source and binary forms, with or without modification, 
	are permitted provided that the following conditions are met:

		* Redistributions of source code must retain the above copyright notice, 
		  this list of conditions and the following disclaimer.
		* Redistributions in binary form must reproduce the above copyright notice, 
		  this list of conditions and the following disclaimer in the documentation 
		  and/or other materials provided with the distribution.
		* Neither the names of University College Dublin, Washington University in 
		  St. Louis, or University of Marburg nor the names of its contributors may 
		  be used to endorse or promote products derived from this software without 
		  specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
	ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
	WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
	IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
	INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
	BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
	OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED 
	OF THE POSSIBILITY OF SUCH DAMAGE.

"""

__date__  = "2 September 2009"
__author__ = "Todd Dolinsky, Nathan Baker, Jens Nielsen, Paul Czodrowski, Jan Jensen, Samir Unni, Yong Huang"
__version__ = "1.5"


import string
import sys
import getopt
import os
import time
import httplib
import re
import glob
import tempfile
from src import pdb
from src import utilities
from src import structures
from src import routines
from src import protein
#from src import server
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
from StringIO import *
from main import *

def printHeader(pagetitle,have_opal=None,jobid=None):
    """
        Function to print html headers
    """
    if jobid:
        if have_opal:
            print "Location: querystatus.cgi?jobid=%s&typeofjob=opal\n" % (jobid,typeOfJob)
        else:
            print "Location: querystatus.cgi?jobid=%s&typeofjob=local\n" & (jobid,typeOfJob)

    #print "Content-type: text/html\n"
    print "<HTML>"
    print "<HEAD>"
    print "\t<TITLE>%s</TITLE>" % pagetitle
    print "\t<link rel=\"stylesheet\" href=\"%s\" type=\"text/css\">\n" % STYLESHEET
    print "</HEAD>"
    return

def redirector(name):
    """
        Prints a page which redirects the user to querystatus.cgi and writes starting time to file
    """

    starttimefile = open('%s%s%s/pdb2pqr_start_time' % (INSTALLDIR, TMPDIR, name), 'w')
    starttimefile.write(str(time.time()))
    starttimefile.close()

    string = ""
    string+= "<html>\n"
    string+= "\t<head>\n"
    string+= "\t\t<meta http-equiv=\"Refresh\" content=\"0; url=%squerystatus.cgi?jobid=%s&calctype=pdb2pqr\">\n" % (WEBSITE, name)
    string+= "\t</head>\n"
    string+= "</html>\n"
    return string

def mainCGI():
    """
        Main driver for running PDB2PQR from a web page
    """
    print "Content-type: text/html\n"
    import cgi
    import cgitb

    cgitb.enable()
    form = cgi.FieldStorage()
    ff = form["FF"].value 
    input = 0

    apbs_input = form.has_key("INPUT")
    typemap = form.has_key("TYPEMAP")
    neutraln = form.has_key("NEUTRALN")
    neutralc = form.has_key("NEUTRALC")

    if HAVE_PDB2PQR_OPAL=="1":
        have_opal = True
        # Opal-specific import statments
        from AppService_client import AppServiceLocator, getAppMetadataRequest, launchJobRequest, launchJobBlockingRequest, getOutputAsBase64ByNameRequest
        from AppService_types import ns0
        from ZSI.TC import String
    else:
        have_opal = False

    if have_opal:
        options = {}
        options["ff"] = ff
        fffile = None
        namesfile = None
    else:
        options = {"extensions":{}}
 
  
    if form.has_key("DEBUMP"):
        options["debump"] = 1
    else:
        options["debump"] = 0
    if form.has_key("OPT"):
        options["opt"] = 1
    else:
        options["opt"] = 0
    if form.has_key("PROPKA"):
        try:
            ph = float(form["PH"].value)
            if ph < 0.0 or ph > 14.0: raise ValueError
            options["ph"] = ph
        except ValueError:
             text = "The entered pH of %.2f is invalid!  " % form["PH"].value
             text += "Please choose a pH between 0.0 and 14.0."
             #print "Content-type: text/html\n"
             print text
             sys.exit(2)
    if form.has_key("PDBID"):
        pdbfile = getPDBFile(form["PDBID"].value)
        pdbfilename = form["PDBID"].value
    elif form.has_key("PDB"):
        pdbfile = StringIO(form["PDB"].value)
        pdbfilename = form["PDB"].filename
        pdbfilename = re.split(r'[/\\]',pdbfilename)[-1]
    if form.has_key("INPUT"):
        input = 1
        options["apbs"] = 1
    if form.has_key("USERFF"):
        if have_opal:
            ffname = form["USERFF"].filename
            ffname = re.split(r'[/\\]',ffname)[-1]
            if ffname[-4:] == ".DAT":
               ffname = ffname[:-4]
            fffile = StringIO(form["USERFF"].value)
            namesfile = StringIO(form["USERNAMES"].value)
            options["ff"] = ffname
            options["userff"] = fffile
            options["usernames"] = namesfile
        else:
            userff = StringIO(form["USERFF"].value)
            usernames = StringIO(form["USERNAMES"].value)
            options["ff"] = "user-defined"
            options["userff"] = userff
            options["usernames"] = usernames
    if form.has_key("FFOUT"):
        if form["FFOUT"].value != "internal":
            options["ffout"] = form["FFOUT"].value
    if form.has_key("CHAIN"):
        options["chain"] = 1
    if form.has_key("WHITESPACE"):
        options["whitespace"] = 1
    if form.has_key("TYPEMAP"):
        options["typemap"] = 1
    if form.has_key("NEUTRALN"):
        options["neutraln"] = 1
    if form.has_key("NEUTRALC"):
        options["neutralc"] = 1
    if form.has_key("LIGAND"):
        if have_opal:
            ligandfilename=str(form["LIGAND"].filename)
            ligandfilename=re.split(r'[/\\]',ligandfilename)[-1]

        # for Windows-style newline compatibility
        templigandfilename = tempfile.mkstemp()[1]
        templigandfile = open(templigandfilename,'w')
        templigandfile.write(form["LIGAND"].value)
        templigandfile.close()
        templigandfile = open(templigandfilename,'rU')
        if have_opal:
            options["ligand"] = templigandfile.read()
        else:
            templigandstring = templigandfile.read() # this variable is used again later to write this file to output
            options["ligand"] = StringIO(templigandstring)
            
        templigandfile.close()
        
    if not have_opal:
        pdbfilestring = pdbfile.read()
        pdblist, errlist = readPDB(StringIO(pdbfilestring))
        dummydef = Definition()
        dummyprot = Protein(pdblist, dummydef)
        if len(pdblist) == 0 and len(errlist) == 0:
            text = "Unable to find PDB file - Please make sure this is "
            text += "a valid PDB file ID!"
            #print "Content-type: text/html\n"
            print text
            sys.exit(2)
        elif dummyprot.numAtoms() > MAXATOMS and "opt" in options:
            text = "<HTML><HEAD>"
            text += "<TITLE>PDB2PQR Error</title>"
            text += "<link rel=\"stylesheet\" href=\"%s\" type=\"text/css\">" % STYLESHEET
            text += "</HEAD><BODY><H2>PDB2PQR Error</H2><P>"
            text += "Due to server limits, we are currently unable to optimize "
            text += "proteins of greater than MAXATOMS atoms on the server (PDB2PQR "
            text += "found %s atoms in the selected PDB file).  If you " % dummyprot.numAtoms()
            text += "want to forgo optimization please try the server again.<P>"
            text += "Otherwise you may use the standalone version of PDB2PQR that "
            text += "is available from the <a href=\"http://pdb2pqr.sourceforge.net\">"
            text += "PDB2PQR SourceForge project page</a>."
            text += "<script type=\"text/javascript\">"
            text += "var gaJsHost = ((\"https:\" == document.location.protocol) ? \"https://ssl.\" : \"http://www.\");"
            text += "document.write(unescape(\"%3Cscript src=\'\" + gaJsHost + \"google-analytics.com/ga.js\' type=\'text/javascript\'%3E%3C/script%3E\"));"
            text += "</script>"
            text += "<script type=\"text/javascript\">"
            text += "try {"
            text += "var pageTracker = _gat._getTracker(\"UA-11026338-3\");"
            for key in options:
                text += "pageTracker._trackPageview(\"/main_cgi/has_%s_%s.html\");" % (key, options[key])
            text += "pageTracker._trackPageview();"
            text += "} catch(err) {}</script>"
            text += "</BODY></HTML>"
            #print "Content-type: text/html\n"
            print text
            sys.exit(2)

    try:
        if have_opal:
            ligandFile=None
            ffFile=None
            namesFile=None
        #else:
        starttime = time.time()
        name = setID(starttime)

        os.makedirs('%s%s%s' % (INSTALLDIR, TMPDIR, name))
        apbsInputFile = open('%s%s%s/apbs_input' % (INSTALLDIR, TMPDIR, name),'w')
        apbsInputFile.write(str(apbs_input))
        apbsInputFile.close()
        typemapInputFile = open('%s%s%s/typemap' % (INSTALLDIR, TMPDIR, name),'w')
        typemapInputFile.write(str(typemap))
        typemapInputFile.close()

        if have_opal:
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
                    ligandFile = ns0.InputFileType_Def('inputFile')
                    ligandFile._name = val
                    ligandFile._contents = options["ligand"]
                elif key=="apbs":
                    key="apbs-input"
                elif key=="chain":
                    key="chain"
                elif key=="whitespace":
                    key="whitespace"
                elif key=="typemap":
                    key="typemap"
                elif key=="ff":
                    val=options[key]
                    key="ff=%s" % val
                    if fffile:
                      ffFile = ns0.InputFileType_Def('inputFile')
                      ffFile._name = val + ".DAT"
                      ffFileString = fffile.read()
                      ffFile._contents = ffFileString
                    if namesfile:
                      namesFile = ns0.InputFileType_Def('inputFile')
                      namesFile._name = val + ".names"
                      namesFileString = namesfile.read()
                      namesFile._contents = namesFileString
                myopts+="--"+str(key)+" "
            myopts+=str(pdbfilename)+" "
            myopts+="%s.pqr" % str(pdbfilename)
            appLocator = AppServiceLocator()
            appServicePort = appLocator.getAppServicePort(PDB2PQR_OPAL_URL)
            # launch job
            req = launchJobRequest()
            req._argList = myopts
            inputFiles = []
            pdbOpalFile = ns0.InputFileType_Def('inputFile')
            pdbOpalFile._name = pdbfilename
            pdbOpalFile._contents = pdbfile.read()
            pdbfile.close()
            inputFiles.append(pdbOpalFile)
            if ligandFile:
              inputFiles.append(ligandFile)
            if ffFile:
              inputFiles.append(ffFile)
            if namesFile:
              inputFiles.append(namesFile)
            req._inputFile=inputFiles
            try:
                resp=appServicePort.launchJob(req)
            except Exception, e:
                printHeader("PDB2PQR Job Submission - Error")
                print "<BODY>\n<P>"
                print "There was an error with your job submission<br>"
                print "</P>"
                print "<script type=\"text/javascript\">"
                print "var gaJsHost = ((\"https:\" == document.location.protocol) ? \"https://ssl.\" : \"http://www.\");"
                print "document.write(unescape(\"%3Cscript src=\'\" + gaJsHost + \"google-analytics.com/ga.js\' type=\'text/javascript\'%3E%3C/script%3E\"));"
                print "</script>"
                print "<script type=\"text/javascript\">"
                print "try {"
                print "var pageTracker = _gat._getTracker(\"UA-11026338-3\");"
                for key in options:
                    print "pageTracker._trackPageview(\"/main_cgi/has_%s_%s.html\");" % (key, options[key])
                print "pageTracker._trackPageview();"
                print "} catch(err) {}</script>"
                print "</BODY>"
                print "</HTML>"
                sys.exit(2)
            #printHeader("PDB2PQR Job Submission",have_opal,jobid=resp._jobID)
            pdb2pqrOpalJobIDFile = open('%s%s%s/pdb2pqr_opal_job_id' % (INSTALLDIR, TMPDIR, name), 'w')
            pdb2pqrOpalJobIDFile.write(resp._jobID)
            pdb2pqrOpalJobIDFile.close()
            print redirector(name)
            # Recording CGI run information for PDB2PQR Opal
            pdb2pqrOpalLogFile = open('%s%s%s/pdb2pqr_opal_log' % (INSTALLDIR, TMPDIR, name), 'w')
            pdb2pqrOpalLogFile.write(str(options)+'\n'+str(ff)+'\n'+str(os.environ["REMOTE_ADDR"]))
            pdb2pqrOpalLogFile.close()

        else:
            #pqrpath = startServer(name)
            statusfile = open('%s%s%s/pdb2pqr_status' % (INSTALLDIR, TMPDIR, name), 'w')
            statusfile.write('running')
            statusfile.close()


            pid = os.fork()
            if pid:
                print redirector(name)
                sys.exit()
            else:
                currentdir = os.getcwd()
                os.chdir("/")
                os.setsid()
                os.umask(0)
                os.chdir(currentdir)
                os.close(1) # not sure if these
                os.close(2) # two lines are necessary
                pqrpath = '%s%s%s/%s.pqr' % (INSTALLDIR, TMPDIR, name, name)
                options["outname"] = pqrpath
                options["verbose"] = ""
                orig_stdout = sys.stdout
                orig_stderr = sys.stderr
                sys.stdout = open('%s%s%s/pdb2pqr_stdout.txt' % (INSTALLDIR, TMPDIR, name), 'w')
                sys.stderr = open('%s%s%s/pdb2pqr_stderr.txt' % (INSTALLDIR, TMPDIR, name), 'w')
                header, lines, missedligands = runPDB2PQR(pdblist, ff, options)
                sys.stdout.close()
                sys.stderr.close()
                sys.stdout = orig_stdout
                sys.stderr = orig_stderr

                endtimefile = open('%s%s%s/pdb2pqr_end_time' % (INSTALLDIR, TMPDIR, name), 'w')
                endtimefile.write(str(time.time()))
                endtimefile.close()

                pqrfile = open(pqrpath, "w")
                pqrfile.write(header)
                for line in lines:
                    # Adding whitespaces if --whitespace is in the options
                    if "whitespace" in options.keys() and options["whitespace"] == 1: 
                        if line[0:4] == 'ATOM':
                            newline = line[0:16] + ' ' + line[16:38] + ' ' + line[38:46] + ' ' + line[46:]
                            pqrfile.write("%s\n" % string.strip(newline))
                        elif line[0:6] == 'HETATM':
                            newline = line[0:16] + ' ' + line[16:38] + ' ' + line[38:46] + ' ' + line[46:]
                            pqrfile.write("%s\n" % string.strip(newline))
                    else: 
                        pqrfile.write("%s\n" % string.strip(line))
                pqrfile.close()
                        
                if input:
                    from src import inputgen
                    from src import psize
                    method = "mg-auto"
                    size = psize.Psize()
                    size.parseInput(pqrpath)
                    size.runPsize(pqrpath)
                    async = 0 # No async files here!
                    myinput = inputgen.Input(pqrpath, size, method, async)
                    myinput.printInputFiles()
                    myinput.dumpPickle()
                            
                endtime = time.time() - starttime
                #createResults(header, input, name, endtime, missedligands)
                logRun(options, endtime, len(lines), ff, os.environ["REMOTE_ADDR"])
                #printHeader("PDB2PQR Job Submission",have_opal,jobid=name)
                if form.has_key("LIGAND"):
                    outputligandfile = open('%s%s%s/%s.mol2' % (INSTALLDIR,TMPDIR, name, name),'w')
                    outputligandfile.write(templigandstring)
                    outputligandfile.close()
                outputpdbfile = open('%s%s%s/%s.pdb' % (INSTALLDIR,TMPDIR,name,name),'w')
                outputpdbfile.write(pdbfilestring)
                outputpdbfile.close()

                statusfile = open('%s%s%s/pdb2pqr_status' % (INSTALLDIR, TMPDIR, name), 'w')
                statusfile.write('complete\n')
                filelist = glob.glob('%s%s%s/%s*' % (INSTALLDIR, TMPDIR, name, name))
                for filename in filelist:
                    statusfile.write(filename+'\n')
                statusfile.close()


    except StandardError, details:
        print details
        createError(name, details)
