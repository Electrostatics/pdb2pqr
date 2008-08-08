#!@WHICHPYTHON@
"""
    Displays the Jmol input page
"""
print 'Content-type: text/html\n\n'

__date__ = "18 June 2008"
__author__ = "Samir Unni"
__version__ = "0.0.1"

from src.aconf import *
import cgi, cgitb, pickle, urllib, os, glob

cgitb.enable()
form = cgi.FieldStorage()

def initVars():
    #global serviceURL

    if not form.has_key("jobid"):
        pass # add code to redirect to PDB2PQR input page here
    else:
        jobid = form['jobid'].value

        #aspFile = open('./tmp/%s/%s-asp' % (logTime, logTime))
        #appServicePort = pickle.load(aspFile)
        #aspFile.close()

        #jobInfoFile = open('./tmp/%s/%s-jobinfo' % (logTime, logTime))
        #jobID = jobInfoFile.read()
        #jobInfoFile.close()

        aoFile = open('%s%s%s/%s-ao' % (INSTALLDIR, TMPDIR, jobid, jobid))
        apbsOptions = pickle.load(aoFile)
        aoFile.close()

        if APBS_OPAL_URL!="":
            from AppService_client import AppServiceLocator, getOutputsRequest
            apbsOpalJobIDFile = open('%s%s%s/apbs_opal_job_id' % (INSTALLDIR, TMPDIR, jobid))
            apbsOpalJobID = apbsOpalJobIDFile.read()
            apbsOpalJobIDFile.close()
            
            appLocator = AppServiceLocator()
            resp = appLocator.getAppServicePort(APBS_OPAL_URL).getOutputs(getOutputsRequest(apbsOpalJobID))
            if not os.access('%s%s%s' % (INSTALLDIR, TMPDIR, jobid), os.F_OK):
                os.mkdir('%s%s%s' % (INSTALLDIR, TMPDIR, jobid))
            for file in resp._outputFile:
                fileName = file._name
                if fileName!="Standard Output" and fileName!="Standard Error":
                    urllib.urlretrieve(file._url, '%s%s%s/%s' % (INSTALLDIR, TMPDIR, jobid, fileName))

        return apbsOptions


def main(apbsOptions):
    cgiFile = "jmol.cgi"
    cgiName = "thisform"
    defaultVisType = "jmol"
    checkJmolType = True
    cssFile = '../pdb2pqr/pdb2pqr.css'
    jobid = form['jobid'].value

    print '<html>'
    print '\t<head>'
    print '\t\t<title>Visualization Configuration</title>'
    print '\t\t<link rel=\"stylesheet\" href=\"%s\" type=\"text/css\">' % cssFile
    print '\t</head>'
    print '\t<body>'
    print '\t\t<h3>Visualization Configuration <span style="color:red">(EXPERIMENTAL)</span></h3>'
    print '\t\t<form action=\"%s\" method=\"post\" enctype=\"multipart/form-data\" name=\"%s\" id=\"%s\">\n' % (cgiFile, cgiName, cgiName)
    print '\t\tSelect the type of visual representation:'
    print '\t\t<ul>'
    print '\t\t\t<input type=\"radio\" name=\"vistype\" value=\"jmol\"',

    if defaultVisType == "jmol":
        print ' checked=\"checked\"',
    print '/> <a href=\"http://jmol.sourceforge.net/\" target=\"_blank\">Jmol</a> <a href=\"http://jmol.sourceforge.net/docs/\" target=\"_blank\"><font title=\"Jmol documentation\">(<span class=\"tooltip\">?</span>)</font></a>'
    print '\t\t<br />'
    print '\t\t<ul>'
    print '\t\t\t<li>Select the type of display:'
    print '\t\t\t\t<ul>'
    if apbsOptions['writePot']:
        print '\t\t\t\t\t<input type=\"radio\" name=\"jmoltype\" value=\"pot\"',
        if checkJmolType:
            print ' checked=\"checked\"',
            checkJmolType = False
        print '/>Electrostatic Potential'
        print '\t\t\t\t\t<br />'

    if apbsOptions['writeLap']:
        print '\t\t\t\t\t<input type=\"radio\" name=\"jmoltype\" value=\"lap\"',
        if checkJmolType:
            print ' checked=\"checked\"',
            checkJmolType = False
        print '/>Laplacian of the potential'
        print '<br />'
    if apbsOptions['writeEdens']:
        print '\t\t\t\t\t<input type=\"radio\" name=\"jmoltype\" value=\"edens\"',
        if checkJmolType:
            print ' checked=\"checked\"',
            checkJmolType = False
        print '/>Energy density'
        print '<br />'
    if apbsOptions['writeNdens']:
        print '\t\t\t\t\t<input type=\"radio\" name=\"jmoltype\" value=\"ndens\"',
        if checkJmolType:
            print ' checked=\"checked\"',
            checkJmolType = False
        print '/>Mobile ion number density'
        print '<br />'
    if apbsOptions['writeQdens']:
        print '\t\t\t\t\t<input type=\"radio\" name=\"jmoltype\" value=\"qdens\"'
        if checkJmolType:
            print ' checked=\"checked\"',
            checkJmolType = False
        print '/>Mobile charge density'
    print '\t\t\t\t</ul>'
    print '\t\t\t</li>'
    print '\t\t\t<li>'
    print '\t\t\t\tMinimum color range:<input type=\"text\" name=\"jmolcolormin\" maxlength=\"10\" value="-5.0"/>'
    print '\t\t\t</li>'
    print '\t\t\t<li>'
    print '\t\t\t\tMaximum color range:<input type=\"text\" name=\"jmolcolormax\" value="5.0"/>'
    print '\t\t\t</li>'
    print '\t\t</ul>'
    #print '\t\t<input type=\"radio\" name=\"vistype\" value=\"vmd\" disabled/> VMD (currently unavailable)'
    print '\t\t<input type=\"hidden\" name=\"jobid\" value=\"%s\"/>' % jobid
    print '\t\t<input type=\"hidden\" name=\"pqrfilename\" value=\"%s\"/>\n' % apbsOptions['pqrFileName']
    print '\t\t\t\t</ul>'
    print '\t\t\t<input type=\"submit\" value=\"Submit\"/>'
    print '\t\t</form>'

    print '\t</body>'
    print '</html>'


main(initVars())
