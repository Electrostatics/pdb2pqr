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
from src.utilities import getEventTrackingString, getTrackingScriptString

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
    cssFile = 'pdb2pqr.css'
    jobid = form['jobid'].value

    string =  """
<html>
    <head>
        {trackingscript}
        <script type="text/javascript">
            {trackingevents}
        </script>
        <title>Visualization</title>
        <link rel="stylesheet" href="3dmol/css/pdb2pqr_3dmol.css" type="text/css">
        <link rel="stylesheet" href="3dmol/css/foundation.css">
        <link rel="stylesheet" type="text/css" href="3dmol/css/pure-min.css" media="screen" />
        <link rel="stylesheet" type="text/css" href="pdb2pqr.css" media="screen" />

        <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
        <script src="http://3dmol.csb.pitt.edu/build/3Dmol.js"></script>
        <script type="text/JavaScript" src="3dmol/visualize_html.js"></script>
        <script type="text/JavaScript" src="3dmol/3dmol.js"></script>
    </head>
    <body>
        <script type="text/javascript">build_page({jobid})</script>
        <script type="text/javascript">getpqr({jobid})</script>
    </body>
</html>""".format(jobid=jobid,
                  trackingevents=getEventTrackingString(category='apbs',
                                                        action='visualize', 
                                                        label=str(os.environ["REMOTE_ADDR"])),
                  trackingscript=getTrackingScriptString(jobid=jobid))

    print string


main(initVars())
