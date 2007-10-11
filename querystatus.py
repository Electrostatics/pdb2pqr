#!@WHICHPYTHON@
"""
  CGI Module for checking on the status of an OPAL job
"""

__date__   = "27 August 2007"
__author__ = "Wes Goodman"

import sys
import cgi
import cgitb
from src.server import STYLESHEET
from AppService_services import AppServiceLocator, getAppMetadataRequest, launchJobRequestWrapper, \
    launchJobBlockingRequestWrapper, getOutputAsBase64ByNameRequestWrapper
from AppService_services_types import ns1
from ZSI.TC import String

serviceURL="@OPALURL@"
refresh=30

cgitb.enable()
form = cgi.FieldStorage()

def printheader(pagetitle,refresh=None):
  print "Content-type: text/html\n"
  print "<HTML>"
  print "<HEAD>"
  if refresh:
    print "\t<META HTTP-EQUIV=\"Refresh\" CONTENT=\"%s\">" % refresh
  print "\t<TITLE>%s</TITLE>" % pagetitle
  print "\t<link rel=\"stylesheet\" href=\"%s\" type=\"text/css\">\n" % STYLESHEET
  print "</HEAD>"
  return

if not form.has_key("jobid"):
  printheader("PDB2PQR Job Status - Error")
  text="<BODY>\n"
  text+="\t<H2>Missing jobid field</H2>\n"
  text+="\t<P>Your request url is missing the jobid field</P>\n"
  text+="</BODY>\n</HTML>"
  print text
  sys.exit(2)

# Construct SOAP request
appLocator = AppServiceLocator()
appServicePort = appLocator.getAppServicePortType(serviceURL)
jobid=form["jobid"].value
try:
  status=appServicePort.queryStatus(str(jobid))
except Exception, e:
  printheader("PDB2PQR Job Status Page - Error", 0)
  print "<BODY>\n<P>"
  print "There was an error with your query request.  This page will not refresh."
  print "</P>\n</BODY>"
  print "</HTML>"
  sys.exit(2)

if status._code==8 or status._code==4:
  # job is done
  printheader("PDB2PQR Job Status Page")
  print "<BODY>\n<P>"
  print "<h3>Status</h3>"
  print "Code: %s<BR>" % status._code
  print "Message: %s<BR>" % status._message
  print "Output Base URL: <a href=%s>%s</a><BR>" % (status._baseURL, status._baseURL)
  print "</P>"
  print "<HR>"
  print "Click here to see your results: <a href=%s>%s</a>" % (status._baseURL, status._baseURL)
  print "</BODY>"
else:
  # job is not done, refresh in 30 seconds
  printheader("PDB2PQR Job Status Page", refresh)
  print "<BODY>\n<P>"
  print "<h3>Status</h3>"
  print "Code: %s<BR>" % status._code
  print "Message: %s<BR>" % status._message
  print "Output Base URL: <a href=%s>%s</a><BR>" % (status._baseURL, status._baseURL)
  print "Page will refresh in %d seconds" % refresh
  print "</P>"
  print "</BODY>"

print "</HTML>"