#!@WHICHPYTHON@
"""
	Jmol applet generator
"""

__date__ = "5 July 2007"
__author__ = "Samir Unni"
__version__ = "0.0.1"

import string, sys, os, pickle, cgi, cgitb, locale, shutil
from sys import stdout, stderr, stdin
from src.aconf import *

def jmolGen():
	file = stdout
	file.write("Content-type: text/html\n\n")
	cgitb.enable()
	

	form = cgi.FieldStorage()

	visType = form["vistype"].value
	jobid = form["jobid"].value
	pqrFileName = form["pqrfilename"].value
	scriptName = "myscript.spt"
	jmolPage = "jmol.html"
	jmolJsFile = "Jmol.js"
	jmolAppFile = "JmolApplet.jar"
	jmolColorMin = form["jmolcolormin"].value
	jmolColorMax = form["jmolcolormax"].value

	try:
		jmolType = form["jmoltype"].value
	except KeyError:
		pass

	dxFilePath = "%s-%s.dx" % (pqrFileName[:-4], jmolType)

	shutil.copyfile(jmolJsFile, '%s%s%s/%s' % (INSTALLDIR, TMPDIR, jobid, jmolJsFile))
	shutil.copyfile(jmolAppFile, '%s%s%s/%s' % (INSTALLDIR, TMPDIR, jobid, jmolAppFile))

	script = open('%s%s%s/%s' % (INSTALLDIR, TMPDIR, jobid, scriptName), 'w')
	script.write('load %s\n' % pqrFileName)
	script.write('isosurface s1 colorscheme bwr color absolute %s %s sasurface map \"%s\"\n' % ( jmolColorMin, jmolColorMax, dxFilePath) )
	script.write('write isosurface pot.jvxl')
	script.close()

	jmolpage = open('%s%s%s/%s' % (INSTALLDIR, TMPDIR, jobid, jmolPage), 'w')
	jmolpage.write('<html>\n')
	jmolpage.write('<head><script type=\"text/javascript\" src=\"%s\"></script></head>\n' % jmolJsFile)
	jmolpage.write('<body>\n')
	jmolpage.write('<script type=\"text/javascript\">\n')
	#jmolpage.write('jmolInitialize(\"../\");\n');
	jmolpage.write('jmolApplet(600, \"script myscript.spt\");\n')
	jmolpage.write('</script>\n')
	jmolpage.write('</body>\n')
	jmolpage.write('</html>\n')

	jmolpage.close()


	sys.stdout.write('<html> <body>')
	sys.stdout.write('<meta http-equiv=\"refresh\" content=\"1;url=%s\">' % ('tmp/%s/%s' % (jobid, jmolPage)))
	sys.stdout.write('</body> </html>')
	sys.stdout.flush()

jmolGen()

