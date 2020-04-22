#!@WHICHPYTHON@
"""
  CGI Module for checking on the status of an OPAL job
"""

__date__   = "4 January 2010"
__author__ = "Wes Goodman, Samir Unni, Yong Huang"

import sys
import cgi
import cgitb
import os,shutil,glob,string,time,urllib
from datetime import timedelta
from src.server import *
from src.aconf import *
from src.utilities import getTrackingScriptString, getEventTrackingString

cgitb.enable()
form = cgi.FieldStorage()

def printheader(pagetitle,refresh=None,jobid=None):
    str = ""
    str+= "<html>\n"
    str+= "<HEAD>\n"
    str+= '<img src="https://raw.githubusercontent.com/Electrostatics/apbs-pdb2pqr/master/apbs/doc/icons/APBS_128_v2.png" style="float:right; position:relative;right:200px; top: 2px;">'
    if refresh:
        str+= "\t<META HTTP-EQUIV=\"Refresh\" CONTENT=\"%s\">\n" % refresh
    str+= "\t<TITLE>%s</TITLE>\n" % pagetitle
    str+= "\t<link rel=\"stylesheet\" href=\"%s\" type=\"text/css\">\n" % STYLESHEET
    str+= getTrackingScriptString(jobid)
    str+= "</HEAD>\n"
    return str

def createcube(dx_input, pqr_input, output):

    with open(dx_input, 'r') as in_f, open(output, 'w') as out_f, open(pqr_input, 'r') as in_pqr:
        out_f.write("CPMD CUBE FILE.\n"
                    "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")

        #Discard comments at top of file.
        line = in_f.readline()
        newline = in_pqr.readline()
        while line.startswith('#'):
            line = in_f.readline()


        split_line = line.split()
        grid_sizes = [int(x)*-1 for x in split_line[-3:]]

        split_line = in_f.readline().split()

        origin = [float(x) for x in split_line[-3:]]

        parameter_fmt = "{:>4} {:>11.6f} {:>11.6f} {:>11.6f}\n"
        atom_num = 0
        while newline.startswith('REMARK'):
            newline = in_pqr.readline()

        try:
            while newline.startswith('ATOM') or newline.startswith('HETATM'):
                newline =  in_pqr.readline()
                new_split_line = newline.split()
                atom_num = new_split_line[1]
        except IndexError:
            pass
        in_pqr.seek(0)
        newline = in_pqr.readline()
        while newline.startswith('REMARK'):
            newline = in_pqr.readline()

        origin_line = parameter_fmt.format(atom_num, *origin)
        out_f.write(origin_line)


        for x in xrange(3):
            split_line = in_f.readline().split()
            grid_dims = [float(item) for item in split_line[-3:]]

            dim_lin = parameter_fmt.format(grid_sizes[x], *grid_dims)
            out_f.write(dim_lin)

        atoms_parameter_fmt = "{:>4} {:>11.6f} {:>11.6f} {:>11.6f} {:>11.6f}\n"
        a = True
        xreal_center = []
        yreal_center = []
        zreal_center = []
        try:
            while newline.startswith('ATOM') or newline.startswith('HETATM'):
                new_split_line = newline.split()
                radius = new_split_line[-1]
                xyz = new_split_line[-5:-2]
                line_atom_num = new_split_line[1]
                atom_radius = new_split_line[-1]
                pqr_lin = atoms_parameter_fmt.format(int(line_atom_num), float(new_split_line[-2]), float(xyz[0]), float(xyz[1]), float(xyz[2]))
                out_f.write(pqr_lin)
                newline = in_pqr.readline()
                xreal_center.append(float(xyz[0]))
                yreal_center.append(float(xyz[1]))
                zreal_center.append(float(xyz[2]))
        except IndexError:
            pass

        x_avg = sum(xreal_center)/float(atom_num)
        y_avg = sum(yreal_center)/float(atom_num)
        z_avg = sum(zreal_center)/float(atom_num)

        #print origin
        #new_origin = []
        #for item in origin:
        #    newitem = item/0.529177
            #new_new = item/2 + newitem/2
        #    new_origin.append(newitem)
        #print new_origin

        #Consume unneeded object lines.
        in_f.readline()
        in_f.readline()

        ##TODO: put atoms here

        value_format = ["{:< 13.5E}"]
        value_format = ' '.join(value_format * 6) + '\n'
        group = []
        line = in_f.readline()
        while not line.startswith('attribute'):
            values = [float(item) for item in line.split()]
            group.extend(values)

            if len(group) >= 6:
                out_f.write(value_format.format(*group))
                group = []

            line = in_f.readline()

        if group:
            group_strs = ["{:< 13.5E}".format(item) for item in group]
            out_f.write(' '.join(group_strs))

def checkprogress(jobid=None,appServicePort=None,calctype=None):
    """
        Finds out if the job has been completed
    """

    if have_opal:

        # construct soap request
        try:
            status=appServicePort.queryStatus(queryStatusRequest(jobid))
        except Exception, e:
            return ["error"]
        if status._code == 4:
            return ["error"]

        if status._code == 8:
            return ["complete",status]
        else:
            return ["running",status]

    else:
        progress = []
        file = open('%s%s%s/%s_status' % (INSTALLDIR,TMPDIR,jobid, form["calctype"].value))

        for line in file.readlines():
            progress.append(line.strip())
        file.close()
        return progress

def mainCGI():
    """
        Main method for determining the query page output
    """
    logopts = {}
    print "Content-type: text/html\n\n"
    calctype = form["calctype"].value

    # prints version error, if it exists
    if form["jobid"].value == 'False':
        print printheader("%s Job Status Page" % calctype.upper())
        progress = "version_mismatch"
        runtime = 0
    elif form["jobid"].value == 'notenoughmem':
        print printheader("%s Job Status Page" % calctype.upper())
        progress = "not_enough_memory"
        runtime = 0
    else:
        progress = None

    #Check for error html
    errorpath = '%s%s%s.html' % (INSTALLDIR, TMPDIR, form["jobid"].value)
    if os.path.isfile(errorpath):
        string = ""
        string+= "<html>\n"
        string+= "\t<head>\n"
        string+= "\t\t<meta http-equiv=\"Refresh\" content=\"0; url=%s%s%s.html\">\n" % (WEBSITE, TMPDIR, form["jobid"].value)
        string+= "\t</head>\n"
        string+= "</html>\n"
        print string
        return

    # prepares for Opal query, if necessary
    if have_opal:
        if calctype=="pdb2pqr":
            opal_url = PDB2PQR_OPAL_URL
        elif calctype=="apbs":
            opal_url = APBS_OPAL_URL
        appLocator = AppServiceLocator()
        appServicePort = appLocator.getAppServicePort(opal_url)
    else:
        appServicePort = None

    # if PDB2PQR, determines if link to APBS calculation should be shown
    if calctype=="pdb2pqr":
        #if(form["apbsinput"].value=="True"): # change to use a file
        #    apbs_input = True
        #else:
        #    apbs_input = False
        apbsInputFile = open('%s%s%s/apbs_input' % (INSTALLDIR, TMPDIR, form["jobid"].value))
        apbs_input = apbsInputFile.read()
        apbsInputFile.close()
        if apbs_input=="True":
            apbs_input = True
        else:
            apbs_input = False

        typemapInputFile = open('%s%s%s/typemap' % (INSTALLDIR, TMPDIR, form["jobid"].value))
        typemap = typemapInputFile.read()
        typemapInputFile.close()
        if typemap=="True":
            typemap = True
        else:
            typemap = False

    if have_opal and progress == None:
        if form["calctype"].value=="pdb2pqr":
            pdb2pqrJobIDFile = open('%s%s%s/pdb2pqr_opal_job_id' % (INSTALLDIR, TMPDIR, form["jobid"].value))
            jobid = pdb2pqrJobIDFile.read()
            pdb2pqrJobIDFile.close()
        elif form["calctype"].value=="apbs":
            apbsJobIDFile = open('%s%s%s/apbs_opal_job_id' % (INSTALLDIR, TMPDIR, form["jobid"].value))
            jobid = apbsJobIDFile.read()
            apbsJobIDFile.close()
    else:
        jobid = form["jobid"].value

    if progress == None:
        cp = checkprogress(jobid,appServicePort,calctype) # finds out status of job
        progress = cp[0]

    #initialize with bogus value just in case
    starttime = time.time()

    if progress == "running" or progress == "complete":
        timefile = open('%s%s%s/%s_start_time' % (INSTALLDIR, TMPDIR, form["jobid"].value, form["calctype"].value))
        starttime = float(timefile.read())
        timefile.close()

    if progress == "running" or (have_opal and progress not in ("version_mismatch",
                                                                "not_enough_memory",
                                                                "error",
                                                                "complete")):
        runtime = time.time()-starttime
        runtime = int(runtime)

    elif progress == "complete":
        endTimeFileString = '%s%s%s/%s_end_time' % (INSTALLDIR, TMPDIR, form["jobid"].value, form["calctype"].value)
        if have_opal and not os.path.isfile(endTimeFileString):
            runtime = time.time()-starttime
            with open(endTimeFileString, 'w') as endTimeFile:
                endTimeFile.write(str(time.time()))
        else:
            with open(endTimeFileString, 'r') as endTimeFile:
                runtime = float(endTimeFile.read())-starttime
    else:
        runtime = -1

    if progress == "running":
        #if have_opal:
        #    resultsurl = cp[1]._baseURL
        #else:
        if calctype=="pdb2pqr":
            resultsurl = '%squerystatus.cgi?jobid=%s&apbsinput=%s&calctype=pdb2pqr' % (WEBSITE, form["jobid"].value, apbs_input)
        else:
            resultsurl = '%squerystatus.cgi?jobid=%s&calctype=apbs' % (WEBSITE, form["jobid"].value)

    if progress == "complete":
        print printheader("%s Job Status Page" % calctype.upper(), jobid=form["jobid"].value)

    elif progress == "error":
        print printheader("%s Job Status Page - Error" % calctype.upper(),0, jobid=form["jobid"].value)

    elif progress == "running": # job is not complete, refresh in 30 seconds
        print printheader("%s Job Status Page" % calctype.upper(), refresh, jobid=form["jobid"].value)

    print "<BODY>\n<P>"
    print "<p></p>"
    print '<div id="content">'
    print "<h3>Status:"

    color = "FA3434"
    image = WEBSITE+"images/red_x.png"

    if progress == "complete":
        color = "2CDE56"
        image = WEBSITE+"images/green_check.png"
    elif progress == "running":
        color = "ffcc00"
        image = WEBSITE+"images/yellow_exclamation.png"

    print "<strong style=\"color:#%s;\">%s</strong>" % (color, progress)
    print "<img src=\"%s\"><br />" % image
    print "</h3>"
    print "Run time: " + str(timedelta(seconds=round(runtime))) + '<br />'
    print "Current time: %s<br />" % time.asctime()
    print "</P>\n<HR>\n<P>"

    if progress == "complete":
        if calctype=="pdb2pqr":
            nexturl = 'apbs_cgi.cgi?jobid=%s' % form["jobid"].value
        else:
            url_3dmol = 'visualize.cgi?jobid=%s&tool=%s' % (form["jobid"].value,'tool_3dmol')
            url_jmol = 'visualize.cgi?jobid=%s&tool=%s' % (form["jobid"].value,'tool_jmol')


        if have_opal:
            resp = appServicePort.getOutputs(getOutputsRequest(jobid))
            filelist = resp._outputFile

        print "Here are the results:<ul>"
        print "<li>Input files<ul>"

        if calctype=="pdb2pqr":
            # this code should be cleaned up once local PDB2PQR runs output the PDB file with the .pdb extension
            if have_opal:
                for i in range(0,len(filelist)):
                    file_name = filelist[i]._name
                    if ((len(file_name) == 4 and '.' not in file_name) or
                        (file_name.lower().endswith(".pdb") and "pdb2pka_output" not in file_name)):
                        print "<li><a href=%s>%s</a></li>" % (filelist[i]._url, filelist[i]._name)

                    if file_name.lower().endswith((".mol", ".mol2", ".names", ".dat")) and "pdb2pka_output" not in file_name:
                        print "<li><a href=%s>%s</a></li>" % (filelist[i]._url, filelist[i]._name)


            else:
                print "<li><a href=%s%s%s/%s.pdb>%s.pdb</a></li>" % (WEBSITE, TMPDIR, jobid, jobid, jobid)

        elif calctype=="apbs":
            if have_opal:
                for i in range(0,len(filelist)):
                    if filelist[i]._name == "apbsinput.in" or filelist[i]._name[-4:] == ".pqr":
                        print "<li><a href=%s>%s</a></li>" % (filelist[i]._url, filelist[i]._name)
            else:
                print "<li><a href=%s%s%s/apbsinput.in>apbsinput.in</a></li>" % (WEBSITE, TMPDIR, jobid)
                print "<li><a href=%s%s%s/%s.pqr>%s.pqr</a></li>" % (WEBSITE, TMPDIR, jobid, jobid, jobid)

        print "</ul></li>"
        print "<li>Output files<ul>"

        queryString = [str(os.environ["REMOTE_ADDR"])]
        # Getting PDB2PQR run log info
        if os.path.isfile('%s%s%s/pdb2pqr_log' % (INSTALLDIR, TMPDIR, jobid)):
            pdb2pqrLogFile=open('%s%s%s/pdb2pqr_log' % (INSTALLDIR, TMPDIR, jobid), 'r')
            logstr=pdb2pqrLogFile.read().split('\n')
            templogopts = eval(logstr[0])
            pdb2pqrLogFile.close()
            queryString.insert(0, templogopts.get('pdb',''))

        if calctype=="pdb2pqr":
            if have_opal:
                for i in range(0,len(filelist)):
                    if filelist[i]._name.endswith((".propka", "-typemap.html")):
                        print "<li><a href=%s>%s</a></li>" % (filelist[i]._url, filelist[i]._name)

                    if filelist[i]._name.endswith(".in") and "pdb2pka_output" not in filelist[i]._name:
                        print "<li><a href=%s>%s</a></li>" % (filelist[i]._url, filelist[i]._name)

                    if filelist[i]._name.endswith(".pqr") and not filelist[i]._name.endswith("background_input.pqr"):
                        print "<li><a href=%s>%s</a></li>" % (filelist[i]._url, filelist[i]._name)

                    #Get the first line of the summary file.
                    if filelist[i]._name.endswith(".summary"):
                        f=urllib.urlopen(filelist[i]._url)
                        summaryLine = f.readline().strip()
                        #logopts["pdb"]=logopts.get("pdb", "") + "|" + summaryLine
                        queryString.append(summaryLine)
                        f.close()
#                logRun(logopts, runtime, pqrOpalFileLength, logff, REMOTE_ADDR)
            else:
                #Get the first line of the summary file.
                summaryFile = '%s%s%s/%s%s' % (INSTALLDIR, TMPDIR, jobid, jobid, ".summary")
                if os.path.isfile(summaryFile):
                    with open(summaryFile) as f:
                        summaryLine = f.readline().strip()
                        #logopts["pdb"]=logopts.get("pdb", "") + "|" + summaryLine
                        queryString.append(summaryLine)

                outputfilelist = glob.glob('%s%s%s/*.propka' % (INSTALLDIR, TMPDIR, jobid))
                for i in range(0,len(outputfilelist)):
                    outputfilelist[i] = os.path.basename(outputfilelist[i])
                for extension in ["-typemap.html", ".pqr", ".in"]:
                    if extension != ".in" or apbs_input != False:
                        if extension == "-typemap.html" and typemap == False:
                            continue
                        outputfilelist.append('%s%s' % (jobid, extension))
                for outputfile in outputfilelist:
                    print "<li><a href=%s%s%s/%s>%s</a></li>" % (WEBSITE, TMPDIR, jobid, outputfile, outputfile)

            logopts['queryPDB2PQR'] = '|'.join(queryString)

                #for extension in ["-typemap.html", ".pqr", ".in"]:
                #    print "<li><a href=%s%s%s/%s%s>%s%s</a></li>" % (WEBSITE, TMPDIR, jobid, jobid, extension, jobid, extension)
        elif calctype=="apbs":
            if have_opal:
                for i in range(0,len(filelist)):
                    if filelist[i]._name[-3:]==".dx":
                        # compressing APBS OpenDX output files
                        currentpath = os.getcwd()
                        zipjobid = filelist[i]._name.split("-")[0]
                        dxfilename = '%s%s%s/%s' % (INSTALLDIR, TMPDIR, zipjobid, filelist[i]._name)
                        urllib.urlretrieve(filelist[i]._url, dxfilename)
                        os.chdir('%s%s%s' % (INSTALLDIR, TMPDIR, zipjobid))
                        # making both the dx file and the compressed file (.gz) available in the directory
                        syscommand = 'cp %s dxbkupfile' % (filelist[i]._name)
                        os.system(syscommand)
                        syscommand = 'gzip -9 ' + filelist[i]._name
                        os.system(syscommand)
                        syscommand = 'mv dxbkupfile %s' % (filelist[i]._name)
                        os.system(syscommand)
                        outputfilezip = filelist[i]._name + '.gz'

                        pqrfilename = '%s%s%s/%s.pqr' % (INSTALLDIR, TMPDIR, zipjobid, zipjobid)
                        cubefilename = '%s%s%s/%s.cube' % (INSTALLDIR, TMPDIR, zipjobid, zipjobid)

                        # making both the cube file and the compressed file (.gz) available in the directory
                        createcube(dxfilename, pqrfilename, cubefilename)
                        cubefilebasename = os.path.basename(cubefilename)

                        syscommand = 'cp %s cubebkupfile' % cubefilebasename
                        os.system(syscommand)
                        syscommand = 'gzip -9 ' + cubefilebasename
                        os.system(syscommand)
                        syscommand = 'mv cubebkupfile %s' % cubefilebasename
                        os.system(syscommand)
                        os.chdir(currentpath)
                        outputcubefilezip = cubefilebasename+".gz"

                        print "<li><a href=%s%s%s/%s>%s</a></li>" % (WEBSITE, TMPDIR, zipjobid, outputfilezip, outputfilezip)
                        print "<li><a href=%s%s%s/%s>%s</a></li>" % (WEBSITE, TMPDIR, zipjobid, outputcubefilezip, outputcubefilezip)

            else:
                outputfilelist = glob.glob('%s%s%s/%s-*.dx' % (INSTALLDIR, TMPDIR, jobid, jobid))
                for dxfile in outputfilelist:
                    # compressing APBS OpenDX output files
                    currentpath = os.getcwd()
                    workingpath = os.path.dirname(dxfile)
                    os.chdir(workingpath)
                    # making both the dx file and the compressed file (.gz) available in the directory
                    syscommand = 'cp %s dxbkupfile' % (os.path.basename(dxfile))
                    os.system(syscommand)
                    syscommand = 'gzip -9 ' + os.path.basename(dxfile)
                    os.system(syscommand)
                    syscommand = 'mv dxbkupfile %s' % (os.path.basename(dxfile))
                    os.system(syscommand)
                    os.chdir(currentpath)
                    outputfilezip = dxfile+".gz"



                    cubefilename = '%s%s%s/%s.cube' % (INSTALLDIR, TMPDIR, jobid, jobid)
                    pqrfilename = '%s%s%s/%s.pqr' % (INSTALLDIR, TMPDIR, jobid, jobid)


                    createcube(dxfile, pqrfilename, cubefilename)

                    print "<li><a href=%s%s%s/%s>%s</a></li>" % (WEBSITE, TMPDIR, jobid, os.path.basename(outputfilezip), os.path.basename(outputfilezip))

                outputcubefilelist = glob.glob('%s%s%s/%s.cube' % (INSTALLDIR, TMPDIR, jobid, jobid))
                for cubefile in outputcubefilelist:
                    # compressing cube output file
                    currentpath = os.getcwd()
                    os.chdir(workingpath)
                    # making both the cube file and the compressed file (.gz) available in the directory
                    syscommand = 'cp %s cubebkupfile' % (os.path.basename(cubefile))
                    os.system(syscommand)
                    syscommand = 'gzip -9 ' + os.path.basename(cubefile)
                    os.system(syscommand)
                    syscommand = 'mv cubebkupfile %s' % (os.path.basename(cubefile))
                    os.system(syscommand)
                    os.chdir(currentpath)
                    outputcubefilezip = cubefile+".gz"

                    print "<li><a href=%s%s%s/%s>%s</a></li>" % (WEBSITE, TMPDIR, jobid, os.path.basename(outputcubefilezip), os.path.basename(outputcubefilezip))

            logopts['queryAPBS'] = '|'.join(queryString)

        if calctype=="pdb2pqr":
            if have_opal:
                outputfilelist = []
                for i in range(0,len(filelist)):
                    if filelist[i]._name.endswith((".DAT", ".txt")):
                        outputfilelist.append((filelist[i]._url, os.path.basename(filelist[i]._name)))
                        #print "<li><a href=%s>%s</a></li>" % (filelist[i]._url, filelist[i]._name)
                if outputfilelist:
                    print "</ul></li>"
                    print "<li>PDB2PKA files<ul>"
                    for outputfile in outputfilelist:
                        print "<li><a href=%s>%s</a></li>" % (outputfile[0], outputfile[1])
            else:
                outputfilelist = glob.glob('%s%s%s/pdb2pka_output/*.DAT' % (INSTALLDIR, TMPDIR, jobid))
                outputfilelist.extend(glob.glob('%s%s%s/pdb2pka_output/*.txt' % (INSTALLDIR, TMPDIR, jobid)))
                outputfilelist = [os.path.basename(outputfile) for outputfile in outputfilelist]
                if outputfilelist:
                    print "</ul></li>"
                    print "<li>PDB2PKA files<ul>"
                    for outputfile in outputfilelist:
                        print "<li><a href=%s%s%s/pdb2pka_output/%s>%s</a></li>" % (WEBSITE, TMPDIR, jobid, outputfile, outputfile)

        print "</ul></li>"
        print "<li>Runtime and debugging information<ul>"

        if have_opal:
            stdouturl = resp._stdOut
            stderrurl = resp._stdErr
        else:
            stdouturl = "%s%s%s/%s_stdout.txt" % (WEBSITE, TMPDIR, jobid, calctype)
            stderrurl = "%s%s%s/%s_stderr.txt" % (WEBSITE, TMPDIR, jobid, calctype)

        print "<li><a href=%s>Program output (stdout)</a></li>" % stdouturl
        print "<li><a href=%s>Program errors and warnings (stderr)</a></li>" % stderrurl


        print "</ul></li></ul>"


        #if have_opal:
        #    resp = appServicePort.getOutputs(getOutputsRequest(jobid))
        #    for opalfile in resp._outputFile:
        #        if opalfile._name[-8:]!="-input.p":
        #            print "<li><a href=%s>%s</a></li>" % (opalfile._url, opalfile._name)
        #    print "<li><a href=%s>Standard output</a></li>" % (resp._stdOut)
        #    print "<li><a href=%s>Standard error</a></li>" % (resp._stdErr)
        #else:
        #    for line in cp[1:]:
        #        line = os.path.basename(line)
        #        if line[-8:]!="-input.p":
        #            if line[-11:]=="_stdout.txt":
        #                printname = "Standard output"
        #            elif line[-11:]=="_stderr.txt":
        #                printname = "Standard error"
        #            else:
        #                printname = line
        #            print "<li><a href=%s>%s</a></li>" % (WEBSITE+TMPDIR+jobid+"/"+line,printname)

        if calctype=="pdb2pqr" and apbs_input and HAVE_APBS:
            print "</ul></p><hr><p><b><a href=%s>Click here</a> to run APBS with your results.</b></p>" % nexturl
        elif calctype=="apbs":
            #print "</ul></p><hr><p><b><a href=%s>Click here</a> to visualize your results in Jmol.</b></p>" % nexturl
            print "</ul></p><hr><p><b>Visualize your results online:"
            print "<ul> <li><a href=%s>3Dmol</a></li><li><a href=%s>Jmol</a></li></ul>" % (url_3dmol, url_jmol)

    elif progress == "error":
        print "There was an error with your query request. This page will not refresh."

        print "</ul></li>"
        print "<li>Runtime and debugging information<ul>"

        if have_opal:
            resp = appServicePort.getOutputs(getOutputsRequest(jobid))
            stdouturl = resp._stdOut
            stderrurl = resp._stdErr

        else:
            stdouturl = "%s%s%s/%s_stdout.txt" % (WEBSITE, TMPDIR, jobid, calctype)
            stderrurl = "%s%s%s/%s_stderr.txt" % (WEBSITE, TMPDIR, jobid, calctype)

        print "<li><a href=%s>Program output (stdout)</a></li>" % stdouturl
        print "<li><a href=%s>Program errors and warnings (stderr)</a></li>" % stderrurl

        print "</ul></li></ul>"

        if have_opal:
            print " <br />If your job has been running for a prolonged period of time and failed with no reason listed in the standard out or standard error, then the job probably timed out and was terminated by the system.<br />"

        print '<br />If you are having trouble running PDB2PQR on the webserver, please download the <a href="http://www.poissonboltzmann.org/docs/downloads/">command line version of PDB2PQR</a> and run the job from there.'



    elif progress == "running":
        print "Page will refresh in %d seconds<br />" % refresh
        print "<HR>"

        if not have_opal:
            print "</ul></li>"
            print "<li>Runtime and debugging information<ul>"
            stdouturl = "%s%s%s/%s_stdout.txt" % (WEBSITE, TMPDIR, jobid, calctype)
            stderrurl = "%s%s%s/%s_stderr.txt" % (WEBSITE, TMPDIR, jobid, calctype)
            print "<li><a href=%s>Program output (stdout)</a></li>" % stdouturl
            print "<li><a href=%s>Program errors and warnings (stderr)</a></li>" % stderrurl
            print "</ul></li></ul>"

        print "<small>Your results will appear at <a href=%s>this page</a>. If you want, you can bookmark it and come back later (note: results are only stored for approximately 12-24 hours).</small>" % resultsurl
    elif progress == "version_mismatch":
        print "The versions of APBS on the local server and on the Opal server do not match, so the calculation could not be completed"

    print "</P>"
    print "<script type=\"text/javascript\">"
    for key in logopts:
        print getEventTrackingString('queryData', key, logopts[key]),
        #print "_gaq.push(['_trackPageview', '/main_cgi/has_%s_%s.html']);" % (key, logopts[key])
    print "</script>"
    print "</div> <!--end content div-->"
    print "</BODY>"
    print "</HTML>"

if __name__ == "__main__" and os.environ.has_key("REQUEST_METHOD"):
    """ Determine if called from command line or CGI """
    refresh=30

    if not form.has_key("jobid") and form["calctype"].value=="pdb2pqr":
        print "Content-type: text/html\n\n"
        print printheader("PDB2PQR Job Status - Error")
        text="<BODY>\n"
        text+="\t<H2>Missing jobid field</H2>\n"
        text+="\t<P>Your request url is missing the jobid field</P>\n"
        text += "<script type=\"text/javascript\">"
        text += "var gaJsHost = ((\"https:\" == document.location.protocol) ? \"https://ssl.\" : \"http://www.\");"
        text += "document.write(unescape(\"%3Cscript src=\'\" + gaJsHost + \"google-analytics.com/ga.js\' type=\'text/javascript\'%3E%3C/script%3E\"));"
        text += "</script>"
        text += "<script type=\"text/javascript\">"
        text += "try {"
        text += "var pageTracker = _gat._getTracker(\"UA-11026338-3\");"
        text += "pageTracker._trackPageview();"
        text += "} catch(err) {}</script>"
        text+="</BODY>\n</HTML>"
        print text
        sys.exit(2)


    if (form["calctype"].value=="pdb2pqr" and HAVE_PDB2PQR_OPAL) or (form["calctype"].value=="apbs" and APBS_OPAL_URL!=""):
        have_opal = True
        from AppService_client import AppServiceLocator, queryStatusRequest, getOutputsRequest
    else:
        have_opal = False

    mainCGI()
