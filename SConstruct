from os.path import isfile
import os

#Are we running from the installed version?
if isfile('pdb2pqr.py.in'):
    SConscript('SConscript-main.py')
else:
    #Allow test targets to run in installed version.
    env = Environment()
    if os.name == 'nt':
        env.Append(ENV={'PATH' : os.environ['PATH']})
    Export('env')  
    pdb2pqr = [File('pdb2pqr.py')]
    Export('pdb2pqr')   
    SConscript('tests/SConscript')
