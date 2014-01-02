import distutils.sysconfig
import os
import sys
import platform
from datetime import datetime

defaultPrefix = os.path.expanduser('~/pdb2pqr/')
defaultURL = 'http://' + platform.node() + '/pdb2pqr/'
defaultMaxAtoms = 10000
codePath = os.getcwd()

pythonBin = sys.executable

buildTime = datetime.today()
productVersion = 'Monson-Dev BRANCH '# + buildTime.strftime('%Y%m%d%H%M')

defaultOpalURL = 'http://nbcr-222.ucsd.edu/opal2/services/pdb2pqr_1.8'
defaultAPBSOpalURL = 'http://nbcr-222.ucsd.edu/opal2/services/apbs_1.3'