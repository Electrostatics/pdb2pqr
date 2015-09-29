from fabric.api import *
import tarfile
import os
import zipfile

env.hosts = ['yourhost.blah.org']

import sys, os
sys.path.append("site_scons")
from defaults import productVersion

pv = productVersion.replace(' ', '_')

class TarWrapper():
    def __init__(self, name, prefix=None):
        self.prefix = prefix
        self.tar = tarfile.open(name, 'w:gz')

    def add(self, name, arcname=None):
        if self.prefix is not None:
            if arcname is not None:
                arcname = os.path.join(self.prefix, arcname)
            else:
                arcname = os.path.join(self.prefix, name)

        self.tar.add(name, arcname)
        print('Packing file: ' + name if arcname is None else arcname)

    def close(self):
        self.tar.close()


def pack():
    tar = start_src_tar()
    if '127.0.0.1' in env.host or 'rocks-86' in env.host:
        tar.add('apbs_libs/linux/apbslib.py','pdb2pka/apbslib.py')
        tar.add('apbs_libs/linux/_apbslib.so','pdb2pka/_apbslib.so')
    if 'PT24098' in env.host:
        tar.add('apbs_libs/osx/apbslib.py','pdb2pka/apbslib.py')
        tar.add('apbs_libs/osx/_apbslib.so','pdb2pka/_apbslib.so')
    tar.close()

def pack_for_nadya():
    tar = start_src_tar('pdb2pqr-src-'+pv+'.tar.gz', 'pdb2pqr-src-'+pv)
    tar.add('apbs_libs/linux/apbslib.py','pdb2pka/apbslib.py')
    tar.add('apbs_libs/linux/_apbslib.so','pdb2pka/_apbslib.so')
    tar.close()

def start_src_tar(name='pdb2pqr.tgz', prefix=None):
    file_list = local('git ls-tree -r --name-only HEAD', capture=True).split('\n')
    tar = TarWrapper(name, prefix)
    for f in file_list:
        tar.add(f)
    return tar

@runs_once
def misc():
    with settings(warn_only=True):
        local('mkdir dist_files')
    tar = start_src_tar('pdb2pqr-src-'+pv+'.tar.gz', 'pdb2pqr-src-'+pv)
    tar.close()
    local('copy Changelog.md "dist_files\PDB2PQR-' + pv + '-ReleaseNotes.txt"')
    local("move pdb2pqr-src-"+pv+'.tar.gz dist_files/')

def deploy(complete_test='F'):
    with settings(warn_only=True):
        if 'rocks-86' in env.host:
            python = '/opt/python/bin/python'
        else:
            python = 'python2.7'

    put('pdb2pqr.tgz', '~/')
    with settings(warn_only=True):
        run('rm -rf tmp')
        run('mkdir tmp')

    with cd('~/tmp/'):
        with settings(warn_only=True):
            run('rm -rf *')
        run('tar -zxvf ~/pdb2pqr.tgz')
        run(python+' scons/scons.py')

        if complete_test.lower().startswith('t'):
            run(python+' scons/scons.py pdb2pka-test')
            run(python+' scons/scons.py -j 4 complete-test')


def install_on_deployed(opal='F'):
    with settings(warn_only=True):
        if 'rocks-86' in env.host:
            python = '/opt/python/bin/python'
            with cd('/share/apps/pdb2pqr_1.9'):
                run('rm -rf *')
        else:
            python = 'python2.7'
            if opal.lower().startswith('t'):
                with cd('~/www/pdb2pqr_opal'):
                    run('rm -rf *')
            else:
                with cd('~/www/pdb2pqr'):
                    run('rm -rf *')

    with cd('~/tmp/'):
        configopts = ''
        if 'PT24098' in env.host:
            configopts += ' APBS=/Users/d3y382/bin/apbs'
            if opal.lower().startswith('t'):
                configopts += ' URL=http://PT24098/d3k084/pdb2pqr_opal'
                configopts += ' PREFIX=/Users/d3k084/www/pdb2pqr_opal/'
                configopts += ' APBS_OPAL=http://nbcr-222.ucsd.edu/opal2/services/apbs_1.3'
                configopts += ' OPAL=http://nbcr-222.ucsd.edu/opal2/services/pdb2pqr_2.0.0'
            else:
                configopts += ' PREFIX=/Users/d3k084/www/pdb2pqr/'
                configopts += ' URL=http://PT24098/d3k084/pdb2pqr'

        elif 'rocks-86' in env.host:
            configopts += ' URL=http://rocks-86.sdsc.edu/pdb2pqr'
            configopts += ' PREFIX=/share/apps/pdb2pqr_1.9/'
            configopts += ' APBS=/opt/apbs/bin/apbs'
            configopts += ' APBS_OPAL=http://nbcr-222.ucsd.edu/opal2/services/apbs_1.3'
            configopts += ' OPAL=http://rocks-86.sdsc.edu/opal2/services/pdb2pqr_2.0.0'

        run(python+' scons/scons.py ' + configopts)
        run(python+' scons/scons.py install ' + configopts)

def build_binary_from_deploy(complete_test='F'):
    with cd('~/tmp/'):
        os_string = 'NOT_SET_FIX_ME'
        if '127.0.0.1' in env.host:
            os_string = 'linux'

        if 'PT24098' in env.host:
            os_string = 'osx'

        run('pyinstaller pdb2pqr.spec')

        name = 'pdb2pqr-' + os_string + '-bin64-' + pv
        run('mv dist/pdb2pqr dist/' + name)
        with cd('dist/'):
            run('tar -zcvf ' + name + '.tar.gz ' + name)
            if complete_test.lower().startswith('t'):
                with cd(name):
                    run('./pdb2pqr --ff=parse --verbose --ligand=examples/ligands/LIG_1ABF.mol2 1ABF 1ABF.pqr')
                    run('./pdb2pqr --with-ph=7.0 --ph-calc-method=pdb2pka --ff=parse --verbose 1a1p 1a1p.pqr')

        get("~/tmp/dist/*.tar.gz","dist_files/")

def linux_test():
    '''
    Push the linux bin to a host and test it.
    '''
    os_string = 'linux'
    name = 'pdb2pqr-' + os_string + '-bin64-' + pv
    put('dist_files/' + name + '.tar.gz', '~/')
    run('tar -zxvf ~/' + name + '.tar.gz')
    with cd(name):
        run('./pdb2pqr --ff=parse --verbose --ligand=examples/ligands/LIG_1ABF.mol2 1ABF 1ABF.pqr')
        run('./pdb2pqr --with-ph=7.0 --ph-calc-method=pdb2pka --ff=parse --verbose 1a1p 1a1p.pqr')

    run('rm -rf '+name)

@runs_once
def build_windows(complete_test = 'F'):
    local('cp apbs_libs/windows/* pdb2pka/')
    local('python27-64 scons/scons.py -c')
    local('python27-64 scons/scons.py')
    if complete_test.lower().startswith('t'):
        local('python27-64 scons/scons.py pdb2pka-test')
        local('python27-64 scons/scons.py -j 7 complete-test')

    build_windows_binary()
    if complete_test.lower().startswith('t'):
        test_windows_binary()

@runs_once
def test_windows_binary():
    name = 'pdb2pqr-windows-bin64-' + pv
    with lcd('dist/' + name):
        local('pdb2pqr --ff=parse --verbose --ligand=examples/ligands/LIG_1ABF.mol2 1ABF 1ABF.pqr')
        local('pdb2pqr --with-ph=7.0 --ph-calc-method=pdb2pka --ff=parse --verbose 1a1p 1a1p.pqr')

@runs_once
def build_windows_binary():
    with settings(warn_only=True):
        local('rm -rf dist')

    local('c:\Python27_64\Scripts\pyinstaller pdb2pqr.spec')
    name = 'pdb2pqr-windows-bin64-' + pv
    local('mv dist/pdb2pqr dist/' + name)

    zip_file = zipfile.ZipFile(name + '.zip', 'w', zipfile.ZIP_DEFLATED)
    for root, _, files in os.walk('dist/'+name):
        new_root = root.split('/', 1)[-1]
        print root
        for f in files:
            zip_file_path = os.path.join(new_root,f)
            real_file_path = os.path.join(root,f)
            zip_file.write(real_file_path, zip_file_path)
            print 'Zipping ' + zip_file_path
    zip_file.close()
    local('mv ' + name + '.zip' + ' dist_files/' + name + '.zip')


def doall(complete_test='F'):
    #misc()
    pack()
    deploy(complete_test)
    #windows(complete_test)
    #zip_windows_binary()

