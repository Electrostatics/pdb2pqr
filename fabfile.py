from fabric.api import *
import tarfile
import os
import zipfile

import fabfile_settings

env.hosts = []

if hasattr(fabfile_settings, 'osx_host'):
    osx_host = fabfile_settings.osx_host
    env.hosts.append(osx_host)
    if '@' in osx_host:
        osx_host = osx_host.split('@')[1]
    if ':' in osx_host:
        osx_host = osx_host.split(':')[0]
else:
    osx_host = ''

if hasattr(fabfile_settings, 'linux_host'):
    linux_host = fabfile_settings.linux_host
    env.hosts.append(linux_host)
    if '@' in linux_host:
        linux_host = linux_host.split('@')[1]
    if ':' in linux_host:
        linux_host = linux_host.split(':')[0]
else:
    linux_host = ''

if hasattr(fabfile_settings, 'run_tests'):
    run_tests = fabfile_settings.run_tests
else:
    run_tests=False

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
    if env.host == linux_host:
        tar.add('apbs_libs/linux/apbslib.py','pdb2pka/apbslib.py')
        tar.add('apbs_libs/linux/_apbslib.so','pdb2pka/_apbslib.so')
    elif env.host == osx_host:
        tar.add('apbs_libs/osx/apbslib.py','pdb2pka/apbslib.py')
        tar.add('apbs_libs/osx/_apbslib.so','pdb2pka/_apbslib.so')
    tar.close()

@runs_once
def pack_for_nbcr():
    create_dist_folder()
    tar = start_src_tar('pdb2pqr-src-nbcr-'+pv+'.tar.gz', 'pdb2pqr-src-'+pv)
    tar.add('apbs_libs/linux/apbslib.py','pdb2pka/apbslib.py')
    tar.add('apbs_libs/linux/_apbslib.so','pdb2pka/_apbslib.so')
    tar.close()
    local("move pdb2pqr-src-nbcr-"+pv+'.tar.gz dist_files\\')

def start_src_tar(name='pdb2pqr.tgz', prefix=None):
    file_list = local('git ls-tree -r --name-only HEAD', capture=True).split('\n')
    tar = TarWrapper(name, prefix)
    for f in file_list:
        tar.add(f)
    return tar

@runs_once
def create_dist_folder():
    with settings(warn_only=True):
        local('mkdir dist_files')

@runs_once
def pack_for_ditro():
    create_dist_folder()
    tar = start_src_tar('pdb2pqr-src-'+pv+'.tar.gz', 'pdb2pqr-src-'+pv)
    tar.close()
    local('copy Changelog.md "dist_files\PDB2PQR-' + pv + '-ReleaseNotes.txt"')
    local("move pdb2pqr-src-"+pv+'.tar.gz dist_files\\')

def deploy():
    python = 'python2.7'

    put('pdb2pqr.tgz', '~/')
    with settings(warn_only=True):
        run('rm -rf tmp')
        run('mkdir tmp')

    with cd('~/tmp/'):
        run('tar -zxvf ~/pdb2pqr.tgz')
        run(python+' scons/scons.py')

        if run_tests:
            run(python+' scons/scons.py -j 4 pdb2pka-test')
            run(python+' scons/scons.py -j 4 complete-test')


def install_on_deployed():
    python = 'python2.7'

    with settings(warn_only=True):
        with cd('~/www/pdb2pqr'):
            run('rm -rf *')

    with cd('~/tmp/'):
        configopts = ''

        if hasattr(fabfile_settings, 'APBS'):
            configopts += ' APBS='+fabfile_settings.APBS

        if hasattr(fabfile_settings, 'URL'):
            configopts += ' URL='+fabfile_settings.URL

        if hasattr(fabfile_settings, 'PREFIX'):
            configopts += ' PREFIX='+fabfile_settings.PREFIX

        if hasattr(fabfile_settings, 'OPAL'):
            configopts += ' OPAL='+fabfile_settings.OPAL

        if hasattr(fabfile_settings, 'APBS_OPAL'):
            configopts += ' APBS_OPAL='+fabfile_settings.APBS_OPAL

#         if True:
#             configopts += ' URL=http://PT24098/d3k084/pdb2pqr_opal'
#             configopts += ' PREFIX=/Users/d3k084/www/pdb2pqr_opal/'
#             configopts += ' APBS_OPAL=http://nbcr-222.ucsd.edu/opal2/services/apbs_1.3'
#             configopts += ' OPAL=http://nbcr-222.ucsd.edu/opal2/services/pdb2pqr_2.0.0'
#         else:
#             configopts += ' PREFIX=/Users/d3k084/www/pdb2pqr/'
#             configopts += ' URL=http://PT24098/d3k084/pdb2pqr'

        run(python+' scons/scons.py ' + configopts)
        run(python+' scons/scons.py install ' + configopts)

def build_binary_from_deploy():
    create_dist_folder()
    with cd('~/tmp/'):
        os_string = 'NOT_SET_FIX_ME'
        if env.host == linux_host:
            os_string = 'linux'

        if env.host == osx_host:
            os_string = 'osx'

        run('pyinstaller pdb2pqr.spec')

        name = 'pdb2pqr-' + os_string + '-bin64-' + pv
        run('mv dist/pdb2pqr dist/' + name)
        with cd('dist/'):
            run('tar -zcvf ' + name + '.tar.gz ' + name)
            if run_tests:
                with cd(name):
                    run('./pdb2pqr --ff=parse --verbose --ligand=examples/ligands/LIG_1ABF.mol2 1ABF 1ABF.pqr')
                    run('./pdb2pqr --with-ph=7.0 --ph-calc-method=pdb2pka --ff=parse --verbose 1a1p 1a1p.pqr')

        get("~/tmp/dist/*.tar.gz","dist_files/")

def linux_bin_cross_platform_test():
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
def build_windows():
    local(r'copy apbs_libs\windows\* pdb2pka\\')
    local(r'python scons\scons.py -c')
    local(r'python scons\scons.py')
    if run_tests:
        local(r'python scons\scons.py -j 4 pdb2pka-test')
        local(r'python scons\scons.py -j 7 complete-test')

    build_windows_binary()
    if run_tests:
        test_windows_binary()

@runs_once
def test_windows_binary():
    name = 'pdb2pqr-windows-bin64-' + pv
    with lcd(r'dist\\' + name):
        local(r'pdb2pqr --ff=parse --verbose --ligand=examples\ligands\LIG_1ABF.mol2 1ABF 1ABF.pqr')
        local('pdb2pqr --with-ph=7.0 --ph-calc-method=pdb2pka --ff=parse --verbose 1a1p 1a1p.pqr')

@runs_once
def build_windows_binary():
    with settings(warn_only=True):
        local('del /Q /F dist')

    local('pyinstaller pdb2pqr.spec')
    name = 'pdb2pqr-windows-bin64-' + pv
    local(r'move /Y dist\pdb2pqr dist\\' + name)

    zip_file = zipfile.ZipFile(name + '.zip', 'w', zipfile.ZIP_DEFLATED)
    for root, _, files in os.walk(r'dist\\'+name):
        new_root = root.split('/', 1)[-1]
        print root
        for f in files:
            zip_file_path = os.path.join(new_root,f)
            real_file_path = os.path.join(root,f)
            zip_file.write(real_file_path, zip_file_path)
            print 'Zipping ' + zip_file_path
    zip_file.close()
    create_dist_folder()
    local('move /Y ' + name + '.zip' + r' dist_files\\' + name + '.zip')


def build_all_tarballs():
    pack_for_ditro()
    pack_for_nbcr()

def build_all_binaries():
    pack()
    deploy()
    build_binary_from_deploy()
    build_windows()

def deploy_and_install():
    pack()
    deploy()
    install_on_deployed()

