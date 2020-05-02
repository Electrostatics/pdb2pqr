#!/usr/bin/env python

#
# Take a mol2 file as input, convert to SMILES, search pka_lig_tool database and produce
# mol2 files of the reference state (all neutral) + mol2 files (with charges + radii) for each charged state + pKa values for the transitions
#
# This script relies on PDB2PQR and OpenBABEL
#
# Copyright (c) Jens Erik Nielsen, University College Dublin, 2011
#

# TODO - eliminate sys and os imports
import sys.exc_info, argv, sys_path
import os
try:
    file_name=__file__
    if file_name[:2]=='./':
        scriptpath=os.getcwd()
    else:
        scriptpath=os.path.join(os.getcwd(),os.path.split(file_name)[0])
        if scriptpath[-1] == "/":
            scriptpath=scriptpath[:-1]
except:
    scriptpath=os.path.split(argv[0])[0]
    if scriptpath=='.':
        scriptpath=os.getcwd()
import logging


_LOGGER = logging.getLogger(__name__)


#
# Path to pKa_lig_tool server
#
server = 'enzyme.ucd.ie'
url = 'http://%s/cgi-bin/pka_lig_tool/download_pka.cgi?%s'
#
# -----
#
OPENBABEL='/software/bin/babel'
import os
if os.environ['OSTYPE']=='darwin':
    OPENBABEL='/sw/bin/babel'
#
if not os.path.isfile(OPENBABEL):
    raise Exception('Configure OPENBABEL in preapre_pKa_ligand.py')
#
# Add to import path
#
pdb2pqr_path=os.path.split(scriptpath)[0]
sys_path.append(pdb2pqr_path)
        
class ligand_pKa:

    def __init__(self,mol2lines):
        #
        # Manage the process
        #
        mol2obj=self.read_mol2(mol2lines)
        #
        SMILES=self.mol2_2_SMILES(mol2lines)
        #
        self.pKa_objects=self.search_pka_ligtool(SMILES)
        #
        #ref_state=self.get_refstate(mol2obj)
        #
        #states=[['reference',ref_state]]
        #for pKa_object in pKa_objects:
        #    states.append(self.getstate(pKa_object,ref_state))
        #
        return 
        
    #
    # -----
    #
        
    def read_mol2(self,mol2lines):
        """Parse the mol2 lines"""
        import StringIO, string
        mol2fileobj=StringIO.StringIO(string.join(mol2lines))
        #
        # Use the mol2 parser in pdb2pqr
        #
        import src.pdb
        mol2object=src.pdb.Mol2Molecule()
        mol2object.read(mol2fileobj)
        return mol2object
        
    #
    # -----
    #
        
    def mol2_2_SMILES(self,mol2lines):
        """Convert the mol2 file to a SMILES string"""
        import tempfile
        dirname=tempfile.mkdtemp()
        mol2file=os.path.join(dirname,'ligand.mol2')
        fd=open(mol2file,'w')
        for line in mol2lines:
            fd.write(line)
        fd.close()
        # 
        # Do the conversion using Babel
        #
        hydfile=os.path.join(dirname,'noHs.mol2')
        smi_file = os.path.join(dirname, 'ligand.smi')
        command='%s -imol2 %s -d -omol2 %s' %(OPENBABEL,mol2file,hydfile)
        os.system(command)
        #
        com = '%s -imol2 %s -osmi %s >/dev/null'%(OPENBABEL,hydfile, smi_file)
        os.system(com)
        #
        fd=open(smi_file,'r')
        lines=fd.readlines()
        SMILES =lines[0].split()[0]
        #
        # Delete the tempdir
        #
        import shutil
        shutil.rmtree(dirname)
        #
        return SMILES

    #
    # -----
    #

    def search_pka_ligtool(self,smiles):
        """Search the pka_lig_tool database for a ligand match"""
        import StringIO, urllib
    
        #
        # Get the XML data from the server
        #
        args=urllib.urlencode({'smiles':smiles})
        thisurl= url %(server,args)
        f=urllib.urlopen(thisurl)
        text=f.read()
        output = StringIO.StringIO(text)
        _LOGGER.info(text)
        _LOGGER.info(smiles)
        return text
        #
        # Parse the XML
        #
        from xml.dom import minidom
        xmldoc = minidom.parse(output)

        ligands = xmldoc.firstChild
        _LOGGER.info('Search-type was: %s'%ligands.attributes['Type'].value)
        for ligand in ligands.childNodes:
            atoms = ligand.getElementsByTagName('Atoms')[0]
            mol2 = ligand.getElementsByTagName('mol2')[0]

            _LOGGER.info('*'*50)
            _LOGGER.info('Found ligand \"%s\"'%ligand.attributes['Name'].value)
            _LOGGER.info("")
            _LOGGER.info('mol2 file')
            _LOGGER.info('-'*50)
            _LOGGER.info(mol2.firstChild.data)
            _LOGGER.info("")
            _LOGGER.info(' Atoms and associated pKa values')
            _LOGGER.info('-'*50)
            for atom in atoms.childNodes:
                _LOGGER.info('%4s %4s %-6s'%(atom.attributes['Name'].value,
                                     atom.attributes['Number'].value,
                                     atom.attributes['Type'].value))
                pkas = atom.firstChild
                if pkas:
                    for pka in pkas.childNodes:
                        _LOGGER.info('          Value:            ',pka.attributes['Value'].value)
                        _LOGGER.info('          Temperature:      ',pka.attributes['Temperature'].value)
                        _LOGGER.info('          pH:               ',pka.attributes['pH'].value)
                        _LOGGER.info('          Solvent:          ',pka.attributes['Solvent'].value)
                        _LOGGER.info('          Salt type:        ',pka.attributes['Salt_type'].value)
                        _LOGGER.info('          Salt conc.:       ',pka.attributes['Salt_conc.'].value)
                        _LOGGER.info('          Titratable group: ',pka.attributes['Titratable_group'].value)
                        _LOGGER.info('          Most bio. rel.:   ',pka.attributes['Most_bio._relavent'].value)
                        _LOGGER.info('          Reference:        ',pka.attributes['Reference'].value)
                        _LOGGER.info('          Comment:          ',pka.attributes['Comment'].value)

            _LOGGER.info('*'*50)
        return

    def get_allhyd_state(self):
    
        return
        
    def get_state(self,allhyd_state,remove_hydrogen,add_hydrogen):
        
        return
    


if __name__=='__main__':
    _LOGGER.info("")
    _LOGGER.info('Get pKa values and structures of protonation states for a ligand')
    _LOGGER.info('Chresten Soendergaard, Paul Czodrowski, Jens Erik Nielsen 2006-2010')
    _LOGGER.info("")
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <file>',version='%prog 1.0')
    parser.add_option('-m','--mol2',dest='mol2',action='store',type='string',
                        help='Mol2 file holding the ligand coordinates. Default: %s', default='5fiv_ligand.mol2')
    parser.add_option('--testall',dest='testall',action='store_true',default=False,
                        help='Test SMILES searcher with all mol2 files in docking set. Default: %default')
    (options, args) = parser.parse_args()
    if not options.testall:
        #
        # Read the mol2 lines
        #
        fd=open(options.mol2)
        mol2lines=fd.readlines()
        fd.close()
        #
        # Call the class
        #
        I=ligand_pKa(mol2lines)
    else:
        mol2dir='/Users/nielsen/lib/development/Experimental_data/dGBind/'
        import os
        dirs=os.listdir(mol2dir)
        failed=[]
        ok=[]
        for dirname in dirs:
            rdir=os.path.join(mol2dir,dirname)
            if os.path.isdir(rdir):
                # This dir contains a mol2file
                files=os.listdir(rdir)
                mol2file=None
                for file in files:
                    if file.split('.')[-1]=='mol2':
                        mol2file=os.path.join(rdir,file)
                        break
                #
                if mol2file:
                    _LOGGER.info(mol2file)
                    fd=open(mol2file)
                    mol2lines=fd.readlines()
                    fd.close()
                    #
                    # Call the class
                    #
                    try:
                        I=ligand_pKa(mol2lines)
                        ok.append()
                    except:
                        failed.append([mol2file,sys.exc_info()[0]])
                        _LOGGER.error('FAILED')
                        _LOGGER.error(sys.exc_info()[0])
        _LOGGER.warning(failed)
        _LOGGER.info('OK',len(ok))
        _LOGGER.warning('FAILED',len(failed))
                        

    
    