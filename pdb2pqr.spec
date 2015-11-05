# -*- mode: python -*-

from os.path import dirname

a = Analysis(['pdb2pqr.py'],
             pathex=None,
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)

# remove lib2to3
targets = []
for module in a.pure:
    if module[0].startswith('lib2to3'):
        targets.append(module)

for target in targets:
    a.pure.remove(target)

pyz = PYZ(a.pure)

import os
if os.name == 'nt':
    exe_name = 'pdb2pqr.exe'
else:
    exe_name = 'pdb2pqr'

exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name=exe_name,
          debug=False,
          strip=None,
          upx=True,
          console=True )

#Include lib2to3 by hand. Needed for networkx 1.10.
import lib2to3
lib23_path = dirname(lib2to3.__file__)

coll = COLLECT(exe, Tree('dat/', 'dat'),
			   Tree('doc/', 'doc'),
			   Tree('examples/', 'examples'),
			   Tree('extensions/', 'extensions'),
			   Tree(lib23_path, 'lib2to3'),
			   [('propka30/Source/protein_bonds.dat', 'propka30/Source/protein_bonds.dat', 'DATA'),
			    ('propka30/Source/ions.list', 'propka30/Source/ions.list', 'DATA'),
			    ('pdb2pka/TITRATION.DAT', 'pdb2pka/TITRATION.DAT', 'DATA'),
				('BINARY_README.md', 'BINARY_README.md', 'DATA')],
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=True,
               name='pdb2pqr')
