# -*- mode: python -*-
a = Analysis(['pdb2pqr.py'],
             pathex=None,
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
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
coll = COLLECT(exe, Tree('dat/', 'dat'),
			   Tree('doc/', 'doc'),
			   Tree('examples/', 'examples'),
			   Tree('extensions/', 'extensions'),
			   [('propka30/Source/protein_bonds.dat', 'propka30/Source/protein_bonds.dat', 'DATA'),
			    ('propka30/Source/ions.list', 'propka30/Source/ions.list', 'DATA'),
				('BINARY_README.md', 'BINARY_README.md', 'DATA')],
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=True,
               name='pdb2pqr')
