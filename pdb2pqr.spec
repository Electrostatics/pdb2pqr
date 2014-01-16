# -*- mode: python -*-
a = Analysis(['pdb2pqr.py'],
             pathex=['C:\\MinGW\\msys\\1.0\\home\\D3Y382\\workspaces\\pdb2pqr\\branches\\monson-dev'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='pdb2pqr.exe',
          debug=False,
          strip=None,
          upx=True,
          console=True )
coll = COLLECT(exe, Tree('dat/', 'dat'),
			   Tree('doc/', 'doc'),
			   Tree('examples/', 'examples'),
			   Tree('extensions/', 'extensions'),
			   [('propka30/Source/protein_bonds.dat', 'propka30/Source/protein_bonds.dat', 'DATA'),
			    ('propka30/Source/ions.list', 'propka30/Source/ions.list', 'DATA')],
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=True,
               name='pdb2pqr')
