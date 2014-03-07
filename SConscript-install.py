Import('env')
Import('compile_targets')

from os.path import dirname, basename, join, isfile
from os import walk
from glob import iglob

compile_target_files = [str(x[0]) for x in compile_targets]

def getAllFiles(root_dir, python_only=False, recursive=True):
    if recursive:
        for root, dirnames, filenames in walk(root_dir):
            for filename in filenames:
                if not python_only or filename.endswith('.py'):
                	yield join(root, filename)
                
    else:
        for filename in iglob(join(root_dir,'*')):
            if not isfile(filename):
                continue
			
            if not python_only or filename.endswith('.py'):
                yield filename
	

def installFile(file_name, build_target='install'):	
	#Don't do compile targets
	if file_name in compile_target_files:
		return
	#ignore pyc and pyo. We get them from the py build.
	if file_name.endswith('.pyc') or  file_name.endswith('.pyo'):
		return
		
	if basename(file_name).startswith('.'):
		return
	
	target = File(file_name)
	if file_name.endswith('.py'):
		compiled_name = file_name+'c'
		result = env.Command(compiled_name, file_name, CompilePython('$TARGET', '$SOURCE'))
		if GetOption("clean"):
			env.Default(result)
		Alias(build_target, env.Install(env['PREFIX']+dirname(file_name), result))
	else:
		Alias(build_target, env.Install(env['PREFIX']+dirname(file_name), target))

#Contrib
for file_name in getAllFiles('contrib/ZSI-2.1-a1/ZSI'):
	installFile(file_name)

installFile('contrib/ZSI-2.1-a1/Copyright')
	
#ProPKA	
for file_name in getAllFiles('propka30/'):
	installFile(file_name)

#Whole directories
for dir_name in ('dat/', 'doc/', 'examples/', 'jmol/'):
    for file_name in getAllFiles(dir_name):
        installFile(file_name)
	#dat = Dir(dir_name)
	#Alias('install', env.Install(env['PREFIX'], dat))

#Compiled Targets
for target in compile_targets:
	file_name = str(target[0])
	if file_name.endswith('.py'):
		compiled_name = file_name+'c'
		result = env.Command(compiled_name, file_name, CompilePython('$TARGET', '$SOURCE'))
		if GetOption("clean"):
			env.Default(result)
		Alias('install', env.Install(env['PREFIX']+dirname(file_name), result))
	else:
		Alias('install', env.Install(env['PREFIX']+dirname(file_name), target))

#PDB2PKA	
for file_name in getAllFiles('pdb2pka/'):
	installFile(file_name)
	
Alias('install', env.Install(env['PREFIX']+'pdb2pka/', 'pdb2pka/pka.py'))
	
#Main Program
for dir_name in ('src/',
                 'extensions/'):
	for file_name in getAllFiles(dir_name, python_only=True):
		installFile(file_name)
		
Alias('install', env.Install(env['PREFIX'], 'pdb2pqr.py'))
Alias('install', env.Install(env['PREFIX'], 'pdb2pqr.cgi'))
Alias('install', env.InstallAs(env['PREFIX']+'/index.html', 'html/server.html'))

for file_name in ('AppService_client.py',
				  'AppService_services.py',
				  'AppService_services_types.py',
				  'AppService_types.py',
				  'AUTHORS',
				  'ChangeLog.md',
				  'COPYING',
				  'main.py',
				  'main_cgi.py',
				  'NEWS',
				  'pdb2pqr.css',
				  'README.md'):
	installFile(file_name)
	
Alias('install', env.Command(env['PREFIX']+'tmp', None,
							[Mkdir('$TARGET'),
			 				Chmod('$TARGET', 0777)]))

#Testing stuff to test installed pdb2pqr	
for dir_name in ('tests/',
                 'scons/',
                 'site_scons/'):
	for file_name in getAllFiles(dir_name):
		installFile(file_name, build_target='install-tests')
		
Alias('install-tests', env.Install(env['PREFIX']+'scons/', 'scons/scons.py'))
Alias('install-tests', env.Install(env['PREFIX']+'site_scons/', 'site_scons/site_init.py'))

installFile('SConstruct', build_target='install-tests')

Alias('install-tests', 'install')
	
#Depends('install', DEFAULT_TARGETS)
#Requires('install', DEFAULT_TARGETS)

