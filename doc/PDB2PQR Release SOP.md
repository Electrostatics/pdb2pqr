# PDB2PQR Release SOP

## Preface
This document explains the Standard Operating Procedure (SOP) to release a new version of PDB2PQR.
This document is designed to make the process go more smoothly.
It should NOT be viewed as a strict Methodology but rather as a list of suggestions and procedures that should be done or accounted for before a new release tar ball is posted on SourceForge (or GitHub) and announced to the pdb2pqr-announce mailing list and on http://www.poissonboltzmann.org/.

## Repository

### Branch
When preparing for a new release, the trunk (or tag if a bug fix version needs to be released) should be branched and/or tagged.
It should be named something obvious such as pdb2pqr-<code>new version number</code>.
All testing should be done on a clean checkout of the branch.
Fixes found during testing should be checked into the branch and merged into the trunk if appropriate.
Merging fixes back into a tag is not required and probably not appropriate.

In the file <code>site_scons/defaults.py</code> change the value of <code>productVersion</code> to reflect the release version.

## Fabric
All of the steps in the Build, Binary builds and Automated Testing can be done automatically for all target platforms by adapting fabric_settings.py and using the fabric.py script found in the pdb2pqr directory. Documentation for the python fabric module and how to use the fabric script can be found at http://docs.fabfile.org/

Build the documentation.

Use the command "fab build_all_binaries" on a Windows machine to build and test pdb2pqr on all target platforms.

If all tests pass and all binaries are built then restore fabric_settings.py to it's default state. Run "fab build_all_tarballs"

After these steps all 3 binary tarballs, the source tarball and the source tarball for NBCR will all be in "dist_file" directory.

## Build
The branch should be tested to see if it configures and compiles correctly on all the target operating systems.
Currently those systems are Linux, Mac OSX 10.6 or newer, and Windows 7.
Testing against the GNU compilers for pdb2pka is sufficient.
Remember, these tests should be done on a fresh checkout the branch.

To configure and make branch:
See INSTALL.md for how to build and install pdb2pqr.

## Binary builds
See INSTALL.md for how to create binary builds. Please note the specific needs of each platform.
For release we build the OSX version on 10.6 (64 bit only) and the Linux version on RHEL or CentOS 6 (32 and 64 bit) to ensure compatibility with RHEL and Unbuntu LTS. We do not support RHEL 5 or earlier.

On each supported platform the build should be tested with the following commands from <code>pdb2pqr/dist/pdb2pqr</code>:

	./pdb2pqr.exe --ff=parse --ph-calc_method=propka --verbose --ligand=examples/ligands/LIG_1ABF.mol2 1abf 1abf.pqr

	./pdb2pqr.exe --ff=parse --ph-calc_method=pdb2pka --verbose 1a1p 1a1p.pqr

If the program does not crash in either case then we can release this binary build.

## Testing

### Automated Testing
After a successful build automated tests should be run.
Run <code>python scons/scons.py test</code>, <code>python scons/scons.py advtest</code>, <code>python scons/scons.py -j7 complete-test</code>, and <code>python scons/scons.py -j4 pdb2pka-test</code>. (Change the -j parameters according to the number of cores available in the testing machine.)
This will test out the most basic case, a more complex case with a pH setting and ligand file, and then test every option in various combinations.
All three tests should run without issue unless changes were made that affect the output of the program (Version numbers and some atom ordering differences are ignored when comparing expected and actual results.).

### Web Site Testing
Like command-line testing, every available option on the web front end should be tested.
All options should also be tested with both opal and non-opal configurations.
Builds should be tested with APBS support turned on and off but every permutation of APBS support and Opal is not required.
Edit <code>build_config.py</code> to configure Opal, Non-Opal, and APBS support builds. Instructions are included in the build_config.py file.
A properly installed and compiled copy of APBS is needed for both Opal and non Opal tests.
Incompatible options and invalid files should be tried to ensure that errors are reported correctly by the web interface.

## Documentation

Documentation should be rebuilt and the newly generated files and checked in before creating the final tar ball.
The documentation must be built only after <code>python scons/scons.py -c</code> (to clean up the codebase) has been run.
Any extra files created by the build will cause extra documentation pages for standard python modules to be created. (This already happens, but after configure it is worse.)
Also it will create documentation pages for configuration py files such as aconf.py. This should be avoided.
Documentation is created by running the genpydoc.sh script in the doc directory.
After generating the new documentation the updated files should be checked in to the branch.
This should be the last thing done before packing up a tar ball so as to reflect any changes to the documentation made during testing.

## Build (Again)
If testing revealed any needed code changes, a fresh clone should be made with another build and automated tests.
The main purpose of this is to ensure that you did not forget to check anything in and everything needed for deployment is in the branch.
Build and test the binary version. One platform should be sufficient.

## Deployment

### Tagging
The branch should now be tagged with an appropriate name. (pdb2pqr-<code>new version number</code>).

### Packaging
After testing is complete, documentation updated, and another build is done it is time to package it up.

Use the last git clone created.
Tar and gzip with the following command from the repos base directory:

	git archive --format=tar --prefix=pdb2pqr-<new version number>/ pdb2pqr-<new version number>:pdb2pqr/ | gzip >git-<new version number>.tar.gz

Build the binaries on all supported targets.

Rename the pdb2pqr folder to

	pdb2pqr-<osx|linux|windows>-bin-<version>

From the dist folder create an archive like so:

#### Linux and OSX

	tar zcvf pdb2pqr-<osx|linux>-bin-<version>.tar.gz pdb2pqr-<osx|linux>-bin-<version>

#### Windows

Compress the

	pdb2pqr-windows-bin-<version>

folder into a zip file.

### Release
The tar and zip files should be uploaded to Source Forge and github.

#### Source Forge
A new folder should be created to reflect the version number.
The new tar ball should be uploaded to that folder.
A copy of ChangeLog.md (which you should be keeping up to date!) renamed to <code>PDB2PQR-<version>-ReleaseNotes.txt</code> (example PDB2PQR-1.7.1a-ReleasNotes.txt) should also be uploaded to the same folder.

#### Github
Follow these direction to create a release on Github
https://github.com/blog/1547-release-your-software

### Announcement
New releases should be announced on pdb2pqr-announce@lists.sourceforge.net and http://www.poissonboltzmann.org/ with the ChangeLog text.

