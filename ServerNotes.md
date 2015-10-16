#Deploying and testing the PDB2PQR website

These directions are taylored for installation on OSX. Adjustments will have to be made for installation on Linux. On OSX specifically you must ensure that web sharing is enabled for the user who will be running pdb2pqr.

## On your remote host:

1. Create www directory to install pdb2pqr in. (We will use /Users/kyle/www for all examples)

2. Create a link to this directory from /Library/WebServer/Documents
  
    sudo ln -s /Users/kyle/www/ /Library/WebServer/Documents/kyle  
 

2. edit /etc/apache2/httpd.conf to allow cgi scipts to execute in the installation directory.

  * add paths for links in the /Library directory in the <IfModule alias_module> section 
    
    ScriptAliasMatch ^/cgi-bin/((?!(?i:webobjects)).*$) "/Library/WebServer/Documents/kyle/pdb2pqr/$1"
	
  * add a Directory section
  
    <Directory "/Library/WebServer/Documents/kyle/pdb2pqr">
      AllowOverride None
      Options ExecCGI
      Order allow,deny
      Allow from all
    </Directory>

## On you local machine

1. Install fabric using pip if you haven't already.

1. in your browser, go to remote_ip/user-name/
  * should work: list a directory

2. clone apbs-pdb2pqr

3. cd pdb2pqr

4. Get the APBS libs
  
Currently the APBS libs are not easily compiled on any of the target platforms. These libraries can be obtained from the prebuilt executables from pdb2pqr 2.0 onward.

Internally we pass around a zip file that has the files alreay to go for each target platform.

5. edit the fabric script (fabfile.py) in the pdb2pqr folder to find your files on the remote machine. 

Sepcifically where the scons command line is built in the fucntion "install_on_deployed" update URL and PREFIX.

In our example:

    ' URL=http://my_host_or_ip/kyle/pdb2pqr'
	' PREFIX=/Users/kyle/www/pdb2pqr/'
	
In the future the fabric script will be configurable in a non-terrible way.

The script does the following steps:
  * pack command in script makes tar 
  * deploy command logs into remote machine 
  * install command pushes into the www directory on remote machine

## Deploy the site for testing

1. run "fab doall"

2. run "fab install_on_deployed"

Note:  3dmol stuff: querrystatus.py and 3dmol directory


