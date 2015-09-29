#Deploying and testing the PDB2PQR website

## On your remote host:

1. Get access to Kyle's machine: pt24098
  * create public www directory

2. edit /etc/apache2/httpd.conf
  * add paths for links in the /Library directory

3. create symbolic links in the /Library directory

## On you local machine

1. in your browser, go to pt24098/d3k084/
  * should work: list a directory

2. clone apbs-pdb2pqr

3. cd pdb2pqr

4. Get the APBS libs
  * unpack apbs libs in pdb2pqr dir (creates apbs_libs folder)

5. edit the fabric script (fabfile.py) in the pdb2pqr folder to find your
files on the remote machine.
  * pack command in script makes tar 
  * deploy command logs into remote machine (pt24098)
  * install command pushes into the www directory on remote machine

## Deploy the site for testing

1. run "fab doall"

2. run "fab install_on_deployed"

Note:  3dmol stuff: querrystatus.py and 3dmol directory


