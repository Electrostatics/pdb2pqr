import sys
import os
import getopt
import shutil
import numpy as np

from urllib.request import urlopen
from io import StringIO
from pkaani.pkaani import calculate_pka
from pkaani.prep_pdb import prep_pdb

def usage_pkaani():
    """
    Show how to use this program!
    """
    print("""
Example usages:

  * If PDB file doesnt exist, it is downloaded and prepared for pKa calculations.

      pkaani -i 1BNZ 
      pkaani -i 1BNZ.pdb

  * Multiple files can be given as inputs

      pkaani -i 1BNZ,1E8L

  * If a specific directory is wanted:

      pkaani -i path_to_file/1BNZ
      pkaani -i path_to_file/1BNZ,path_to_file/1E8L


  Arguments: -i: Input files. Inputs can be given with or without 
                 file extension (.pdb). If PDB file is under a 
                 specific directory (or will be downloaded) the path
                 can also be given as path_to_file/PDBFILE. 
                 Multiple PDB files can be given 
                 by using "," as separator (i.e. pkaani -i 1BNZ,1E8L).
""")

def handle_arguments_pkaani():
    inp_file = None
    prep_files = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:", ["help","inp="])
    except getopt.GetoptError:
        usage_pkaani()
        sys.exit(-1)
		
    for opt, arg in opts:
	

        if opt in ('-h', "--help"):
            usage_pkaani()
            sys.exit(-1)
		
        elif opt in ("-i", "--inp"):
            inp_file=[x.strip() for x in arg.split(',')]
			
        else:
            assert False, usage_pkaani()

    if inp_file is None:
        usage_pkaani()
        sys.exit(-1)
    # Input PDB file is mandatory!
    if len(inp_file)==0: # is None:
        print("@> ERROR: A PDB file is mandatory!")
        usage_pkaani()
        sys.exit(-1)

    return inp_file

def main():

    input_files = handle_arguments_pkaani()   
    pdbfiles=np.array(input_files)
  
    #first prepare PDB files for pkaani 
    for inputpdb in pdbfiles:
        pdbid=inputpdb.rsplit('.', 1)[0]
        pdbfile=pdbid+".pdb"
        file_exist=True
        if not os.path.exists(pdbfile):
            file_exist=False
            base=os.path.basename(pdbfile)
            dpdbid=base.rsplit('.', 1)[0]

            print("File %s is not accessible" % pdbfile)
            print("Downloading : http://www.rcsb.org/pdb/files/%s.pdb" % dpdbid)
        
            url = 'http://www.rcsb.org/pdb/files/%s.pdb' % dpdbid
            
            file = urlopen(url)
            contents = file.read().decode('utf-8')
            file.close()
            file = StringIO(contents)

            outfile=pdbfile
            
            with open(outfile, 'w') as f2:
              for line in contents:
                f2.write(line)
                
        prep_pdb(pdbfile)            
        file_exist=True
                            
    #CALCULATER PKA
    pkadict=calculate_pka(pdbfiles,writefile=True)

    #RENAME FILES PROPERLY
    for inputpdb in pdbfiles:
        pdbid=inputpdb.rsplit('.', 1)[0]
        pdbfile=pdbid+".pdb"

        if os.path.exists(pdbfile):
           oldf=pdbfile
           newf=pdbid+"_pkaani.pdb"
           os.rename(oldf,newf)
           if file_exist:
              oldf=pdbid+"_0.pdb"
              newf=pdbfile
              os.rename(pdbid+"_0.pdb",pdbfile)



if __name__ == "__main__":
    main()
    
    
