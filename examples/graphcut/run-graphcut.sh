#!/bin/bash -x

#
#  Put everything in an output directory
#
PDB_OUT=out
mkdir ${PDB_OUT}
syspython="/usr/local/Cellar/python/2.7.9/Frameworks/Python.framework/Versions/2.7/bin/python2.7"


for pdb in 1P9I 4IL7 1DF4 2GUS 1YZM 4G3O 4CVD 3HRO 3V1A 1KFN 1ZVA 3LAA 2CMP 1ZX6 4ZMK 2REY 1JQ0 1ZLM 3CP1 2IGD 1USE 2IWN 3NGP 1AIL 2GZV 3LLB 1YIB 3DVI 3IDW 4CFI 2P5K 4NPN 1NH9 4Q2Q 3KZD 2XXC 1YQB 1SJV 2VWR 1L2P 2CJJ 2H2C 3QE1 5EE2 3LMO 1DUP 1R5Q 2QBV 2NSN 2CKX 4AXT 1OX3 2FI9 2SFA 2J6B 1Z0P 3ONJ 4PGR 2IVY 4O7Q 1LN4 4I2T 1TIG 4POY 2J9V 2OO2 4N6T 1R7J 3QC7 4S11 1OPC 2PMR 2FWG 3ZSL 2NWD 1RIS 1RJ1 4YPC 2D8E 4GQM 1J2A 3US6 3KT2 3T1S; do
	echo "******** ${pdb} starting"

        #
        #  Run the python graph cut code
        #
	#pkadir=../titration\ curves/pdb2pka/${pdb}
	#graphdir=../titration\ curves/graphcut/${pdb}
        #mkdir "${graphdir}"
	#python main.py "${pkadir}" "${graphdir}" | tee "${graphdir}/stdout.txt"

        #
        #  Run pdb2pqr graph cut code
        #
        PDB2PKA_OUT=${PDB_OUT}/${pdb}
        mkdir ${PDB2PKA_OUT}
	${syspython} "../../pdb2pqr.py" "--verbose" \
                "--titration-state-method=pdb2pka" "--ff=parse" \
                "--pdb2pka-out=${PDB2PKA_OUT}" \
                "${pdb}" "${PDB2PKA_OUT}/$pdb.pqr" | tee "${PDB2PKA_OUT}/stdout.txt"
                #>& "${graphdir}/stdout.txt" \
                # | tee "${graphdir}/stdout.txt"
	echo "******** ${pdb} done"
done
