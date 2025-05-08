#!/bin/bash

# Usage: ./test_pHs_pkaani_propka.sh [input_pdb]

if [ $# -ne 1 ]; then
    echo "Usage: $0 input.pdb"
    exit 1
fi

pdb_file="$1"
base_name=$(basename "$pdb_file" .pdb)

# count residues in input PDB
# only works if there are no ANISOU records, only ATOM records in the PDB
asps=$(grep "CA  ASP" "$pdb_file" | wc -l)
glus=$(grep "CA  GLU" "$pdb_file" | wc -l)
his_s=$(grep "CA  HIS" "$pdb_file" | wc -l)
lys=$(grep "CA  LYS" "$pdb_file" | wc -l)
args=$(grep "CA  ARG" "$pdb_file" | wc -l)
tyrs=$(grep "CA  TYR" "$pdb_file" | wc -l)

# Run PROPKA
pdb2pqr30 --ffout AMBER --titration-state-method propka --with-ph 5.0 "$pdb_file" "${base_name}_propka.pqr"

ashs=$(grep "CA  ASH" "${base_name}_propka.pqr" | wc -l)
glhs=$(grep "CA  GLH" "${base_name}_propka.pqr" | wc -l)
hips=$(grep "CA  HIP" "${base_name}_propka.pqr" | wc -l)
pqr_lys=$(grep "CA  LYS" "${base_name}_propka.pqr" | wc -l)
pqr_arg=$(grep "CA  ARG" "${base_name}_propka.pqr" | wc -l)
pqr_tyr=$(grep "CA  TYR" "${base_name}_propka.pqr" | wc -l)

echo "PROPKA RESULTS ${pdb_file}"
echo "---------------------------------------"
echo "Aspartates in pdb file: (none should be protonated after pdb2pqr)"
echo $asps 
echo "Protonated aspartates in pqr file:"
echo $ashs
echo "Glutamates in pdb file: (none should not be protonated after pdb2pqr)"
echo $glus
echo "Protonated glutamates in pqr file:"
echo $glhs
echo "Histidines in pdb file: (all should be protonated after pdb2pqr)"
echo $his_s
echo "Protonated histidines in pqr file:"
echo $hips
echo "protonated Lysines in pdb file: (all should be protonated after pdb2pqr)"
echo $lys
echo "Lysines in pqr file:"
echo $pqr_lys
echo "protonated Arginines in pdb file: (all should be protonated after pdb2pqr)"
echo $args
echo "Arginines in pqr file:"
echo $pqr_arg
echo "protonated Tyrosines in pdb file: (all should be protonated after pdb2pqr)"
echo $tyrs
echo "Tyrosines in pqr file:"
echo $pqr_tyr
echo "---------------------------------------"
printf "\n"

# Run PKAANI
pdb2pqr30 --ffout AMBER --titration-state-method pkaani --with-ph 5.0 "$pdb_file" "${base_name}_pkaani.pqr"

ashs=$(grep "CA  ASH" "${base_name}_pkaani.pqr" | wc -l)
glhs=$(grep "CA  GLH" "${base_name}_pkaani.pqr" | wc -l)
hips=$(grep "CA  HIP" "${base_name}_pkaani.pqr" | wc -l)
pqr_lys=$(grep "CA  LYS" "${base_name}_pkaani.pqr" | wc -l)
pqr_arg=$(grep "CA  ARG" "${base_name}_pkaani.pqr" | wc -l)
pqr_tyr=$(grep "CA  TYR" "${base_name}_pkaani.pqr" | wc -l)

echo "PKAANI RESULTS ${pdb_file}"
echo "---------------------------------------"
echo "Aspartates in pdb file: (most should not be protonated after pdb2pqr)"
echo $asps 
echo "Protonated aspartates in pqr file:"
echo $ashs
echo "Glutamates in pdb file: (most should not be protonated after pdb2pqr)"
echo $glus
echo "Protonated glutamates in pqr file:"
echo $glhs
echo "Histidines in pdb file: (all should be protonated after pdb2pqr)"
echo $his_s
echo "Protonated histidines in pqr file:"
echo $hips
echo "protonated Lysines in pdb file: (none should be deprotonated after pdb2pqr)"
echo $lys
echo "Lysines in pqr file:"
echo $pqr_lys
echo "protonated Arginines in pdb file: (none should be deprotonated after pdb2pqr)"
echo $args
echo "Arginines in pqr file:"
echo $pqr_arg
echo "protonated Tyrosines in pdb file: (none should be deprotonated after pdb2pqr)"
echo $tyrs
echo "Tyrosines in pqr file:"
echo $pqr_tyr
echo "---------------------------------------"
printf "\n"
