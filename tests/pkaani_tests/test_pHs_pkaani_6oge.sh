#!/bin/bash

# pkaani test code.
# count number of residues at a protonation state before pdb2pqr at a pH of 5.0
# All the Asp, Glu, and His residues should be protonated, none the LYS, ARG, and TYR
# residues should be deprotonated. Ensuring that this is the case

asps=$(cat 6oge_de.pdb | grep "CA  ASP" | wc -l)
glus=$(cat 6oge_de.pdb | grep "CA  GLU" | wc -l)
his_s=$(cat 6oge_de.pdb | grep "CA  HIS" | wc -l)
lys=$(cat 6oge_de.pdb | grep "CA  LYS" | wc -l)
args=$(cat 6oge_de.pdb | grep "CA  ARG" | wc -l)
tyrs=$(cat 6oge_de.pdb | grep "CA  TYR" | wc -l)
pdb2pqr30 --ffout AMBER --titration-state-method propka --with-ph 5.0 6oge_de.pdb 6oge_de_propka.pqr 

ashs=$(cat 6oge_de_propka.pqr | grep "CA  ASH" | wc -l)
glhs=$(cat 6oge_de_propka.pqr | grep "CA  GLH" | wc -l)
hips=$(cat 6oge_de_propka.pqr | grep "CA  HIP" | wc -l)
pqr_lys=$(cat 6oge_de_propka.pqr | grep "CA  LYS" | wc -l)
pqr_arg=$(cat 6oge_de_propka.pqr | grep "CA  ARG" | wc -l)
pqr_tyr=$(cat 6oge_de_propka.pqr | grep "CA  TYR" | wc -l)


echo "PROPKA RESULTS"
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

pdb2pqr30 --ffout AMBER --titration-state-method pkaani --with-ph 5.0 6oge_de.pdb 6oge_de_pkaani.pqr 

ashs=$(cat 6oge_de_pkaani.pqr | grep "CA  ASH" | wc -l)
glhs=$(cat 6oge_de_pkaani.pqr | grep "CA  GLH" | wc -l)
hips=$(cat 6oge_de_pkaani.pqr | grep "CA  HIP" | wc -l)
pqr_lys=$(cat 6oge_de_pkaani.pqr | grep "CA  LYS" | wc -l)
pqr_arg=$(cat 6oge_de_pkaani.pqr | grep "CA  ARG" | wc -l)
pqr_tyr=$(cat 6oge_de_pkaani.pqr | grep "CA  TYR" | wc -l)

echo "PKAANI RESULTS"
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
