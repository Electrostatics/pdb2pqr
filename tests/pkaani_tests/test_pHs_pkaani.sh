#!/bin/bash

# pkaani test code.
# count number of residues at a protonation state before pdb2pqr at a pH of 1.1.
# All the Asp, Glu, and His residues should be protonated, none the LYS, ARG, and TYR
# residues should be deprotonated. Ensuring that this is the case

asps=$(cat 1brs.pdb | grep "CA  ASP" | wc -l)
glus=$(cat 1brs.pdb | grep "CA  GLU" | wc -l)
his_s=$(cat 1brs.pdb | grep "CA  HIS" | wc -l)
lys=$(cat 1brs.pdb | grep "CA  LYS" | wc -l)
args=$(cat 1brs.pdb | grep "CA  ARG" | wc -l)
tyrs=$(cat 1brs.pdb | grep "CA  TYR" | wc -l)
pdb2pqr30 --ffout AMBER --titration-state-method pkaani --with-ph 1.1 1brs.pdb 1brs.pqr > pkaani_out.txt

ashs=$(cat 1brs.pqr | grep "CA  ASH" | wc -l)
glhs=$(cat 1brs.pqr | grep "CA  GLH" | wc -l)
hips=$(cat 1brs.pqr | grep "CA  HIP" | wc -l)
pqr_lys=$(cat 1brs.pqr | grep "CA  LYS" | wc -l)
pqr_arg=$(cat 1brs.pqr | grep "CA  ARG" | wc -l)
pqr_tyr=$(cat 1brs.pqr | grep "CA  TYR" | wc -l)



echo "Aspartates in pdb file: (all should be protonated after pdb2pqr)"
echo $asps 
echo "Protonated aspartates in pqr file:"
echo $ashs

echo "Glutamates in pdb file: (all should be protonated after pdb2pqr)"
echo $glus
echo "Protonated glutamates in pqr file:"
echo $glhs

echo "Histidines in pdb file: (all should be protonated after pdb2pqr)"
echo $his_s
echo "Protonated histidines in pqr file:"
echo $hips

echo "Lysines in pdb file: (none should be deprotonated after pdb2pqr)"
echo $lys
echo "Lysines in pqr file:"
echo $pqr_lys

echo "Arginines in pdb file: (none should be deprotonated after pdb2pqr)"
echo $args
echo "Arginines in pqr file:"
echo $pqr_arg

echo "Tyrosines in pdb file: (none should be deprotonated after pdb2pqr)"
echo $tyrs
echo "Tyrosines in pqr file:"
echo $pqr_tyr


