#!/bin/bash

#**************************************************************
# This script will run a series of the same DROP command,
# setDihedral, to manipulate an extended Villin Headpiece
# protein into its native folded structure. See README.md
# for more information.
#
# Usage: $ bash fold-villin.sh /path/to/compiled/drop
# Ex if drop binary is in parent directory of git repository: 
# $ bash fold-villin.sh ../../
#
# Author: Brian Andrews
# Last Date Modified: 6/28/2023
#**************************************************************

path=$1
declare -a backbone_angles=( phi psi )

touch movie.xyz

#start with residue 2 because gmx rama only outputs residues with both backbone dihedrals defined
residue=2
while read line
do
  i=0
  for angle in ${line}
  do
    echo ${angle} ${backbone_angles[i]}

    ${path}./drop -f setDihedral -i villin-unfolded-conect.pdb -d ${backbone_angles[i]} -a ${angle} -n ${residue} -o output.xyz
    cat output.xyz >> movie.xyz

    i=$((i+1))
  done

residue=$((residue+1))
done < villin_backbone_ramachandran_angles.txt

