#!/bin/bash

#*************************************************************
# This script runs a series of gromacs commands using an
# input pdb/gro/gromacs compatible file to generate conect
# records for DROP to read.
#
# Run: $ bash generate-conect-gromacs.sh input-file output-file
#
# Author: Brian Andrews
# Last Date Modified: 6/28/2023
#*************************************************************

input=$1
output=$2 #use .pdb

#pick any force field and water model
gmx pdb2gmx -f ${input} -o pep.gro -p pep.top -ignh

#need a box to avoid errors
gmx editconf -f pep.gro -o pep_box.gro -bt cubic -box 4

#need to generate a tpr file - maxwarn is just to fend off warnings about charge and energy cutoffs, which don't matter here
gmx grompp -f em.mdp -c pep_box.gro -p pep.top -o output.tpr -maxwarn 10

#print out structure with conect records
#choose your selection - should only really be Protein
gmx trjconv -f pep_box.gro -s output.tpr -o ${output} -conect yes

#add REMARK at the beginning of pdb file
echo "REMARK This file was generated using generate-conect-gromacs.sh from https://github.com/andrewsb8/DROP in scripts/." > /tmp/file.tmp
cat ${output} >> /tmp/file.tmp
mv /tmp/file.tmp ${output}

#remove extra files from gromacs
rm *.tpr *.gro *.top *.itp mdout.mdp
