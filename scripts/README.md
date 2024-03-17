This directory contains scripts which can be useful when using DROP.

It includes:
- ```generate-conect-gromacs.sh```: Adds hydrogens and generates CONECT records for an input pdb file using a series of [GROMACS](https://gitlab.com/gromacs/gromacs) commands. GROMACS must be installed to use this script
  - ```em.mdp```: file required for ```generate-conect-gromacs.sh``` to work. Must be in same directory as ```generate-conect-gromacs.sh``` 
  - ```Villin-Headpiece/``` includes an example input (```villin-unfolded.pdb```) and output (```villin-unfolded-conect.pdb```) of this script
