# DROP
Dihedral Rotation Of Proteins (DROP)

This personal project inspired by my research during my PhD studies at Drexel University with Professor Brigita Urbanc. At Drexel, I study polypeptide/protein dynamics using Molecular Dynamics (MD) software; i.e. GROMACs, VMD, etc. One of the pain points I have recognized during my studies that I have not found a sufficient solution for was the manual manipulation of proteins and molecules without painful manual calculation or file editing to deliberately explore configuration spaces. Of particular interest in my research is dihedral angles, which are directly associated with higher order structures in proteins. So the inception of this project manifested as an ability to computationally manipulate a polypeptide as if it were built as a physical molecular model with a chemistry set from school.

This project is still in its inception phase. But will serve as an exercise in many other facets of software development as well.

- Learn about parallelization (in C, although it should translate well to other languages)
- Create a makefile to make this a useful command line tool
- Learn to implement command line flags (before I have just kept track of command line option orders which is not good practice)
- Use Git more directly to simulate a team environment in software development
- Develop an environment where addition of functions and methods is an streamlined process

Currently (7/20/2021), the main executable is Steric_Rama_MPI.c and includes compilation and running commands at the top as comments. This will be changed. But, for now, the compiled executable will read in a pdb file of a tripeptide, identify the dihedrals, and rotate about the two dihedral angles of the central amino acid residue, and print out a "trajectory" of these rotations in an xyz file which currently needs to be named manually by the user in Steric_Rama_MPI.c before compilation. Example videos of triglycine and trialanine are in the movies directory.
