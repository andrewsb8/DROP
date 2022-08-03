# DROP
Dihedral Rotation Of Proteins (DROP)

This personal project inspired by my research during my PhD studies at Drexel University with Professor Brigita Urbanc. At Drexel, I study polypeptide/protein dynamics using Molecular Dynamics (MD) software; i.e. GROMACs, VMD, etc. One of the pain points I have recognized during my studies that I have not found a sufficient solution for was the manual manipulation of proteins and molecules without painful manual calculation or file editing to deliberately explore configuration spaces. Of particular interest in my research are dihedral angles, which are directly associated with higher order structures in proteins. So the inception of this project manifested as an ability to computationally manipulate a polypeptide as if it were built as a physical molecular model with a chemistry set from school.

This project is still in its inception phase. But will serve as an exercise in many other facets of software development as well. Some examples of these exercises are the following, where 'x' is completed and 'p' is progress:

- [ ] Learn about parallelization (in C, although it should translate well to other languages)
- [p] Learn to make and use a makefile
- [ ] Make this a useful command line tool
- [ ] Learn to implement command line flags (before I have just kept track of command line option orders which is not good practice)
- [x] Use Git more directly to simulate a team environment in software development
- [ ] Develop an environment where addition of functions and methods to this project is a streamlined process
- [x] Relearn writing recursive algorithms... apparently! (See readProtein.c method countCovalentBonds for details)
- [ ] Learn to use package versioning within Git
- [p] Code Profiling

As of 6/1/2022:

 The main executables are in src/run and there are multiple. This is because some of my [academic work](https://www.researchgate.net/profile/Brian-Andrews-11) has inspired investigation into the direct effect of short ranged interactions on the structure of amino acids and polypeptide chains. The executables include compilation and running commands at the top as comments and this will be changed. But, for now, the compiled executables will read in a pdb file with CONECT records of a protein, identify the backbone dihedrals, and rotate about the two dihedral angles of the central amino acid residue, and print out a "trajectory" of these rotations in an xyz file which currently needs to be named manually by the user in Steric_Rama_MPI.c before compilation. Example videos of triglycine, trialanine, butane, and central Ala in trialanine's side chain (some still are only in xyz, which can be opened with VMD) are in the movies directory. These movies are generated by VMD using the output xyz trajectories of the steric executable.
