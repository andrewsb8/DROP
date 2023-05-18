# DROP
Dihedral Rotation Of Proteins (DROP)

This is a work in progress personal project which manipulates protein structures via rotating dihedral angles. A major refactor to make this into a command line application has begun on 5/15/2023.

You can compile DROP with ```make``` in the parent directory. Although, currently the only thing that can be done is read a pdb file with CONECT records (see example files) to populate ```struct protein``` (see ```src/readProtein/readProtein.h```).

See options with ```./drop --help```

See different commands with ```./drop -c```

See options for the different commands with ```./drop -f [command string] --help```

The only current working command is ```./drop -f trial -i [input pdb file]```

NOTE: files in 'results' cannot be currently replicated by compiling from source. The c files used for completing those analyses are in 'src/run' and are currently being refactored into the new code. The code and results are being preserved until testing shows the results are reproducible after refactoring is completed.
