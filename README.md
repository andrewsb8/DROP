# DROP
Dihedral Rotation Of Proteins (DROP)

Copyright 2023- Brian Andrews.

The repository contains a work-in-progress command line tool which can manipulate protein structures by rotating dihedral angles in the molecule. This tool will be used to investigate accessible regions of the high-dimensional conformation space of proteins based on atomic overlaps. More factors may be added in the future.

You can compile DROP with ```make``` in the parent directory. Although it is recommended to make a ```build``` directory and execute ```make -C /path/to/DROP```.

See options with ```./drop --help```

See different commands with ```./drop -c```

See options for the different commands with ```./drop -f [command string] --help```

The only current working command is ```./drop -f setDihedral -i [input pdb file]...```. This command allows the user to change a single dihedral angle based on user input.

NOTE: files in 'results' cannot be currently replicated by compiling from source. The c files used for completing those analyses are in 'src/run' and are currently being refactored into the new code. The code and results are being preserved until testing shows the results are reproducible after refactoring is completed.
