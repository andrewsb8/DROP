# DROP
Dihedral Rotation Of Proteins (DROP)

Copyright 2023- Brian Andrews.

The repository contains a work-in-progress command line tool which can manipulate protein structures by rotating dihedral angles. This tool will be used to investigate accessible regions of the high-dimensional conformation space of proteins based on atomic overlaps. More factors may be added to the analysis in the future.

### Installation

You can compile DROP with ```make``` in the parent directory. (Not working yet) Although it is recommended to make a ```build``` directory and execute ```make -C /path/to/DROP```.

### Usage

The tool is organized as follows:
- There are "parent" options (-c, -f, --help, etc) consisting of the main executable and the first option. Ex: ```./drop -c``` or ```./drop -f ```. Everything after the "parent" option are "child" options for commands.

- See the "parent" options: ```./drop --help```

- See different commands for ```-f```: ```./drop -c```

- See "child" options for the different commands: ```./drop -f [command string from -c] --help```

After compilation, you can export the path to DROP to use the tool without the ```./```.

The only current working command is ```./drop -f setDihedral -i [input pdb file]...```. This command allows the user to change a single dihedral angle based on user input.

NOTE: files in 'results' cannot be currently replicated by compiling from source. The c files used for completing those analyses are in 'src/run' and are currently being refactored into the new code. The code and results are being preserved until testing shows the results are reproducible after refactoring is completed.
