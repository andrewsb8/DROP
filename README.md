# DROP
Dihedral Rotation Of Proteins (DROP)

Copyright 2023- Brian Andrews.

### What this tool is:

The repository contains a work-in-progress command line tool which can manipulate protein structures by rotating dihedral angles. It has tools which will be useful for preparing protein structures for visualizations in presentations, lectures, or papers and preparing protein structures for simulation. This tool will be used to investigate accessible regions of the high-dimensional conformation space of proteins based on atomic overlaps. The scope of this project will continue to expand over time with additional modules for analyzing and manipulating molecular structures.

### What this tool is not:

This is not a pdb (or any structure file) preprocessing tool. It requires a sanitized pdb file and currently require CONECT records for covalent bonds. I personally use [GROMACS](https://www.gromacs.org/) ([gitlab](https://gitlab.com/gromacs/gromacs)) to produce ordered pdb files and CONECT records. See ```scripts/generate-conect/``` for more details. The requirement for CONECT records will be removed for some modules in the future (see Issue #15 for more details).

### Installation

You can compile DROP with ```make``` in the parent directory.

### Usage

The tool is organized as follows:
- There are "parent" options (-c, -f, --help, etc) consisting of the main executable and the first option. Ex: ```./drop -c``` or ```./drop -f ```. Everything after the "parent" option are "child" options for commands.

- See the "parent" options: ```./drop --help```

- See different commands for ```-f```: ```./drop -c```

- See "child" options for the different commands: ```./drop -f [command string from -c] --help```

After compilation, you can export the path to the DROP executable to use the tool without the ```./```.

### Examples

To measure the dihedral angles of a given structure, use ```measureDihedrals```. To manipulate protein structures, you can use ```setDihedral``` or ```setDihedralList```. The former allows for the modification of one dihedral angle and the latter allows the user to provide a list of dihedral angles in a file to change in the structure. First example for ```setDihedral```:

```./drop -f setDihedral -i example_files/ILE_conect_110.pdb -n 2 -d phi -a -60 -o output.pdb```

which changes the dihedral angle phi of residue 2 to -60 degrees. You can visualize ```output.pdb``` with your favorite visualization tool like VMD or PyMol. ```setDihedralList``` has examples in ```examples/README```, but here's an example command:

```./drop -f setDihedralList -i examples/setDihedralList/Polyarginine/poly-R-beta.pdb -d examples/setDihedralList/Polyarginine/beta-to-helix.txt -o poly-R-helix.pdb```

### Known Issues

- "Parent" option ```-f``` or ```-c``` must come first after the executable. See Issue #10 for more details.

NOTE: files in 'src/archive/results' cannot be currently replicated by compiling from source. The c files used for completing those analyses are in 'src/archive/run' and are currently being refactored into the new code. The code and results are being preserved until testing shows the results are reproducible after refactoring is completed.
