# DROP
Dihedral Rotation Of Proteins (DROP)

Copyright 2023- Brian Andrews.

<<<<<<< HEAD
### What this tool is:

The repository contains a work-in-progress command line tool which can manipulate protein structures by rotating dihedral angles. This tool will be used to investigate accessible regions of the high-dimensional conformation space of proteins based on atomic overlaps. More factors may be added to the analysis in the future.
=======
The repository contains a work-in-progress command line tool which can manipulate protein structures by rotating dihedral angles. It has tools which will be useful for preparing protein structures for visualizations in presentations, lectures, or papers and preparing protein structures for simulation. This tool will be used to investigate accessible regions of the high-dimensional conformation space of proteins based on atomic overlaps. More factors may be added to the analysis in the future.
>>>>>>> origin/main

### What this tool is not:

This is not a pdb (or any structure file) preprocessing tool. It requires a sanitized pdb file with CONECT records for covalent bonds. I personally use [GROMACS](https://www.gromacs.org/) ([gitlab](https://gitlab.com/gromacs/gromacs)) to do this. See ```scripts/generate-conect/``` for more details.

### Installation

You can compile DROP with ```make``` in the parent directory. Add drop to your path to use it without specifying a path to the binary after execution.

### Usage

The tool is organized as follows:
- There are "parent" options (-c, -f, --version or --help) consisting of the main executable and the first option. Ex: ```./drop -c```, ```./drop -f ```, or ```./drop -?```. Everything after the "parent" option are "child" options for commands.

- See the "parent" options: ```./drop --help```

- See different commands for ```-f```: ```./drop -c```

- See "child" options for the different commands: ```./drop -f [command string from -c] --help```

<<<<<<< HEAD
=======
After compilation, you can export the path to the DROP executable to use the tool without the ```./```.

### Examples

To measure the dihedral angles of a given structure, use ```measureDihedrals```. To manipulate protein structures, you can use ```setDihedral``` or ```setDihedralList```. The former allows for the modification of one dihedral angle and the latter allows the user to provide a list of dihedral angles in a file to change in the structure. First example for ```setDihedral```:

```./drop -f setDihedral -i example_files/ILE_conect_110.pdb -n 2 -d phi -a -60 -o output.pdb```

And you can visualize ```output.pdb``` with your favorite visualization tool like VMD or PyMol. ```setDihedralList``` has examples in ```examples/README```, but here's an example command:

```./drop -f setDihedralList -i examples/setDihedralList/Polyarginine/poly-R-beta.pdb -d examples/setDihedralList/Polyarginine/beta-to-helix.txt -o poly-R-helix.pdb```

>>>>>>> origin/main
### Known Issues

- The "parent" options have to be first when specifying options due to how ```argp``` was implemented.
- Using this tool with large protein leads to a segfault. See Issue [#6](https://github.com/andrewsb8/DROP/issues/6).
