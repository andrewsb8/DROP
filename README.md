# DROP
Dihedral Rotation Of Proteins (DROP)

Copyright 2023- Brian Andrews.

### What this tool is:

The repository contains a work-in-progress command line tool which can manipulate protein structures by rotating dihedral angles. It has tools which will be useful for preparing protein structures for visualizations in presentations, lectures, or papers and preparing protein structures for simulation. This tool will be used to investigate accessible regions of the high-dimensional conformation space of proteins based on atomic overlaps. The scope of this project will continue to expand over time with additional modules for analyzing and manipulating molecular structures.

### What this tool is not:

This is not a pdb (or any structure file) preprocessing tool. It requires a sanitized pdb file and currently requires CONECT records for covalent bonds. I personally use [GROMACS](https://www.gromacs.org/) ([gitlab](https://gitlab.com/gromacs/gromacs)) to produce ordered pdb files and CONECT records. See ```scripts/generate-conect/``` for more details. The requirement for CONECT records will be removed for some modules in the future (see Issue #15 for more details).

### Installation

You can compile DROP with ```make``` in the parent directory.

### Usage

- See commands: ```./drop --help```

- See command options: ```./drop [command] --help```

After compilation, you can export the path to the DROP executable to use the tool without the ```./```.

### Examples

To measure the dihedral angles of a given structure, use ```measureDihedrals```. To manipulate protein structures, you can use ```setDihedral``` or ```setDihedralList```. The former allows for the modification of one dihedral angle and the latter allows the user to provide a list of dihedral angles in a file to change in the structure. First example for ```setDihedral```:

```./drop setDihedral -i example_files/ILE_conect_110.pdb -n 2 -d phi -a -60 -o output.pdb```

which changes the dihedral angle phi of residue 2 to -60 degrees. You can visualize ```output.pdb``` with your favorite visualization tool like VMD or PyMol. ```setDihedralList``` has examples in ```examples/setDihedralList/README.md```, but here's an example command:

```./drop setDihedralList -i examples/setDihedralList/Polyarginine/poly-R-beta.pdb -d examples/setDihedralList/Polyarginine/beta-to-helix.txt -o poly-R-helix.pdb```

### Papers Featuring DROP

If trying to reproduce results from a specific paper below, you may find past versions, scripts, and raw data in the Releases tab.

- B. Andrews. Amino Acid Residue-Specific Ramachandran Distributions Derived from a Simple Mean Field Potential. Physical Chemistry Au. 2024. 10.1021/acsphyschemau.4c00064.

### Citation

If you use this work for research or presentations, please consider citing this repository.

```
@software{Andrews_DROP_2024,
  author = {Andrews, Brian},
  title = {{Dihedral Rotation of Proteins (DROP)}},
  url = {https://github.com/andrewsb8/DROP/tree/main},
  year = {2024}
}
```

### Formatting

The ```indent``` package is used to format this code with the following command:

```indent -kr -ts4 [input file]```
