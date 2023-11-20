Command Line: ./drop -f stericClashes -i test_files/clashes_test_conect.pdb -lc false

Reading structure file: test_files/clashes_test_conect.pdb

Covalent Bond Matrix:
1 1 1 1 2 2 2 3 3 4 4 5 5 6 6 7 7 7 6 7 7 7 8 8 8 5 6 6 7 7 8 8 8 9 9 
2 2 2 3 3 3 4 4 5 5 6 6 7 7 8 8 8 7 8 8 8 9 9 9 6 7 7 8 8 9 9 9 10 10 
2 2 3 3 3 4 4 5 5 6 6 7 7 8 8 8 7 8 8 8 9 9 9 6 7 7 8 8 9 9 9 10 10 
2 3 3 3 4 4 5 5 6 6 7 7 8 8 8 7 8 8 8 9 9 9 6 7 7 8 8 9 9 9 10 10 
1 1 1 2 2 3 3 4 4 5 5 6 6 6 5 6 6 6 7 7 7 4 5 5 6 6 7 7 7 8 8 
2 2 3 3 4 4 5 5 6 6 7 7 7 6 7 7 7 8 8 8 5 6 6 7 7 8 8 8 9 9 
2 3 3 4 4 5 5 6 6 7 7 7 6 7 7 7 8 8 8 5 6 6 7 7 8 8 8 9 9 
1 1 2 2 3 3 4 4 5 5 5 4 5 5 5 6 6 6 3 4 4 5 5 6 6 6 7 7 
2 3 3 4 4 5 5 6 6 6 5 6 6 6 7 7 7 4 5 5 6 6 7 7 7 8 8 
1 1 2 2 3 3 4 4 4 3 4 4 4 5 5 5 2 3 3 4 4 5 5 5 6 6 
2 3 3 4 4 5 5 5 4 5 5 5 6 6 6 3 4 4 5 5 6 6 6 7 7 
1 1 2 2 3 3 3 2 3 3 3 4 4 4 1 2 2 3 3 4 4 4 5 5 
2 3 3 4 4 4 3 4 4 4 5 5 5 2 3 3 4 4 5 5 5 6 6 
1 1 2 2 2 1 2 2 2 3 3 3 2 3 3 4 4 5 5 5 6 6 
2 3 3 3 2 3 3 3 4 4 4 3 4 4 5 5 6 6 6 7 7 
1 1 1 2 3 3 3 4 4 4 3 4 4 5 5 6 6 6 7 7 
2 2 3 4 4 4 5 5 5 4 5 5 6 6 7 7 7 8 8 
2 3 4 4 4 5 5 5 4 5 5 6 6 7 7 7 8 8 
3 4 4 4 5 5 5 4 5 5 6 6 7 7 7 8 8 
1 1 1 2 2 2 3 4 4 5 5 6 6 6 7 7 
2 2 3 3 3 4 5 5 6 6 7 7 7 8 8 
2 3 3 3 4 5 5 6 6 7 7 7 8 8 
1 1 1 4 5 5 6 6 7 7 7 8 8 
2 2 5 6 6 7 7 8 8 8 9 9 
2 5 6 6 7 7 8 8 8 9 9 
5 6 6 7 7 8 8 8 9 9 
1 1 2 2 3 3 3 4 4 
2 3 3 4 4 4 5 5 
1 1 2 2 2 3 3 
2 3 3 3 4 4 
1 1 1 2 2 
2 2 3 3 
2 3 3 
1 1 
2 



Number of dihedrals identified in structure: 6
Calculating initial dihedral angles.
Columns: Angle, Angle Type (phi, psi, etc), Residue Name, Residue Number, Dihedral Number
-59.601438 phi ILE 2 0
-60.475241 phi GLY 3 1
55.183015 psi GLY 1 2
-149.898428 psi ILE 2 3
-109.185350 .590  19.340  19.520  1.00  0.00ILE ILE 2 4
-77.266600 380  1.00  0.00           H
ATOMILE ILE 2 5

Done reading structure file: test_files/clashes_test_conect.pdb

Clash: (AtomType AtomNumber)-(AtomType AtomNumber) [Allowed Distance, Distance in Structure]
Clash: (O9)-(O28) [2.700000, 2.384387]
Clash: (H11)-(H15) [1.900000, 1.763859]
Clash: (C16)-(N29) [2.800000, 2.247243]
Clash: (C16)-(H30) [2.200000, 1.635268]
Clash: (H19)-(N29) [2.200000, 1.724210]
Clash: (H19)-(H30) [1.900000, 0.835284]
Clash: (O28)-(O36) [2.700000, 2.376426]
Number of steric clashes: 7

