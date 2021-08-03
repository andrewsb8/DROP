struct _atoms
{
  int atom_number;
  char atom_type[4];
  char residue[4];
  int residue_number;
  double coordinates[3];
  char atom_name[1];

  /*
  This section will describe the covalent_bondArray below. For atom i, it is an array of ints from
  j = i+1 to prot->number_of_atoms - 1. Each int is the number of covalent bonds between atom i and
  atom j.

  Since each atom i has one and j must be great than i, this is a visualization of the arrays

  i(atom number) array
  1 [covalent bonds between 1 and 2, ... , covalent bonds between 1 and last atom]
  2 [covalent bonds between 2 and 3, ... , covalent bonds between 2 and last atom]
  .
  .
  .
  last atom -1 [covalent bonds between 2 and last atom]
  */

  int* covalent_bondArray; //stores number of covalent bonds between atoms i and j where *!*j > i*!* which cuts memory required for this in half!!!
  int len_covalent_bondArray; //stores the length of the above array so that it is convenient to retrieve
};

struct _bonds
{
  int bond_atomNumbers[2];
};

struct _dihedrals
{
  int dihedral_atomNumbers[4];
  double dihedral_angle;
  char *phi_or_psi[3]; //stores which type of dihedral it is
  char *dihedral_identifier[4]; // stores residue type and number i.e. GLY1. NOTE: not implemented yet
};

struct protein
{
  struct _atoms *atoms;
  struct _bonds *bonds;
  struct _dihedrals *dihedrals;
  int number_of_residues;
  int number_of_atoms;
  int number_of_bonds;
  int number_of_dihedrals; //count of dihedrals **only backbone dihedrals for now**
  int expected_num_dihedrals; //in an unblocked polypeptide, this should be equal to number_of_dihedrals
};

void readPDB(struct protein *prot,char *filename);
void readPDBbonds(struct protein *prot, char *filename);
void makeBondMatrix(struct protein *prot);
void countCovalentBonds(struct protein *prot);
int recursivePairSearch(struct protein *prot, int atom1, int atom2, int found, int *covalentBondCount);
void identifyDihedrals(struct protein *prot);
void printXYZ(struct protein *prot);
void writeXYZ(struct protein *prot, char *filename, char *comment, char type, int frame, int rank);
void writeXYZmultiframe(struct protein *prot,char *filename, char *comment, int frame, int rank);
void writeXYZsingleframe(struct protein *prot,char *filename, char *comment, int rank);
void writePDB(struct protein *prot,char *filename,char type,int frame);
void writePDBmultiframe(struct protein *prot,char *filename,int frame);
void writePDBsingleframe(struct protein *prot,char *filename);
