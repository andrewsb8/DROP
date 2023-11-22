struct _atoms
{
  int atom_number;
  char atom_type[5]; //5 instead of 4 for the termination character
  char residue[5]; //5 instead of 4 for the termination character
  int residue_number;
  double coordinates[3];
  char atom_name[2]; //2 instead of 1 for the termination character

  /*
  This section will describe the covalent_bondArray below. For atom i, it is an array of ints from
  j = i+1 to prot->number_of_atoms - 1 (N). Each int is the number of covalent bonds between atom i and
  atom j.

  This is a visualization of the arrays, where there is an array for each atom:

  i(atom number) array
  1             [covalent bonds between atoms 1 and 2, covalent bonds between atoms 1 and 3, ... , covalent bonds between atom 1 and last atom] (length prot->number_of_atoms - 2, atom 1 not included)
  2             [covalent bonds between atoms 2 and 3, covalent bonds between atoms 2 and 4, ... , covalent bonds between atom 2 and last atom]
  .
  .
  .
  N             [covalent bonds between 2nd to last atom and last atom] (this has length 1 by definition)

  For a single molecule, there should be no zeroes in these arrays!
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
  char *dihedral_angType[4]; //stores type of dihedral. e.g. phi, psi, chi1, chi2, omega (omega not implemented yet). Lower Case!!!!!
  char *dihedral_resName[4]; // stores residue name
  int dihedral_resNum; //stores residue number
};

struct _residues
{
  int backbone_atoms[8]; //backbone atom numbers for a residue. 8 is the most backbone atoms possible (N- or C-terminal (COOH) glycine with 8 atoms)
  int sidechain_atoms[18]; //side chain atom numbers for a residue. 18 is the most side chain atoms possible (Arg or Trp)
};

struct protein
{
  struct _atoms *atoms;
  struct _bonds *bonds;
  struct _dihedrals *dihedrals;
  struct _residues *residues;
  int number_of_residues;
  int number_of_atoms;
  int number_of_bonds;
  int number_of_dihedrals; //count of dihedrals
  int expected_num_dihedrals; //in an unblocked polypeptide, this should be equal to number_of_dihedrals
};

//list of backbone atoms including NH3 and COOH termini atoms
static char *backbone_atom_list[15] = { "N" ,"H1", "H2", "H3", "HN", "HA", "HA1", "HA2", "CA", "C", "O", "OT", "OT1", "OT2", "HT2" };
static int size_bb_atom_list = sizeof(backbone_atom_list) / sizeof(backbone_atom_list)[0];

//dihedral definitions for use in identifyDihedrals
static int numberDihedralTypes = 2;
static char *backboneDihedralDefinitions[2][4] = { //can't use int to set this array size?
  //backbone
  {"C", "N", "CA", "C"}, //phi
  {"N", "CA", "C", "N"},  //psi
};

static char *chi1DihedralDefinitions[3][4] = {
  {"N", "CA", "CB", "HB3"}, //Ala "chi 1" (use HB3 b/c no other amino acid has this atom type)
  {"N", "CA", "CB", "CG1"}, //Ile, Val chi 1
  {"N", "CA", "CB", "CG"} //Leu chi 1
};

static char *chi2DihedralDefinitions[2][4] = {
  {"CA", "CB", "CG1", "CD"}, //Ile chi 2
  {"CA", "CB", "CG", "CD1"} //Leu chi 2
};

static char *chi3DihedralDefinitions[1][4];

static char *chi4DihedralDefinitions[1][4];

static char *chi5DihedralDefinitions[1][4];

static char *customDihedralDefintions[1][4] = {
  {"C", "C", "C", "C"} //dihedral for butane
};


void readPDB(struct protein *prot,char *filename, FILE *log_file);
bool isBackbone(char *atomtype);
char * substr(char * s, int x, int y);
char * removeSpaces(char *string);
void readPDBbonds(struct protein *prot, char *filename, FILE *log_file);
void makeBondMatrix(struct protein *prot);
void countCovalentBonds(struct protein *prot, FILE *log_file);
int recursivePairSearch(struct protein *prot, int previousAtom, int atom1, int atom2, int found, int *covalentBondCount);
void printCovalentBondMatrix(struct protein *prot, FILE *log_file);
void identifyDihedrals(struct protein *prot);
void printXYZ(struct protein *prot);
void writeXYZ(struct protein *prot, char *filename, char *comment, char type, int frame, int rank);
void writeXYZmultiframe(struct protein *prot,char *filename, char *comment, int frame, int rank);
void writeXYZsingleframe(struct protein *prot,FILE *fp, char *comment, int rank);
void writePDB(struct protein *prot,char *filename,char type,int frame, bool conect);
void writePDBmultiframe(struct protein *prot, FILE *fp, int frame);
void writePDBsingleframe(struct protein *prot, FILE *fp);
void writePDBConect(struct protein *prot, FILE *fp);
