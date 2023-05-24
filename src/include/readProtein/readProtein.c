#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "readProtein.h"
#include "../dihedralRotation/dihedralRotation.h"

void readPDB(struct protein *prot, char *filename, FILE *log_file)
{
  FILE *fp;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  char * line_split;
  int count;
  const int numberOfEntries = 11; //this is standard for PDB formats. 11 entries per line with delimiters between.
  int line_number = 0; //keep track of what line number I am on while reading the file
  char *ptr; //pointer to use for strtod (string to double) function later in this function

  int lens;
  int check = 1;

  fp = fopen(filename, "r");

  //allocate memory for the atoms struct to store information
  size_t size = sizeof(struct _atoms);
  prot->atoms = (struct _atoms*) malloc(size);

  //read file line by line until EOF
  while((read = getline(&line,&len,fp)) != -1)
  {
    //Quantities to get for ATOM entries: Atom Number, Atom name, Atom type, residue number, residue (all in a struct), Atom poitions (own array/table within the struct)
    if(strcmp(substr(line,0,4), "ATOM") == 0)
    {
      //reallocate memory dynamically which allows for a flexible number of atoms from entry pdb
      if(line_number > 0)
      {
        prot->atoms = (struct _atoms*) realloc(prot->atoms, size*(line_number+1));
      }

      //since pdb files have a standard format, assignments are handled manually as opposed to using if/else if or switch cases to be concise
      //structure of pdb file can be found here: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
      char *atomNum = substr(line,6,12);
      prot->atoms[line_number].atom_number = atoi(atomNum);
      free(atomNum);

      char *atomType = substr(line,12,16);
      strcpy(prot->atoms[line_number].atom_type, removeSpaces(atomType));
      free(atomType);

      char *residueName = substr(line,17,20);
      strcpy(prot->atoms[line_number].residue, removeSpaces(residueName));
      free(residueName);

      char *residueNumber = substr(line,22,26);
      prot->atoms[line_number].residue_number = atoi(residueNumber);
      free(residueNumber);

      char *xPos = substr(line,30,38);
      char *yPos = substr(line,38,46);
      char *zPos = substr(line,46,54);
      prot->atoms[line_number].coordinates[0] = strtod(xPos, &ptr);
      prot->atoms[line_number].coordinates[1] = strtod(yPos, &ptr);
      prot->atoms[line_number].coordinates[2] = strtod(zPos, &ptr);
      free(xPos);
      free(yPos);
      free(zPos);

      char *atomName = substr(line,77,78);
      strcpy(prot->atoms[line_number].atom_name, removeSpaces(atomName));
      free(atomName);

      line_number++;

    }

  }

  fclose(fp);

  if (line_number == 0)
  {
    fprintf(stderr, "ERROR: No ATOM entries in the input file. Exiting.\n");
    fprintf(log_file, "ERROR: No ATOM entries in the input file. Exiting\n");
    exit(1);
  }

  prot->number_of_atoms = line_number;
  prot->number_of_residues = prot->atoms[line_number-1].residue_number;

  readPDBbonds(prot, filename, log_file);
  identifyDihedrals(prot);

  //log initial dihedral angle values
  fprintf(log_file, "Number of dihedrals identified in structure: %d\n", prot->number_of_dihedrals);
  fprintf(log_file, "Calculating initial dihedral angles.\nColumns: Angle, Angle Type (phi, psi, etc), Residue Name, Residue Number\n");

  for(int i = 0; i < prot->number_of_dihedrals; i++)
  {
    prot->dihedrals[i].dihedral_angle = calculateDihedral(prot, i);
    fprintf(log_file, "%f %s %s %d\n", prot->dihedrals[i].dihedral_angle, prot->dihedrals[i].dihedral_angType, prot->dihedrals[i].dihedral_resName, prot->dihedrals[i].dihedral_resNum);
  }
  fprintf(log_file, "\n");

}

//function to return a portion of a string with user defined indices
//taken from here: https://stackoverflow.com/a/10375855
//need to brush up on how pointers work *shrug*
char * substr(char * s, int x, int y)
{
  size_t size_str = y-x;
  char * ret = malloc(size_str+1);
  char * p = ret;
  char * q = &s[x];
  assert(ret != NULL);
  while(x < y)
  {
    *p++ = *q++;
    x++;
  }
  *p++ = '\0';
  return ret;
}

//remove spaces from strings
//taken from here: https://www.educative.io/edpresso/how-to-remove-spaces-from-a-string-in-c
char * removeSpaces(char *string)
{
  int space_count = 0;
  for(int i = 0; string[i] != '\0'; i++)
  {
    if(string[i] != ' ')
    {
      string[space_count] = string[i];
      space_count++;
    }
  }
  string[space_count] = '\0';
  return string;
}

void readPDBbonds(struct protein *prot, char *filename, FILE *log_file)
{
  FILE *fp;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  char * line_split;
  int count;
  const int numberOfEntries = 3; //this is standard for PDB formats. 3 entries for the CONECT lines.
  int line_number = 0; //keep track of what line number I am on while reading the file

  int lens;
  int check = 1;

  fp = fopen(filename, "r");

  //allocate memory for the bonds struct to store information
  size_t size = sizeof(struct _bonds);
  prot->bonds = (struct _bonds*) malloc(size);

  //read file line by line until EOF
  while((read = getline(&line,&len,fp)) != -1)
  {
    //keep count of number of "words" or "entries" in a line of a pdb
    count = 0;

    //Get the first word of each line by splitting the string
    line_split = strtok(line, " \t");
    lens = strlen(line_split);

    //create string array for whole line of pdb file
    char stringT[numberOfEntries][lens];

    //For this function, I only want lines that start with CONECT.
    if(strcmp(line_split, "CONECT") == 0)
    {
      while(line_split != NULL)
      {
        strcpy(stringT[count], line_split); //see if strings are the same
        count++;
        line_split = strtok(NULL, " \t");
      }

      //reallocate memory dynamically which allows for a flexible number of atoms from entry pdb
      if(line_number > 0)
      {
        prot->bonds = (struct _bonds*) realloc(prot->bonds, size*(line_number+1));
      }

      //the second and third entries of the line are the atom numbers in the bond
      for(int k = 1; k < 3; k++)
      {
        prot->bonds[line_number].bond_atomNumbers[k-1] = atoi(stringT[k]);
      }
      line_number++;

    }

  }

  fclose(fp);

  if (line_number == 0)
  {
    fprintf(stderr, "ERROR: No CONECT entries in the input file. No bonds were read or inferred from structure. Exiting.\n");
    fprintf(log_file, "ERROR: No CONECT entries in the input file. No bonds were read or inferred from structure. Exiting.\n");
    exit(1);
  }

  prot->number_of_bonds = line_number;

  makeBondMatrix(prot);
  countCovalentBonds(prot, log_file);

  printf("\n\n");
  return;
}

//initiate covalent bond arrays for each atom. place 1s where atom i is covalently bonded to atom j.
//the way this is written may not catch all covalent bonds on first pass because this only checks first column of CONECT record
//those left will be caught in the recursive search.
void makeBondMatrix(struct protein *prot)
{
  //last atom will not be considered to avoid weird memory artifacts. It should not have any unique bonds anyway
  for(int t = 0; t < prot->number_of_atoms-1; t++)
  {
    prot->atoms[t].len_covalent_bondArray = prot->number_of_atoms - (t+1);
    prot->atoms[t].covalent_bondArray = (int*) calloc(prot->atoms[t].len_covalent_bondArray, sizeof(int*));

    for(int u = 0; u < prot->number_of_bonds; u++)
    {
      //if first atom in bond is the atom currently looking at, put a 1 in the location of the second atom in bond
      if(prot->bonds[u].bond_atomNumbers[0] == t+1)
      {
        prot->atoms[t].covalent_bondArray[prot->bonds[u].bond_atomNumbers[1] - (t+2)] = 1; //see readProtein.h for confusing indexing of this data structure
      }
    }
  }

  return;
}

//count covalent bonds between atoms i and j via recursive search using the CONECT records
void countCovalentBonds(struct protein *prot, FILE *log_file)
{
  int warning = 0;
  static int covalentBondCount = 0;
  for(int h = 0; h < prot->number_of_atoms; h++)
  {
    for(int p = h+1; p < prot->number_of_atoms; p++)
    {
      //only want to call search function for pairs of atoms that aren't directly covalently bonded to one another
      if(prot->atoms[h].covalent_bondArray[p-h-1] != 1)
      {
        recursivePairSearch(prot, 0, h+1, p+1, 0, &covalentBondCount);
        if (covalentBondCount == 0)
        {
          warning = warning + 1;
          fprintf(stderr, "WARNING %d: There is a zero in your covalent bond matrix. "
          "This means you either:\n\n1. Have at least one noncovalently bonded atyom "
          "in your structure.\n2. Are missing or have incorrect CONECT records.\n\nPlease "
          "check your pdb file as the results will not be accurate with this structure. "
          "Dihedral angles may be unrecognized because bond information is missing. See "
          "your log file for more details and your covalent bond matrix.\n\n", warning);
          fprintf(log_file, "WARNING %d: There is a zero in your covalent bond matrix. "
          "This means you either:\n\n1. Have at least one noncovalently bonded atom "
          "in your structure.\n2. Are missing or have incorrect CONECT records.\n\nPlease "
          "check your pdb file as the results will not be accurate with this structure. "
          "Dihedral angles may be unrecognized because bond information is missing.\n\n", warning);
        }
        prot->atoms[h].covalent_bondArray[p-h-1] = covalentBondCount;
        covalentBondCount = 0;
      }
    }
  }

  printCovalentBondMatrix(prot, log_file);
  return;
}

//recursive search through CONECT records to find path from atom i to atom j
//previousAtom is used to prevent infinite recursion from traveling backwards through the path of connected atoms
//CONECT records have two columns and some atom numbers may only exist in one of the two columns so this is accounted
//for by actively checking both columns for the desired atom number.
int recursivePairSearch(struct protein *prot, int previousAtom, int atom1, int atom2, int found, int *covalentBondCount)
{
  for(int z = 0; z < prot->number_of_bonds; z++)
  {
    if(prot->bonds[z].bond_atomNumbers[0] == atom1) //atom is in first column of bonds list
    {
      if(prot->bonds[z].bond_atomNumbers[1] == previousAtom) //don't go backwards
      {
        continue;
      }
      if(prot->bonds[z].bond_atomNumbers[1] == atom2) //condition for killing recursion
      {
        *covalentBondCount+=1;
        found = 1;
        return found;
      }
      found = recursivePairSearch(prot, atom1, prot->bonds[z].bond_atomNumbers[1], atom2, found, covalentBondCount);

      if(found == 1) //stops loop if search is completed
      {
        break;
      }
    }

    else if(prot->bonds[z].bond_atomNumbers[1] == atom1) // atom is in second column of bond list
    {
      if(prot->bonds[z].bond_atomNumbers[0] == previousAtom) //don't go backwards
      {
        continue;
      }
      if(prot->bonds[z].bond_atomNumbers[0] == atom2) //condition for killing recursion
      {
        *covalentBondCount+=1;
        found = 1;
        return found;
      }
      found = recursivePairSearch(prot, atom1, prot->bonds[z].bond_atomNumbers[0], atom2, found, covalentBondCount);

      if(found == 1)
      {
        break;
      }
    }
  }

  if(found == 1)
  {
    *covalentBondCount+=1;
  }

  return found;

}

//prints the covalent bond arrays, aka matrix. typically used only for log file after pdb processing.
void printCovalentBondMatrix(struct protein *prot, FILE *log_file)
{
  fprintf(log_file, "Covalent Bond Matrix:\n");
  for(int u = 0; u < prot->number_of_atoms; u++)
  {
    for(int v = 0; v < prot->atoms[u].len_covalent_bondArray; v++)
    {
      fprintf(log_file, "%d ", prot->atoms[u].covalent_bondArray[v]);
    }
    fprintf(log_file, "\n");
  }

  fprintf(log_file, "\n\n");
  return;
}

/*
This function has an issue where the CONECT section of a pdb (or whatever
is read into the bonds struct) needs to be in order.

Ex: The following will read as a dihedral (if the atom types are defined as a
dihedral in the program)

CONECT 1 2
CONECT 2 3
CONECT 3 4

But the following won't

CONECT 2 1
CONECT 2 3
CONECT 3 4

This could produce problems. Especially in cases where the atom numbering gets
weird or potentially out of order if this is to be expanded.

May copy recursive strategy used for the bond matrix above
5/18/2023: OR.. just use the covalent bond matrix already formed. Each atom
included in the dihedral has to be bonded and are identified by a '1' in
the bond matrix......
*/
void identifyDihedrals(struct protein *prot)
{
  //this number needs to be updated to account for amino acid type and side
  //chain dihedrals. then it should be used to produce potential warnings or
  //errors if one hasn't been thrown due to covalent bond matrix error.
  //move to it's own function?
  prot->expected_num_dihedrals = (prot->number_of_residues*2)-2;

  //allocate memory for the dihedrals struct to store information
  size_t size = sizeof(struct _dihedrals);
  prot->dihedrals = (struct _dihedrals*) malloc(size);

  //check variable to compare strings
  int check;
  //variables to save index of bond struct which contains pairs of atoms within a dihedral
  int pairOne_index;
  int pairTwo_index;

  prot->number_of_dihedrals = 0;

  for(int m = 0; m < numberDihedralTypes; m++)
  {
    for(int n = 0; n < prot->number_of_bonds; n++)
    {
      //if two atoms bonded are the same as the first two of the dihedral definition, search for the second pair
      if(strcmp( dihedralDefinitions[m][0], prot->atoms[prot->bonds[n].bond_atomNumbers[0]-1].atom_type ) == 0 && strcmp( dihedralDefinitions[m][1], prot->atoms[prot->bonds[n].bond_atomNumbers[1]-1].atom_type ) == 0)
      {
        pairOne_index = n;
        for(int p = 0; p < prot->number_of_bonds; p++)
        {
          //search for second pair
          if(strcmp( dihedralDefinitions[m][2], prot->atoms[prot->bonds[p].bond_atomNumbers[0]-1].atom_type ) == 0 && strcmp( dihedralDefinitions[m][3], prot->atoms[prot->bonds[p].bond_atomNumbers[1]-1].atom_type ) == 0)
          {
            pairTwo_index = p;
            //multiple other bonds satisfy the above condition. Need to check for bond between the two pairs of bonds
            for(int r = 0; r < prot->number_of_bonds; r++)
            {
              if(prot->bonds[pairOne_index].bond_atomNumbers[1] == prot->bonds[r].bond_atomNumbers[0] && prot->bonds[pairTwo_index].bond_atomNumbers[0] == prot->bonds[r].bond_atomNumbers[1])
              {
                prot->dihedrals[prot->number_of_dihedrals].dihedral_atomNumbers[0] = prot->bonds[n].bond_atomNumbers[0];
                prot->dihedrals[prot->number_of_dihedrals].dihedral_atomNumbers[1] = prot->bonds[n].bond_atomNumbers[1];
                prot->dihedrals[prot->number_of_dihedrals].dihedral_atomNumbers[2] = prot->bonds[p].bond_atomNumbers[0];
                prot->dihedrals[prot->number_of_dihedrals].dihedral_atomNumbers[3] = prot->bonds[p].bond_atomNumbers[1];
                //m==0 is manually assigned to be phi in dihedralDefinitions in readProtein.h
                if(m==0)
                {
                  strcpy(prot->dihedrals[prot->number_of_dihedrals].dihedral_angType, "Phi");
                }
                //m==1 is manually assigned to be psi in dihedralDefinitions in readProtein.h
                if(m==1)
                {
                  strcpy(prot->dihedrals[prot->number_of_dihedrals].dihedral_angType, "Psi");
                }
                //add more conditions for chi angles. Output is random bytes for some amino acids.
                //but may need more systematic ordering of dihedralDefinitions

                //residue is identified by third atom in dihedral because that will always be in the ith residue. Could have used second atom also.
                //this will not work if omega is to be incorporated. but currently no plans to do so.
                strcpy(prot->dihedrals[prot->number_of_dihedrals].dihedral_resName, prot->atoms[ prot->dihedrals[prot->number_of_dihedrals].dihedral_atomNumbers[2] -1 ].residue);
                prot->dihedrals[prot->number_of_dihedrals].dihedral_resNum = prot->atoms[ prot->dihedrals[prot->number_of_dihedrals].dihedral_atomNumbers[2] -1 ].residue_number;

                prot->number_of_dihedrals += 1;
                if(prot->number_of_dihedrals > 0)
                {
                  prot->dihedrals = (struct _dihedrals*) realloc(prot->dihedrals, size*(prot->number_of_dihedrals + 1));
                }
              }
            }

          }
        }
      }
    }
  }
  return;
}

void printXYZ(struct protein *prot)
{
  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    size_t bufsz = snprintf(NULL, 0, "%s %f %f %f", prot->atoms[i].atom_type, prot->atoms[i].coordinates[0], prot->atoms[i].coordinates[1], prot->atoms[i].coordinates[2]);
    char* line = malloc(bufsz+1);
    sprintf(line, "%s %f %f %f", prot->atoms[i].atom_type, prot->atoms[i].coordinates[0], prot->atoms[i].coordinates[1], prot->atoms[i].coordinates[2]);

    printf("%s\n", line);
    free(line);
  }
  printf("\n\n");
  return;
}

void writeXYZ(struct protein *prot,char *filename,char *comment,char type,int frame, int rank)
{
  //use type variable to decide if writing a multiframe or single frame xyz
  switch(type){
    case 'm' :
      //printf("Writing multi-frame xyz file.\n");
      writeXYZmultiframe(prot, filename, comment, frame, rank);
      break;
    case 's' :
      printf("Writing single-structure xyz file.\n");
      writeXYZsingleframe(prot, filename, comment, rank);
      break;
    default :
      printf("Please specify single-structure (s) or multi-frame (m) pdb option.\n");
      break;
  }
  return;
}

void writeXYZsingleframe(struct protein *prot,char *filename, char *comment, int rank)
{
  FILE *fp;
  fp = fopen(filename, "w");

  fprintf(fp, "%d\n", prot->number_of_atoms);
  fprintf(fp, "Comment: %s\n", comment);
  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    size_t bufsz = snprintf(NULL, 0, "%s %f %f %f", prot->atoms[i].atom_type, prot->atoms[i].coordinates[0], prot->atoms[i].coordinates[1], prot->atoms[i].coordinates[2]);
    char* line = malloc(bufsz+1);
    sprintf(line, "%s %f %f %f", prot->atoms[i].atom_type, prot->atoms[i].coordinates[0], prot->atoms[i].coordinates[1], prot->atoms[i].coordinates[2]);

    fprintf(fp, "%s\n", line);
    free(line);
  }
  fclose(fp);
  printf("Completed XYZ structure generation. Filename: %s.\n", filename);
  return;
}

void writeXYZmultiframe(struct protein *prot,char *filename, char *comment, int frame, int rank)
{
  FILE *fp;
  if(frame > 0)
  {
    fp = fopen(filename, "a+");
  }
  else
  {
    fp = fopen(filename, "w");
  }


  fprintf(fp, "%d\n", prot->number_of_atoms);
  fprintf(fp, "%s\n", comment);
  char line[100];
  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    size_t bufsz = snprintf(NULL, 0, "%s %f %f %f", prot->atoms[i].atom_type, prot->atoms[i].coordinates[0], prot->atoms[i].coordinates[1], prot->atoms[i].coordinates[2]);
    char* line = malloc(bufsz+1);
    sprintf(line, "%s %f %f %f", prot->atoms[i].atom_type, prot->atoms[i].coordinates[0], prot->atoms[i].coordinates[1], prot->atoms[i].coordinates[2]);

    fprintf(fp, "%s\n", line);
    fflush(fp);
    free(line);
  }
  fclose(fp);
  return;
}

void writePDB(struct protein *prot,char *filename,char type,int frame)
{
  //use type variable to decide if writing a multiframe or single frame pdb
  switch(type){
    case 'm' :
      printf("Writing multi-frame pdb file.\n");
      writePDBmultiframe(prot, filename, frame);
      break;
    case 's' :
      printf("Writing single-structure pdb file.\n");
      writePDBsingleframe(prot, filename);
      break;
    default :
      printf("Please specify single-structure (s) or multi-frame (m) pdb option.\n");
      break;
  }
  return;
}

void writePDBmultiframe(struct protein *prot,char *filename, int frame)
{
  printf("placeholder\n");
  return;
}

void writePDBsingleframe(struct protein *prot,char *filename)
{
  FILE *fp;
  fp = fopen(filename, "w");

  fprintf(fp, "MODEL\t1\n");
  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    char line[40];
    sprintf(line, "%4s   %5d  %s %3s %d %f %f %f %s", "ATOM", prot->atoms[i].atom_number, prot->atoms[i].atom_type, prot->atoms[i].residue, prot->atoms[i].residue_number, prot->atoms[i].coordinates[0], prot->atoms[i].coordinates[1], prot->atoms[i].coordinates[2], prot->atoms[i].atom_name);

    printf("%s", line);
    fprintf(fp, "%s", line);
  }
  fprintf(fp, "ENDMDL\n");
  return;
}
