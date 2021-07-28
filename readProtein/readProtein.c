#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "readProtein.h"
#include "../dihedralRotation/dihedralRotation.h"

void readPDB(struct protein *prot, char *filename)
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
    //keep count of number of "words" or "entries" in a line of a pdb
    count = 0;

    //Get the first word of each line by splitting the string
    line_split = strtok(line, " \t");
    lens = strlen(line_split);

    //create string array for whole line of pdb file
    char stringT[numberOfEntries][lens];

    //For this function, I only want lines that start with ATOM.
    //Quantities to get: Atom Number, Atom name, Atom type, residue number, residue (all in a struct)
    //Atom poitions (own array/table within the struct)
    if(strcmp(line_split, "ATOM") == 0)
    {
      while(line_split != NULL)
      {
        //printf("%s\n", line_split);
        strcpy(stringT[count], line_split); //see if strings are the same
        count++;
        line_split = strtok(NULL, " \t");
      }

      //reallocate memory dynamically which allows for a flexible number of atoms from entry pdb
      if(line_number > 0)
      {
        prot->atoms = (struct _atoms*) realloc(prot->atoms, size*(line_number+1));
      }

      //since pdb files have a standard format, assignments are handled manually as opposed to using if/else if or switch cases to be concise
      prot->atoms[line_number].atom_number = atoi(stringT[1]);
      strcpy(prot->atoms[line_number].atom_type, stringT[2]);
      strcpy(prot->atoms[line_number].residue, stringT[3]);
      prot->atoms[line_number].residue_number = atoi(stringT[4]);
      for(int k = 5; k < 8; k++)
      {
        prot->atoms[line_number].coordinates[k-5] = strtod(stringT[k], &ptr);
      }
      strcpy(prot->atoms[line_number].atom_name, stringT[10]);
      line_number++;

    }

  }

  prot->number_of_atoms = line_number;
  prot->number_of_residues = prot->atoms[line_number-1].residue_number;

}



void readPDBbonds(struct protein *prot, char *filename)
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

    //For this function, I only want lines that start with ATOM.
    //Quantities to get: Atom Number, Atom name, Atom type, residue number, residue (all in a struct)
    //Then all of the atoms basically become an array or list of struct data types
    if(strcmp(line_split, "CONECT") == 0)
    {
      while(line_split != NULL)
      {
        //printf("%s\n", line_split);
        strcpy(stringT[count], line_split); //see if strings are the same
        count++;
        line_split = strtok(NULL, " \t");
      }

      //reallocate memory dynamically which allows for a flexible number of atoms from entry pdb
      if(line_number > 0)
      {
        prot->bonds = (struct _bonds*) realloc(prot->bonds, size*(line_number+1));
      }

      //since pdb files have a standard format, assignments are handled manually as opposed to using if/else if or switch cases to be concise
      for(int k = 1; k < 3; k++)
      {
        prot->bonds[line_number].bond_atomNumbers[k-1] = atoi(stringT[k]);
      }
      line_number++;

    }

  }

  prot->number_of_bonds = line_number;

  makeBondMatrix(prot);

}

void makeBondMatrix(struct protein *prot)
{
  for(int t = 0; t < prot->number_of_atoms; t++)
  {
    prot->bonds[t].len_covalent_bondArray = prot->number_of_atoms - (t+1);
    prot->bonds[t].covalent_bondArray = (int*) calloc(prot->number_of_atoms, sizeof(int*));

    //now need to search through bonds to count the number of covalent bonds between two atom numbers
    for(int u = t+1; u < prot->number_of_atoms; u++)
    {
      prot->bonds[t].covalent_bondArray[u] = countCovalentBonds(&prot, prot->atoms[t].atom_number, prot->atoms[u].atom_number)
    }
  }
}

int countCovalentBonds(struct protein *prot, int atom1, int atom2)
{
  //I think this function is going to have to be recursive....... and oh my god I didn't even consider proline....

  return 0;
}


void identifyDihedrals(struct protein *prot)
{
  //the number of dihedrals a protein has is 2*number of residues - 2 (N and C terminal only typically have 1 dihedral angle each in unblocked polypeptides)
  prot->expected_num_dihedrals = (prot->number_of_residues*2)-2;

  //dihedral definitions
  const int numberDihedralTypes = 2;
  char *dihedralDefinitions[3][4] = {
    {"C", "N", "CA", "C"}, //phi
    {"N", "CA", "C", "N"}  //psi
  };

  //allocate memory for the dihedrals struct to store information
  size_t size = sizeof(struct _dihedrals);
  prot->dihedrals = (struct _dihedrals*) malloc(size);

  //check variable to compare strings
  int check;
  //variables to save index of bond struct which contains pairs of
  //atoms within a dihedral
  int pairOne_index;
  int pairTwo_index;

  prot->number_of_dihedrals = 0;

  for(int m = 0; m < numberDihedralTypes; m++)
  {
    //printf("%s\n", dihedralDefinitions[m][3]);
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
                if(m==0)
                {
                  prot->dihedrals[prot->number_of_dihedrals].phi_or_psi[3] = "Phi";
                }
                else
                {
                  prot->dihedrals[prot->number_of_dihedrals].phi_or_psi[3] = "Psi";
                }

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
}

void writeXYZ(struct protein *prot,char *filename,char *comment,char type,int frame, int rank)
{
  //use type variable to decide if writing a multiframe or single frame pdb
  //there are slight syntax differences in the file formats
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

    //printf("%s", line);
    fprintf(fp, "%s\n", line);
    fflush(fp);
    free(line);
  }
  fclose(fp);
}

void writePDB(struct protein *prot,char *filename,char type,int frame)
{
  //use type variable to decide if writing a multiframe or single frame pdb
  //there are slight syntax differences in the file formats
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
}

void writePDBmultiframe(struct protein *prot,char *filename, int frame)
{
  printf("placeholder\n");
}

void writePDBsingleframe(struct protein *prot,char *filename)
{
  FILE *fp;
  fp = fopen(filename, "w");

  fprintf(fp, "MODEL\t1\n");
  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    char line[0];
    sprintf(line, "%4s   %5d  %s %3s %d %f %f %f %s", "ATOM", prot->atoms[i].atom_number, prot->atoms[i].atom_type, prot->atoms[i].residue, prot->atoms[i].residue_number, prot->atoms[i].coordinates[0], prot->atoms[i].coordinates[1], prot->atoms[i].coordinates[2], prot->atoms[i].atom_name);

    printf("%s", line);
    fprintf(fp, "%s", line);
  }
  fprintf(fp, "ENDMDL\n");
}
