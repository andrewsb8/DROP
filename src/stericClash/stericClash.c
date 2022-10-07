#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../readProtein/readProtein.h"
#include "stericClash.h"
#include "../vectorCalculus/vectorCalculus.h"

//normally allowed via Ramachandran JMB 1963 in Angstroms
//Values:           C-C, C-O, C-N, C-H, O-O, O-N, O-H, N-N, N-H, H-H
//struct VDW radii = {3.2, 2.8, 2.9, 2.4, 2.8, 2.7, 2.4, 2.7, 2.4, 2.0};
//inner limit allowed via Ramachandran JMB 1963 in Angstroms
struct VDW radii = {3.0, 2.7, 2.8, 2.2, 2.7, 2.6, 2.2, 2.6, 2.2, 1.9};

//outer limit via Ramachandran JMB 1963
//struct VDW radii = {};

/*
Inside this function, the output of strcmp should be equal to zero. But, for whatever reason, it needs to be 10 here.
I can't explain this as strcmp is used normally elsewhere. But this problem is minor and not worth figuring out right now.
This note is just acknowledge my inevitable future confusion and say that "I know, it should be 0, but it has to be 10".
*/
int countClashes(struct protein *prot)
{
  double distance;
  double min_distance_allowed;
  int count = 0;
  //loop through all atom combinations i and j where j > i
  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    for(int j = i+1; j < prot->number_of_atoms; j++)
    {
      //check that at least x number of covalent bonds are between the atoms being compared
      //printf("%d %s %d %s %d\n", prot->atoms[i].atom_number, prot->atoms[i].atom_type, prot->atoms[j].atom_number, prot->atoms[j].atom_type, prot->atoms[i].covalent_bondArray[j-i-1]);
      if(prot->atoms[i].covalent_bondArray[j-i-1] > 3)
      {
        //if this condition is satisfied, check distance between atoms, and the difference of van der waals distances
        double *bond_vector = vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates);
        distance = vectorMagnitude(bond_vector);
        free(bond_vector);
        min_distance_allowed = getVDWRadii(&radii, prot->atoms[i].atom_name, prot->atoms[j].atom_name);
        //printf("%d %s %d %s %f %f\n", prot->atoms[i].atom_number, prot->atoms[i].atom_type, prot->atoms[j].atom_number, prot->atoms[j].atom_type, distance, min_distance_allowed);
        if(distance < min_distance_allowed)
        {
          //printf("%d %s %d %s %f %f\n", prot->atoms[i].atom_number, prot->atoms[i].atom_type, prot->atoms[j].atom_number, prot->atoms[j].atom_type, distance, min_distance_allowed);
          //return 1; //clash found
          count += 1;
        }
      }
    }
  }

  return count;
}

double getVDWRadii(struct VDW *radii, char *atom_name, char *atom_name_2)
{
  switch(*atom_name)
  {
    case 'N' :
      switch(*atom_name_2)
      {
        case 'N':
          return radii->NN;
        case 'H':
          return radii->NH;
        case 'C':
          return radii->CN;
        case 'O':
          return radii->ON;
      }
    case 'C' :
      switch(*atom_name_2)
      {
        case 'N':
          return radii->CN;
        case 'H':
          return radii->CH;
        case 'C':
          return radii->CC;
        case 'O':
          return radii->CO;
      }
    case 'O' :
      switch(*atom_name_2)
      {
        case 'N':
          return radii->ON;
        case 'H':
          return radii->OH;
        case 'C':
          return radii->CO;
        case 'O':
          return radii->OO;
      }
    case 'H' :
      switch(*atom_name_2)
      {
        case 'N':
          return radii->NH;
        case 'H':
          return radii->HH;
        case 'C':
          return radii->CH;
        case 'O':
          return radii->OH;
      }
    default :
      printf("No atom name defined for %s.\n", *atom_name);
      break;
  }
}

/*
  Below this comment is the old function from before I utilized the covalent bond matrix (see readProtein.c)
*/
/*int checkClashes(struct protein *prot)
{
  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    for(int j = i+6; j < prot->number_of_atoms; j++)
    {
      //below is a check print statement that I will leave in to debug in future
      //printf("%d %d %s %d %s %d\n",i,j, prot->atoms[i].atom_name,strcmp(prot->atoms[i].atom_name,"N"),prot->atoms[j].atom_name, strcmp(prot->atoms[j].atom_name,"C"));
      if( strcmp(prot->atoms[i].atom_name,"N") == 10  && strcmp(prot->atoms[j].atom_name,"N") == 10 && (prot->atoms[i].residue_number != prot->atoms[j].residue_number) )
      {
        //printf("No conflict! N-N: %d %d %f %f\n",i,j, vectorMagnitude(vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates)), (radii.N + radii.N));
        if( vectorMagnitude(vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates)) < (radii.N + radii.N) )
        {
          printf("N-N: %d %d %f %f\n",i,j, vectorMagnitude(vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates)), (radii.N + radii.N));
        }
      }
      else if( strcmp(prot->atoms[i].atom_name,"C") == 10 && strcmp(prot->atoms[j].atom_name,"C") == 10 && (prot->atoms[i].residue_number != prot->atoms[j].residue_number) )
      {
        //printf("No conflict! C-C: %d %d %f %f\n",i,j, vectorMagnitude(vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates)), (radii.C + radii.C));
        if( vectorMagnitude(vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates)) < (radii.C + radii.C) )
        {
          printf("C-C: %d %d %f %f\n",i,j, vectorMagnitude(vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates)), (radii.C + radii.C));
        }
      }
      else if( strcmp(prot->atoms[i].atom_name,"N") == 10 && strcmp(prot->atoms[j].atom_name,"C") == 10 && (prot->atoms[i].residue_number != prot->atoms[j].residue_number) )
      {
        //printf("No conflict! N-C: %d %d %f %f\n",i,j, vectorMagnitude(vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates)), (radii.C + radii.N));
        if( vectorMagnitude(vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates)) < (radii.N + radii.N) )
        {
          printf("N-C: %d %d %f %f\n",i,j, vectorMagnitude(vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates)), (radii.C + radii.N));
        }
      }
      else if( strcmp(prot->atoms[i].atom_name,"C") == 10 && strcmp(prot->atoms[j].atom_name,"N") == 10 && (prot->atoms[i].residue_number != prot->atoms[j].residue_number) )
      {
        //printf("No conflict! C-N: %d %d %f %f\n",i,j, vectorMagnitude(vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates)), (radii.C + radii.N));
        if( vectorMagnitude(vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates)) < (radii.N + radii.N) )
        {
          printf("C-N: %d %d %f %f\n",i,j, vectorMagnitude(vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates)), (radii.C + radii.N));
        }
      }
    }
  }
}*/
