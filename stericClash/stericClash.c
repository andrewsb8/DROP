#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../readProtein/readProtein.h"
#include "stericClash.h"
#include "../vectorCalculus/vectorCalculus.h"

struct VDW radii = {1, 1.5, 1.35, 1.35};

/*
Inside this function, the output of strcmp should be equal to zero. But, for whatever reason, it needs to be 10 here.
I can't explain this as strcmp is used normally elsewhere. But this problem is minor and not worth figuring out right now.
This note is just acknowledge my inevitable future confusion and say that "I know, it should be 0, but it has to be 10".
*/
int checkClashes(struct protein *prot)
{
  double distance;
  double min_distance_allowed;
  //loop through all atom combinations i and j where j > i
  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    for(int j = i+1; j < prot->number_of_atoms; j++)
    {
      //check that at least x number of covalent bonds are between the atoms being compared
      if(prot->atoms[i].covalent_bondArray[j-i-1] > 2)
      {
        //if this condition is satisfied, check distance between atoms, and the difference of van der waals distances
        distance = vectorMagnitude(vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates));
        min_distance_allowed = getVDWRadii(&radii, prot->atoms[i].atom_name) + getVDWRadii(&radii, prot->atoms[j].atom_name);
        //printf("%s %s %f %f\n", prot->atoms[i].atom_name, prot->atoms[j].atom_name, distance, min_distance_allowed);
        if(distance < min_distance_allowed)
        {
          //printf("%d %s %d %s %f %f\n", prot->atoms[i].atom_number, prot->atoms[i].atom_name, prot->atoms[j].atom_number, prot->atoms[j].atom_name, distance, min_distance_allowed);
          return 1; //clash found
        }
      }
    }
  }

  return 0;
}

double getVDWRadii(struct VDW *radii, char *atom_name)
{
  switch(*atom_name)
  {
    case 'N' :
      return radii->N;
    case 'C' :
      return radii->C;
    case 'O' :
      return radii->O;
    case 'H' :
      return radii->H;
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
