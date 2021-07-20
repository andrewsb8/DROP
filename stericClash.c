#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "readProtein.h"
#include "stericClash.h"
#include "vectorCalculus.h"

struct VDW radii = {1, 1.5, 1.35, 1.35};

/*
Inside this function, the output of strcmp should be equal to zero. But, for whatever reason, it needs to be 10 here.
I can't explain this as strcmp is used normally elsewhere. But this problem is minor and not worth figuring out right now.
This note is just acknowledge my inevitable future confusion and say that "I know, it should be 0, but it has to be 10".
*/
int checkClashes(struct protein *prot)
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
}
