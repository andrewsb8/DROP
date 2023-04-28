#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../readProtein/readProtein.h"
#include "../dihedralRotation/dihedralRotation.h"
#include "../vectorCalculus/vectorCalculus.h"

long double PI = 3.14159265358979323846264338327950288419716939937510L;

double calculateDihedral(struct protein *prot, int dihedralNumber)
{
  double angle;
  double angle_in_rads;
  double otherangle;
  double dot;
  double norm;
  double sign;

  double *vec1 = vectorSubtract(prot->atoms[prot->dihedrals[dihedralNumber].dihedral_atomNumbers[0] - 1].coordinates, prot->atoms[prot->dihedrals[dihedralNumber].dihedral_atomNumbers[1] - 1].coordinates);
  double *vec2 = vectorSubtract(prot->atoms[prot->dihedrals[dihedralNumber].dihedral_atomNumbers[1] - 1].coordinates, prot->atoms[prot->dihedrals[dihedralNumber].dihedral_atomNumbers[2] - 1].coordinates);
  double *vec3 = vectorSubtract(prot->atoms[prot->dihedrals[dihedralNumber].dihedral_atomNumbers[2] - 1].coordinates, prot->atoms[prot->dihedrals[dihedralNumber].dihedral_atomNumbers[3] - 1].coordinates);

  double *cross1 = crossProduct(vec1,vec2);
  double *cross2 = crossProduct(vec2,vec3);

  double *doublecross = crossProduct(cross1, cross2);

  norm = vectorMagnitude(cross1)*vectorMagnitude(cross2);
  dot = dotProduct(cross1, cross2);
  angle_in_rads = acosl( dot/norm );
  angle = (180/PI)*angle_in_rads;

  double ip = dotProduct(vec1, doublecross);
  if (ip < 0)
  {
    sign = -1;
  }
  else
  {
    sign = 1;
  }
  otherangle = sign*(180/PI)*atan2( vectorMagnitude(doublecross), dot );

  //printf("angle and otherangle: %f %f\n", angle, otherangle);

  free(vec1);
  free(vec2);
  free(vec3);

  free(cross1);
  free(cross2);

  return angle;
}

//temp function to be removed, returns sign of angle between vectors
double determineSign(struct protein *prot, int dihedralNumber)
{
  double *vec1 = vectorSubtract(prot->atoms[prot->dihedrals[dihedralNumber].dihedral_atomNumbers[0] - 1].coordinates, prot->atoms[prot->dihedrals[dihedralNumber].dihedral_atomNumbers[1] - 1].coordinates);
  double *vec2 = vectorSubtract(prot->atoms[prot->dihedrals[dihedralNumber].dihedral_atomNumbers[1] - 1].coordinates, prot->atoms[prot->dihedrals[dihedralNumber].dihedral_atomNumbers[2] - 1].coordinates);
  double *vec3 = vectorSubtract(prot->atoms[prot->dihedrals[dihedralNumber].dihedral_atomNumbers[2] - 1].coordinates, prot->atoms[prot->dihedrals[dihedralNumber].dihedral_atomNumbers[3] - 1].coordinates);

  double *cross1 = crossProduct(vec1,vec2);
  double *cross2 = crossProduct(vec2,vec3);

  double *doublecross = crossProduct(cross1, cross2);

  double ip = dotProduct(vec1, doublecross);
  if (ip < 0)
  {
    return -1;
  }
  else
  {
    return 1;
  }

}

void updatePositions(struct protein *prot, double newPositions[3], int atomNumber)
{
  for(int j = 0; j < 3; j++)
  {
    prot->atoms[atomNumber].coordinates[j] = newPositions[j];
  }
}

double rotateDihedral(struct protein *prot, int dihedralNumber, double dihedralAngle, double dihedralAngleChange, int bb_or_sc, int chi)
{
  //translate all atoms such that the second atom of the dihedral of interest is at the origin
  int atom_to_origin = prot->dihedrals[dihedralNumber].dihedral_atomNumbers[1];
  int atom_rotation_index = prot->dihedrals[dihedralNumber].dihedral_atomNumbers[2];

  double translation[3]; //copy of translation vector
  for(int i = 0; i < 3; i++)
  {
    translation[i] = prot->atoms[atom_to_origin-1].coordinates[i];
  }

  //translate all atoms
  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    double *tmp = vectorSubtract(translation, prot->atoms[i].coordinates);
    //equate tmp values to prot->atoms[i].coordinates using a loop i guess
    updatePositions(prot, tmp, i);
    free(tmp);
  }

  //rotate about z axis to x-z plane
  double angleToXZPlane = atanl((prot->atoms[atom_rotation_index-1].coordinates[1])/(prot->atoms[atom_rotation_index-1].coordinates[0]));
  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    double *tmp = vectorRotate(prot->atoms[i].coordinates,2,-angleToXZPlane);
    updatePositions(prot, tmp, i);
    free(tmp);
  }

  //rotate about y axis so that the axis between the middle two atoms of the dihedral is aligned with the z axis
  double angleToZAxis = atanl((prot->atoms[atom_rotation_index-1].coordinates[0])/(prot->atoms[atom_rotation_index-1].coordinates[2]));
  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    double *tmp = vectorRotate(prot->atoms[i].coordinates,1,-angleToZAxis);
    updatePositions(prot, tmp, i);
    free(tmp);
  }

  //printXYZ(prot);

  //rotate all atoms after bond in question about the z axis by the desired change here
  if(bb_or_sc == 1) //1 indicates rotating backbone
  {
    for(int i = atom_rotation_index-1; i < prot->number_of_atoms; i++)
    {
      double *tmp = vectorRotate(prot->atoms[i].coordinates,2,(PI/180.0)*dihedralAngleChange);
      updatePositions(prot, tmp, i);
      free(tmp);
    }
  }
  else //rotate side chain instead
  {
    if(chi == 1) //which chi angle is being rotated, chi 1 or chi 2?
    {
      //int ala2_sidechain_temp[4] = {17,18,19,20}; //temporary to test sidechain rotation method. Need better, general implementation
      int ile2_chi1[13] = {7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 ,19};
      //int leu2_chi1[13] = {7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 ,19};
      //int val2_chi1[10] = {7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
      for(int i = 0; i < 13; i++)
      {
        double *tmp = vectorRotate(prot->atoms[ile2_chi1[i]-1].coordinates,2,(PI/180.0)*dihedralAngleChange);
        updatePositions(prot, tmp, ile2_chi1[i]-1);
        free(tmp);
      }
    }
    else if(chi == 2)
    {
      int ile2_chi2[6] = {14, 15, 16, 17, 18, 19};
      //int leu2_chi2[9] = {11, 12, 13, 14, 15, 16, 17, 18 ,19};
      for(int i = 0; i < 6; i++)
      {
        double *tmp = vectorRotate(prot->atoms[ile2_chi2[i]-1].coordinates,2,(PI/180.0)*dihedralAngleChange);
        updatePositions(prot, tmp, ile2_chi2[i]-1);
        free(tmp);
      }
    }
  }

  //undo the first three steps in reverse order
  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    double *tmp = vectorRotate(prot->atoms[i].coordinates,1,angleToZAxis);
    updatePositions(prot, tmp, i);
    free(tmp);
  }

  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    double *tmp = vectorRotate(prot->atoms[i].coordinates,2,angleToXZPlane);
    updatePositions(prot, tmp, i);
    free(tmp);
  }

  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    double *tmp = vectorAdd(translation, prot->atoms[i].coordinates);
    updatePositions(prot, tmp, i);
    free(tmp);
  }
}

//This method will be used to make successive rotations about the same dihedral
//For now, the rotation methods in this file translate protein back to position in
//box. But that is ok for now.
double rotateDihedral_noTranslate(struct protein *prot, int dihedralNumber, double dihedralAngle, double dihedralAngleChange)
{
  int atom_to_origin = prot->dihedrals[dihedralNumber].dihedral_atomNumbers[1] - 1;
  //rotate all atoms about the z axis by the desired change here
  for(int i = atom_to_origin+1; i < prot->number_of_atoms; i++)
  {
    double *tmp = vectorRotate(prot->atoms[i].coordinates,2,(PI/180.0)*dihedralAngleChange);
    updatePositions(prot, tmp, i);
    free(tmp);
  }
}
