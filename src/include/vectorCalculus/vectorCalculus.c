#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vectorCalculus.h"

//Calculates the magnitude of a vector
double vectorMagnitude(double vector1[3])
{
  double result = 0;
  for(int i = 0; i < 3; i++)
  {
    result += vector1[i]*vector1[i];
  }
  return sqrt(result);
}

//calculates the dot product component-wise of two vectors
double dotProduct(double vector1[3], double vector2[3])
{
  double result = 0;
  for(int i = 0; i < 3; i++)
  {
    result += vector1[i]*vector2[i];
  }
  return result;
}

//calculate the sum of two vectors
double * vectorAdd(double vector1[3], double vector2[3])
{
  double *result = malloc(sizeof(vector1)); //NOTE: this line throws a warning but does not break the code
  for(int i = 0; i < 3; i++)
  {
    result[i] = vector2[i] + vector1[i];
  }
  return result;
}

//calculate the difference of two vectors
double * vectorSubtract(double vector1[3], double vector2[3])
{
  double *result = malloc(sizeof(vector1)); //NOTE: this line throws a warning but does not break the code
  for(int i = 0; i < 3; i++)
  {
    result[i] = vector2[i] - vector1[i];
  }
  return result;
}

//calculates the cross product component wise of 2 vectors
//returns a double array!
double * crossProduct(double vector1[3], double vector2[3])
{
  double *result = malloc(sizeof(vector1)); //NOTE: this line throws a warning but does not break the code
  result[0] = vector1[1]*vector2[2] - vector2[1]*vector1[2];
  result[1] = vector1[2]*vector2[0] - vector2[2]*vector1[0];
  result[2] = vector1[0]*vector2[1] - vector2[0]*vector1[1];
  return result;
}

double * matrixVectorMult(double matrix[3][3], double vector[3])
{
  double *result = malloc(sizeof(vector));

  for(int n = 0; n < 3; n++)
  {
    result[n] = (matrix[n][0]*vector[0]) + (matrix[n][1]*vector[1]) + (matrix[n][2]*vector[2]);
  }

  return result;
}

double * vectorRotate(double vector1[3], int axisOfRotation, double angle)
{
  double mat[3][3][3] = {
    {{1,0,0},{0,cosl(angle),-sinl(angle)},{0,sinl(angle),cosl(angle)}}, //x
    {{cosl(angle),0,sinl(angle)},{0,1,0},{-sinl(angle),0,cosl(angle)}}, //y
    {{cosl(angle),-sinl(angle),0},{sinl(angle),cosl(angle),0},{0,0,1}} //z
  };

  return matrixVectorMult(mat[axisOfRotation], vector1);
}

double * cartesianToPolar(double vector[3])
{
  double *result = malloc(sizeof(vector));
  result[0] = vectorMagnitude(vector);
  result[1] = atan(vector[1]/vector[0]);
  result[2] = atan( sqrt( (vector[0]*vector[0]) + (vector[1]*vector[1]) )/vector[2] );

  return result;
}

double * polarToCartesian(double vector[3])
{
  double *result = malloc(sizeof(vector));
  result[0] = vector[0]*sin(vector[2])*cos(vector[1]);
  result[1] = vector[0]*sin(vector[2])*sin(vector[1]);
  result[2] = vector[0]*cos(vector[2]);

  return result;
}
