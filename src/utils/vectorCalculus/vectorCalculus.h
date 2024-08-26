#ifndef VECTORCALCULUS_H_
#define VECTORCALCULUS_H_

double vectorMagnitude(double vector[3]);
double dotProduct(double vector1[3], double vector2[3]);
double *vectorSubtract(double vector1[3], double vector2[3]);
double *vectorAdd(double vector1[3], double vector2[3]);
double *crossProduct(double vector1[3], double vector2[3]);
double *matrixVectorMult(double matrix[3][3], double vector[3]);
double *vectorRotate(double vector1[3], int axisOfRotation, double angle);
double *cartesianToPolar(double vector[3]);
double *polarToCartesian(double vector[3]);

#endif
