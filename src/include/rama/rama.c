#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void writeRamaDistribution(FILE *file, float phi, float psi, float value)
{
  fprintf(file, "%f %f %f\n", phi, psi, value);
}
