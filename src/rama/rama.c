#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "../readProtein/readProtein.h"

void writeRamaDistribution(float phi, float psi, float value)
{
  FILE *fp;
  fp = fopen("alanine_SC_RAMA.txt", "a+");
  fprintf(fp, "%f %f %f\n", phi, psi, value);
  fflush(fp);
  fclose(fp);
}
