#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "../readProtein/readProtein.h"

void writeRamaDistribution(char *filename, int frame, float phi, float psi, float value)
{
  FILE *fp;
  if(frame > 1)
  {
    fp = fopen(filename, "a+");
  }
  else
  {
    fp = fopen(filename, "w");
  }
  fprintf(fp, "%f %f %f\n", phi, psi, value);
  fflush(fp);
  fclose(fp);
}
