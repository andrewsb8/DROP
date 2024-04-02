#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fatal.h"

void drop_fatal(FILE *log, char *message)
{
  fprintf(log, "%s", message);
  fprintf(stderr, "%s", message);
  exit(1);
}
