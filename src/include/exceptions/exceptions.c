#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "exceptions.h"

void drop_fatal(FILE *log, char *message)
{
  fprintf(log, "%s", message);
  fprintf(stderr, "%s", message);
  exit(1);
}

void drop_warning(FILE *log, char *message)
{
  fprintf(log, "%s", message);
  fprintf(stderr, "%s", message);
}
