#include <stdio.h>
#include <unistd.h>
#include "fileHandling.h"

int fileExists(char *filename)
{
  return access(filename, F_OK);
}
