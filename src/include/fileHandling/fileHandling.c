#include <stdio.h>
#include <unistd.h>
#include "fileHandling.h"

//does the file exist? returns -1 if not
int fileExists(char *filename)
{
  return access(filename, F_OK);
}
