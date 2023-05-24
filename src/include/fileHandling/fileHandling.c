#include <stdio.h>
#include <unistd.h>
#include "fileHandling.h"

//does the file exist?
int fileExists(char *filename)
{
  return access(filename, F_OK);
}
