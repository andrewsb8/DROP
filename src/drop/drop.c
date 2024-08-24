#include <stdio.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "commands.h"

const char *program_bug_address =
  "https://github.com/andrewsb8/DROP/issues";
const char *program_version = "DROP Version 0.0.1";
static char doc[] =
  "DROP (Dihedral Rotation Of Proteins) -- A command line tool to investigate protein structures via direct manipulation of dihedral angles.";

int
main (int argc, char *argv[])
{
  fprintf(stderr, "%s\n", argv[1]);
  if (strcmp(argv[1], "-?") == 0 || strcmp(argv[1], "--help") == 0)
      {
          fprintf(stderr, "this is the help message\n");
      }
  else if (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0)
      {
          fprintf(stderr, "%s\n", program_version);
      }
  else
  {
      if (!findCommand (argv[1], argc, argv))
          {
              fprintf(stderr, "Command %s not recognized.\n", argv[1]);
              //print usage
          }
  }
  return 0;
}
