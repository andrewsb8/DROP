#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>

#include "../logging/logging.h"
#include "../readProtein/readProtein.h"

//does the file exist? returns -1 if not
int
fileExists (char *filename)
{
  return access (filename, F_OK);
}

//makes a string of argv to print command line to log
char *
makeStringArgv (int argc, char **argv)
{
  char *strng = argv[0];
  for (int i = 1; i < argc; i++)
	{
	  strcat (strng, " ");
	  strcat (strng, argv[i]);
	}
  return strng;
}

void
processInput (struct protein *prot, char *input_file, FILE * log,
			  bool calc_bond_matrix, bool print_bond_matrix, int argc, char **argv)
{
  //make a string of argv arguments for the log file output
  char *stringArgv = makeStringArgv (argc, argv);
  fprintf (log, "Command Line: %s\n\n", stringArgv);
  if (fileExists (input_file) == -1)
	{
	  drop_fatal (log, "ERROR: Input file does not exist. Exiting.\n");
	}

  fprintf (log, "Reading structure file: %s\n\n", input_file);

  readPDB (prot, input_file, log, calc_bond_matrix, print_bond_matrix);

  fprintf (log, "Done reading structure file: %s\n\n", input_file);
  return;
}


void
writeFileLine (FILE * file, char *message)
{
  fprintf (file, "%s", message);
}
