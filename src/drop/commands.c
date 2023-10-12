#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <argp.h>

#include "commands.h"
#include "../dropanalysis/setDihedral.h"
#include "../dropanalysis/countClashes.h"

const char *commandList[][2]={
  {"setDihedral", "Change a single user-specified dihedral angle for a given residue."},
  {"countClashes", "Counts the number of atomic overlaps according to atomic radii used by Ramachandran."}
};
const int commandListLen = sizeof(commandList)/sizeof(commandList[0]);

void printCommandList()
{

  fprintf(stderr, "These are the commands available:\n\n");
  fprintf(stderr, "EXAMPLE > #: function - description -\n\n");

  for (int i = 0; i < commandListLen; i++)
  {
    fprintf(stderr, "%d: ", i+1);
    for (int j = 0; j < 2; j++) //2 accounts for the command and the description
    {
      fprintf(stderr, "%s - ", commandList[i][j]);
    }
    fprintf(stderr, "\n\n");
  }
}

//assigns -f option and argument to null terminator for the
//child function which won't be able to parse an -f option
void stripFArgv(int argc, char **argv)
{
  for(int i = 0; i < argc; i++)
  {
    if( (strcmp(argv[i], "-f") == 0) ||  (strcmp(argv[i], "--func") == 0))
    {
      argv[i] = "\0";
      argv[i++] = "\0";
    }
  }
}

//assigns all of argv to the null terminator
void stripAllArgv(int argc, char **argv)
{
  for(int i = 0; i < argc; i++)
  {
    argv[i] = "\0";
  }
}

//makes a string of argv to print command line to log
char * makeStringArgv(int argc, char **argv)
{
  char * strng;
  strcpy(strng, argv[0]);
  for(int i = 1; i < argc; i++)
  {
    strcat(strng, " ");
    strcat(strng, argv[i]);
  }
  return strng;
}

bool findCommand(char *arg, int argc, char **argv)
{
  bool found = false;
  //make a string of argv arguments for the log file output
  char *stringArgv = makeStringArgv(argc, argv);

  stripFArgv(argc, argv); //get rid of '-f FUNCTION' to avoid argp errors

  //search available commands or functions
  if ( strcmp(arg, commandList[0][0]) == 0 )
  {
    setDihedral(argc, argv, stringArgv);
    found = true;
  }
  else if ( strcmp(arg, commandList[1][0]) == 0 )
  {
    countClashes(argc, argv, stringArgv);
    found = true;
  }

  if(found)
  {
    stripAllArgv(argc, argv); //all argv is made the null terminator so main will not throw an error
  }

  return found;
}
