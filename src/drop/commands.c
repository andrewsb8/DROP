#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <argp.h>

#include "commands.h"
#include "../dropanalysis/setDihedral.h"
#include "../dropanalysis/measureDihedrals.h"
#include "../dropanalysis/stericClashes.h"

const char *commandList[][2]={
  {"measureDihedrals", "Parses structure and provides log with structure and dihedral information."},
  {"setDihedral", "Change a single user-specified dihedral angle for a given residue."},
  {"setDihedralList", "Change several dihedrals from a user-provided input file with a list of dihedrals and angles."},
  {"stericClashes", "Counts the number of atomic overlaps according to atomic radii used by Ramachandran."}
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

bool findCommand(char *func, int argc, char **argv)
{
  bool found = false;
  //make a string of argv arguments for the log file output
  char *stringArgv = makeStringArgv(argc, argv);

  //search available commands or functions
  if ( strcmp(func, commandList[0][0]) == 0 )
  {
    measureDihedrals(argc, argv, stringArgv);
    found = true;
  }
  else if ( strcmp(func, commandList[1][0]) == 0 )
  {
    setDihedral(argc, argv, stringArgv);
    found = true;
  }
  else if ( strcmp(func, commandList[2][0]) == 0 )
  {
    setDihedralList(argc, argv, stringArgv);
    found = true;
  }
  else if ( strcmp(func, commandList[3][0]) == 0 )
  {
    stericClashes(argc, argv, stringArgv);
    found = true;
  }

  if(found)
  {
    stripAllArgv(argc, argv); //all argv is made the null terminator so parser will not throw an error passing back to main
  }

  return found;
}
