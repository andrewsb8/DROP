#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "measureDihedrals.h"
#include "../include/readProtein/readProtein.h"
#include "../include/fileHandling/fileHandling.h"

struct arguments
{
  char *input_file;
  char *output_file;
  char *log_file;
  bool bond_matrix;
};

static int measureDihedralsParse(int key, char *arg, struct argp_state *state)
{
  struct arguments *a = state->input;
  switch(key)
  {
      case 'i':
      {
        a->input_file = arg;
        break;
      }
      case 'o':
      {
        a->output_file = arg;
        break;
      }
      case 'l':
      {
        a->log_file = arg;
        break;
      }
      case 'b':
      {
        a->bond_matrix = atoi(arg);
        break;
      }
      case 'f':
      {
        break;
      }

  }
  return 0;
}

void measureDihedrals(int argc, char **argv, char *stringArgv)
{
  struct argp_option measureDihedralsOptions[] =
  {
    { 0, 0, 0, 0, "./drop -f measureDihedral Options:\n" },
    { "input", 'i', "[Input File]", 0, "Input pdb file" },
    { "output", 'o', "[Output File]", 0, "Output file. Options: see -e for options." },
    { "log", 'l', "[Log File]", 0, "Output log file" },
    { "bond_matrix", 'b', "[Boolean]", 0, "Choose whether or not to print bond matrix to log file. Default: true" },
    { "", 'f', "", OPTION_HIDDEN, "" }, //gets rid of error for -f flag
    { 0 }
  };

  //DEFAULTS
  struct arguments args = {NULL, "output.pdb", "drop.log", 1, NULL};
  //parse options
  struct argp measureDihedralsArgp = { measureDihedralsOptions, measureDihedralsParse, 0, 0 };
  argp_parse(&measureDihedralsArgp, argc, argv, 0, 0, &args);

  struct protein prot;
  FILE *log = fopen(args.log_file, "w");
  processInput(&prot, args.input_file, log, 0, 0, stringArgv);

  fclose(log);
  return;
}
