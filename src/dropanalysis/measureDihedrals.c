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
    { 0, 0, 0, 0, "./drop -f setDihedral Options:\n" },
    { "input", 'i', "[Input File]", 0, "Input pdb file" },
    { "output", 'o', "[Output File]", 0, "Output file. Options: see -e for options." },
    { "log", 'l', "[Log File]", 0, "Output log file" },
    { "", 'f', "", OPTION_HIDDEN, "" }, //gets rid of error for -f flag
    { 0 }
  };

  //DEFAULTS
  struct arguments args = {NULL, "output.pdb", "drop.log", NULL};
  //parse options
  struct argp measureDihedralsArgp = { measureDihedralsOptions, measureDihedralsParse, 0, 0 };
  argp_parse(&measureDihedralsArgp, argc, argv, 0, 0, &args);

  if (fileExists(args.input_file) == -1)
  {
    fprintf(stderr, "ERROR: Input file does not exist. Exiting.\n");
    exit(1);
  }

  struct protein prot;
  FILE *log = fopen(args.log_file, "w");
  processInput(&prot, args.input_file, log, 0, args.bond_matrix, stringArgv);

  fclose(log);
  return;
}
