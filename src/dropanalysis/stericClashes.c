#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "stericClashes.h"
#include "../utils/readProtein/readProtein.h"
#include "../utils/stericClash/stericClash.h"
#include "../utils/fileHandling/fileHandling.h"

struct arguments
{
  char *input_file;
  char *log_file;
  bool list_clashes;
  bool bond_matrix;
};

static int stericClashesParse(int key, char *arg, struct argp_state *state)
{
  struct arguments *a = state->input;
  switch(key)
  {
      case 'i':
      {
        a->input_file = arg;
        break;
      }
      case 'l':
      {
        a->log_file = arg;
        break;
      }
      case 'c':
      {
        a->list_clashes = atoi(arg);
        break;
      }
      case 'b':
      {
        a->bond_matrix = atoi(arg);
      }

  }
  return 0;
}

void stericClashes(int argc, char **argv, char *stringArgv)
{
  struct argp_option stericClashesOptions[] =
  {
    { 0, 0, 0, 0, "./drop -f stericClashes Options:\n" },
    { "input", 'i', "[Input File]", 0, "Input pdb file" },
    { "log", 'l', "[Log File]", 0, "Output log file" },
    { "bond_matrix", 'b', "[Boolean]", 0, "Choose whether or not to print bond matrix. Default: true" },
    { "list_clashes", 'c', "[Boolean]", 0, "Choose to list atomic clashes in log file. 0 does not print list. Default: 1." },
    { "", 'f', "", OPTION_HIDDEN, "" }, //gets rid of error for -f flag
    { 0 }
  };

  //DEFAULTS
  struct arguments args = {NULL, "drop.log", 1, 1};
  //parse options
  struct argp stericClashesArgp = { stericClashesOptions, stericClashesParse, 0, 0, NULL };
  argp_parse(&stericClashesArgp, argc, argv, 0, 0, &args);

  struct protein prot;
  FILE *log = fopen(args.log_file, "w");
  processInput(&prot, args.input_file, log, 1, args.bond_matrix, stringArgv);

  if(!args.list_clashes)
  {
    fprintf(log, "Individual atomic clashes will not be reported.\n");
  }
  else
  {
    fprintf(log, "Individual atomic clashes will be reported.\n");
  }

  int num_clashes = countClashes(&prot, log, args.list_clashes);

  fprintf(log, "Number of steric clashes: %d\n\n", num_clashes);

  fclose(log);
  return;
}
