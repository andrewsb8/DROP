#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "stericClashes.h"
#include "../include/readProtein/readProtein.h"
#include "../include/stericClash/stericClash.h"
#include "../include/fileHandling/fileHandling.h"

struct arguments
{
  char *input_file;
  char *log_file;
  bool list_clashes;
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

  }
  return 0;
}

void stericClashes(int argc, char **argv, char *stringArgv)
{
  struct argp_option stericClashesOptions[] =
  {
    { 0, 0, 0, 0, "./drop -f setDihedral Options:\n" },
    { "input", 'i', "[Input File]", 0, "Input pdb file" },
    { "log", 'l', "[Log File]", 0, "Output log file" },
    { "list_clashes", 'c', "[Boolean]", 0, "Choose to list atomic clashes in log file. 0 does not print list. Default: 1." },
    { 0 }
  };

  //DEFAULTS
  struct arguments args = {NULL, "drop.log", 1};
  //parse options
  struct argp stericClashesArgp = { stericClashesOptions, stericClashesParse, 0, 0 };
  argp_parse(&stericClashesArgp, argc, argv, 0, 0, &args);

  if (fileExists(args.input_file) == -1)
  {
    fprintf(stderr, "ERROR: Input file does not exist. Exiting.\n");
    exit(1);
  }

  //log command line inputs
  FILE *log = fopen(args.log_file, "w");
  fprintf(log, "Command Line: %s\n\n", stringArgv);

  //initialize protein struct and begin analysis
  fprintf(log, "Reading structure file: %s\n\n", args.input_file);
  struct protein prot;
  readPDB(&prot, args.input_file, log);

  fprintf(log, "Done reading structure file: %s\n\n", args.input_file);

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
