#include <stdio.h>
#include <argp.h>
#include <string.h>

#include "trial.h"
#include "../include/readProtein/readProtein.h"

struct arguments
{
  char *input_file;
  char *log_file;
};

static int trial_parse(int key, char *arg, struct argp_state *state)
{
  struct arguments *a = state->input;
  switch(key)
  {
      //require this option to be present, not just argument
      case 'i':
      {
        //TO DO: check if file exists first and terminate if it does not
        a->input_file = arg;
        break;
      }

      case 'l':
      {
        a->log_file = arg;
        break;
      }

  }
  return 0;
}

void trial(int argc, char **argv)
{
  struct argp_option trial_options[] =
  {
    { 0, 0, 0, 0, "./drop -f trial Options:\n" },
    { "input", 'i', "[Input File]", 0, "Input pdb file" },
    { "log", 'l', "[Log File]", 0, "Output log file" },
    { 0 }
  };

  //DEFAULTS
  struct arguments args = {NULL, "drop.log"};

  struct argp trial_argp = { trial_options, trial_parse, 0, 0 };
  argp_parse(&trial_argp, argc, argv, 0, 0, &args);

  //log inputs
  FILE *log = fopen(args.log_file, "w");
  fprintf(log, "Command Line: drop -f trial ");
  for(int k = 1; k < argc-2; k++)
  {
    fprintf(log, "%s ", argv[k]);
  }
  fprintf(log, "\n\n");

  //initialize protein struct and begin analysis
  struct protein prot;
  readPDB(&prot, args.input_file, log);

  return;
}
