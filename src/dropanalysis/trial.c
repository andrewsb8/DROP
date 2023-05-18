#include <stdio.h>
#include <argp.h>

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
      case 'i': //require this option to be present, not just argument
      {
        a->input_file = arg;
        break;
      }

      case 'l':
      {
        a->log_file = arg; //not overriding defaults
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
  struct arguments args = {NULL, "log.txt"};
  struct argp trial_argp = { trial_options, trial_parse, 0, 0 };
  argp_parse(&trial_argp, argc, argv, 0, 0, &args);
  printf("logfile name: %s\n", args.log_file);

  //initialize protein struct and begin analysis
  struct protein prot;
  readPDB(&prot, args.input_file);

  return;
}
