#include <stdio.h>
#include <argp.h>

#include "trial.h"
#include "../include/readProtein/readProtein.h"

char *filename; //change to an arguments struct and keep global

static int trial_parse(int key, char *arg, struct argp_state *state)
{
  switch(key)
  {
      case 'i': //require this option to be present, not just argument
      {
        fprintf(stderr, "this is the i flag in trial.c.\n");
        filename = arg;
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
    { "input", 'i', "[Input File]", 0, "input pdb file" },
    { 0 }
  };

  struct argp trial_argp = { trial_options, trial_parse, 0, 0 };
  argp_parse(&trial_argp, argc, argv, 0, 0, 0);

  printf("filename: %s\n", filename);
  struct protein prot;
  readPDB(&prot, filename);

  return;
}
