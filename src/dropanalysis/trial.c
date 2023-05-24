#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>

#include "trial.h"
#include "../include/readProtein/readProtein.h"
#include "../include/dihedralRotation/dihedralRotation.h"
#include "../include/fileHandling/fileHandling.h"

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

  }
  return 0;
}

void trial(int argc, char **argv, char *stringArgv)
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
  //parse options
  struct argp trial_argp = { trial_options, trial_parse, 0, 0 };
  argp_parse(&trial_argp, argc, argv, 0, 0, &args);

  //if input file does not exist, exit
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

  //current tests - going to remove and replace with whatever function trial becomes
  rotateDihedral(&prot, 0, prot.dihedrals[0].dihedral_angle, 2, 1, 0);
  rotateDihedral(&prot, 1, prot.dihedrals[1].dihedral_angle, 2, 1, 0);
  rotateDihedral(&prot, 2, prot.dihedrals[2].dihedral_angle, 2, 0, 1);
  rotateDihedral(&prot, 3, prot.dihedrals[3].dihedral_angle, 2, 0, 2);

  for(int i = 0; i < prot.number_of_dihedrals; i++)
  {
    prot.dihedrals[i].dihedral_angle = calculateDihedral(&prot, i);
    printf("%f\n", prot.dihedrals[i].dihedral_angle);
  }
  printf("\n");

  fclose(log);
  return;
}
