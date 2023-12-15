#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "stericScan.h"
#include "../include/readProtein/readProtein.h"
#include "../include/dihedralRotation/dihedralRotation.h"
#include "../include/fileHandling/fileHandling.h"

struct arguments
{
  char *input_file;
  char *output_file;
  char *log_file;
  int res_number;
  char dih_type[5];
  double angle;
  char *extension;
  bool conect;
  bool bond_matrix;
};

static int stericScanParse(int key, char *arg, struct argp_state *state)
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
      case 'n':
      {
        a->res_number = atoi(arg);
        break;
      }
      case 'd':
      {
        strcpy(a->dih_type, arg);
        break;
      }
      case 'a':
      {
        a->angle = atof(arg);
        break;
      }
      case 'e':
      {
        a->extension = arg;
        break;
      }
      case 'c':
      {
        a->conect = atoi(arg);
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

void stericScan(int argc, char **argv, char *stringArgv)
{
  struct argp_option stericScanOptions[] =
  {
    { 0, 0, 0, 0, "./drop -f setDihedral Options:\n" },
    { "input", 'i', "[Input File]", 0, "Input pdb file" },
    { "output", 'o', "[Output File]", 0, "Output file. Options: see -e for options." },
    { "log", 'l', "[Log File]", 0, "Output log file" },
    { "resnum", 'n', "INT", 0, "Residue Number" },
    { "dihtype", 'd', "[Dihedral Type]", 0, "phi, psi" },
    { "dihangle", 'a', "FLOAT", 0, "Target dihedral angle in degrees" },
    { "extension", 'e', "[Output File Extension]", 0, "Options: pdb, xyz" },
    { "conect", 'c', "BOOL", 0, "Include CONECT records in PDB. 0 does not print conect. Default: 0." },
    { "bond_matrix", 'b', "[Boolean]", 0, "Choose whether or not to print bond matrix to log file. Default: true" },
    { "", 'f', "", OPTION_HIDDEN, "" }, //gets rid of error for -f flag
    { 0 }
  };

  //DEFAULTS
  struct arguments args = {NULL, "output.pdb", "drop.log", 1, "phi", 0, "pdb", 0, 1, NULL};
  //parse options
  struct argp stericScanArgp = { stericScanOptions, stericScanParse, 0, 0 };
  argp_parse(&stericScanArgp, argc, argv, 0, 0, &args);

  struct protein prot;
  FILE *log = fopen(args.log_file, "w");
  inputInfo(&prot, args.input_file, log, args.bond_matrix, stringArgv);

  fclose(log);
  return;
}
