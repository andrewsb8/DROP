#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "setDihedralList.h"
#include "../include/readProtein/readProtein.h"
#include "../include/dihedralRotation/dihedralRotation.h"
#include "../include/fileHandling/fileHandling.h"

struct arguments
{
  char *input_file;
  char *input_dih_list;
  char *output_file;
  char *log_file;
  char *extension;
  bool conect;
  bool bond_matrix;
};

static int setDihedralListParse(int key, char *arg, struct argp_state *state)
{
  struct arguments *a = state->input;
  switch(key)
  {
      case 'i':
      {
        a->input_file = arg;
        break;
      }
      case 'id':
      {
        a->input_dih_list = arg;
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

void setDihedralList(int argc, char **argv, char *stringArgv)
{
  struct argp_option setDihedralListOptions[] =
  {
    { 0, 0, 0, 0, "./drop -f setDihedral Options:\n" },
    { "input", 'i', "[Input File]", 0, "Input pdb file" },
    { "input_dih_list", 'id', "[Input Dihedral List File]", 0, "Input pdb file" },
    { "output", 'o', "[Output File]", 0, "Output file. Options: see -e for options." },
    { "log", 'l', "[Log File]", 0, "Output log file" },
    { "extension", 'e', "[Output File Extension]", 0, "Options: pdb, xyz" },
    { "conect", 'c', "BOOL", 0, "Include CONECT records in PDB. 0 does not print conect. Default: 0." },
    { "bond_matrix", 'b', "[Boolean]", 0, "Choose whether or not to print bond matrix to log file. Default: true" },
    { "", 'f', "", OPTION_HIDDEN, "" }, //gets rid of error for -f flag
    { 0 }
  };

  //DEFAULTS
  struct arguments args = {NULL, NULL, "output.pdb", "drop.log", 1, "phi", 0, "pdb", 0, 1, NULL};
  //parse options
  struct argp setDihedralListArgp = { setDihedralListOptions, setDihedralParse, 0, 0 };
  argp_parse(&setDihedralListArgp, argc, argv, 0, 0, &args);

  if (fileExists(args.input_file) == -1)
  {
    fprintf(stderr, "ERROR: Input file does not exist. Exiting.\n");
    exit(1);
  }
  if (fileExists(args.input_dih_list) == -1)
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
  readPDB(&prot, args.input_file, log, args.bond_matrix);

  fprintf(log, "Done reading structure file: %s\n\n", args.input_file);

  fprintf(log, "Starting structure manipulation from dihedral list: %s\n", args.input_dih_list);
  FILE *dih_list;
  dih_list = fopen(args.input_dih_list, "r");
  while((read = getline(&line,&len,fp)) != -1)
  {
    //parse line - residue number, dihedral type, dihedral angle

    //find dihedral to change based on user input
    int index = findDihedral(&prot, args.res_number, args.dih_type);
    if (index == -1)
    {
      fprintf(log, "Error: dihedral angle %s in residue number %d was not found.\n\n", args.dih_type, args.res_number);
      fprintf(stderr, "Error: dihedral angle %s in residue number %d was not found.\n\n", args.dih_type, args.res_number);
      exit(1);
    }
    else
    {
      fprintf(log, "Found dihedral number: %d\n\n", index);
    }

    //is the angle being changed the backbone or side chain?
    bool backbone;
    if( strcmp(args.dih_type, "phi") == 0 || strcmp(args.dih_type, "psi") == 0 )
    {
      backbone = true;
    }
    else
    {
      backbone = false;
    }

    //calculate the angle change based on current angle and angle defined by command line
    double dih_angle_change = args.angle - prot.dihedrals[index].dihedral_angle ;
    fprintf(log, "Changing dihedral angle %s in residue number %d by %f degrees.\n\n", args.dih_type, args.res_number, dih_angle_change);

    //rotate the dihedral
    rotateDihedral(&prot, index, dih_angle_change, backbone);
    prot.dihedrals[index].dihedral_angle = calculateDihedral(&prot, index);

    fprintf(log, "Rotation complete. Please check the accuracy of the operation.\nUser input angle: %f\nAngle after rotation: %f\n\n", args.angle, prot.dihedrals[index].dihedral_angle);

  }
  //print out structure after rotation
  if(strcmp(args.extension, "pdb") == 0)
  {
    writePDB(&prot, args.output_file, 's', 0, args.conect);
  }
  else if(strcmp(args.extension, "xyz") == 0)
  {
    writeXYZ(&prot, args.output_file, "Frame 1", 's', 0, 0);
  }
  else
  {
    fprintf(log, "Error: File extension for output not recognized.\n");
    fprintf(stderr, "Error: File extension for output not recognized.\n");
    exit(1);
  }

  fclose(log);
  return;
}
