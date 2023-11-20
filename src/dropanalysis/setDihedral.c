#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "setDihedral.h"
#include "../include/readProtein/readProtein.h"
#include "../include/dihedralRotation/dihedralRotation.h"
#include "../include/fileHandling/fileHandling.h"

struct arguments
{
  char *input_file;
  char *output_file;
  char *log_file;
  int res_number;
  char dih_type[4];
  double angle;
  char *extension;
};

static int setDihedralParse(int key, char *arg, struct argp_state *state)
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
      }
      case 'd':
      {
        strcpy(a->dih_type, arg);
      }
      case 'a':
      {
        a->angle = atoi(arg);
      }
      case 'e':
      {
        a->extension = arg;
      }

  }
  return 0;
}

void setDihedral(int argc, char **argv, char *stringArgv)
{
  struct argp_option setDihedralOptions[] =
  {
    { 0, 0, 0, 0, "./drop -f setDihedral Options:\n" },
    { "input", 'i', "[Input File]", 0, "Input pdb file" },
    { "output", 'o', "[Output File]", 0, "Output xyz file" },
    { "log", 'l', "[Log File]", 0, "Output log file" },
    { "resnum", 'n', "INT", 0, "Residue Number" },
    { "dihtype", 'd', "[Dihedral Type]", 0, "phi, psi" },
    { "dihangle", 'a', "DOUBLE", 0, "Target dihedral angle in degrees" },
    { "extension", 'e', "[Output File Extension]", 0, "Options: pdb, xyz" },
    { 0 }
  };

  //DEFAULTS
  struct arguments args = {NULL, "output.pdb", "drop.log", 1, "phi", 0, "pdb"};
  //parse options
  struct argp setDihedralArgp = { setDihedralOptions, setDihedralParse, 0, 0 };
  argp_parse(&setDihedralArgp, argc, argv, 0, 0, &args);

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

  //stores number associated with the chi or side chain torsion
  //default is zero, but will be read by user input at some point
  int chi_num = 0;
  if(!backbone)
  {
    fprintf(stderr, "Need to do something here for the chi angles.\n\n");
  }

  //calculate the angle change based on current angle and angle defined by command line
  double dih_angle_change = args.angle - prot.dihedrals[index].dihedral_angle ;
  fprintf(log, "Changing dihedral angle %s in residue number %d by %f degrees.\n\n", args.dih_type, args.res_number, dih_angle_change);

  //rotate the dihedral
  rotateDihedral(&prot, index, dih_angle_change, backbone, chi_num);
  prot.dihedrals[index].dihedral_angle = calculateDihedral(&prot, index);

  fprintf(log, "Rotation complete. Please check the accuracy of the operation.\nUser input angle: %f\nAngle after rotation: %f\n\n", args.angle, prot.dihedrals[index].dihedral_angle);

  //print out structure after rotation
  if(strcmp(args.extension, "pdb") == 0)
  {
    writePDB(&prot, args.output_file, 's', 0);
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
