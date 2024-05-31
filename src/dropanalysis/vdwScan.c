#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "vdwScan.h"
#include "../include/readProtein/readProtein.h"
#include "../include/dihedralRotation/dihedralRotation.h"
#include "../include/vdwEnergy/vdwEnergy.h"
#include "../include/fileHandling/fileHandling.h"

struct arguments
{
  char *input_file;
  char *output_file;
  char *log_file;
  int res_number;
  double resolution;
  bool bond_matrix;
};

static int vdwScanParse(int key, char *arg, struct argp_state *state)
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
      case 'r':
      {
        a->resolution = atof(arg);
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

void vdwScan(int argc, char **argv, char *stringArgv)
{
  struct argp_option vdwScanOptions[] =
  {
    { 0, 0, 0, 0, "./drop -f stericScan Options:\n" },
    { "input", 'i', "[Input File]", 0, "Input pdb file" },
    { "output", 'o', "[Output File]", 0, "Output .txt file with three columns: phi, psi, average number of clashes." },
    { "log", 'l', "[Log File]", 0, "Output log file" },
    { "resnum", 'n', "INT", 0, "Residue Number for analysis. Default: 2 (first amino acid will typically not have both backbone angles defined)." },
    { "resolution", 'r', "DOUBLE", 0, "Resolution of Ramachandran space (and therefore dihedral rotation magnitude). Default: 2 deg" },
    { "bond_matrix", 'b', "[Boolean]", 0, "Choose whether or not to print bond matrix to log file. Default: true" },
    { "", 'f', "", OPTION_HIDDEN, "" }, //gets rid of error for -f flag
    { 0 }
  };

  //DEFAULTS
  struct arguments args = {NULL, "rama.txt", "drop.log", 2, 2, 1};
  //parse options
  struct argp vdwScanArgp = { vdwScanOptions, vdwScanParse, 0, 0 };
  argp_parse(&vdwScanArgp, argc, argv, 0, 0, &args);

  struct protein prot;
  FILE *log = fopen(args.log_file, "w");
  processInput(&prot, args.input_file, log, 1, args.bond_matrix, stringArgv);

  /*TO DO
  how to have loops depend on the number of chi angles?
  rotate "highest" chi angle and have nested lists do the rest. check stericClashes at each step, then continue rotating
  */

  //set backbone dihedral angles to top left of Ramachandran distribution
  //might want a conditional here to avoid unnecessary processing...
  int phi_index = findDihedral(&prot, args.res_number, "phi", log);
  double phi_change = -179 - prot.dihedrals[phi_index].dihedral_angle ;
  fprintf(log, "Changing dihedral angle %s in residue number %d by %f degrees.\n\n", "phi", args.res_number, phi_change);
  rotateDihedral(&prot, phi_index, phi_change, 1);
  prot.dihedrals[phi_index].dihedral_angle = calculateDihedral(&prot, phi_index);

  int psi_index = findDihedral(&prot, args.res_number, "psi", log);
  double psi_change = 179 - prot.dihedrals[psi_index].dihedral_angle ;
  fprintf(log, "Changing dihedral angle %s in residue number %d by %f degrees.\n\n", "psi", args.res_number, psi_change);
  rotateDihedral(&prot, psi_index, psi_change, 1);
  prot.dihedrals[psi_index].dihedral_angle = calculateDihedral(&prot, psi_index);

  int chi1_index = findDihedral(&prot, args.res_number, "chi1", log);
  int chi2_index = findDihedral(&prot, args.res_number, "chi2", log);

  //for now, just hard code loops for chi1, phi, and psi to do alanine and valine
  FILE *output = fopen(args.output_file, "w+");
  double range = 360 / args.resolution;
  double energy_sum = 0;
  char *message[40];

  /*//loop for AAs with chi2
  for(int i = 0; i < range; i++)
  {
    //psi loop
    for(int j = 0; j < range; j++)
    {
      double clashes = 0;
      //chi 1 loop
      for(int k = 0; k < range; k++)
      {
        //chi 2 loop
        for(int m = 0; m < range; m++)
        {

          energy_sum += calculateVDWEnergy(&prot, log);
          rotateDihedral(&prot, chi2_index, args.resolution, 0);
          prot.dihedrals[chi2_index].dihedral_angle = calculateDihedral(&prot, chi2_index);

        }

        rotateDihedral(&prot, chi2_index, args.resolution, 0);
        prot.dihedrals[chi2_index].dihedral_angle = calculateDihedral(&prot, chi2_index);

        rotateDihedral(&prot, chi1_index, args.resolution, 0);
        prot.dihedrals[chi1_index].dihedral_angle = calculateDihedral(&prot, chi1_index);

      }

      sprintf(message, "%f %f %f\n", prot.dihedrals[phi_index].dihedral_angle, prot.dihedrals[psi_index].dihedral_angle, energy_sum/(range*range));
      printf("%s", message);
      writeFileLine(output, message);

      rotateDihedral(&prot, psi_index, -args.resolution, 1);
      prot.dihedrals[psi_index].dihedral_angle = calculateDihedral(&prot, psi_index);

    }

    writeFileLine(output, "\n");

    rotateDihedral(&prot, phi_index, args.resolution, 1);
    prot.dihedrals[phi_index].dihedral_angle = calculateDihedral(&prot, phi_index);

  }*/

  //loop structure for only amino acids with a chi1
  //phi loop
  for(int i = 0; i < range; i++)
  {
    //psi loop
    for(int j = 0; j < range; j++)
    {
      double clashes = 0;
      //chi 1 loop
      for (int k = 0; k < range; k++)
      {
        energy_sum += calculateVDWEnergy(&prot, log);
        rotateDihedral(&prot, chi1_index, args.resolution, 0);
        prot.dihedrals[chi1_index].dihedral_angle = calculateDihedral(&prot, chi1_index);
      }

      sprintf(message, "%f %f %f\n", prot.dihedrals[phi_index].dihedral_angle, prot.dihedrals[psi_index].dihedral_angle, energy_sum/(range));
      printf("%s", message);
      //writeFileLine(output, message);

      rotateDihedral(&prot, chi1_index, args.resolution, 0);
      prot.dihedrals[chi1_index].dihedral_angle = calculateDihedral(&prot, chi1_index);

      rotateDihedral(&prot, psi_index, -args.resolution, 1);
      prot.dihedrals[psi_index].dihedral_angle = calculateDihedral(&prot, psi_index);

    }

    //writeFileLine(output, "\n");

    rotateDihedral(&prot, phi_index, args.resolution, 1);
    prot.dihedrals[phi_index].dihedral_angle = calculateDihedral(&prot, phi_index);

  }

  fclose(output);
  fclose(log);
  return;
}
