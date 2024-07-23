#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <math.h>
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
  double gamma;
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
      case 'g':
      {
        a->gamma = atof(arg);
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
    { "gamma", 'g', "[Boolean]", 0, "Factor for cutoff for allowed states based on steric cutoffs. Default: 1.0" },
    { "", 'f', "", OPTION_HIDDEN, "" }, //gets rid of error for -f flag
    { 0 }
  };

  //DEFAULTS
  struct arguments args = {NULL, "rama.txt", "drop.log", 2, 2, 1, 1.0};
  //parse options
  struct argp vdwScanArgp = { vdwScanOptions, vdwScanParse, 0, 0 };
  argp_parse(&vdwScanArgp, argc, argv, 0, 0, &args);

  struct protein prot;
  FILE *log = fopen(args.log_file, "w");
  processInput(&prot, args.input_file, log, 1, args.bond_matrix, stringArgv);

  //array -> [phi index, psi index, chi1 index, ..., chi5 index]
  //value is -1 for any dihedral not detected
  int *dihedral_indices = findDihedrals(&prot, args.res_number, log);

  //set backbone dihedral angles to top left of Ramachandran distribution
  double phi_change = -179 - prot.dihedrals[dihedral_indices[0]].dihedral_angle ;
  fprintf(log, "Changing dihedral angle %s in residue number %d by %f degrees.\n\n", "phi", args.res_number, phi_change);
  rotateDihedral(&prot, dihedral_indices[0], phi_change, 1);
  prot.dihedrals[dihedral_indices[0]].dihedral_angle = calculateDihedral(&prot, dihedral_indices[0]);

  double psi_change = 179 - prot.dihedrals[dihedral_indices[1]].dihedral_angle ;
  fprintf(log, "Changing dihedral angle %s in residue number %d by %f degrees.\n\n", "psi", args.res_number, psi_change);
  rotateDihedral(&prot, dihedral_indices[1], psi_change, 1);
  prot.dihedrals[dihedral_indices[1]].dihedral_angle = calculateDihedral(&prot, dihedral_indices[1]);

  //for now, just hard code loops for chi1, phi, and psi to do alanine and valine
  /*FILE *output = fopen(args.output_file, "w+");
  double range = 360 / args.resolution;
  double energy_sum = 0;
  double energy = 0;
  int clashes = 0;
  char *message[40];

  //loop for AAs with chi2
  for(int i = 0; i < range; i++)
  {
    //psi loop
    for(int j = 0; j < range; j++)
    {
      double energy_sum = 0;
      double norm = range*range; //normalization for AAs with 2 chi angles
      //chi 1 loop
      for(int k = 0; k < range; k++)
      {
        //chi 2 loop
        for(int m = 0; m < range; m++)
        {
          energy = calculateVDWEnergy(&prot, args.gamma, log);
          if(!isnan(energy))
          {
            energy_sum += energy;
          }
          else
          {
            norm = norm - 1; //subtract from normalization constant to reflect correct number of structures considered
          }
          rotateDihedral(&prot, chi2_index, args.resolution, 0);
          prot.dihedrals[chi2_index].dihedral_angle = calculateDihedral(&prot, chi2_index);

        }

        rotateDihedral(&prot, chi2_index, args.resolution, 0);
        prot.dihedrals[chi2_index].dihedral_angle = calculateDihedral(&prot, chi2_index);

        rotateDihedral(&prot, chi1_index, args.resolution, 0);
        prot.dihedrals[chi1_index].dihedral_angle = calculateDihedral(&prot, chi1_index);

      }

      // if no side chain configurations exist without a steric clash, assign high energy
      if(energy_sum == 0 && norm == 0){
        energy_sum = 500000;
        norm = 1;
      }

      sprintf(message, "%f %f %f\n", prot.dihedrals[phi_index].dihedral_angle, prot.dihedrals[psi_index].dihedral_angle, energy_sum/norm);
      printf("%s", message);
      writeFileLine(output, message);

      rotateDihedral(&prot, psi_index, -args.resolution, 1);
      prot.dihedrals[psi_index].dihedral_angle = calculateDihedral(&prot, psi_index);

    }

    writeFileLine(output, "\n");

    rotateDihedral(&prot, phi_index, args.resolution, 1);
    prot.dihedrals[phi_index].dihedral_angle = calculateDihedral(&prot, phi_index);

  }*/

  //fclose(output);
  fclose(log);
  return;
}
