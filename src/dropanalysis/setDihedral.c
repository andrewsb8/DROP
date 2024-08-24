#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "setDihedral.h"
#include "../utils/readProtein/readProtein.h"
#include "../utils/dihedralRotation/dihedralRotation.h"
#include "../utils/fileHandling/fileHandling.h"
#include "../utils/logging/logging.h"

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

static int
setDihedralParse (int key, char *arg, struct argp_state *state)
{
  struct arguments *a = state->input;
  switch (key)
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
		a->res_number = atoi (arg);
		break;
	  }
	case 'd':
	  {
		strcpy (a->dih_type, arg);
		break;
	  }
	case 'a':
	  {
		a->angle = atof (arg);
		break;
	  }
	case 'e':
	  {
		a->extension = arg;
		break;
	  }
	case 'c':
	  {
		a->conect = atoi (arg);
		break;
	  }
	case 'b':
	  {
		a->bond_matrix = atoi (arg);
		break;
	  }
	case 'f':
	  {
		break;
	  }

	}
  return 0;
}

void
setDihedral (int argc, char **argv)
{
  struct argp_option setDihedralOptions[] = {
	{0, 0, 0, 0, "./drop -f setDihedral Options:\n"},
	{"input", 'i', "[Input File]", 0, "Input pdb file"},
	{"output", 'o', "[Output File]", 0,
	 "Output file. Options: see -e for options."},
	{"log", 'l', "[Log File]", 0, "Output log file"},
	{"resnum", 'n', "INT", 0, "Residue Number"},
	{"dihtype", 'd', "[Dihedral Type]", 0,
	 "phi, psi, chi1, chi2, ... , chi5"},
	{"dihangle", 'a', "FLOAT", 0, "Target dihedral angle in degrees"},
	{"extension", 'e', "[Output File Extension]", 0, "Options: pdb, xyz"},
	{"conect", 'c', "BOOL", 0,
	 "Include CONECT records in PDB. 0 does not print conect. Default: 0."},
	{"bond_matrix", 'b', "[Boolean]", 0,
	 "Choose whether or not to print bond matrix to log file. Default: true"},
	{"", 'f', "", OPTION_HIDDEN, ""},	//gets rid of error for -f flag
	{0}
  };

  //DEFAULTS
  struct arguments args =
	{ NULL, "output.pdb", "drop.log", 1, "phi", 0, "pdb", 0, 1 };
  //parse options
  struct argp setDihedralArgp =
	{ setDihedralOptions, setDihedralParse, 0, 0 };
  argp_parse (&setDihedralArgp, argc, argv, 0, 0, &args);

  struct protein prot;
  FILE *log = fopen (args.log_file, "w");
  processInput (&prot, args.input_file, log, 0, 0, argc, argv);

  //find dihedral to change based on user input
  int index = findDihedral (&prot, args.res_number, args.dih_type);
  if (index == -1)
	{
	  char message[70];
	  sprintf (message,
			   "ERROR: Dihedral type %s in residue %d not found. Exiting.\n",
			   args.dih_type, args.res_number);
	  drop_fatal (log, message);
	}

  //is the angle being changed the backbone or side chain?
  bool backbone;
  if (strcmp (args.dih_type, "phi") == 0
	  || strcmp (args.dih_type, "psi") == 0)
	{
	  backbone = true;
	}
  else
	{
	  backbone = false;
	}

  //calculate the angle change based on current angle and angle defined by command line
  double dih_angle_change = args.angle - prot.dihedrals[index].dihedral_angle;
  fprintf (log,
		   "Changing dihedral angle %s in residue number %d by %f degrees.\n\n",
		   args.dih_type, args.res_number, dih_angle_change);

  //rotate the dihedral
  rotateDihedral (&prot, index, dih_angle_change, backbone);
  prot.dihedrals[index].dihedral_angle = calculateDihedral (&prot, index);

  fprintf (log,
		   "Rotation complete. Please check the accuracy of the operation.\nUser input angle: %f\nAngle after rotation: %f\n\n",
		   args.angle, prot.dihedrals[index].dihedral_angle);

  //print out structure after rotation
  if (strcmp (args.extension, "pdb") == 0)
	{
	  writePDB (&prot, args.output_file, 's', 0, args.conect);
	}
  else if (strcmp (args.extension, "xyz") == 0)
	{
	  writeXYZ (&prot, args.output_file, "Frame 1", 's', 0, 0);
	}
  else
	{
	  drop_fatal (log, "Error: File extension for output not recognized.\n");
	  exit (1);
	}

  fclose (log);
  return;
}
