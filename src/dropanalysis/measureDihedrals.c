#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "measureDihedrals.h"
#include "../utils/readProtein/readProtein.h"
#include "../utils/fileHandling/fileHandling.h"

struct arguments
{
  char *input_file;
  char *log_file;
  char *dihedral_list;
  bool bond_matrix;
};

static int
measureDihedralsParse (int key, char *arg, struct argp_state *state)
{
  struct arguments *a = state->input;
  switch (key)
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
	case 'b':
	  {
		a->bond_matrix = atoi (arg);
		break;
	  }
	case 'd':
	  {
		a->dihedral_list = arg;
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
measureDihedrals (int argc, char **argv)
{
  struct argp_option measureDihedralsOptions[] = {
	{0, 0, 0, 0, "./drop -f measureDihedrals Options:\n"},
	{"input", 'i', "[Input File]", 0, "Input pdb file"},
	{"log", 'l', "[Log File]", 0, "Output log file"},
	{"bond_matrix", 'b', "[Boolean]", 0,
	 "Choose whether or not to print bond matrix to log file. Default: true"},
	{"dihedral_list", 'd', "[Output File]", 0,
	 "Optional output file name to print the dihedral information to a file that can be modified and used as an input file for setDihedralList."},
	{"", 'f', "", OPTION_HIDDEN, ""},	//gets rid of error for -f flag
	{0}
  };

  //DEFAULTS
  struct arguments args = { NULL, "drop.log", NULL, 1 };
  //parse options
  struct argp measureDihedralsArgp =
	{ measureDihedralsOptions, measureDihedralsParse, 0, 0 };
  argp_parse (&measureDihedralsArgp, argc, argv, 0, 0, &args);

  fprintf(stderr, "%s\n", args.input_file);

  struct protein prot;
  FILE *log = fopen (args.log_file, "w");
  processInput (&prot, args.input_file, log, 0, 0, argc, argv);

  if (args.dihedral_list)
	{
	  FILE *dih_list = fopen (args.dihedral_list, "w");
	  fprintf (dih_list, "#resnum residue dihedral-angle angle\n");
	  for (int i = 0; i < prot.number_of_dihedrals; i++)
		{
		  fprintf (dih_list, "%d %s %s %f\n",
				   prot.dihedrals[i].dihedral_resNum,
				   *prot.dihedrals[i].dihedral_resName,
				   *prot.dihedrals[i].dihedral_angType,
				   prot.dihedrals[i].dihedral_angle);
		}
	  fclose (dih_list);
	}

  fclose (log);
  return;
}
