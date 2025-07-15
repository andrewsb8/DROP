#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "measureDihedral.h"
#include "../utils/readProtein/readProtein.h"
#include "../utils/fileHandling/fileHandling.h"
#include "../utils/dihedralRotation/dihedralRotation.h"
#include "../utils/logging/logging.h"

struct arguments {
	char *input_file;
	char *log_file;
	int res_number;
	char dih_type[5];
};

static int measureDihedralParse(int key, char *arg,
								 struct argp_state *state)
{
	struct arguments *a = state->input;
	switch (key) {
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
	case 'f':
		{
			break;
		}

	}
	return 0;
}

void measureDihedral(int argc, char **argv)
{
	struct argp_option measureDihedralOptions[] = {
		{ 0, 0, 0, 0, "./drop measureDihedrals Options:\n" },
		{ "input", 'i', "[Input File]", 0, "Input pdb file" },
		{ "log", 'l', "[Log File]", 0, "Output log file" },
		{ "resnum", 'n', "INT", 0, "Residue Number" },
		{ "dihtype", 'd', "[Dihedral Type]", 0,
		 "phi, psi, chi1, chi2, ... , chi5" },
		{ 0 }
	};

	//DEFAULTS
	struct arguments args = { NULL, "drop.log", 2, "phi" };
	//parse options
	struct argp measureDihedralArgp =
		{ measureDihedralOptions, measureDihedralParse, 0, 0 };
	argp_parse(&measureDihedralArgp, argc, argv, 0, 0, &args);

	struct protein prot;
	FILE *log = fopen(args.log_file, "w");
	processInput(&prot, args.input_file, log, 0, 0, argc, argv);

	//find dihedral to change based on user input
	int index = findDihedral(&prot, args.res_number, args.dih_type);
	if (index == -1) {
		char message[70];
		sprintf(message,
				"ERROR: Dihedral type %s in residue %d not found. Exiting.\n",
				args.dih_type, args.res_number);
		drop_fatal(log, message);
	} else {
        printf("%f\n", prot.dihedrals[index].dihedral_angle);
	}

	fclose(log);
	return;
}
