#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "setDihedralList.h"
#include "../utils/readProtein/readProtein.h"
#include "../utils/dihedralRotation/dihedralRotation.h"
#include "../utils/fileHandling/fileHandling.h"
#include "../utils/logging/logging.h"

struct arguments {
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
	switch (key) {
	case 'i':
		{
			a->input_file = arg;
			break;
		}
	case 'd':
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

void setDihedralList(int argc, char **argv)
{
	struct argp_option setDihedralListOptions[] = {
		{0, 0, 0, 0, "./drop setDihedralList Options:\n"},
		{"input", 'i', "[Input File]", 0, "Input pdb file"},
		{"input_dih_list", 'd', "[Input Dihedral List File]", 0,
		 "Input pdb file"},
		{"output", 'o', "[Output File]", 0,
		 "Output file. Options: see -e for options."},
		{"log", 'l', "[Log File]", 0, "Output log file"},
		{"extension", 'e', "[Output File Extension]", 0,
		 "Options: pdb, xyz"},
		{"conect", 'c', "BOOL", 0,
		 "Include CONECT records in PDB. 0 does not print conect. Default: 0."},
		{"bond_matrix", 'b', "[Boolean]", 0,
		 "Choose whether or not to print bond matrix to log file. Default: true"},
		{0}
	};

	//DEFAULTS
	struct arguments args =
	    { NULL, NULL, "output.pdb", "drop.log", "pdb", 0, 1 };
	//parse options
	struct argp setDihedralListArgp =
	    { setDihedralListOptions, setDihedralListParse, 0, 0 };
	argp_parse(&setDihedralListArgp, argc, argv, 0, 0, &args);

	if (fileExists(args.input_dih_list) == -1) {
		fprintf(stderr,
			"ERROR: Input dihedral list does not exist. Exiting.\n");
		exit(1);
	}

	struct protein prot;
	FILE *log = fopen(args.log_file, "w");
	processInput(&prot, args.input_file, log, 0, 0, argc, argv);

	fprintf(log, "Starting structure manipulation from dihedral list: %s\n",
		args.input_dih_list);
	FILE *dih_list;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	dih_list = fopen(args.input_dih_list, "r");
	while ((read = getline(&line, &len, dih_list)) != -1) {
		if (line[0] != '#') {
			//parse space separated line - residue number dihedral type dihedral angle
			char *line_split = strtok(line, " \t");
			char stringT[4][10];
			int count = 0;
			while (line_split != NULL) {
				strcpy(stringT[count], line_split);
				count++;
				line_split = strtok(NULL, " \t");
			}

			if (count != 4) {
				char message[120];
				sprintf(message,
					"ERROR: Line in dihedral list has more or fewer than 4 elements in line %s of the input dihedral list. Exiting.\n",
					line);
				drop_fatal(log, message);
			}

			int res_number = atoi(stringT[0]);
			char *residue = stringT[1];
			char *dih_type = stringT[2];
			float angle = atof(stringT[3]);

			//find dihedral to change based on user input
			int index = findDihedral(&prot, res_number, dih_type);
			if (index == -1) {
				char message[70];
				sprintf(message,
					"ERROR: Dihedral type %s in residue %d on line %s not found. Exiting.\n",
					dih_type, res_number, line);
				drop_fatal(log, message);
			}
			//is the angle being changed the backbone or side chain?
			bool backbone;
			if (strcmp(dih_type, "phi") == 0
			    || strcmp(dih_type, "psi") == 0) {
				backbone = true;
			} else {
				backbone = false;
			}

			//calculate the angle change based on current angle and angle defined by command line
			double dih_angle_change =
			    angle - prot.dihedrals[index].dihedral_angle;
			fprintf(log,
				"Changing dihedral angle %s in residue number %d by %f degrees.\n\n",
				dih_type, res_number, dih_angle_change);

			//rotate the dihedral
			rotateDihedral(&prot, index, dih_angle_change,
				       backbone);
			prot.dihedrals[index].dihedral_angle =
			    calculateDihedral(&prot, index);

			fprintf(log,
				"Rotation complete. Please check the accuracy of the operation.\nUser input angle: %f\nAngle after rotation: %f\n\n",
				angle, prot.dihedrals[index].dihedral_angle);
		}
	}
	//print out structure after rotation
	if (strcmp(args.extension, "pdb") == 0) {
		writePDB(&prot, args.output_file, 's', 0, args.conect);
	} else if (strcmp(args.extension, "xyz") == 0) {
		writeXYZ(&prot, args.output_file, "Frame 1", 's', 0, 0);
	} else {
		drop_fatal(log,
			   "Error: File extension for output not recognized.\n");
	}

	fclose(log);
	return;
}
