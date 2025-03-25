#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "vdwScanSC.h"
#include "../utils/readProtein/readProtein.h"
#include "../utils/dihedralRotation/dihedralRotation.h"
#include "../utils/vdwEnergy/vdwEnergy.h"
#include "../utils/fileHandling/fileHandling.h"
#include "../utils/logging/logging.h"

struct arguments {
	char *input_file;
	char *output_file;
	char *log_file;
	int res_number;
	double resolution;
	bool bond_matrix;
	double gamma;
};

static int vdwScanSCParse(int key, char *arg, struct argp_state *state)
{
	struct arguments *a = state->input;
	switch (key) {
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

void vdwScanSC(int argc, char **argv)
{
	struct argp_option vdwScanOptions[] = {
		{ 0, 0, 0, 0, "./drop vdwScanSC Options:\n" },
		{ "input", 'i', "[Input File]", 0, "Input pdb file" },
		{ "output", 'o', "[Output File]", 0,
		 "Output .txt file with three columns: phi, psi, average Lennard-Jones energy in kJ/mol"
		 },
		{ "log", 'l', "[Log File]", 0, "Output log file" },
		{ "resnum", 'n', "INT", 0,
		 "Residue Number for analysis. Default: 2 (first amino acid will typically not have both backbone angles defined)."
		 },
		{ "resolution", 'r', "DOUBLE", 0,
		 "Resolution of Ramachandran space (and therefore dihedral rotation magnitude). Default: 2 deg"
		 },
		{ "bond_matrix", 'b', "[Boolean]", 0,
		 "Choose whether or not to print bond matrix to log file. Default: true"
		 },
		{ "gamma", 'g', "[Float]", 0,
		 "Factor for cutoff for allowed states based on steric cutoffs. Default: 1.0"
		 },
		{ 0 }
	};

	//DEFAULTS
	struct arguments args = { NULL, "rama-sc.txt", "drop.log", 2, 2, 1, 1.0 };
	//parse options
	struct argp vdwScanArgp = { vdwScanOptions, vdwScanSCParse, 0, 0 };
	argp_parse(&vdwScanArgp, argc, argv, 0, 0, &args);

	struct protein prot;
	FILE *log = fopen(args.log_file, "w");
	processInput(&prot, args.input_file, log, 1, args.bond_matrix, argc,
				 argv);

	//array -> [phi index, psi index, chi1 index, ..., chi5 index]
	//value is -1 for any dihedral not detected
	int *dihedral_indices = findDihedrals(&prot, args.res_number, log);
	if (dihedral_indices[0] == -1 || dihedral_indices[1] == -1) {
		char message[141];
		sprintf(message,
				"ERROR: One or both backbone dihedral angles not found in residue %d not found. Need to detect valid backbone to validate structure. Exiting.\n",
				args.res_number);
		drop_fatal(log, message);
	}

	//for now, just hard code loops for chi1, phi, and psi to do alanine and valine
	FILE *output = fopen(args.output_file, "w+");
	double range = 360 / args.resolution;
	double norm_factor = 1;		//normalization factor for averaging

	int index_start = -1;
	//if chi2 exists, want to normalize over states of chi3 through chi5, for as many as exist
	//if chi2 does not exist, then norm_factor will be 1, and each energy will be calculated individually
	if (dihedral_indices[3] != -1) {
		// start normalizing at chi3 (h=4)
        for (int h = 4; h < sizeDihedralList; h++) {
      		if (dihedral_indices[h] != -1)	//if index is -1, dihedral not detected and not included in normalization
      		{
     			norm_factor *= range;
      		}
       	}
	}
	double energy_sum = 0;
	double energy = 0;
	char message[50];
	int norm;

	//loop through dihedrals - TO DO: maybe do recursion here?
	//chi 1 loop
	for (int k = 0; k < range; k++) {
	    // if chi2 does not exist, reset norm factor here
	    if (dihedral_indices[3] == -1) {
			norm = norm_factor;
			energy_sum = 0;
		}
		//chi 2 loop
		for (int m = 0; m < range; m++) {
            if (dihedral_indices[3] != -1) {
                norm = norm_factor;
                energy_sum = 0;
            }
			//chi3 loop
			for (int n = 0; n < range; n++) {
				//chi4 loop
				for (int p = 0; p < range; p++) {
					//chi5 loop
					for (int q = 0; q < range; q++) {
						energy =
							calculateVDWEnergy(&prot, args.gamma);
						if (!isnan(energy)) {
							energy_sum += energy;
						} else {
							norm = norm - 1;	//subtract from normalization constant to reflect correct number of structures considered
						}
						if (dihedral_indices[6] != -1) {
							rotateDihedral
								(&prot,
								 dihedral_indices
								 [6], args.resolution, 0);
							prot.dihedrals
								[dihedral_indices
								 [6]].dihedral_angle
								=
								calculateDihedral
								(&prot, dihedral_indices[6]);
						} else {
							break;
						}

					}

					if (dihedral_indices[5]
						!= -1) {
						rotateDihedral
							(&prot,
							 dihedral_indices
							 [5], args.resolution, 0);
						prot.dihedrals
							[dihedral_indices
							 [5]].dihedral_angle
							=
							calculateDihedral
							(&prot, dihedral_indices[5]);
					} else {
						break;
					}

				}

				if (dihedral_indices[4] != -1) {
					rotateDihedral(&prot,
								   dihedral_indices
								   [4], args.resolution, 0);
					prot.dihedrals
						[dihedral_indices
						 [4]].dihedral_angle =
						calculateDihedral
						(&prot, dihedral_indices[4]);
				} else {
					break;
				}
			}

			if (dihedral_indices[3] != -1) {
				rotateDihedral(&prot,
							   dihedral_indices
							   [3], args.resolution, 0);
				prot.dihedrals[dihedral_indices
							   [3]].dihedral_angle =
					calculateDihedral(&prot, dihedral_indices[3]);
			} else {
				break;
			}
		}

		// if no side chain configurations exist without a steric clash, assign high energy
		if (energy_sum == 0 && norm == 0) {
			//energy_sum = 500000;
			energy_sum = 10;
			norm = 1;
		}

		// if no chi2, only plot chi1 angle and energy
		if (dihedral_indices[3] == -1) {
            sprintf(message, "%f %f\n",
			prot.dihedrals[dihedral_indices[2]].dihedral_angle,
			energy_sum / norm);
		} else { // plot chi1, chi2, energy
            sprintf(message, "%f %f %f\n",
			prot.dihedrals[dihedral_indices[2]].dihedral_angle,
			prot.dihedrals[dihedral_indices[3]].dihedral_angle,
			energy_sum / norm);
		}
		printf("%s", message);
		writeFileLine(output, message);

		if (dihedral_indices[2] != -1) {
			rotateDihedral(&prot,
						   dihedral_indices[2],
						   args.resolution, 0);
			prot.dihedrals[dihedral_indices[2]].dihedral_angle =
				calculateDihedral(&prot, dihedral_indices[2]);
		} else {
			break;
		}
	}

	fclose(output);
	fclose(log);
	return;
}
