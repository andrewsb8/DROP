#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <argp.h>

#include "commands.h"
#include "../dropanalysis/measureDihedrals.h"
#include "../dropanalysis/setDihedral.h"
#include "../dropanalysis/setDihedralList.h"
#include "../dropanalysis/stericClashes.h"
#include "../dropanalysis/stericScan.h"
#include "../dropanalysis/vdwScan.h"

const char *commandList[][2] = {
	{ "measureDihedrals",
	 "Parses structure and provides log with structure and dihedral information."
	 },
	{ "setDihedral",
	 "Change a single user-specified dihedral angle for a given residue."
	 },
	{ "setDihedralList",
	 "Change several dihedrals from a user-provided input file with a list of dihedrals and angles."
	 },
	{ "stericClashes",
	 "Counts the number of atomic overlaps according to atomic radii used by Ramachandran."
	 },
	{ "stericScan",
	 "Calculates the average number of steric clashes in amino acid structures in Ramachandran Space."
	 },
	{ "vdwScan",
	 "Calculates the average Lennard-Jones energy in amino acid structures in Ramachandran Space." }
};

const int commandListLen = sizeof(commandList) / sizeof(commandList[0]);

void printCommandList()
{

	fprintf(stderr, "These are the commands available:\n");
	fprintf(stderr, "EXAMPLE > #: function - description -\n");

	for (int i = 0; i < commandListLen; i++) {
		fprintf(stderr, "%d: ", i + 1);
		for (int j = 0; j < 2; j++)	//2 accounts for the command and the description
		{
			fprintf(stderr, "%s - ", commandList[i][j]);
		}
		fprintf(stderr, "\n");
	}
}

bool findCommand(int argc, char **argv)
{
	bool found = true;

	//search available commands or functions
	if (strcmp(argv[1], commandList[0][0]) == 0) {
		measureDihedrals(argc, argv);
	} else if (strcmp(argv[1], commandList[1][0]) == 0) {
		setDihedral(argc, argv);
	} else if (strcmp(argv[1], commandList[2][0]) == 0) {
		setDihedralList(argc, argv);
	} else if (strcmp(argv[1], commandList[3][0]) == 0) {
		stericClashes(argc, argv);
	} else if (strcmp(argv[1], commandList[4][0]) == 0) {
		stericScan(argc, argv);
	} else if (strcmp(argv[1], commandList[5][0]) == 0) {
		vdwScan(argc, argv);
	} else {
		found = false;
	}

	return found;

}
