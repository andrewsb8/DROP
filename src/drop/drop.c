#include <stdio.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "commands.h"

const char *usage_messg =
	"USAGE: drop command [OPTIONS]. Use drop -? or drop --help for more information.\n";
const char *program_bug_address =
	"https://github.com/andrewsb8/DROP/issues";
const char *program_version = "DROP Version 2024.1";
const char *doc =
	"DROP (Dihedral Rotation Of Proteins) -- A command line tool to investigate protein structures via direct manipulation of dihedral angles.";


int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "%s", usage_messg);
	} else if (strcmp(argv[1], "-?") == 0
			   || strcmp(argv[1], "--help") == 0) {
		fprintf(stderr, "%s\n", program_version);
		fprintf(stderr, "%s\n\n", doc);
		printCommandList();
		fprintf(stderr, "For more information: drop [command] -?\n");
		fprintf(stderr, "Report bugs to: %s\n", program_bug_address);
	} else if (strcmp(argv[1], "-v") == 0
			   || strcmp(argv[1], "--version") == 0) {
		fprintf(stderr, "%s\n", program_version);
	} else {
		if (!findCommand(argc, argv)) {
			fprintf(stderr, "Command %s not recognized.\n", argv[1]);
			fprintf(stderr, "%s", usage_messg);
		}
	}
	return 0;
}
