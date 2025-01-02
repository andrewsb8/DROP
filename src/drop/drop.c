#include <stdio.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "commands.h"
#include "../dropinfo/dropinfo.h"

int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "%s", usage_messg);
	} else if (strcmp(argv[1], "-?") == 0
			   || strcmp(argv[1], "--help") == 0) {
		fprintf(stderr, "%s - %s\n", program_name, program_version);
		fprintf(stderr, "%s\n\n", program_desc);
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
