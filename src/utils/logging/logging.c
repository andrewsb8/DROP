#include <argp.h>
#include <stdio.h>
#include <stdlib.h>

#include "logging.h"

void logWelcome(FILE * log)
{
    fprintf(log, "%s\n", argp_program_version);
    //fprintf(log, "%s\n", doc);
    fprintf(log, "%s\n", argp_program_bug_address);
    fflush(log);
}

void logArgv(FILE * log, int argc, char **argv)
{
	//argp puts command at end of argv after parsing
	fprintf(log, "Command line: %s %s ", argv[0], argv[argc - 1]);
	for (int i = 1; i < argc - 1; i++) {
		fprintf(log, "%s ", argv[i]);
	}
	fprintf(log, "\n\n");
	fflush(log);
	return;
}

void drop_fatal(FILE * log, char *message)
{
	fprintf(log, "%s", message);
	fprintf(stderr, "%s", message);
	exit(1);
}

void drop_warning(FILE * log, char *message)
{
	fprintf(log, "%s", message);
	fprintf(stderr, "%s", message);
	fflush(log);
}

void drop_info(FILE * log, char *message)
{
	fprintf(log, "%s", message);
	fprintf(stderr, "%s", message);
	fflush(log);
}
