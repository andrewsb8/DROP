#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>

#include "../logging/logging.h"
#include "../readProtein/readProtein.h"

//does the file exist? returns -1 if not
int fileExists(char *filename)
{
	return access(filename, F_OK);
}

void
processInput(struct protein *prot, char *input_file, FILE * log,
			 bool calc_bond_matrix, bool print_bond_matrix, int argc,
			 char **argv)
{
	logWelcome(log);
	//make a string of argv arguments for the log file output
	logArgv(log, argc, argv);
	if (fileExists(input_file) == -1) {
		drop_fatal(log, "ERROR: Input file does not exist. Exiting.\n");
	}

	fprintf(log, "Reading structure file: %s\n\n", input_file);

	readPDB(prot, input_file, log, calc_bond_matrix, print_bond_matrix);

	fprintf(log, "Done reading structure file: %s\n\n", input_file);
	fflush(log);
	return;
}

void writeFileLine(FILE * file, char *message)
{
	fprintf(file, "%s", message);
	fflush(file);
}
