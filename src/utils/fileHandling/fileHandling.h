#ifndef FILEHANDLING_H_
#define FILEHANDLING_H_

int fileExists (char *filename);
void processInput (struct protein *prot, char *input_file, FILE * log,
				   bool calc_bond_matrix, bool print_bond_matrix,
				   int argc, char **argv);
char *makeStringArgv (int argc, char **argv);
void writeFileLine (FILE * file, char *message);

#endif
