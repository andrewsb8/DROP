#ifndef FILEHANDLING_H_
#define FILEHANDLING_H_

int fileExists(char *filename);
void inputInfo(struct protein *prot, char *input_file, FILE *log, bool print_bond_matrix, char *stringArgv);

#endif
