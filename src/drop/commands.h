#ifndef COMMANDS_H_
#define COMMANDS_H_

extern const char *commandList[][2];	/* list of commands and descriptions */
extern const int commandListLen;	/* number of commands */

void printCommandList ();
bool findCommand (char *arg, int argc, char **argv);
void stripFArgv (int argc, char **argv);
void stripAllArgv (int argc, char **argv);
char *makeStringArgv (int argc, char **argv);

#endif
