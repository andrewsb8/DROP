#ifndef COMMANDS_H_
#define COMMANDS_H_

extern const char *commandList[][2];	/* list of commands and descriptions */
extern const int commandListLen;	/* number of commands */

void printCommandList();
bool findCommand(int argc, char **argv);

#endif
