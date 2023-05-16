#ifndef COMMANDS_H_
#define COMMANDS_H_

extern const char *commandList[][2];
extern const int commandListLen;

void printCommandList();
bool findCommand(char *arg, int argc, char **argv);
char ** stripFArgv(int argc, char **argv);

#endif
