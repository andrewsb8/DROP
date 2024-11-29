#ifndef LOGGING_H_
#define LOGGING_H_

void logWelcome(FILE * log);
void logArgv(FILE * log, int argc, char **argv);
void drop_fatal(FILE * log, char *message);
void drop_warning(FILE * log, char *message);
void drop_info(FILE * log, char *message);

#endif
