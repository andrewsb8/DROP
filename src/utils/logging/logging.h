#ifndef LOGGING_H_
#define LOGGING_H_

void drop_fatal(FILE * log, char *message);
void drop_warning(FILE * log, char *message);
void drop_info(FILE * log, char *message);

#endif
