#include <stdio.h>
#include <stdlib.h>

#include "logging.h"

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
