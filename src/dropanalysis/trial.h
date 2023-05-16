#ifndef TRIAL_H_
#define TRIAL_H_

extern struct argp readpdb_argp;

static int trial_parse(int key, char *arg, struct argp_state *state);
void trial(int argc, char **argv);

#endif
