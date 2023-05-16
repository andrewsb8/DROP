#include <stdio.h>
#include <argp.h>
#include <string.h>
#include <stdbool.h>

#include "commands.h"

const char *argp_program_bug_address = "https://github.com/andrewsb8/DROP/issues";
const char *argp_program_version = "DROP Version 0.0.1";
static char doc[] = "DROP (Dihedral Rotation Of Proteins) -- A tool to investigate protein structures.";

int num_master_options_read = 0;

static int parse_opt(int key, char *arg, struct argp_state *state)
{
  //only allow one "parent" option to be executed at a time
  if (num_master_options_read == 0)
  {
    switch(key)
    {
        // function definition case
        case 'f':
        {
          //NULL case is handled by requiring -f in struct argp_option
          if (!findCommand(arg, state->argc, state->argv))
          {
            fprintf(stderr, "This is not an option for -f. See argp -c for details.\n");
            argp_usage(state);
          }
          num_master_options_read++;
          break;
        }

        //print list of available commands
        case 'c':
        {
          if ( (arg == NULL) )
          {
            printCommandList();
          }
          else //block currently does not execute
          {
            fprintf(stderr, "This option does not accept an argument.\n");
            argp_usage(state);
          }
          num_master_options_read++;
          break;
        }

        case ARGP_KEY_END:
        {
          // Not enough arguments.
          if (state->arg_num < 2)
          {
            fprintf(stderr, "No arguments provided.\n");
            argp_usage(state);
          }
    			break;
        }
    }
  }
  return 0;
}

int main(int argc, char *argv[])
{
  static struct argp_option options[] =
  {
    { "func", 'f', "[Function String]", 0, "Identify Function" },
    { "commands", 'c', 0,  OPTION_ARG_OPTIONAL, "See options for -f" },
    { 0, 0, 0, 0, "To see options for a specific function: ./drop -f [FUNC] --help"},
    { 0, 0, 0, 0, "Informational Options:" },
    { 0 }
  };

  struct argp argp = { options, parse_opt, 0, 0 };

  return argp_parse(&argp, argc, argv, 0, 0, 0);
}
