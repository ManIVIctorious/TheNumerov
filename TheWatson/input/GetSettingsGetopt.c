
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>

#include "data.h"
#include "settings.h"
#include "gitversion.h"
#include "ConvertString.h"

void usage(void);

void GetSettingsGetopt(settings* set, data* data, int argc, char** argv){

// optstring contains a list of all short option indices
//  if followed by one colon an argument is required,
//  if followed by two colons an optional argument can be set
    const char * optstring = "hi:m:M:o:a:Pt:X:Y:Z:V";

// create an array of option structs
/* struct description
 *{{{
 *  const char * option string
 *  int          has_arg (no_argument (0), required_argument (1), optional_argument (2))
 *  int        * flag    if !NULL : points to a variable set to the value given by val
 *  int          val     return value if flat == NULL
 *}}}*/
    const struct option longopts[] = {
        {"help",              no_argument, NULL, 'h'},
        {"version",           no_argument, NULL, 'V'},
        {"pipe-input",        no_argument, NULL, 'P'},
    // double values:
        {"threshold",   required_argument, NULL, 't'},
    // string values:
        {"input-file",  required_argument, NULL, 'i'},
        {"mode-file",   required_argument, NULL, 'm'},
        {"masses-file", required_argument, NULL, 'M'},
        {"output-file", required_argument, NULL, 'o'},
        {"append-file", required_argument, NULL, 'a'},
    // string array values:
        {"zeta-x",      required_argument, NULL, 'X'},
        {"zeta-y",      required_argument, NULL, 'Y'},
        {"zeta-z",      required_argument, NULL, 'Z'},
    // requires zero termination
        { NULL , 0 , NULL , 0 }
    };

// optind (provided by <getopt.h>)
//  should be initialized with 1 prior to every parsing of the command line
    optind = 1;
    int * longindex = NULL;
    while(optind < argc){

        int control = getopt_long(argc, argv, optstring, longopts, longindex);

    // getopt_long returns the integer "val" set in longopts array
    //  or -1 when end of options is reached
    //  => compare return values with "val" to determine the desired action
        if(control == -1){ break; }

        switch(control){

        // print version information
            case 'V':
                printf("\t%s, version: %s\n", program_invocation_short_name, gitversion);
                exit(EXIT_SUCCESS);

        // print help message
            case 'h':
                usage();
                exit(EXIT_SUCCESS);

        // set threshold
            case 't':
                set->threshold = convertstring_to_double(optarg, "threshold", NULL);
                break;

        // read input from pipe instead of input file
            case 'P':
                set->input_coordinates = "/dev/stdin";
                break;

        // set input file
            case 'i':
                set->input_coordinates = optarg;
                break;

        // set zeta_{x,y,z} colon separated string arrays
            case 'X':
                set->zeta_x = optarg;
                break;

            case 'Y':
                set->zeta_y = optarg;
                break;

            case 'Z':
                set->zeta_z = optarg;
                break;

        // add optarg to mode file list
        //  and increment dimension
            case 'm':
                data->dimension++;
                set->modelist = realloc( set->modelist, data->dimension * sizeof(char*) );
                set->modelist[data->dimension - 1] = optarg;
                break;

            case 'M':
                set->masses_file = optarg;
                set->masses_from_modes = 0;
                break;

        // set output file
            case 'o':
                set->output_file = optarg;
                set->output_fmode = "w";
                break;

        // set output file in append mode
            case 'a':
                set->output_file = optarg;
                set->output_fmode = "a";
                break;

        }
    }
}
