
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "settings.h"
#include "ConvertString.h"

void usage(void);

void GetSettingsGetopt(settings* prefs, int argc, char** argv){

// optstring contains a list of all short option indices
//  if followed by one colon an argument is required,
//  if followed by two colons an optional argument can be set
    const char * optstring = "hi:m:M:o:a:Pt:";

// create an array of option structs
/* struct description
 *{{{
 *  const char * option string
 *  int          has_arg (no_argument (0), required_argument (1), optional_argument (2))
 *  int        * flag    if !NULL : points to a variable set to the value given by val
 *  int          val     return value if flat == NULL
 *}}}*/
    const struct option longopts[] = {
        {"help",              no_argument, 0, 'h'},
        {"pipe-input",        no_argument, 0, 'P'},
    // double values:
        {"threshold",   required_argument, 0, 't'},
    // string values:
        {"input-file",  required_argument, 0, 'i'},
        {"mode-file",   required_argument, 0, 'M'},
        {"output-file", required_argument, 0, 'o'},
        {"append-file", required_argument, 0, 'a'},
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

        // print help message
            case 'h':
                usage();
                exit(EXIT_SUCCESS);

        // set threshold
            case 't':
                prefs->threshold = convertstring_to_double(optarg, "threshold", NULL);
                break;

        // read input from pipe instead of input file
            case 'P':
                prefs->input_coordinates = "/dev/stdin";
                break;

        // set input file
            case 'i':
                prefs->input_coordinates = optarg;
                break;

        // add optarg to mode file list
        //  and increment dimension
            case 'M':
            case 'm':
                prefs->dimension++;
                prefs->modelist = realloc( prefs->modelist, prefs->dimension * sizeof(char*) );
                prefs->modelist[prefs->dimension - 1] = optarg;
                break;

        // set output file
            case 'o':
                prefs->output_file = optarg;
                prefs->output_fmode = "w";
                break;

        // set output file in append mode
            case 'a':
                prefs->output_file = optarg;
                prefs->output_fmode = "a";
                break;

        }
    }
}
