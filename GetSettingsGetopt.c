
#include <stdlib.h>
#include <getopt.h>

#include "settings.h"

void Help();

settings GetSettingsGetopt(settings prefs, int argc, char** argv){

// optstring contains a list of all short option indices
//  if followed by one colon an argument is required,
//  if followed by two colons an optional argument can be set
    const char * optstring = "hi:m:M:o:Pt:";
    const struct option longopts[] = {
    //  *name:      option name,
    //  has_arg:    if option requires argument,
    //  *flag:      if set to NULL getopt_long() returns val,
    //              else it returns 0 and flag points to a variable set to val
    //  val:        value to return
        {"help",              no_argument, 0, 'h'},
        {"threshold",   required_argument, 0, 't'},
        {"pipe-input",        no_argument, 0, 'P'},
        {"input-file",  required_argument, 0, 'i'},
        {"mode-file",   required_argument, 0, 'M'},
        {"output-file", required_argument, 0, 'o'},
    // requires zero termination
        { 0 , 0 , 0 , 0 }
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
                Help();
                exit(EXIT_SUCCESS);

        // set threshold
            case 't':
                prefs.threshold = atof(optarg);
                break;

        // read input from pipe instead of input file
            case 'P':
                prefs.input_coordinates = "/dev/stdin";
                prefs.input_coordinates_set++;
                break;

        // set input file
            case 'i':
                prefs.input_coordinates = optarg;
                prefs.input_coordinates_set++;
                break;

        // add optarg to mode file list
        //  and increment dimension
            case 'M':
            case 'm':
                prefs.dimension++;
                prefs.modelist = realloc( prefs.modelist, prefs.dimension * sizeof(char*) );
                prefs.modelist[prefs.dimension - 1] = optarg;
                break;

        // set output file
            case 'o':
                prefs.output_file = optarg;
                prefs.output_file_set++;
                break;

        }
    }



    return prefs;
}
