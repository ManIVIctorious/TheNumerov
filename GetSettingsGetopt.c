
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "typedefinitions.h"

// Dependencies
int Help(char* filename, settings defaults);


// provided prototypes
settings GetSettingsGetopt(settings defaults, int argc, char** argv);

settings GetSettingsGetopt(settings defaults, int argc, char** argv){

    int control = 0;
    int * longindex = NULL;
    settings preferences = defaults;

//---------------------------------------------------------------------------------
//  FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS
//---------------------------------------------------------------------------------
    if(argc == 1){ exit(Help(argv[0], defaults)); }

    // optstring contains a list of all short option indices,
    //  indices followed by a colon are options requiring an argument.
    const char         * optstring = "hm:k:v:n:l:u:N:s:ac:i:dPo:t:TM:D:f:C:";
    const struct option longopts[] = {
    //  *name:      option name,
    //  has_arg:    if option requires argument,
    //  *flag:      if set to NULL getopt_long() returns val,
    //              else it returns 0 and flag points to a variable set to val
    //  val:        value to return
        {"help",                   no_argument, 0, 'h'},
        {"pipe",                   no_argument, 0, 'P'},
    // Boolian values:
        {"analyze",                no_argument, 0, 'a'},
        {"dipole",                 no_argument, 0, 'd'},
        {"no-spacing-check",       no_argument, 0, 'T'},
    // integer values:
        {"dimension",        required_argument, 0, 'D'},
        {"n-stencil",        required_argument, 0, 'n'},
        {"spline",           required_argument, 0, 's'},
        {"n-out",            required_argument, 0, 'N'},
    // double values:
        {"fkin",             required_argument, 0, 'k'},
        {"fpot",             required_argument, 0, 'v'},
        {"fdipole",          required_argument, 0, 'f'},
        {"fmu",              required_argument, 0, 'M'},
        {"dq-threshold",     required_argument, 0, 't'},
        {"lower-bound",      required_argument, 0, 'l'},
        {"upper-bound",      required_argument, 0, 'u'},
    // string values:
        {"masses",           required_argument, 0, 'm'},
        {"input-file",       required_argument, 0, 'i'},
        {"output-file",      required_argument, 0, 'o'},
        {"coriolis-input",   required_argument, 0, 'c'},
    // flags:
        {"mkl",                    no_argument, &preferences.Eigensolver, 1},
        {"armadillo",              no_argument, &preferences.Eigensolver, 2},
    // requires zero termination
        { 0 , 0 , 0 , 0 }
    };

    optind = 1; // option index starting by 1, provided by <getopt.h>
    while(optind < argc){

    // control is the integer representation of the corresponding option, e.g. 'x' = 120
    //  control = -1 corresponds to the end of the options
        control = getopt_long(argc, argv, optstring, longopts, longindex);

    // iterate over options control
        switch(control){

        // print help messages
            case 'h':
                Help(argv[0], defaults);
                exit(0);

        // Pipe: read input from stdin
            case 'P':
                strncpy(preferences.input_file, "/dev/stdin", _MaxSettingsStringLength_);
                break;


        // Boolian values
            case 'a':
                preferences.analyze = 1;
                break;

            case 'd':
                preferences.dipole = 1;
                break;

            case 'T':
                preferences.check_spacing = 0;
                break;


        // integer values
            case 'D':
                preferences.dimension   = atoi(optarg);
                break;

            case 'n':
                preferences.n_stencil   = atoi(optarg);
                break;

            case 's':
                preferences.n_spline    = atoi(optarg);
                break;

            case 'N':
                preferences.n_out       = atoi(optarg);
                break;


        // double values
            case 'k':
                preferences.ekin_factor = atof(optarg);
                break;

            case 'v':
                preferences.epot_factor = atof(optarg);
                break;

            case 'f':
                preferences.DipToAsm    = atof(optarg);
                break;

            case 'M':
                preferences.mu_factor   = atof(optarg);
                break;

            case 't':
                preferences.threshold   = atof(optarg);
                break;

            case 'l':
                preferences.e_min       = atof(optarg);
                break;

            case 'u':
                preferences.e_max       = atof(optarg);
                break;


        // string values
            case 'm':
            // copy optarg to string and ensure zero termination
                strncpy(preferences.masses_string, optarg, _MaxSettingsStringLength_);
                preferences.masses_string[_MaxSettingsStringLength_ - 1] = '\0';
                break;

            case 'i':
            // copy optarg to string and ensure zero termination
                strncpy(preferences.input_file, optarg, _MaxSettingsStringLength_);
                preferences.input_file[_MaxSettingsStringLength_ - 1] = '\0';
                break;

            case 'c':
            // copy optarg to string and ensure zero termination
                strncpy(preferences.coriolis_file, optarg, _MaxSettingsStringLength_);
                preferences.coriolis_file[_MaxSettingsStringLength_ - 1] = '\0';
                break;

            case 'o':
            // copy optarg to string and ensure zero termination
                strncpy(preferences.output_file, optarg, _MaxSettingsStringLength_);
                preferences.output_file[_MaxSettingsStringLength_ - 1] = '\0';
                break;


        }
    }

// return new settings struct "preferences"
    return preferences;
}
