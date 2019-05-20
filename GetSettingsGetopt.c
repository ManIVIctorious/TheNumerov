
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "typedefinitions.h"

// Dependencies
int Help();


// provided prototypes
settings GetSettingsGetopt(settings preferences, int argc, char** argv);

settings GetSettingsGetopt(settings preferences, int argc, char** argv){

    int control = 0;
    int * longindex = NULL;

//---------------------------------------------------------------------------------
//  FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS
//---------------------------------------------------------------------------------
    // optstring contains a list of all short option indices,
    //  indices followed by a colon are options requiring an argument.
    const char         * optstring = "a::c:d::e:f:hi:k:l:m:n:o:s:t:u:v:C:D:M:N:PT::";
    const struct option longopts[] = {
    //  *name:      option name,
    //  has_arg:    if option requires argument,
    //  *flag:      if set to NULL getopt_long() returns val,
    //              else it returns 0 and flag points to a variable set to val
    //  val:        value to return
        {"help",                   no_argument, 0, 'h'},
        {"pipe",                   no_argument, 0, 'P'},
    // Boolian values:
        {"analyze",          optional_argument, 0, 'a'},
        {"dipole",           optional_argument, 0, 'd'},
        {"no-spacing-check", optional_argument, 0, 'T'},
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
        {"ext-dipole-file",  required_argument, 0, 'e'},
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
                Help();
                exit(EXIT_SUCCESS);

        // Pipe: read input from stdin
            case 'P':
                preferences.input_file = "/dev/stdin";
                ++preferences.input_file_set;
                break;


        // Boolian values
            case 'a':
                if(optarg == NULL){ preferences.analyze = 1; }
                else{
                    preferences.analyze     = atoi(optarg);
                    if(preferences.analyze == 0){
                        if(strncasecmp("true", optarg, 4) == 0) { preferences.analyze = 1; }
                        else                                    { preferences.analyze = 0; }
                    }
                }
                break;

            case 'd':
                if(optarg == NULL){ preferences.dipole = 1; }
                else{
                    preferences.dipole      = atoi(optarg);
                    if(preferences.dipole == 0){
                        if(strncasecmp("true", optarg, 4) == 0) { preferences.dipole = 1; }
                        else                                    { preferences.dipole = 0; }
                    }
                }
                break;

            case 'T':
                if(optarg == NULL){ preferences.check_spacing = 0; }
                else{
                    preferences.check_spacing = atoi(optarg);
                    if(preferences.check_spacing == 0){
                        if(strncasecmp("true", optarg, 4) == 0) { preferences.check_spacing = 0; }
                        else                                    { preferences.check_spacing = 1; }
                    }
                }
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
                if(preferences.masses_string_set){ free(preferences.masses_string); }
                preferences.masses_string = optarg;
                ++preferences.masses_string_set;
                break;

            case 'i':
            // copy optarg to string and ensure zero termination
                if(preferences.input_file_set){ free(preferences.input_file); }
                preferences.input_file = optarg;
                ++preferences.input_file_set;
                break;

            case 'e':
            // copy optarg to string and ensure zero termination
                if(preferences.ext_dip_file){ free(preferences.ext_dip_file); }
                preferences.ext_dip_file = optarg;
                ++preferences.ext_dip_file_set;
                preferences.dipole = 1;
                break;

            case 'c':
            // copy optarg to string and ensure zero termination
                if(preferences.coriolis_file_set){ free(preferences.coriolis_file); }
                preferences.coriolis_file = optarg;
                ++preferences.coriolis_file_set;
                break;

            case 'o':
            // copy optarg to string and ensure zero termination
                if(preferences.output_file_set){ free(preferences.output_file); }
                preferences.output_file = optarg;
                ++preferences.output_file_set;
                break;


        }
    }

// return new settings struct "preferences"
    return preferences;
}
