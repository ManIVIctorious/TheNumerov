
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "typedefinitions.h"
#include "controlfile.h"

// Provided prototypes
void GetSettingsControlFile(char* inputfile, settings* preferences);


void GetSettingsControlFile(char* inputfile, settings* preferences){

    int i;
    size_t wordlength;
    char * optarg;

// create keyword list:
//  const char * keyword;
//  int          set;
//  const int    identifier;
//  char         value[_MaxSettingsStringLength_];
    struct keywords keywordlist[] = {
    // Boolian values:
        {"Analyze",                 0,  'a',  "" },
        {"Dipole",                  0,  'd',  "" },
        {"Check_Spacing",           0,  'T',  "" },
    // integer values:
        {"Dimensionality",          0,  'D',  "" },
        {"Stencil_Size",            0,  'n',  "" },
        {"Interpolation_points",    0,  's',  "" },
        {"N_Eigenstates",           0,  'N',  "" },
    // double values:
        {"Kin_E_Factor",            0,  'k',  "" },
        {"Pot_E_Factor",            0,  'v',  "" },
        {"Dipole_Factor",           0,  'f',  "" },
        {"IMOI_Factor",             0,  'M',  "" },
        {"Spacing_Threshold",       0,  't',  "" },
        {"Lower_Bound",             0,  'l',  "" },
        {"Upper_Bound",             0,  'u',  "" },
    // string values:
        {"Reduced_Masses",          0,  'm',  "" },
        {"Input_File",              0,  'i',  "" },
        {"External_Dipole_File",    0,  'e',  "" },
        {"Output_File",             0,  'o',  "" },
        {"Coriolis_File",           0,  'c',  "" },
    // other:
        {"Eigensolver",             0,  'E',  "" },
    // requires zero termination
        { NULL , 0 , 0 , "" }
    };

// get keyword values
    ControlFileParser(inputfile, keywordlist, 0);

// assign values to variables
    i = 0;
    while(keywordlist[i].keyword != NULL){

    // only iterate over the values set by the parsing function
        if(!keywordlist[i].set){
            ++i;
            continue;
        }
        optarg = keywordlist[i].value;

        switch(keywordlist[i].identifier){

        // Boolian values
            case 'a':
                preferences->analyze = atoi(optarg);
                if(preferences->analyze == 0){
                    if(strncasecmp("true", optarg, 4) == 0) { preferences->analyze = 1; }
                    else                                    { preferences->analyze = 0; }
                }
                break;

            case 'd':
                preferences->dipole = atoi(optarg);
                if(preferences->dipole == 0){
                    if(strncasecmp("true", optarg, 4) == 0) { preferences->dipole = 1; }
                    else                                    { preferences->dipole = 0; }
                }
                break;

            case 'T':
                preferences->check_spacing = atoi(optarg);
                if(preferences->check_spacing == 0){
                    if(strncasecmp("true", optarg, 4) == 0) { preferences->check_spacing = 1; }
                    else                                    { preferences->check_spacing = 0; }
                }
                break;


        // integer values
            case 'D':
                preferences->dimension   = atoi(optarg);
                break;

            case 'n':
                preferences->n_stencil   = atoi(optarg);
                break;

            case 's':
                preferences->n_spline    = atoi(optarg);
                break;

            case 'N':
                preferences->n_out       = atoi(optarg);
                break;


        // double values
            case 'k':
                preferences->ekin_factor = atof(optarg);
                break;

            case 'v':
                preferences->epot_factor = atof(optarg);
                break;

            case 'f':
                preferences->DipToAsm    = atof(optarg);
                break;

            case 'M':
                preferences->mu_factor   = atof(optarg);
                break;

            case 't':
                preferences->threshold   = atof(optarg);
                break;

            case 'l':
                preferences->e_min       = atof(optarg);
                break;

            case 'u':
                preferences->e_max       = atof(optarg);
                break;


        // string values
            case 'm':
            // memory allocation
                if( (wordlength = strlen(optarg)) < _MaxSettingsStringLength_){
                    preferences->masses_string = malloc( (wordlength+1) * sizeof(char) );
                    if(preferences->masses_string == NULL){ perror("preferences->masses_string"); exit(EXIT_FAILURE); }
                }else{
                    fprintf(stderr, "\n (-) Error: Keyword length for \"%s\" too long"
                                    "\n     Aborting...\n\n", keywordlist [i].keyword
                    );
                    exit(EXIT_FAILURE);
                }
            // copy optarg to string and ensure zero termination
                strncpy(preferences->masses_string, optarg, wordlength);
                preferences->masses_string[wordlength] = '\0';
                ++preferences->masses_string_set;
                break;

            case 'i':
            // memory allocation
                if( (wordlength = strlen(optarg)) < _MaxSettingsStringLength_){
                    preferences->input_file = malloc( (wordlength+1) * sizeof(char));
                    if(preferences->input_file == NULL){ perror("preferences->input_file"); exit(EXIT_FAILURE); }
                }else{
                    fprintf(stderr, "\n (-) Error: Keyword length for \"%s\" too long"
                                    "\n     Aborting...\n\n", keywordlist [i].keyword
                    );
                    exit(EXIT_FAILURE);
                }
            // copy optarg to string and ensure zero termination
                strncpy(preferences->input_file, optarg, wordlength);
                preferences->input_file[wordlength] = '\0';
                ++preferences->input_file_set;
                break;

            case 'e':
            // memory allocation
                if( (wordlength = strlen(optarg)) < _MaxSettingsStringLength_ ){
                    preferences->ext_dip_file = malloc( (wordlength+1) * sizeof(char) );
                    if(preferences->ext_dip_file == NULL){ perror("preferences->ext_dip_file"); exit(EXIT_FAILURE); }
                }else{
                    fprintf(stderr, "\n (-) Error: Keyword length for \"%s\" too long"
                                    "\n     Aborting...\n\n", keywordlist [i].keyword
                    );
                    exit(EXIT_FAILURE);
                }
            // copy optarg to string and ensure zero termination
                strncpy(preferences->ext_dip_file, optarg, wordlength);
                preferences->ext_dip_file[wordlength] = '\0';
                ++preferences->ext_dip_file_set;
                preferences->dipole = 1;
                break;

            case 'c':
                if( (wordlength = strlen(optarg)) < _MaxSettingsStringLength_){
                    preferences->coriolis_file = malloc( (wordlength+1) * sizeof(char));
                    if(preferences->coriolis_file == NULL){ perror("preferences->coriolis_file"); exit(EXIT_FAILURE); }
                }else{
                    fprintf(stderr, "\n (-) Error: Keyword length for \"%s\" too long"
                                    "\n     Aborting...\n\n", keywordlist [i].keyword
                    );
                    exit(EXIT_FAILURE);
                }
            // copy optarg to string and ensure zero termination
                strncpy(preferences->coriolis_file, optarg, wordlength);
                preferences->coriolis_file[wordlength] = '\0';
                ++preferences->coriolis_file_set;
                break;

            case 'o':
                if( (wordlength = strlen(optarg)) < _MaxSettingsStringLength_){
                    preferences->output_file = malloc( (wordlength+1) * sizeof(char));
                    if(preferences->output_file == NULL){ perror("preferences->output_file"); exit(EXIT_FAILURE); }
                }else{
                    fprintf(stderr, "\n (-) Error: Keyword length for \"%s\" too long"
                                    "\n     Aborting...\n\n", keywordlist [i].keyword
                    );
                    exit(EXIT_FAILURE);
                }
            // copy optarg to string and ensure zero termination
                strncpy(preferences->output_file, optarg, wordlength);
                preferences->output_file[wordlength] = '\0';
                ++preferences->output_file_set;
                break;


        // other
            case 'E':
                preferences->Eigensolver = atoi(optarg);
                if(preferences->Eigensolver == 0){

                    if(strncasecmp(optarg, "Intel_MKL_FEAST", strlen("Intel_MKL_FEAST")) == 0){
                        preferences->Eigensolver = 1;
                    }
                    else if(strncasecmp(optarg, "ARMADILLO_ARPACK", strlen("ARMADILLO_ARPACK")) == 0){
                        preferences->Eigensolver = 2;
                    }

                }
                break;
        }
        ++i;
    }

}
