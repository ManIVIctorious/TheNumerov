
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "ConvertString.h"
#include "controlfile.h"
#include "settings.h"

// Provided prototypes
void GetSettingsControlFile(char* inputfile, settings* set);
static char* set_string_values(char* source);
void free_set_string_values(void);


void GetSettingsControlFile(char* inputfile, settings* set){

// create keyword list:
/* struct description
 *{{{
 *  char      * keyword;
 *  int         set;
 *  const int   identifier;
 *  char        value[_MaxSettingsStringLength_];
 *}}}*/
    struct keywords keywordlist[] = {
    // Boolean values (i.e. char):
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
    // NULL keyword as end condition
        { NULL , 0 , 0 , "" }
    };

// get keyword values
    ControlFileParser(inputfile, keywordlist, 0);

// assign values to variables in settings struct
    for(int i = 0; keywordlist[i].keyword; ++i){

    // only iterate over the values set by the parsing function
        if(!keywordlist[i].set){ continue; }

    // create helper optarg (the string array)
        char * optarg = keywordlist[i].value;

    // switch over identifiers
        switch(keywordlist[i].identifier){

        // Boolean values
            case 'a':
                set->analyze = (char)convertstring_to_bool(optarg, keywordlist[i].keyword, NULL);
                break;

            case 'd':
                set->dipole = (char)convertstring_to_bool(optarg, keywordlist[i].keyword, NULL);
                break;

            case 'T':
                set->check_spacing = (char)convertstring_to_bool(optarg, keywordlist[i].keyword, NULL);
                break;


        // integer values
            case 'D':
                set->dimension   = (int)convertstring_to_long(optarg, keywordlist[i].keyword, NULL);
                break;

            case 'n':
                set->n_stencil   = (int)convertstring_to_long(optarg, keywordlist[i].keyword, NULL);
                break;

            case 's':
                if( strcasecmp(optarg, "none") == 0 ){ set->n_spline = 0; }
                else{
                    set->n_spline    = (int)convertstring_to_long(optarg, keywordlist[i].keyword, NULL);
                }
                break;

            case 'N':
                set->n_out       = (int)convertstring_to_long(optarg, keywordlist[i].keyword, NULL);
                break;


        // double values
            case 'k':
                set->ekin_factor = convertstring_to_double(optarg, keywordlist[i].keyword, NULL);
                break;

            case 'v':
                set->epot_factor = convertstring_to_double(optarg, keywordlist[i].keyword, NULL);
                break;

            case 'f':
                set->DipToAsm    = convertstring_to_double(optarg, keywordlist[i].keyword, NULL);
                break;

            case 'M':
                set->mu_factor   = convertstring_to_double(optarg, keywordlist[i].keyword, NULL);
                break;

            case 't':
                set->threshold   = convertstring_to_double(optarg, keywordlist[i].keyword, NULL);
                break;

            case 'l':
                set->e_min       = convertstring_to_double(optarg, keywordlist[i].keyword, NULL);
                break;

            case 'u':
                set->e_max       = convertstring_to_double(optarg, keywordlist[i].keyword, NULL);
                break;


        // string values:
        // strings set here must be freed afterwards (free_set_string_values())
            case 'm':
                set->masses_string = set_string_values(optarg);
                set->masses_string_set++;
                break;

            case 'i':
                set->input_file = set_string_values(optarg);
                break;

            case 'e':
                set->ext_dip_file = set_string_values(optarg);
                set->dipole = 1;
                break;

            case 'c':
                set->coriolis_file = set_string_values(optarg);
                break;

            case 'o':
                set->output_file = set_string_values(optarg);
                set->output_file_set++;
                break;


        // other
            case 'E':
                if     ( strcasecmp(optarg, "Intel_MKL_FEAST")  == 0 ){ set->Eigensolver = 1; }
                else if( strcasecmp(optarg, "ARMADILLO_ARPACK") == 0 ){ set->Eigensolver = 2; }
                else{ set->Eigensolver = (int)convertstring_to_long(optarg, keywordlist[i].keyword, NULL); }
                break;
        }
    }
}


static int     set_string_values_counter      = 0;
static char ** set_string_values_pointer_list = NULL;

static char* set_string_values(char* source){

    ++set_string_values_counter;
    char ** pointer = set_string_values_pointer_list;

// reallocate memory for pointer list
    pointer = realloc(pointer, set_string_values_counter * sizeof(char*));
    if(pointer == NULL){ perror("set_string_values"); exit(errno); }

// get size of source and allocate memory
    int control = strlen(source);
    pointer[set_string_values_counter-1] = malloc( (control+1) * sizeof(char) );
    if(pointer[set_string_values_counter-1] == NULL){ perror(source), exit(errno); }

// copy source to target string and ensure zero termination
    strcpy(pointer[set_string_values_counter-1], source);
    pointer[set_string_values_counter-1][control] = '\0';

// update set_string_values_pointer and return
    set_string_values_pointer_list = pointer;
    return pointer[set_string_values_counter-1];
}


void free_set_string_values(void){

    for(int i = 0; i < set_string_values_counter; ++i){
        free(set_string_values_pointer_list[i]);
        set_string_values_pointer_list[i] = NULL;
    }
    free(set_string_values_pointer_list);
    set_string_values_pointer_list = NULL;
}
