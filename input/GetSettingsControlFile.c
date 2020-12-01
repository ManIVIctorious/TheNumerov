
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "ConvertString.h"
#include "controlfile.h"
#include "settings.h"

// Provided prototypes
void GetSettingsControlFile(char* inputfile, settings* set);
void free_set_string_values(void);

static char* set_string_values(char* source);


void GetSettingsControlFile(char* inputfile, settings* set){

// create keyword list:
/* struct description
 *{{{
 *  char      * keyword;
 *  int         set;
 *  const int   identifier;
 *  char     ** value;
 *}}}*/
    struct keywords keywordlist[] = {
    // Boolean values (i.e. char):
        {"Analyze",                 0,  'a', NULL },
        {"Dipole",                  0,  'd', NULL },
        {"Check_Spacing",           0,  'T', NULL },
    // integer values:
        {"Dimensionality",          0,  'D', NULL },
        {"Stencil_Size",            0,  'n', NULL },
        {"Interpolation_points",    0,  's', NULL },
        {"N_Eigenstates",           0,  'N', NULL },
    // double values:
        {"Kin_E_Factor",            0,  'k', NULL },
        {"Pot_E_Factor",            0,  'v', NULL },
        {"Dipole_Factor",           0,  'f', NULL },
        {"IMOI_Factor",             0,  'M', NULL },
        {"Spacing_Threshold",       0,  't', NULL },
        {"Lower_Bound",             0,  'l', NULL },
        {"Upper_Bound",             0,  'u', NULL },
    // string values:
        {"Reduced_Masses",          0,  'm', NULL },
        {"Input_File",              0,  'i', NULL },
        {"External_Dipole_File",    0,  'e', NULL },
        {"Output_File",             0,  'o', NULL },
        {"Coriolis_File",           0,  'c', NULL },
    // other:
        {"Eigensolver",             0,  'E', NULL },
    // NULL keyword as end condition
        { NULL , 0 , 0 , NULL }
    };

// get keyword values
    TokeniseControlFile(inputfile, keywordlist, 0);

// assign values to variables in settings struct
    for(int i = 0; keywordlist[i].keyword; ++i){

    // only iterate over the values set by the parsing function
        if(!keywordlist[i].set){ continue; }

    // print a warning if keyword was set multiple times
        if( keywordlist[i].set > 1 ){
            fprintf(stderr,
                " (-) Warning: Keyword \"%s\" was set multiple times, only last value will be considered\n"
                , keywordlist[i].keyword
            );
        }

    // create helper optarg (the string array)
        char ** optarg = keywordlist[i].value;
        char *  optarglast = optarg[keywordlist[i].set - 1];

    // switch over identifiers
        switch(keywordlist[i].identifier){

        // Boolean values
            case 'a':
                set->analyze = (char)convertstring_to_bool(optarglast, keywordlist[i].keyword, NULL);
                break;

            case 'd':
                set->dipole = (char)convertstring_to_bool(optarglast, keywordlist[i].keyword, NULL);
                break;

            case 'T':
                set->check_spacing = (char)convertstring_to_bool(optarglast, keywordlist[i].keyword, NULL);
                break;


        // integer values
            case 'D':
                set->dimension   = (int)convertstring_to_long(optarglast, keywordlist[i].keyword, NULL);
                break;

            case 'n':
                set->n_stencil   = (int)convertstring_to_long(optarglast, keywordlist[i].keyword, NULL);
                break;

            case 's':
                if( strcasecmp(optarglast, "none") == 0 ){ set->n_spline = 0; }
                else{
                    set->n_spline    = (int)convertstring_to_long(optarglast, keywordlist[i].keyword, NULL);
                }
                break;

            case 'N':
                set->n_out       = (int)convertstring_to_long(optarglast, keywordlist[i].keyword, NULL);
                break;


        // double values
            case 'k':
                set->ekin_factor = convertstring_to_double(optarglast, keywordlist[i].keyword, NULL);
                break;

            case 'v':
                set->epot_factor = convertstring_to_double(optarglast, keywordlist[i].keyword, NULL);
                break;

            case 'f':
                set->DipToAsm    = convertstring_to_double(optarglast, keywordlist[i].keyword, NULL);
                break;

            case 'M':
                set->mu_factor   = convertstring_to_double(optarglast, keywordlist[i].keyword, NULL);
                break;

            case 't':
                set->threshold   = convertstring_to_double(optarglast, keywordlist[i].keyword, NULL);
                break;

            case 'l':
                set->e_min       = convertstring_to_double(optarglast, keywordlist[i].keyword, NULL);
                break;

            case 'u':
                set->e_max       = convertstring_to_double(optarglast, keywordlist[i].keyword, NULL);
                break;


        // string values:
        // strings set here must be freed afterwards (free_set_string_values())
            case 'm':
                set->masses_string = set_string_values(optarglast);
                break;

            case 'i':
                set->input_file = set_string_values(optarglast);
                break;

            case 'e':
                set->ext_dip_file = set_string_values(optarglast);
                set->dipole = 1;
                break;

            case 'c':
                set->coriolis_file = set_string_values(optarglast);
                break;

            case 'o':
                set->output_file = set_string_values(optarglast);
                set->output_file_set++;
                break;


        // other
            case 'E':
                if     ( strcasecmp(optarglast, "Intel_MKL_FEAST")  == 0 ){ set->Eigensolver = 1; }
                else if( strcasecmp(optarglast, "ARMADILLO_ARPACK") == 0 ){ set->Eigensolver = 2; }
                else{ set->Eigensolver = (int)convertstring_to_long(optarglast, keywordlist[i].keyword, NULL); }
                break;
        }
    }

    free_keywordlistvalues(keywordlist);
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
