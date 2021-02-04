
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>

#include "ConvertString.h"
#include "settings.h"
#include "gitversion.h"

// Dependencies
void usage(void);


// provided prototypes
void GetSettingsGetopt(int argc, char** argv, settings* preferences);

void GetSettingsGetopt(int argc, char** argv, settings* preferences){

/* the following flags are already handled in the calling function
 * and therefore can't be used here:
 *
 * -h, --help           print the usage() function and exit
 * -C, --control-file   read control file to populate settings
 */

// optstring contains a list of all short option indices,
//  indices followed by a colon are options requiring an argument.
    const char * optstring = "a::c:d::e:f:i:k:l:m:n:o:p::s:t:u:v:w:C:D:F::M:N:PT::V";

// create an array of option structs
/* struct description
 *{{{
 *  const char * option string
 *  int          has_arg (no_argument (0), required_argument (1), optional_argument (2))
 *  int        * flag    if !NULL : points to a variable set to the value given by val
 *  int          val     return value if flat == NULL
 *}}}*/
    const struct option longopts[] = {
        {"version",                no_argument, NULL, 'V'},
        {"pipe",                   no_argument, NULL, 'P'},
    // Boolean values:
        {"frequencies",      optional_argument, NULL, 'F'},
        {"analyse",          optional_argument, NULL, 'a'},
        {"dipole",           optional_argument, NULL, 'd'},
        {"no-spacing-check", optional_argument, NULL, 'T'},
        {"periodic",         optional_argument, NULL, 'p'},
    // integer values:
        {"dimension",        required_argument, NULL, 'D'},
        {"n-stencil",        required_argument, NULL, 'n'},
        {"spline",           required_argument, NULL, 's'},
        {"n-out",            required_argument, NULL, 'N'},
    // double values:
        {"fkin",             required_argument, NULL, 'k'},
        {"fpot",             required_argument, NULL, 'v'},
        {"fdipole",          required_argument, NULL, 'f'},
        {"fmu",              required_argument, NULL, 'M'},
        {"dq-threshold",     required_argument, NULL, 't'},
        {"watson-threshold", required_argument, NULL, 'w'},
        {"lower-bound",      required_argument, NULL, 'l'},
        {"upper-bound",      required_argument, NULL, 'u'},
    // string values:
        {"masses",           required_argument, NULL, 'm'},
        {"input-file",       required_argument, NULL, 'i'},
        {"ext-dipole-file",  required_argument, NULL, 'e'},
        {"output-file",      required_argument, NULL, 'o'},
        {"coriolis-input",   required_argument, NULL, 'c'},
    // flags:
      // Eigensolver bit-mask
      // 1  Intel MKL FEAST
      // 2  Armadillo ARPACK
      // 4  Not implemented yet
        {"mkl",                    no_argument, (int*)&preferences->Eigensolver, 1},
        {"armadillo",              no_argument, (int*)&preferences->Eigensolver, 2},
    // requires zero termination
        { NULL , 0 , NULL , 0 }
    };


    optind = 1; // option index starting by 1, provided by <getopt.h>
    while(optind < argc){

    // control is the integer representation of the corresponding option, e.g. 'x' = 120
    //  control = -1 corresponds to the end of the options
        int control = getopt_long(argc, argv, optstring, longopts, NULL);

    // compare control with "val" from longopts to determine the desired action
        if(control == -1){ break; }
        switch( control ){

        // print version information
            case 'V':
                printf("\t%s, version: %s\n", program_invocation_short_name, gitversion);
                exit(EXIT_SUCCESS);

        // Pipe: read input from stdin
            case 'P':
                preferences->input_file = "/dev/stdin";
                break;


        // Boolean values
            case 'F':
                if(optarg == NULL){ preferences->frequencies = 1; }
                else{
                    preferences->frequencies = (char)convertstring_to_bool(optarg, "Frequencies", NULL);
                }
                break;

            case 'a':
                if(optarg == NULL){ preferences->analyse = 1; }
                else{
                    preferences->analyse = (char)convertstring_to_bool(optarg, "Analyse", NULL);
                }
                break;

            case 'd':
                if(optarg == NULL){ preferences->dipole = 1; }
                else{
                    preferences->dipole = (char)convertstring_to_bool(optarg, "Dipole", NULL);
                }
                break;

            case 'p':
                if(optarg == NULL){ preferences->periodic = 1; }
                else{
                    preferences->periodic = (char)convertstring_to_bool(optarg, "Periodic", NULL);
                }
                break;

            case 'T':
                if(optarg == NULL){ preferences->check_spacing = 0; }
                else{
                    preferences->dipole = (char)convertstring_to_bool(optarg, "Spacing check", NULL);
                }
                break;


        // integer values
            case 'D':
                preferences->dimension = (int)convertstring_to_long(optarg, "Dimension", NULL);
                break;

            case 'n':
                preferences->n_stencil = (int)convertstring_to_long(optarg, "Stencil size", NULL);
                break;

            case 's':
                preferences->n_spline  = (int)convertstring_to_long(optarg, "Number of spline points", NULL);
                break;

            case 'N':
                preferences->n_out     = (int)convertstring_to_long(optarg, "Number of Eigenstates", NULL);
                break;


        // double values
            case 'k':
                preferences->kJpermol_to_oue = convertstring_to_double(optarg, "ekin factor", NULL);
                break;

            case 'v':
                preferences->epot_to_oue = convertstring_to_double(optarg, "epot factor", NULL);
                break;

            case 'f':
                preferences->dip_to_Asm  = convertstring_to_double(optarg, "convert dip to Asm", NULL);
                break;

            case 'M':
                preferences->InvInertia_to_molpergAasq = convertstring_to_double(optarg, "ERIT to kJ/mol", NULL);
                break;

            case 't':
                preferences->threshold   = convertstring_to_double(optarg, "Threshold", NULL);
                break;

            case 'w':
                preferences->InvInertiaThreshold = convertstring_to_double(optarg, "Watson Threshold", NULL);
                break;

            case 'l':
                preferences->e_min       = convertstring_to_double(optarg, "Minimum eigenvalue", NULL);
                break;

            case 'u':
                preferences->e_max       = convertstring_to_double(optarg, "Maximum eigenvalue", NULL);
                break;


        // string values
            case 'm':
            // point to optarg
                preferences->masses_string = optarg;
                break;

            case 'i':
            // point to optarg
                preferences->input_file = optarg;
                break;

            case 'e':
            // point to optarg
                preferences->ext_dip_file = optarg;
                preferences->dipole = 1;
                break;

            case 'c':
            // point to optarg
                preferences->coriolis_file = optarg;
                break;

            case 'o':
            // point to optarg
                preferences->output_file = optarg;
                break;

        }
    }
}
