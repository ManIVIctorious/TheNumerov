
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

int InputNormalMode(char *inputfile, int start, double **modedisplacement, double **mass);
int CoriolisCoefficients(int n_atoms, double *mode1, double *mode2, double *zeta_x, double *zeta_y, double *zeta_z);
int Help(char *app_name);


int main(int argc, char **argv){
//------------------------------------------------------------------------------------------------------------------
//  Default values  Default values  Default values  Default values  Default values  Default values  Default values
//------------------------------------------------------------------------------------------------------------------
    int    verbose   = 0;       // set level of verbosity
    int    dimension = 0;       // number of included modes
    double threshold = 1E-10;   // threshold for number comparison

// files
    char  * operation = "w";                    // whether to write or append to output-file
    char  * verbout  = NULL;                    // destination of verbosity output
    char  * outfile  = NULL;                    // standard output file
    char ** modelist = malloc(sizeof(char*));   // list of input modefiles
    if(modelist == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for modelist");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(2);
    }

    int i = 0;
    int j = 0;
    int k = 0;
//------------------------------------------------------------------------------------------------------------------
//  FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS
//------------------------------------------------------------------------------------------------------------------
    if(argc == 1){ exit(Help(argv[0])); }
    // optstring contains a list of all short option indices,
    //  indices followed by a colon are options requiring an argument.
    const char         * optstring = "havV:t:m:o:";
    const struct option longopts[] = {
    //  *name:      option name,
    //  has_arg:    if option requires argument,
    //  *flag:      if set to NULL getopt_long() returns val,
    //              else it returns 0 and flag points to a variable set to val
    //  val:        value to return
        {"help",                  no_argument, NULL, 'h'},
        {"append",                no_argument, NULL, 'a'},
        {"verbose",               no_argument, NULL, 'v'},
        {"verb-to-file",    required_argument, NULL, 'V'},
        {"threshold",       required_argument, NULL, 't'},
        {"modefile",        required_argument, NULL, 'm'},
        {"outputfile",      required_argument, NULL, 'o'},
    };

    optind = 1; // option index starting by 1, provided by <getopt.h>
    while(optind < argc){

    // i is the integer representation of the corresponding option, e.g. x = 120
    //  i = -1 corresponds to the end of the options
        i = getopt_long(argc, argv, optstring, longopts, &j);

    // iterate over options i
        switch(i){
            case 'h':
                exit(Help(argv[0]));
                break;

            case 'v':
                ++verbose;
                break;

            case 'V':
                ++verbose;
                verbout = optarg;
                break;

            case 'm':
            // increment dimension by 1
                ++dimension;
            // set number of entries in modelist to dimension
                modelist              = realloc(modelist, dimension * sizeof(char*));
            // allocate memory for the new entry
                modelist[dimension-1] =  malloc(strlen(optarg) * sizeof(char));
            // copy the content of optarg into the newly allocated char array
                strncpy(modelist[dimension-1], optarg, strlen(optarg));
                break;

            case 'o':
                outfile = optarg;
                break;

            case 'a':
                operation = "a";
                break;

            case 't':
                threshold = atof(optarg);
                break;


            default:
                exit(Help(argv[0]));
        }
    }

//------------------------------------------------------------------------------------------------------------------
//  Error handling  Error handling  Error handling  Error handling  Error handling  Error handling  Error handling
//------------------------------------------------------------------------------------------------------------------

// Error code cheat sheet:
//      1: Wrong input (program handling)
//      2: Memory allocation problem
//      3: Problem with input function
//      4: File inconsistencies

// check if dimension is at least 2
    if(dimension < 2){
        fprintf(stderr, "\n (-) At least two modes have to be considered");
        fprintf(stderr, "\n     Please check your input");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }

// check input argument
    for(i = 0; i < dimension; ++i){
        if(modelist[i] == NULL){
            fprintf(stderr, "\n (-) Please specify valid coordinates and mode files");
            fprintf(stderr, "\n     Aborting...\n\n");
            exit(1);
        }
    }


//------------------------------------------------------------------------------------------------------------------
//   Deklaration   Deklaration   Deklaration   Deklaration   Deklaration   Deklaration   Deklaration   Deklaration
//------------------------------------------------------------------------------------------------------------------
    int m;       // for all m: m in {x,y,z}
    int control  = 0;

// coordinates and modes input
    int n_atoms = 0;
    double * modes  = NULL;  // freed
    double * masses = NULL;  // freed

// correction of the moment of inertia
    double     * zeta = NULL;                   // freed

// auxiliar arrays
    double * auxmode1 = NULL;   // freed
    double * auxmode2 = NULL;   // freed


    FILE * fdverb = NULL;   // verbosity output file handler
    FILE * fdout  = NULL;   // output file handler

// set level of verbosity
    if(verbose == 0){
        fdverb = fopen("/dev/null", "w");
    }
    else if (verbose > 0){
        if(verbout == NULL){
            fdverb = stderr;
        }else{
            fdverb = fopen(verbout, "w");
        }
    }

// open output-file
    if(outfile == NULL){
        fdout = stdout;
    }else{
        fdout = fopen(outfile, operation);
    }


//------------------------------------------------------------------------------------------------------------------
// Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input
//------------------------------------------------------------------------------------------------------------------
// input of modes
    modes  = malloc(sizeof(double));
    masses = malloc(sizeof(double));
    if(modes == NULL || masses == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for mode1 or mass");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(2);
    }

    n_atoms = InputNormalMode(modelist[0], 0, &modes, &masses);
    for(i = 1; i < dimension; ++i){
        control  = InputNormalMode(modelist[i], i*n_atoms, &modes, &masses);

    // check if all files contain the same number of atoms
        if(control/(i+1) != n_atoms){
            fprintf(stderr, "\n (-) Error in reading input files.");
            fprintf(stderr, "\n     Number of atoms doesn't match.");
            fprintf(stderr, "\n     Aborting...\n\n");
            exit(4);
        }
    }

// check for same masses in mode files (atom order/type)
    for(i = 0; i < dimension; ++i){
        for(j = 0; j < n_atoms; ++j){
            if( labs(masses[j] - masses[j + i*n_atoms]) > threshold ){
                fprintf(stderr, "\n (-) Error in reading input files.");
                fprintf(stderr, "\n     Atom masses differ between modes.\n");
                fprintf(stderr, "\n     Aborting...\n\n");
                exit(4);
            }
        }
    }


// output input for control
    fprintf(fdverb, "\nNumber of Atoms: %4d", n_atoms);
    for(i = 0; i < dimension; ++i){
        fprintf(fdverb, "\n    Mode Number: %4d\n", i);
        fprintf(fdverb, "\t  dx%d         ", i);
        fprintf(fdverb, "\t  dy%d         ", i);
        fprintf(fdverb, "\t  dz%d         ", i);
        fprintf(fdverb, "\t  masses       ");
        fprintf(fdverb, "\n");
        for(j = 0; j < n_atoms; ++j){
            fprintf(fdverb, "\t% .8le\t% .8le\t% .8le",
                            modes[j*3     + i*n_atoms*3],
                            modes[j*3 + 1 + i*n_atoms*3],
                            modes[j*3 + 2 + i*n_atoms*3]
                   );
            fprintf(fdverb, "\t% .8le", masses[j]);
            fprintf(fdverb, "\n");
        }
    }


//------------------------------------------------------------------------------------------------------------------
// Coriolis coefficients  Coriolis coefficients  Coriolis coefficients  Coriolis coefficients  Coriolis coefficients
//------------------------------------------------------------------------------------------------------------------

    zeta = calloc(dimension * dimension * 3, sizeof(double));
    if(zeta == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for zeta");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(2);
    }

// fill auxiliary mode arrays:
    auxmode1 = calloc(3*n_atoms, sizeof(double));
    auxmode2 = calloc(3*n_atoms, sizeof(double));

// calculate upper triangle of the zeta tensors (upper triangle of dim x dim)
//  [i][j] values are -[j][i] values and main diagonal remains zero
    for(i = 0; i < dimension; ++i){
        for(j = i+1; j < dimension; ++j){
            for(k = 0; k < n_atoms*3; ++k){
                auxmode1[k] = modes[i*n_atoms*3 + k];
                auxmode2[k] = modes[j*n_atoms*3 + k];
            }

            fprintf(fdverb, "\ndim = %d|%d\n", i, j);
            for(k = 0; k < n_atoms*3; ++k){
                fprintf(fdverb, "\t% le", auxmode1[k]);
                if(k%3 == 2) fprintf(fdverb, "\n");
            }
            fprintf(fdverb, "\n");
            for(k = 0; k < n_atoms*3; ++k){
                fprintf(fdverb, "\t% le", auxmode2[k]);
                if(k%3 == 2) fprintf(fdverb, "\n");
            }

            CoriolisCoefficients( n_atoms,
                                  auxmode1,
                                  auxmode2,
                                  &zeta[(i*dimension + j)*3 + 0],
                                  &zeta[(i*dimension + j)*3 + 1],
                                  &zeta[(i*dimension + j)*3 + 2]
                            );

            zeta[(j*dimension + i)*3 + 0] = -zeta[(i*dimension + j)*3 + 0];
            zeta[(j*dimension + i)*3 + 1] = -zeta[(i*dimension + j)*3 + 1];
            zeta[(j*dimension + i)*3 + 2] = -zeta[(i*dimension + j)*3 + 2];
        }
    }
    free(modes); modes = NULL;
    free(auxmode1); auxmode1 = NULL;
    free(auxmode2); auxmode2 = NULL;

    fprintf(fdverb, "\n");
// output Coriolis coefficient tensor (dimension x dimension x 3)
    fprintf(fdout, "#  Coriolis coefficients zeta^a_ij (a in {x,y,z}, i,j in mode_{0,...,n})\n");
    fprintf(fdout, "#\tModeA\tModeB\t         x         \t         y         \t         z\n");
    for(i = 0; i < dimension; ++i){
        for(j = 0; j < dimension; ++j){
            fprintf(fdout, "\t  %d\t  %d", i, j);
            for(m = 0; m < 3; ++m){
                fprintf(fdout, "\t% .12le", zeta[(i*dimension + j)*3 + m]);
            }
            fprintf(fdout, "\n");
        }
    }


//------------------------------------------------------------------------------------------------------------------
//      Close file descriptors    Close file descriptors    Close file descriptors    Close file descriptors
//------------------------------------------------------------------------------------------------------------------
    fclose(fdverb); fdverb = NULL;
    fclose(fdout);  fdout  = NULL;
    return 0;
}

int CoriolisCoefficients(int n_atoms, double *mode1, double *mode2, double *zeta_x, double *zeta_y, double *zeta_z){

    int i;

    for(i = 0; i < n_atoms; ++i){
        *zeta_x += mode1[i*3 + 1]*mode2[i*3 + 2] - mode1[i*3 + 2]*mode2[i*3 + 1]; // dy1*dz2 - dz1*dy2
        *zeta_y += mode1[i*3 + 2]*mode2[i*3    ] - mode1[i*3    ]*mode2[i*3 + 2]; // dz1*dx2 - dx1*dz2
        *zeta_z += mode1[i*3    ]*mode2[i*3 + 1] - mode1[i*3 + 1]*mode2[i*3    ]; // dx1*dy2 - dy1*dx2
    }

    return 0;
}

int Help(char *app_name){

    printf("\nAvailable options for %s:", app_name);

    printf("\n\nFlags not requiring arguments:");
    printf("\n\t-h|--help           Print this help dialogue");
    printf("\n\t-a|--append         Append to file instead of overwriting it");
    printf("\n\t-v|--verbose        Increase verbosity of program (default to stderr)");

    printf("\n\nFlags which require an argument:");

    printf("\n\t-m|--modefile       Name of file to get mode displacement coordinates,");
    printf("\n\t                      can be called multiple times (at least twice)");

    printf("\n\t-o|--outputfile     Name of outputfile");
    printf("\n\t-t|--threshold      Threshold for number comparison (default 1E-10)");
    printf("\n\t-V|--verb-to-file   Write verbose output to file instead of stderr");

    printf("\n\n");

    return 0;
}
