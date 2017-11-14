
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <getopt.h>

int InputComFile(char *inputfile, double **x, double **y, double **z);
int InputMasses(char *inputfile, double **m);
int InputCoriolisCoefficients(char *inputfile, double ****zeta, int dimension);
int InvertMatrix(gsl_matrix *Matrix, gsl_matrix *InvMatrix, int dimension);
int Help(char *app_name);


int main(int argc, char **argv){
//------------------------------------------------------------------------------------------------------------------
//  Default values  Default values  Default values  Default values  Default values  Default values  Default values
//------------------------------------------------------------------------------------------------------------------
    int    verbose   = 0;       // set level of verbosity
    int    legend    = 0;       // whether to output header with column description
    int    dimension = 0;       // number of included modes
    double threshold = 1E-10;   // threshold for number comparison

// files
    char * operation = "w";     // whether to write or append to output-file
    char * verbout  = NULL;     // destination of verbosity output
    char * outfile  = NULL;     // standard output file
    char * comfile  = NULL;     // input comfile
    char * massfile = NULL;     // file containing the respective atom masses
    char * zetafile = NULL;     // file containing the Coriolis coefficients

// deviation from equilibrium position by i^th mode
    double  * deviation = malloc(sizeof(double));   // freed
    if(deviation == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for deviation");
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
    const char         * optstring = "halLvV:t:c:o:d:m:M:z:";
    const struct option longopts[] = {
    //  *name:      option name,
    //  has_arg:    if option requires argument,
    //  *flag:      if set to NULL getopt_long() returns val,
    //              else it returns 0 and flag points to a variable set to val
    //  val:        value to return
        {"help",                  no_argument, NULL, 'h'},
        {"legend",                no_argument, NULL, 'l'},
        {"legend-only",           no_argument, NULL, 'L'},
        {"append",                no_argument, NULL, 'a'},
        {"verbose",               no_argument, NULL, 'v'},
        {"verb-to-file",    required_argument, NULL, 'V'},
        {"threshold",       required_argument, NULL, 't'},
        {"coordinates",     required_argument, NULL, 'c'},
        {"zetafile",        required_argument, NULL, 'z'},
        {"masses",          required_argument, NULL, 'm'},
        {"deviation",       required_argument, NULL, 'd'},
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

            case 'l':
                legend = 1;
                break;

            case 'L':
                legend = 2;
                break;

            case 'V':
                ++verbose;
                verbout = optarg;
                break;

            case 'c':
                comfile = optarg;
                break;

            case 'm':
            case 'M':
                massfile = optarg;
                break;

            case 'z':
                zetafile = optarg;
                break;

            case 'd':
            // increment dimension by 1
                ++dimension;
            // set number of entries in deviation array to dimension
                deviation              = realloc(deviation, dimension*sizeof(double));
                deviation[dimension-1] = atof(optarg);
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
    if(comfile == NULL){
        fprintf(stderr, "\n (-) Please specify valid coordinates");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }
    if(massfile == NULL){
        fprintf(stderr, "\n (-) Please specify a valid masses file");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }


//------------------------------------------------------------------------------------------------------------------
//   Deklaration   Deklaration   Deklaration   Deklaration   Deklaration   Deklaration   Deklaration   Deklaration
//------------------------------------------------------------------------------------------------------------------
    int m, n;       // for all m,n: m,n in {x,y,z}
    int control  = 0;
    double entry = 0.0;

// coordinates and modes input
    int     n_atoms = 0;
    double   * x    = NULL;   // freed
    double   * y    = NULL;   // freed
    double   * z    = NULL;   // freed
    double   * mass = NULL;   // freed
    double *** zeta = NULL;

// center of mass coordinates and moment of inertia tensor
    double x0, y0, z0, tot_mass;
    double MomentOfInertia[9];

// correction of the moment of inertia
    gsl_matrix * CorrMomentOfInertia = NULL;    // freed
    gsl_matrix * mu = NULL;                     // freed


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
//  Output legend   Output legend   Output legend   Output legend   Output legend   Output legend   Output legend
    if(legend > 0){
        fprintf(fdout, "#");
        for(i = 0; i < dimension; ++i){
            fprintf(fdout, "\t deviation[%d]    ", i);
        }
        for(m = 0; m < 3; ++m){
            for(i = 0; i < dimension; ++i){
                for(j = i+1; j < dimension; ++j){
                    fprintf(fdout, "\t zeta^%c_%d%d       ", "xyz"[m], i, j);
                }
            }
        }
        for(m = 0; m < 3; ++m){
            for(n = m; n < 3; ++n){
                fprintf(fdout, "\t mu_%c%c           ", "xyz"[m], "xyz"[n]);
            }
        }
        fprintf(fdout, "\n");
    }

    if(legend == 2){ exit(0); }
//  Output legend   Output legend   Output legend   Output legend   Output legend   Output legend   Output legend
//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
// Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input
//------------------------------------------------------------------------------------------------------------------
// allocate memory for x, y, z and mass
    x = malloc(sizeof(double));
    y = malloc(sizeof(double));
    z = malloc(sizeof(double));
    mass = malloc(sizeof(double));
    if(x == NULL || y == NULL || z == NULL || mass == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for x, y, z or mass");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(2);
    }

// input com or coords file
    n_atoms = InputComFile(comfile, &x, &y, &z);
    if(n_atoms < 0){
        fprintf(stderr, "\n (-) Error in reading input file \"%s\"", comfile);
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(3);
    }

// input masses file
    control = InputMasses(massfile, &mass);
    if(control != n_atoms){
        fprintf(stderr, "\n (-) Error in input file \"%s\":", massfile);
        fprintf(stderr, "\n     Number of atoms (%d) doesn't match number of masses (%d)", n_atoms, control);
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(3);
    }

// input zeta coefficients
// allocate memory for a 3 dimensional double array of the form
//  dimension*dimension*3 (where dimension is the normal mode's number
//  and 3 refers to the {x,y,z} coordinates respectively)
    zeta = malloc(dimension * sizeof(double **));
    for(i = 0; i < dimension; ++i){
        zeta[i] = malloc(dimension * sizeof(double *));
        for(j = 0; j < dimension; ++j){
            zeta[i][j] = malloc(3 * sizeof(double));
        }
    }

    control = InputCoriolisCoefficients(zetafile, &zeta, dimension);
    if(control*3 != dimension*dimension*3){
        fprintf(stderr, "\n (-) Error in input file \"%s\"", zetafile);
        fprintf(stderr, "\n     Number of entries (%d) does not match expectations (%d)", control, dimension*dimension*3);
        fprintf(stderr, "\n     Please check your input (e.g. do dimensionalities match?)");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit (3);
    }

// output input for control
    fprintf(fdverb, "\nNumber of Atoms:\t%d", n_atoms);
    fprintf(fdverb, "\nInput coordinates:\n");
    fprintf(fdverb, "\t         x         ");
    fprintf(fdverb, "\t         y         ");
    fprintf(fdverb, "\t         z         ");
    fprintf(fdverb, "\t    atomic mass    ");
    fprintf(fdverb, "\n");
    for(i = 0; i < n_atoms; ++i){
        fprintf(fdverb, "\t% .12le\t% .12le\t% .12le\t% .12le\n", x[i], y[i], z[i], mass[i]);
    }

    fprintf(fdverb, "\nInput Coriolis coefficients\n");
    fprintf(fdverb, "\tModeA\tModeB");
    for(m = 0; m < 3; ++m){
        fprintf(fdverb, "\tzeta[ModeA][ModeB][%c]", "xyz"[m]);
    }
    fprintf(fdverb, "\n");
    for(i = 0; i < dimension; ++i){
        for(j = 0; j < dimension; ++j){
            fprintf(fdverb, "\t  %d\t  %d", i, j);
            for(m = 0; m < 3; ++m){
                fprintf(fdverb, "\t% .12le", zeta[i][j][m]);
            }
            fprintf(fdverb, "\n");
        }
    }


//------------------------------------------------------------------------------------------------------------------
// Moment of inertia  Moment of inertia  Moment of inertia  Moment of inertia  Moment of inertia  Moment of inertia
//------------------------------------------------------------------------------------------------------------------
// calculate center of mass
    x0 = y0 = z0 = tot_mass = 0.0;
    for(i = 0; i < n_atoms; ++i){
        x0 += x[i] * mass[i];
        y0 += y[i] * mass[i];
        z0 += z[i] * mass[i];

        tot_mass  += mass[i];
    }
    x0 /= tot_mass;
    y0 /= tot_mass;
    z0 /= tot_mass;

// output center of mass
    fprintf(fdverb, "\nCenter of mass\n");
    fprintf(fdverb, "\t% .12le", x0);
    fprintf(fdverb, "\t% .12le", y0);
    fprintf(fdverb, "\t% .12le", z0);
    fprintf(fdverb, "\n");

// translate molecule center of mass to origin
    for(i = 0; i < n_atoms; ++i){
        x[i] -= x0;
        y[i] -= y0;
        z[i] -= z0;
    }
    x0 = y0 = z0 = 0.0;

// calculate moment of inertia tensor
    for(i = 0; i < n_atoms; ++i){
    // main diagonal
        MomentOfInertia[0] += mass[i] * (y[i]*y[i] + z[i]*z[i]);  // 11 xx
        MomentOfInertia[4] += mass[i] * (x[i]*x[i] + z[i]*z[i]);  // 22 yy
        MomentOfInertia[8] += mass[i] * (x[i]*x[i] + y[i]*y[i]);  // 33 zz
    // upper triangle
        MomentOfInertia[1] -= mass[i] * x[i] * y[i];              // 12 xy
        MomentOfInertia[2] -= mass[i] * x[i] * z[i];              // 13 xz
        MomentOfInertia[5] -= mass[i] * y[i] * z[i];              // 23 yz
    }
    // lower triangle
    MomentOfInertia[3] = MomentOfInertia[1];    // 21 yx
    MomentOfInertia[6] = MomentOfInertia[2];    // 31 zx
    MomentOfInertia[7] = MomentOfInertia[5];    // 32 zy

    free(x);       x = NULL;
    free(y);       y = NULL;
    free(z);       z = NULL;
    free(mass); mass = NULL;

// output moment of inertia tensor
    fprintf(fdverb, "\nMoment of inertia tensor\n");
    for(m = 0; m < 3; ++m){
        for(n = 0; n < 3; ++n){
            fprintf(fdverb, "\t% .12le", MomentOfInertia[m*3 + n]);
        }
        fprintf(fdverb, "\n");
    }


//------------------------------------------------------------------------------------------------------------------
// Coriolis corrected moment of inertia  Coriolis corrected moment of inertia  Coriolis corrected moment of inertia
//------------------------------------------------------------------------------------------------------------------
    CorrMomentOfInertia = gsl_matrix_calloc(3, 3);
    if(CorrMomentOfInertia == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation of GSL_Matrix \"CorrMomentOfInertia\"");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(2);
    }

    fprintf(fdverb, "\nCoriolis correction of moment of inertia (a,b in {x,y,z}, i,j,k in mode_{0,...,n})\n");
    fprintf(fdverb, "\ta|b\ti|k\tj|k\tzeta^a_ik * zeta^b_jk * Q_i * Q_j\n");
    fprintf(fdverb, "\t-------------------------------------------------------------------------\n");

// m,n in {x,y,z}, i,j,k in mode_{0,...,n}
    for(m = 0; m < 3; ++m){
        for(n = 0; n < 3; ++n){

            for(i = 0; i < dimension; ++i){
                for(j = i+1; j < dimension; ++j){

                    for(k = 0; k < j; ++k){
                        entry = zeta[i][j][m] * zeta[k][j][n] * deviation[i] * deviation[k];
                        gsl_matrix_set(CorrMomentOfInertia, m, n, gsl_matrix_get(CorrMomentOfInertia, m, n) + entry);

                        fprintf(fdverb, "\t%c|%c\t%d|%d\t%d|%d", "xyz"[m], "xyz"[n], i, j, k, j);
                        fprintf(fdverb, "\tzeta^%c_%d,%d * zeta^%c_%d,%d * Q_%d * Q_%d =", "xyz"[m], i, j, "xyz"[n], k, j, i, k);
                        fprintf(fdverb, "\t% .12le\n", entry);
                    }
                    for(k = j+1; k < dimension; ++k){
                        entry = zeta[i][j][m] * zeta[k][j][n] * deviation[i] * deviation[k];
                        gsl_matrix_set(CorrMomentOfInertia, m, n, gsl_matrix_get(CorrMomentOfInertia, m, n) + entry);

                        fprintf(fdverb, "\t%c|%c\t%d|%d\t%d|%d", "xyz"[m], "xyz"[n], i, j, k, j);
                        fprintf(fdverb, "\tzeta^%c_%d,%d * zeta^%c_%d,%d * Q_%d * Q_%d =", "xyz"[m], i, j, "xyz"[n], k, j, i, k);
                        fprintf(fdverb, "\t% .12le\n", entry);
                    }

                // switch indices
                    for(k = 0; k < i; ++k){
                        entry = zeta[j][i][m] * zeta[k][i][n] * deviation[j] * deviation[k];
                        gsl_matrix_set(CorrMomentOfInertia, m, n, gsl_matrix_get(CorrMomentOfInertia, m, n) + entry);

                        fprintf(fdverb, "\t%c|%c\t%d|%d\t%d|%d", "xyz"[m], "xyz"[n], j, i, k, i);
                        fprintf(fdverb, "\tzeta^%c_%d,%d * zeta^%c_%d,%d * Q_%d * Q_%d =", "xyz"[m], j, i, "xyz"[n], k, i, j, k);
                        fprintf(fdverb, "\t% .12le\n", entry);
                    }
                    for(k = i+1; k < dimension; ++k){
                        entry = zeta[j][i][m] * zeta[k][i][n] * deviation[j] * deviation[k];
                        gsl_matrix_set(CorrMomentOfInertia, m, n, gsl_matrix_get(CorrMomentOfInertia, m, n) + entry);

                        fprintf(fdverb, "\t%c|%c\t%d|%d\t%d|%d", "xyz"[m], "xyz"[n], j, i, k, i);
                        fprintf(fdverb, "\tzeta^%c_%d,%d * zeta^%c_%d,%d * Q_%d * Q_%d =", "xyz"[m], j, i, "xyz"[n], k, i, j, k);
                        fprintf(fdverb, "\t% .12le\n", entry);
                    }

                }
            }

        }
    }

// output subtrahend for the correction of the moment of inertia
    fprintf(fdverb, "\nSubtrahend for correction of the moment of inertia\n");
    for(m = 0; m < 3; ++m){
        for(n = 0; n < 3; ++n){
            fprintf(fdverb, "\t% .12le", gsl_matrix_get(CorrMomentOfInertia, m, n));
        }
        fprintf(fdverb, "\n");
    }

// subtract subtrahend from MomentOfInertia and save in gsl_matrix CorrMomentOfInertia
    for(m = 0; m < 3; ++m){
        for(n = 0; n < 3; ++n){
            gsl_matrix_set(CorrMomentOfInertia, m, n, MomentOfInertia[m*3 + n] - gsl_matrix_get(CorrMomentOfInertia, m, n));
        }
    }

// output corrected Moment of Inertia
    fprintf(fdverb, "\nCorrected moment of inertia\n");
    for(m = 0; m < 3; ++m){
        for(n = 0; n < 3; ++n){
            fprintf(fdverb, "\t% .12le", gsl_matrix_get(CorrMomentOfInertia, m, n));
        }
        fprintf(fdverb, "\n");
    }


//------------------------------------------------------------------------------------------------------------------
//Output mode deviations Output mode deviations Output mode deviations Output mode deviations Output mode deviations
    for(i = 0; i < dimension; ++i){
        fprintf(fdout, "\t% 16.12le", deviation[i]);
    }
//Output mode deviations Output mode deviations Output mode deviations Output mode deviations Output mode deviations
//------------------------------------------------------------------------------------------------------------------
//  Output upper triangle of zeta (without main diagonal)    Output upper triangle of zeta (without main diagonal)
    for(m = 0; m < 3; ++m){
        for(i = 0; i < dimension; ++i){
            for(j = i+1; j < dimension; ++j){
                fprintf(fdout, "\t% 16.12le", zeta[i][j][m]);
            }
        }
    }
//  Output upper triangle of zeta (without main diagonal)    Output upper triangle of zeta (without main diagonal)
//------------------------------------------------------------------------------------------------------------------
    free(deviation); deviation = NULL;
    free(zeta); zeta = NULL;


//------------------------------------------------------------------------------------------------------------------
//  Effective Reciprocal Inertia Tensor   Effective Reciprocal Inertia Tensor   Effective Reciprocal Inertia Tensor
//------------------------------------------------------------------------------------------------------------------
// allocate memory for ERIT
    mu = gsl_matrix_calloc(3, 3);
    if(mu == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation of GSL_Matrix \"mu\"");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(2);
    }

// calculate ERIT by inverting corrected moment of inertia
    InvertMatrix(CorrMomentOfInertia, mu, 3);
    gsl_matrix_free(CorrMomentOfInertia); CorrMomentOfInertia = NULL;

// output inverse of corrected moment of inertia (mu)
    fprintf(fdverb, "\nEffective Reciprocal Inertia Tensor (mu)\n");
    for(m = 0; m < 3; ++m){
        for(n = 0; n < 3; ++n){
            fprintf(fdverb, "\t% .12le",gsl_matrix_get(mu, m, n));
        }
        fprintf(fdverb, "\n");
    }


//------------------------------------------------------------------------------------------------------------------
//Output upper triangle of mu  Output upper triangle of mu  Output upper triangle of mu  Output upper triangle of mu
    for(m = 0; m < 3; ++m){
        for(n = m; n < 3; ++n){
            fprintf(fdout, "\t% .12le", gsl_matrix_get(mu, m, n));
        }
    }
//Output upper triangle of mu  Output upper triangle of mu  Output upper triangle of mu  Output upper triangle of mu
//------------------------------------------------------------------------------------------------------------------
    gsl_matrix_free(mu); mu = NULL;


    fprintf(fdout,  "\n");
    fprintf(fdverb, "\n");
    fclose(fdout);  fdout  = NULL;
    fclose(fdverb); fdverb = NULL;
    return 0;
}

int Help(char *app_name){

    printf("\nAvailable options for %s:", app_name);

    printf("\n\nFlags not requiring arguments:");
    printf("\n\t-h|--help           Print this help dialogue");
    printf("\n\t-a|--append         Append to file instead of overwriting it");
    printf("\n\t-l|--legend         Precede output with a header describing each column");
    printf("\n\t-l|--legend-only    Like -l|--legend but quit after header output");
    printf("\n\t-v|--verbose        Increase verbosity of program (default to stderr)");

    printf("\n\nFlags which require an argument:");
    printf("\n\t-c|--coordinates    Name of file to get coordinates");

    printf("\n\t-z|--zetafile       File containing the Coriolis coefficients.");
    printf("\n\t-m|-M|--modefile    Masses file containing atomic masses.");
    printf("\n\t-d|--deviation      Actual deviation from coordinates by mode,");
    printf("\n\t                      can be called multiple times, at least twice.");

    printf("\n\t-o|--outputfile     Name of outputfile");
    printf("\n\t-t|--threshold      Threshold for number comparison (default 1E-10)");
    printf("\n\t-V|--verb-to-file   Write verbose output to file instead of stderr");

    printf("\n\n");

    return 0;
}
