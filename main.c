
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <getopt.h>

int InputComFile(char *inputfile, double **x, double **y, double **z);
int InputNormalMode(char *inputfile, int start, double **modedisplacement, double **mass);
int InvertMatrix(gsl_matrix *Matrix, gsl_matrix *InvMatrix, int dimension);
int CoriolisCoefficients(int n_atoms, double *mode1, double *mode2, double *zeta_x, double *zeta_y, double *zeta_z);
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
    char  * operation = "w";                    // whether to write or append to output-file
    char  * verbout  = NULL;                    // destination of verbosity output
    char  * outfile  = NULL;                    // standard output file
    char  * comfile  = NULL;                    // input comfile
    char ** modelist = malloc(sizeof(char*));   // list of input modefiles
    if(modelist == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for modelist");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(2);
    }

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
    const char         * optstring = "halvV:t:c:m:o:d:";
    const struct option longopts[] = {
    //  *name:      option name,
    //  has_arg:    if option requires argument,
    //  *flag:      if set to NULL getopt_long() returns val,
    //              else it returns 0 and flag points to a variable set to val
    //  val:        value to return
        {"help",                  no_argument, NULL, 'h'},
        {"legend",                no_argument, NULL, 'l'},
        {"append",                no_argument, NULL, 'a'},
        {"verbose",               no_argument, NULL, 'v'},
        {"verb-to-file",    required_argument, NULL, 'V'},
        {"threshold",       required_argument, NULL, 't'},
        {"coordinates",     required_argument, NULL, 'c'},
        {"modefile",        required_argument, NULL, 'm'},
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

            case 'V':
                ++verbose;
                verbout = optarg;
                break;

            case 'c':
                comfile = optarg;
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

            case 'd':
            // increment k by 1
                ++k;
            // set number of entries is deviation array to k
                deviation = realloc(deviation, k*sizeof(double));
                deviation[k-1] = atof(optarg);
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

// check if there are as many entries in deviation as in modelist
    if(k != dimension){
        fprintf(stderr, "\n (-) Number of deviation entries and number of modes differ");
        fprintf(stderr, "\n     Please check your input");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }

// check input argument
    if(comfile == NULL){
        fprintf(stderr, "\n (-) Please specify valid coordinates and mode files");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }
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
    int m, n;       // for all m,n: m,n in {x,y,z}
    int control  = 0;
    double entry = 0.0;

// coordinates and modes input
    int n_atoms = 0;
    double * x  = NULL;      // freed
    double * y  = NULL;      // freed
    double * z  = NULL;      // freed
    double * modes  = NULL;  // freed
    double * masses = NULL;  // freed

// center of mass coordinates and moment of inertia tensor
    double x0, y0, z0, tot_mass;
    double MomentOfInertia[9];

// correction of the moment of inertia
    double     * zeta = NULL;                   // freed
    gsl_matrix * CorrMomentOfInertia = NULL;    // freed
    gsl_matrix * mu = NULL;                     // freed

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

// open outputfile
    if(outfile == NULL){
        fdout = stdout;
    }else{
        fdout = fopen(outfile, operation);
    }

//------------------------------------------------------------------------------------------------------------------
// Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input
//------------------------------------------------------------------------------------------------------------------
// input of com file
    x = malloc(sizeof(double));
    y = malloc(sizeof(double));
    z = malloc(sizeof(double));
    if(x == NULL || y == NULL || z == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for x, y or z");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(2);
    }
    n_atoms = InputComFile(comfile, &x, &y, &z);
    if(n_atoms < 0){
        fprintf(stderr, "\n (-) Error in reading input file \"%s\"", comfile);
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(3);
    }

// input of modes
    modes  = malloc(sizeof(double));
    masses = malloc(sizeof(double));
    if(modes == NULL || masses == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for mode1 or mass");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(2);
    }

    for(i = 0; i < dimension; ++i){
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
    fprintf(fdverb, "\nNumber of Atoms:\t%d", n_atoms);
    fprintf(fdverb, "\nInput coordinates:");
    for(i = 0; i < dimension; ++i){
        fprintf(fdverb, "\n   Mode No: %d", i);
        fprintf(fdverb, "\n   Mode deviation: % 20.14lf\n", deviation[i]);
        fprintf(fdverb, "\t x             ");
        fprintf(fdverb, "\t y             ");
        fprintf(fdverb, "\t z             ");
        fprintf(fdverb, "\t dx%d          ", i+1);
        fprintf(fdverb, "\t dy%d          ", i+1);
        fprintf(fdverb, "\t dz%d          ", i+1);
        fprintf(fdverb, "\t mass          ");
        fprintf(fdverb, "\n");
        for(j = 0; j < n_atoms; ++j){
            fprintf(fdverb, "\t% .8le\t% .8le\t% .8le", x[j], y[j], z[j]);
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
// Moment of inertia  Moment of inertia  Moment of inertia  Moment of inertia  Moment of inertia  Moment of inertia
//------------------------------------------------------------------------------------------------------------------

// calculate center of mass
    x0 = y0 = z0 = tot_mass = 0.0;
    for(i = 0; i < n_atoms; ++i){
        x0 += x[i] * masses[i];
        y0 += y[i] * masses[i];
        z0 += z[i] * masses[i];

        tot_mass  += masses[i];
    }
    x0 /= tot_mass;
    y0 /= tot_mass;
    z0 /= tot_mass;

// output center of mass
    fprintf(fdverb, "\nCenter of mass\n");
    fprintf(fdverb, "\t% 15.8lf", x0);
    fprintf(fdverb, "\t% 15.8lf", y0);
    fprintf(fdverb, "\t% 15.8lf", z0);
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
        MomentOfInertia[0] += masses[i] * (y[i]*y[i] + z[i]*z[i]);  // 11 xx
        MomentOfInertia[4] += masses[i] * (x[i]*x[i] + z[i]*z[i]);  // 22 yy
        MomentOfInertia[8] += masses[i] * (x[i]*x[i] + y[i]*y[i]);  // 33 zz
    // upper triangle
        MomentOfInertia[1] -= masses[i] * x[i] * y[i];              // 12 xy
        MomentOfInertia[2] -= masses[i] * x[i] * z[i];              // 13 xz
        MomentOfInertia[5] -= masses[i] * y[i] * z[i];              // 23 yz
    }
    // lower triangle
    MomentOfInertia[3] = MomentOfInertia[1];    // 21 yx
    MomentOfInertia[6] = MomentOfInertia[2];    // 31 zx
    MomentOfInertia[7] = MomentOfInertia[5];    // 32 zy

    free(x); x = NULL;
    free(y); y = NULL;
    free(z); z = NULL;

// output moment of inertia tensor
    fprintf(fdverb, "\nMoment of inertia tensor\n");
    for(m = 0; m < 3; ++m){
        for(n = 0; n < 3; ++n){
            fprintf(fdverb, "\t% 20.14lf", MomentOfInertia[m*3 + n]);
        }
        fprintf(fdverb, "\n");
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

            fprintf(fdverb, "\ndim = %d%d\n", i, j);
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

// output Coriolis coefficient tensor (dimension x dimension x 3)
    fprintf(fdverb, "\nCoriolis coefficients zeta^a_ij (a in {x,y,z}, i,j in mode_{0,...,n})\n");
    fprintf(fdverb, "\tModeA\tModeB\t\tx\t\ty\t\tz\n");
    for(i = 0; i < dimension; ++i){
        for(j = 0; j < dimension; ++j){
            fprintf(fdverb, "\t%d\t%d", i, j);
            for(m = 0; m < 3; ++m){
                fprintf(fdverb, "\t% le", zeta[(i*dimension + j)*3 + m]);
            }
            fprintf(fdverb, "\n");
        }
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
    fprintf(fdverb, "\tab\tik\tjk\tzeta^a_ik * zeta^b_jk * Q_i * Q_j\n");
    fprintf(fdverb, "\t-------------------------------------------------------------------------\n");
// m,n in {x,y,z}, i,j,k in mode_{0,...,n}
    for(m = 0; m < 3; ++m){
        for(n = 0; n < 3; ++n){

            for(i = 0; i < dimension; ++i){
                for(j = i+1; j < dimension; ++j){
                    for(k = 0; k < j; ++k){
                    // corresponds to zeta[i][j][m] * zeta[k][j][n] * deviation[i] * deviation[k]
                        entry = zeta[(i*dimension + j)*3 + m] * zeta[(k*dimension + j)*3 + n] * deviation[i] * deviation[k];
                        gsl_matrix_set(CorrMomentOfInertia, m, n, gsl_matrix_get(CorrMomentOfInertia, m, n) + entry);

                        fprintf(fdverb, "\t%d%d", m,n);
                        fprintf(fdverb, "\t%d%d", i,j);
                        fprintf(fdverb, "\t%d%d", k,j);
                        fprintf(fdverb, "\tzeta^%d_%d%d ", m, i, j);
                        fprintf(fdverb, "* zeta^%d_%d%d ", n, k, j);
                        fprintf(fdverb, "* Q_%d * Q_%d ", i, k);
                        fprintf(fdverb, "= % le", entry);
                        fprintf(fdverb, "\n");
                    }

                    for(k = j+1; k < dimension; ++k){
                    // corresponds to zeta[i][j][m] * zeta[k][j][n] * deviation[i] * deviation[k]
                        entry = zeta[(i*dimension + j)*3 + m] * zeta[(k*dimension + j)*3 + n] * deviation[i] * deviation[k];
                        gsl_matrix_set(CorrMomentOfInertia, m, n, gsl_matrix_get(CorrMomentOfInertia, m, n) + entry);

                        fprintf(fdverb, "\t%d%d", m,n);
                        fprintf(fdverb, "\t%d%d", i,j);
                        fprintf(fdverb, "\t%d%d", k,j);
                        fprintf(fdverb, "\tzeta^%d_%d%d ", m, i, j);
                        fprintf(fdverb, "* zeta^%d_%d%d ", n, k, j);
                        fprintf(fdverb, "* Q_%d * Q_%d ", i, k);
                        fprintf(fdverb, "= % le", entry);
                        fprintf(fdverb, "\n");
                    }

                    for(k = 0; k < i; ++k){
                    // corresponds to zeta[j][i][m] * zeta[k][i][n] * deviation[j] * deviation[k]
                        entry = zeta[(j*dimension + i)*3 + m] * zeta[(k*dimension + i)*3 + n] * deviation[j] * deviation[k];
                        gsl_matrix_set(CorrMomentOfInertia, m, n, gsl_matrix_get(CorrMomentOfInertia, m, n) + entry);

                        fprintf(fdverb, "\t%d%d", m,n);
                        fprintf(fdverb, "\t%d%d", j,i);
                        fprintf(fdverb, "\t%d%d", k,i);
                        fprintf(fdverb, "\tzeta^%d_%d%d ", m, j, i);
                        fprintf(fdverb, "* zeta^%d_%d%d ", n, k, i);
                        fprintf(fdverb, "* Q_%d * Q_%d ", j, k);
                        fprintf(fdverb, " = % le", entry);
                        fprintf(fdverb, "\n");
                    }
                    for(k = i+1; k < dimension; ++k){
                    // corresponds to zeta[j][i][m] * zeta[k][i][n] * deviation[j] * deviation[k]
                        entry = zeta[(j*dimension + i)*3 + m] * zeta[(k*dimension + i)*3 + n] * deviation[j] * deviation[k];
                        gsl_matrix_set(CorrMomentOfInertia, m, n, gsl_matrix_get(CorrMomentOfInertia, m, n) + entry);

                        fprintf(fdverb, "\t%d%d", m,n);
                        fprintf(fdverb, "\t%d%d", j,i);
                        fprintf(fdverb, "\t%d%d", k,i);
                        fprintf(fdverb, "\tzeta^%d_%d%d ", m, j, i);
                        fprintf(fdverb, "* zeta^%d_%d%d ", n, k, i);
                        fprintf(fdverb, "* Q_%d * Q_%d ", j, k);
                        fprintf(fdverb, "= % le", entry);
                        fprintf(fdverb, "\n");
                    }
                }
            }
        }
    }
//------------------------------------------------------------------------------------------------------------------
//  Output legend   Output legend   Output legend   Output legend   Output legend   Output legend   Output legend
    if(legend == 1){
        fprintf(fdout, "#");
        for(i = 0; i < dimension; ++i){
            fprintf(fdout, "\t deviation[%d]    ", i);
        }
    for(m = 0; m < 3; ++m){
        for(i = 0; i < dimension; ++i){
            for(j = i+1; j < dimension; ++j){
                fprintf(fdout, "\t zeta^%d_%d%d       ", m, i, j);
            }
        }
    }
    for(m = 0; m < 3; ++m){
        for(n = 0; n < 3; ++n){
            fprintf(fdout, "\t mu_%d%d           ", m, n);
        }
    }
    fprintf(fdout, "\n");
    }
//  Output legend   Output legend   Output legend   Output legend   Output legend   Output legend   Output legend
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
                fprintf(fdout, "\t% 16.12le", zeta[(i*dimension + j)*3 + m]);
            }
        }
    }
//  Output upper triangle of zeta (without main diagonal)    Output upper triangle of zeta (without main diagonal)
//------------------------------------------------------------------------------------------------------------------
    free(deviation); deviation = NULL;
    free(zeta); zeta = NULL;

// output subtrahend for the correction of the moment of inertia
    fprintf(fdverb, "\nSubtrahend for correction of the moment of inertia\n");
    for(m = 0; m < 3; ++m){
        for(n = 0; n < 3; ++n){
            fprintf(fdverb, "\t% 20.14lf",gsl_matrix_get(CorrMomentOfInertia, m, n));
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
            fprintf(fdverb, "\t% 20.14lf",gsl_matrix_get(CorrMomentOfInertia, m, n));
        }
        fprintf(fdverb, "\n");
    }


    mu = gsl_matrix_calloc(3, 3);
    if(mu == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation of GSL_Matrix \"mu\"");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(2);
    }
    InvertMatrix(CorrMomentOfInertia, mu, 3);
    gsl_matrix_free(CorrMomentOfInertia); CorrMomentOfInertia = NULL;

// output inverse of corrected moment of inertia (mu)
    fprintf(fdverb, "\nInverse of corrected moment of inertia (mu)\n");
    for(m = 0; m < 3; ++m){
        for(n = 0; n < 3; ++n){
            fprintf(fdverb, "\t% 20.14lf",gsl_matrix_get(mu, m, n));
        }
        fprintf(fdverb, "\n");
    }
//------------------------------------------------------------------------------------------------------------------
//Output upper triangle of mu  Output upper triangle of mu  Output upper triangle of mu  Output upper triangle of mu
    for(m = 0; m < 3; ++m){
        for(n = 0; n < 3; ++n){
            fprintf(fdout, "\t% .12le", gsl_matrix_get(mu, m, n));
        }
    }
//Output upper triangle of mu  Output upper triangle of mu  Output upper triangle of mu  Output upper triangle of mu
//------------------------------------------------------------------------------------------------------------------
    gsl_matrix_free(mu); mu = NULL;


    fprintf(fdout, "\n");
    fclose(fdout);  fdout  = NULL;
    fclose(fdverb); fdverb = NULL;
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
    printf("\n\t-l|--legend         Precede output with a header describing each column");
    printf("\n\t-v|--verbose        Increase verbosity of program (default to stderr)");

    printf("\n\nFlags which require an argument:");
    printf("\n\t-c|--coordinates    Name of file to get coordinates");

    printf("\n\t-m|--modefile       Name of file to get mode displacement coordinates,");
    printf("\n\t                      can be called multiple times (at least twice)");
    printf("\n\t-d|--deviation      Actual deviation from coordinates by mode, can be");
    printf("\n\t                      called multiple times, same number as coordinates");

    printf("\n\t-o|--outputfile     Name of outputfile");
    printf("\n\t-t|--threshold      Threshold for number comparison (default 1E-10)");
    printf("\n\t-V|--verb-to-file   Write verbose output to file instead of stderr");

    printf("\n\n");


    return 0;
}
