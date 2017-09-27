
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>

int InputComFile(char *inputfile, double **x, double **y, double **z);
int InputNormalMode(char *inputfile, int start, double **modedisplacement, double **mass);
int InvertMatrix(gsl_matrix *Matrix, gsl_matrix *InvMatrix, int dimension);

int CoriolisCoefficients(int n_atoms, double *mode1, double *mode2, double *zeta);

int main(int argc, char **argv){

//------------------------------------------------------------------------------------------------------------------
//   Deklaration   Deklaration   Deklaration   Deklaration   Deklaration   Deklaration   Deklaration   Deklaration
//------------------------------------------------------------------------------------------------------------------
    int m, n;       // m,n always in {x,y,z}
    int i, j, k;
    int control = 0;

    double entry = 0;

// input files
    char * comfile = NULL;
    //char modefiles[2][1024] = {"mode_35", "mode_36"};

// coordinates and modes input
    int n_atoms = 0;
    double * x  = NULL;         // freed
    double * y  = NULL;         // freed
    double * z  = NULL;         // freed
    double * modes     = NULL;  // freed
    double * masses    = NULL;  // freed
    double * deviation = NULL;  // freed

// auxiliar arrays
    double * auxmode1 = NULL;   // freed
    double * auxmode2 = NULL;   // freed


//------------------------------------------------------------------------------------------------------------------
//  Default values  Default values  Default values  Default values  Default values  Default values  Default values
//------------------------------------------------------------------------------------------------------------------

    int    verbose   = 1;       // set level of verbosity
    int    dimension = 2;       // number of included modes
    double threshold = 1E-10;   // threshold for number comparison






    double auxzeta[3];
    double zeta[dimension][dimension][3];

// center of mass coordinates and moment of inertia tensor
    double x0, y0, z0, tot_mass;
    double MomentOfInertia[9];
    gsl_matrix * CorrMomentOfInertia = NULL; // freed
    gsl_matrix * mu = NULL;     // freed


    deviation = malloc(dimension * sizeof(double));
    deviation[0] = -0.05162332776183;
    deviation[1] = -0.15486998328547;

// set verbose level
    FILE *fdverb = fopen("/dev/null", "w");
    if(verbose == 1){
        fclose(fdverb);
        fdverb = stderr;
    }

    char modefiles[2][1024] = {"mode_35", "mode_36"};
    comfile  = "2D_scan_resorcinol-A_b3lyp_ccl4_modes_35_36_dr=-0.05162332776183_-0.15486998328547.com";

//------------------------------------------------------------------------------------------------------------------
// Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input
//------------------------------------------------------------------------------------------------------------------

// Error code cheat sheet:
//      1: Wrong input (program handling)
//      2: Memory allocation problem
//      3: Problem with input function
//      4: File inconsistencies

// check input argument
    if(comfile == NULL){
        fprintf(stderr, "\n (-) Please specify valid coordinates and mode files");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }
    for(i = 0; i < dimension; ++i){
        if(modefiles[i] == NULL){
            fprintf(stderr, "\n (-) Please specify valid coordinates and mode files");
            fprintf(stderr, "\n     Aborting...\n\n");
            exit(1);
        }
    }

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
        control  = InputNormalMode(modefiles[i], i*n_atoms, &modes, &masses);

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
    fprintf(fdverb, "\nInput coordinates:", n_atoms);
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

// fill auxiliary mode arrays:
    auxmode1 = malloc(3*n_atoms*sizeof(double));
    auxmode2 = malloc(3*n_atoms*sizeof(double));

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

            auxzeta[0] = auxzeta[1] = auxzeta[2] = 0.0;
            CoriolisCoefficients(n_atoms, auxmode1, auxmode2, auxzeta);

            zeta[i][j][0] =  auxzeta[0];
            zeta[i][j][1] =  auxzeta[1];
            zeta[i][j][2] =  auxzeta[2];
            zeta[j][i][0] = -auxzeta[0];
            zeta[j][i][1] = -auxzeta[1];
            zeta[j][i][2] = -auxzeta[2];
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
            fprintf(fdverb, "\t%d\t%d\t% le\t% le\t% le\n", i, j, zeta[i][j][0], zeta[i][j][1], zeta[i][j][2]);
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
    fprintf(fdverb, "\t---------------------------------------------------------\n");
// m,n in {x,y,z}, i,j,k in mode_{0,...,n}
    for(m = 0; m < 3; ++m){
        for(n = 0; n < 3; ++n){

            for(i = 0; i < dimension; ++i){
                for(j = i+1; j < dimension; ++j){
                    for(k = 0; k < j; ++k){

                        entry = zeta[i][j][m] * zeta[k][j][n] * deviation[i] * deviation[k];
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

                        entry = zeta[i][j][m] * zeta[k][j][n] * deviation[i] * deviation[k];
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

                        entry = zeta[j][i][m] * zeta[k][i][n] * deviation[j] * deviation[k];
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

                        entry = zeta[j][i][m] * zeta[k][i][n] * deviation[j] * deviation[k];
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
    free(deviation); deviation = NULL;

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
    gsl_matrix_free(mu); mu = NULL;


    return 0;
}

int CoriolisCoefficients(int n_atoms, double *mode1, double *mode2, double *zeta){

    int i;

    for(i = 0; i < n_atoms; ++i){
        zeta[0] += mode1[i*3 + 1]*mode2[i*3 + 2] - mode1[i*3 + 2]*mode2[i*3 + 1]; // dy1*dz2 - dz1*dy2
        zeta[1] += mode1[i*3 + 2]*mode2[i*3    ] - mode1[i*3    ]*mode2[i*3 + 2]; // dz1*dx2 - dx1*dz2
        zeta[2] += mode1[i*3    ]*mode2[i*3 + 1] - mode1[i*3 + 1]*mode2[i*3    ]; // dx1*dy2 - dy1*dx2
    }

    return 0;
}
