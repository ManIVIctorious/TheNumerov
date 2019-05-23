
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <errno.h>
#include "settings.h"

settings SetDefaultSettings();
settings GetSettingsGetopt(settings prefs, int argc, char** argv);

int InputNormalMode(char* inputfile, double* *mode, double* *mass);
int   InputFileList(char* inputfile, int dimension, double** *Q, char** *coordslist);
int    InputComFile(char* inputfile, double* x, double* y, double* z, int max_lines);

void MomentOfInertia(int dim, int n_atoms, double tot_mass, double* mass, double* x, double* y, double* z, double** zeta, double** Q, double* I, int state_index);
void InvertMatrix(gsl_matrix *Matrix, gsl_matrix *InvMatrix, int dimension);


int main(int argc, char **argv){

    int     i, j, k;
    int     control;
    int     n_atoms;
    double tot_mass;

    FILE * fdout;

    settings prefs = SetDefaultSettings();
             prefs = GetSettingsGetopt(prefs, argc, argv);


//--------------------------------------------------------------------------------------
//  Coriolis Coefficients Zeta  Coriolis Coefficients Zeta  Coriolis Coefficients Zeta
//--------------------------------------------------------------------------------------
//      Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input
//--------------------------------------------------------------------------------------
// input mode files:
//------------------------------
    if(prefs.dimension < 2){
        fprintf(stderr,
                "\n (-) Error: At least two modes have to be considered"
                "\n     Please check your input, aborting...\n\n"
            );
        exit(EXIT_FAILURE);
    }
// save data of each mode in modes array and mass in mass array
//  mode[D][entry_row * 3     ] = x_displacement
//  mode[D][entry_row * 3 + 1 ] = y_displacement
//  mode[D][entry_row * 3 + 2 ] = z_displacement
    double ** mode = malloc(prefs.dimension * sizeof(double*));
    if(mode == NULL){ perror("mode"); exit(errno); }
    for(i = 0; i < prefs.dimension; ++i){ mode[i] = NULL; }

    double ** mass = malloc(prefs.dimension * sizeof(double*));
    if(mass == NULL){ perror("mass"); exit(errno); }
    for(i = 0; i < prefs.dimension; ++i){ mass[i] = NULL; }

    n_atoms = InputNormalMode( prefs.modelist[0], &(mode[0]), &(mass[0]) );
    for(i = 1; i < prefs.dimension; ++i){
        control = InputNormalMode( prefs.modelist[i], &(mode[i]), &(mass[i]) );

    // check if all files contain the same number of atoms
        if( control != n_atoms ){
            fprintf(stderr,
                "\n (-) Error in reading mode files. Number of atoms in \"%s\" (\"%d\") does not"
                "\n     match number of atoms in \"%s\" (\"%d\"). Aborting...\n\n"
                , prefs.modelist[0], n_atoms, prefs.modelist[i], control
            );
            exit(EXIT_FAILURE);
        }
    }

// check for same masses in mode files (atom order/type)
    for(i = 1; i < prefs.dimension; ++i){
        for(j = 0; j < n_atoms; ++j){

            if( (mass[0][j] - mass[i][j])*(mass[0][j] - mass[i][j]) > prefs.threshold*prefs.threshold ){
                fprintf(stderr,
                        "\n (-) Error in reading mode files. Atom masses differ considerably"
                        "\n     between modes in entry number %d."
                        "\n             %s\t%s"
                        "\n     Value:  %le\t%le"
                        "\n     Aborting...\n\n"
                        , j, prefs.modelist[i-1], prefs.modelist[i], mass[i-1][j], mass[i][j]
                    );
                exit(EXIT_FAILURE);
            }
        }
    }

// if masses check succeeds free up memory:
//  so that only mass[0][entries] remains
    for(i = 1; i < prefs.dimension; ++i){
        free( mass[i] ); mass[i] = NULL;
    }

// calculate the system's total mass
    tot_mass = 0.0;
    for(i = 0; i < n_atoms; ++i){
        tot_mass += mass[0][i];
    }


// Norm mode files:
//------------------------------
    for(i = 0; i < prefs.dimension; ++i){
        double norm = 0.0;
        for(j = 0; j < n_atoms; ++j){
            norm += mode[i][j*3    ]*mode[i][j*3    ]; // x*x
            norm += mode[i][j*3 + 1]*mode[i][j*3 + 1]; // y*y
            norm += mode[i][j*3 + 2]*mode[i][j*3 + 2]; // z*z
        }

        norm = sqrt(norm);

        #ifdef debug_normalmode
            fprintf(stderr, "\t||%s|| = % lf\n", prefs.modelist[i], norm);
        #endif

        for(j = 0; j < n_atoms; ++j){
            mode[i][j*3    ] /= norm;
            mode[i][j*3 + 1] /= norm;
            mode[i][j*3 + 2] /= norm;
        }
    }


    #ifdef debug_normalmode
    // print mode files
        fprintf(stderr, "Number of Atoms: %d\n", n_atoms);
        for(i = 0; i < prefs.dimension; ++i){
            fprintf(stderr, "#%s:\n#", prefs.modelist[i]);
            fprintf(stderr, "\t    %s      ", "dx");
            fprintf(stderr, "\t    %s      ", "dy");
            fprintf(stderr, "\t    %s      ", "dz");
            fprintf(stderr, "\t    %s", "masses");
            fprintf(stderr, "\n");

            for(j = 0; j < n_atoms; ++j){
                for(k = 0; k < 3; ++k){
                    fprintf(stderr, "\t% .8le", mode[i][j*3 + k]);
                }
                fprintf(stderr, "\t% .8le\n", mass[0][j]);
            }
            fprintf(stderr, "\n");
        }
    #endif


//------------------------------------------------------------------------------
// Calculate Coriolis Coefficients Zeta    Calculate Coriolis Coefficients Zeta
//------------------------------------------------------------------------------
    double ** zeta = malloc( 3 * sizeof(double*) );
    if(zeta == NULL){ perror("zeta"); exit(errno); }
    for(i = 0; i < 3; ++i){
        zeta[i] = calloc(prefs.dimension * prefs.dimension, sizeof(double));
        if(zeta[i] == NULL){ perror("zeta[i]"); exit(errno); }
    }

// The Coriolis Coefficients are calculated by the following formulas:
//      zeta_x = dy1*dz2 - dz1*dy2
//      zeta_y = dz1*dx2 - dx1*dz2
//      zeta_z = dx1*dy2 - dy1*dx2
// where the numbers denote the respective modes
// Per definition therefore, the coefficients of a mode with itself is zero
//  and the coefficient of mode 1 with mode 2 is -(mode 2 with mode 1).
//  => Only the upper triangle has to be calculated!
    for(i = 0; i < prefs.dimension; ++i){
        for(j = i+1; j < prefs.dimension; ++j){

            for(k = 0; k < n_atoms; ++k){
                zeta[0][i*prefs.dimension + j] += mode[i][k*3 + 1]*mode[j][k*3 + 2] - mode[i][k*3 + 2]*mode[j][k*3 + 1];    // dy1*dz2 - dz1*dy2
                zeta[1][i*prefs.dimension + j] += mode[i][k*3 + 2]*mode[j][k*3    ] - mode[i][k*3    ]*mode[j][k*3 + 2];    // dz1*dx2 - dx1*dz2
                zeta[2][i*prefs.dimension + j] += mode[i][k*3    ]*mode[j][k*3 + 1] - mode[i][k*3 + 1]*mode[j][k*3    ];    // dx1*dy2 - dy1*dx2
            }

        }
    }

// free up memory
    for(i = 0; i < prefs.dimension; ++i){ free(mode[i]); mode[i] = NULL; }
    free(mode); mode = NULL;


// fill the lower triangle
    for(i = 0; i < prefs.dimension; ++i){
        for(j = i+1; j < prefs.dimension; ++j){

            zeta[0][j * prefs.dimension + i] = -zeta[0][i * prefs.dimension + j];
            zeta[1][j * prefs.dimension + i] = -zeta[1][i * prefs.dimension + j];
            zeta[2][j * prefs.dimension + i] = -zeta[2][i * prefs.dimension + j];

        }
    }

//--------------------------------------------------------------------------------------
//    Output  Output  Output  Output  Output  Output  Output  Output  Output  Output
//--------------------------------------------------------------------------------------
// output Coriolis coefficients
//--------------------------------------------------
    fdout = fopen(prefs.output_file, "w");

    fprintf(fdout, "#Modes:");
    for(i = 0; i < prefs.dimension; ++i){
        for(j = i+1; j < prefs.dimension; ++j){
            fprintf(fdout, "\t %s|%s", prefs.modelist[i], prefs.modelist[j]);
        }
    }
    fprintf(fdout, "\n");

    for(k = 0; k < 3; ++k){
        fprintf(fdout, "Zeta_%c:", "xyz"[k]);
        for(i = 0; i < prefs.dimension; ++i){
            for(j = i+1; j < prefs.dimension; ++j){
                fprintf(fdout, "\t% .12le", zeta[k][prefs.dimension * i + j]);
            }
        }
        fprintf(fdout, "\n");
    }

    fclose(fdout);


//--------------------------------------------------------------------------------------
//   Moment of inertia    Moment of inertia    Moment of inertia    Moment of inertia
//--------------------------------------------------------------------------------------
// get list of displacements Q and files to be read
//--------------------------------------------------
    if(!prefs.input_coordinates_set){ return 0; }

    char   ** coordslist = NULL;
    double ** Q = malloc(prefs.dimension * sizeof(double*));
    if( Q == NULL ){ perror("Q"); exit(errno); }
    for(i = 0; i < prefs.dimension; ++i){ Q[i] = NULL; }

    int coordslist_length = InputFileList(prefs.input_coordinates, prefs.dimension, &Q, &coordslist);


// print header
//--------------------------------------------------
    fdout = fopen(prefs.output_file, "a");

    fprintf(fdout, "#");
    for(i = 0; i < prefs.dimension; ++i){
        fprintf(fdout, "\t Q[%2d]              ", i);
    }
    for(i = 0; i < 3; ++i){
        for(j = i; j < 3; ++j){
            fprintf(fdout, "\t mu_%c%c             ", "xyz"[i], "xyz"[j]);
        }
    }
    fprintf(fdout, "\n");


// read coordinate files and start processing
//--------------------------------------------------
    for(i = 0; i < coordslist_length; ++i){
    // allocate memory for x, y and z coordinate
        double I[9];
        double * x = malloc(n_atoms * sizeof(double));
        double * y = malloc(n_atoms * sizeof(double));
        double * z = malloc(n_atoms * sizeof(double));
        if( x == NULL ){ perror("x"); exit(errno); }
        if( y == NULL ){ perror("y"); exit(errno); }
        if( z == NULL ){ perror("z"); exit(errno); }

        control = InputComFile(coordslist[i], x, y, z, n_atoms);

        if(control != n_atoms){
            fprintf(stderr,
                    "\n (-) Error in reading coordinate file \"%s\""
                    "\n     Number of coordinates (%d) does not match number of atoms (%d)"
                    "\n     Aborting...\n\n"
                    , coordslist[i], control, n_atoms
                );
            exit(EXIT_FAILURE);
        }

        #ifdef debug_coords
        // output coordinates for control
            fprintf(stderr, "\n%d\t%s\n", i+1, coordslist[i]);
            fprintf(stderr,
                    "\nNumber of Atoms:\t%d"
                    "\nInput coordinates:\n"
                    "\t         x         "
                    "\t         y         "
                    "\t         z         "
                    "\t    atomic mass    "
                    "\n"
                    , n_atoms
                );

            for(j = 0; j < n_atoms; ++j){
                fprintf(stderr, "%2d\t% .12le\t% .12le\t% .12le\t% .12le\n", j+1, x[j], y[j], z[j], mass[0][j]);
            }
            fprintf(stderr, "\n");
        #endif

    // calculate moment of inertia tensor and directly apply the Watson corrections
    //--------------------------------------------------------------------------------
        #ifdef debug_moment_of_inertia
            fprintf(stderr, "\n\n# Filename:\t%s\n", coordslist[i]);
        #endif
        MomentOfInertia(prefs.dimension, n_atoms, tot_mass, mass[0], x, y, z, zeta, Q, I, i);
        free(x); x = NULL;
        free(y); y = NULL;
        free(z); z = NULL;

    // calculate Effective Reciprocal Inertia Tensor mu
    //--------------------------------------------------------------------------------
        gsl_matrix * CorrMomentOfInertia = gsl_matrix_calloc(3, 3);
        gsl_matrix * mu =                  gsl_matrix_calloc(3, 3);

        for(j = 0; j < 3; ++j){
            for(k = 0; k < 3; ++k){
                gsl_matrix_set(CorrMomentOfInertia, j, k, I[j*3 + k]);
            }
        }

        InvertMatrix(CorrMomentOfInertia, mu, 3);
        gsl_matrix_free(CorrMomentOfInertia); CorrMomentOfInertia = NULL;

        #ifdef debug_mu
            fprintf(stderr, "\nEffective Reciprocal Inertia Tensor (mu)\n");
            for(j = 0; j < 3; ++j){
                for(k = 0; k < 3; ++k){
                    fprintf(stderr, "\t% .12le",gsl_matrix_get(mu, j, k));
                }
                fprintf(stderr, "\n");
            }
        #endif

    /* Output unique values of mu (upper triangle with main diagonal)
    //{{{

        Since I is symmetric (I = Iᵀ) and its correction is symmetric => I' = (I')ᵀ
        its inverse mu = (I')⁻¹ must be symmetric as well:

            A = Aᵀ      A * A⁻¹ = 𝟙     𝟙 = 𝟙ᵀ      (AB)ᵀ = BᵀAᵀ

            A * A⁻¹ = 𝟙 = 𝟙ᵀ = (A * A⁻¹)ᵀ = (A⁻¹)ᵀ * Aᵀ

                  A * A⁻¹ = (A⁻¹)ᵀ * Aᵀ       | A = Aᵀ
                  A * A⁻¹ = (A⁻¹)ᵀ * A        | A * A⁻¹ = 𝟙 = A⁻¹ * A
                  A⁻¹ * A = (A⁻¹)ᵀ * A        | * A⁻¹
            A⁻¹ * A * A⁻¹ = (A⁻¹)ᵀ * A * A⁻¹

                      A⁻¹ = (A⁻¹)ᵀ

        If A symmetric its inverse must be symmetric as well, q.e.d

    //}}}*/

    // output Q
        for(j = 0; j < prefs.dimension; ++j){
            fprintf(fdout, "\t% .12le", Q[j][i]);
        }
    // output mu:
        for(j = 0; j < 3; ++j){
            for(k = j; k < 3; ++k){
                fprintf(fdout, "\t% .12le", gsl_matrix_get(mu, j, k));
            }
        }
        fprintf(fdout, "\n");

        gsl_matrix_free(mu); mu = NULL;
    }

    fclose(fdout); fdout = NULL;

// free unused memory
    for(i = 0; i < coordslist_length; ++i){ free(coordslist[i]); coordslist[i] = NULL; }
    free(coordslist); coordslist = NULL;
    for(i = 0; i < prefs.dimension; ++i){ free(Q[i]); Q[i] = NULL; }
    free(Q); Q = NULL;
    for(i = 0; i < 3; ++i){ free(zeta[i]); zeta[i] = NULL; } free(zeta);
    zeta = NULL;
    free(mass[0]); mass[0] = NULL;
    free(mass);    mass    = NULL;

    return 0;
}
