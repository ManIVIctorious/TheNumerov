
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <errno.h>
#include "settings.h"

// dependencies
settings SetDefaultSettings();
settings GetSettingsGetopt(settings prefs, int argc, char** argv);

int InputNormalMode(char* inputfile, double* *mode, double* *mass);
int ProcessFileList(settings prefs);


int main(int argc, char **argv){

    int     i, j, k;
    int     control;

    settings prefs = SetDefaultSettings();
             prefs = GetSettingsGetopt(prefs, argc, argv);


//--------------------------------------------------------------------------------------
//  Coriolis Coefficients Zeta  Coriolis Coefficients Zeta  Coriolis Coefficients Zeta
//--------------------------------------------------------------------------------------
// input mode files:
//--------------------------------------------------
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

    prefs.n_atoms = InputNormalMode( prefs.modelist[0], &(mode[0]), &(mass[0]) );
    for(i = 1; i < prefs.dimension; ++i){
        control = InputNormalMode( prefs.modelist[i], &(mode[i]), &(mass[i]) );

    // check if all files contain the same number of atoms
        if( control != prefs.n_atoms ){
            fprintf(stderr,
                "\n (-) Error in reading mode files. Number of atoms in \"%s\" (\"%d\") does not"
                "\n     match number of atoms in \"%s\" (\"%d\"). Aborting...\n\n"
                , prefs.modelist[0], prefs.n_atoms, prefs.modelist[i], control
            );
            exit(EXIT_FAILURE);
        }
    }

// check for same masses in mode files (atom order/type)
    for(i = 1; i < prefs.dimension; ++i){
        for(j = 0; j < prefs.n_atoms; ++j){

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
// point prefs.atom_masses to the position of masses[0]
    prefs.atom_masses = mass[0];

// calculate the system's total mass
    prefs.tot_mass = 0.0;
    for(i = 0; i < prefs.n_atoms; ++i){
        prefs.tot_mass += prefs.atom_masses[i];
    }


// Norm mode files:
//------------------------------
    for(i = 0; i < prefs.dimension; ++i){
        double norm = 0.0;
        for(j = 0; j < prefs.n_atoms; ++j){
            norm += mode[i][j*3    ]*mode[i][j*3    ]; // x*x
            norm += mode[i][j*3 + 1]*mode[i][j*3 + 1]; // y*y
            norm += mode[i][j*3 + 2]*mode[i][j*3 + 2]; // z*z
        }

        norm = sqrt(norm);

    #ifdef debug_normalmode
        fprintf(stderr, "\t||%s|| = % lf\n", prefs.modelist[i], norm);
    #endif

        for(j = 0; j < prefs.n_atoms; ++j){
            mode[i][j*3    ] /= norm;
            mode[i][j*3 + 1] /= norm;
            mode[i][j*3 + 2] /= norm;
        }
    }


#ifdef debug_normalmode
//{{{ print mode files
    fprintf(stderr, "Number of Atoms: %d\n", prefs.n_atoms);
    for(i = 0; i < prefs.dimension; ++i){
        fprintf(stderr, "#%s:\n#", prefs.modelist[i]);
        fprintf(stderr, "\t    %s      ", "dx");
        fprintf(stderr, "\t    %s      ", "dy");
        fprintf(stderr, "\t    %s      ", "dz");
        fprintf(stderr, "\t    %s", "masses");
        fprintf(stderr, "\n");

        for(j = 0; j < prefs.n_atoms; ++j){
            for(k = 0; k < 3; ++k){
                fprintf(stderr, "\t% .8le", mode[i][j*3 + k]);
            }
            fprintf(stderr, "\t% .8le\n", prefs.atom_masses[j]);
        }
        fprintf(stderr, "\n");
    }
//}}}
#endif


//------------------------------------------------------------------------------
// Calculate Coriolis Coefficients Zeta    Calculate Coriolis Coefficients Zeta
//------------------------------------------------------------------------------
    prefs.zeta = malloc( 3 * sizeof(double*) );
    if(prefs.zeta == NULL){ perror("prefs.zeta"); exit(errno); }
    for(i = 0; i < 3; ++i){
        prefs.zeta[i] = calloc(prefs.dimension * prefs.dimension, sizeof(double));
        if(prefs.zeta[i] == NULL){ perror("prefs.zeta[i]"); exit(errno); }
    }

/* Calculate Coriolis Coefficients
//--------------------------------------------------
//{{{
    The coefficients are calculated for every possible combination of two modes
        zeta_x = dy1*dz2 - dz1*dy2
        zeta_y = dz1*dx2 - dx1*dz2
        zeta_z = dx1*dy2 - dy1*dx2
    where the numbers denote the first and second mode of the respective mode-pair

    Above formula leads to the following findings:
        -> the coefficients of a mode with itself is zero
                f(1,1) = f(2,2) = 0
        -> the coefficients of mode 1 with mode 2 have the same absolute value
            as the ones of mode 2 with mode 1 but with inverted parity:
                f(1,2) = -f(2,1)

    Leading to the conclusion, that only the upper triangle has to be calculated
    and the lower one can be based on the upper triangle results

//}}}*/
    for(i = 0; i < prefs.dimension; ++i){
        for(j = i+1; j < prefs.dimension; ++j){

            for(k = 0; k < prefs.n_atoms; ++k){
                prefs.zeta[0][i*prefs.dimension + j] += mode[i][k*3 + 1]*mode[j][k*3 + 2] - mode[i][k*3 + 2]*mode[j][k*3 + 1];    // dy1*dz2 - dz1*dy2
                prefs.zeta[1][i*prefs.dimension + j] += mode[i][k*3 + 2]*mode[j][k*3    ] - mode[i][k*3    ]*mode[j][k*3 + 2];    // dz1*dx2 - dx1*dz2
                prefs.zeta[2][i*prefs.dimension + j] += mode[i][k*3    ]*mode[j][k*3 + 1] - mode[i][k*3 + 1]*mode[j][k*3    ];    // dx1*dy2 - dy1*dx2
            }

        }
    }

// modes are not required anymore, free up memory
    for(i = 0; i < prefs.dimension; ++i){ free(mode[i]); mode[i] = NULL; }
    free(mode); mode = NULL;

// fill the lower triangle
    for(i = 0; i < prefs.dimension; ++i){
        for(j = i+1; j < prefs.dimension; ++j){

            prefs.zeta[0][j * prefs.dimension + i] = -prefs.zeta[0][i * prefs.dimension + j];
            prefs.zeta[1][j * prefs.dimension + i] = -prefs.zeta[1][i * prefs.dimension + j];
            prefs.zeta[2][j * prefs.dimension + i] = -prefs.zeta[2][i * prefs.dimension + j];

        }
    }


// output Coriolis coefficients
//--------------------------------------------------
    if(prefs.output_file_set){
        prefs.fdout = fopen(prefs.output_file, "w");
    }
    if( prefs.fdout == NULL ){ perror(prefs.output_file); exit(errno); }

    fprintf(prefs.fdout, "#Modes:");
    for(i = 0; i < prefs.dimension; ++i){
        for(j = i+1; j < prefs.dimension; ++j){
            fprintf(prefs.fdout, "\t %s|%s", prefs.modelist[i], prefs.modelist[j]);
        }
    }
    fprintf(prefs.fdout, "\n");

    for(k = 0; k < 3; ++k){
        fprintf(prefs.fdout, "Zeta_%c:", "xyz"[k]);
        for(i = 0; i < prefs.dimension; ++i){
            for(j = i+1; j < prefs.dimension; ++j){
                fprintf(prefs.fdout, "\t% .12le", prefs.zeta[k][prefs.dimension * i + j]);
            }
        }
        fprintf(prefs.fdout, "\n");
    }


// if no coordinate files are provided return early
    if(!prefs.input_coordinates_set){
        fclose(prefs.fdout); prefs.fdout = NULL;
        for(i = 0; i < 3; ++i){ free(prefs.zeta[i]); prefs.zeta[i] = NULL; } free(prefs.zeta);
        prefs.zeta = NULL;

        return 0;
    }


//--------------------------------------------------------------------------------------
//   Moment of inertia    Moment of inertia    Moment of inertia    Moment of inertia
//--------------------------------------------------------------------------------------
// print header
//--------------------------------------------------
    fprintf(prefs.fdout, "#");
    for(i = 0; i < prefs.dimension; ++i){
        fprintf(prefs.fdout,
                    "\t q[%2d]              "
                    , i
        );
    }
    for(i = 0; i < 3; ++i){
        for(j = i; j < 3; ++j){
            fprintf(prefs.fdout,
                    "\t mu_%c%c             "
                    , "xyz"[i], "xyz"[j]
            );
        }
    }
    fprintf(prefs.fdout, "\n");

/* read and process input coordinates list:
//--------------------------------------------------
//{{{
    This file is of the following structure:
        <dimension> columns representing the displacement
            1       column  representing the name of a file
                    containing the system's coordinates

    For each coordinate file in the list the Effective Reciprocal Moment of Inertia Tensor
    is evaluated and printed to the output file.

//}}}*/
    ProcessFileList(prefs);

// close files and free unused memory
    fclose(prefs.fdout); prefs.fdout = NULL;
    for(i = 0; i < 3; ++i){ free(prefs.zeta[i]); prefs.zeta[i] = NULL; } free(prefs.zeta);
    prefs.zeta = NULL;
    free(mass[0]); mass[0] = NULL;
    free(mass);    mass    = NULL;

    return 0;
}
