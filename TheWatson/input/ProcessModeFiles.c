
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

#include "data.h"
#include "settings.h"

// dependencies
int InputNormalMode(char* inputfile, double* *mode, char read_masses, double* *mass);


void ProcessModeFiles(settings* set, data* data){

// save data of each mode in modes array and mass in mass array
//  mode[D][entry_row * 3 + 0] = x_displacement
//  mode[D][entry_row * 3 + 1] = y_displacement
//  mode[D][entry_row * 3 + 2] = z_displacement
    double ** mode = malloc(data->dimension * sizeof(double*));
    if(mode == NULL){ perror("mode"); exit(errno); }
    for(int i = 0; i < data->dimension; ++i){ mode[i] = NULL; }

    double ** mass = malloc(data->dimension * sizeof(double*));
    if(mass == NULL){ perror("mass"); exit(errno); }
    for(int i = 0; i < data->dimension; ++i){ mass[i] = NULL; }

    data->n_atoms = InputNormalMode( set->modelist[0], &(mode[0]), set->masses_from_modes, &(mass[0]) );
    for(int i = 1; i < data->dimension; ++i){
        int control = InputNormalMode( set->modelist[i], &(mode[i]), set->masses_from_modes, &(mass[i]) );

    // check if all files contain the same number of atoms
        if( control != data->n_atoms ){
            fprintf(stderr,
                "\n (-) Error in reading mode files. Number of atoms in \"%s\" (\"%d\") does not"
                "\n     match number of atoms in \"%s\" (\"%d\"). Aborting...\n\n"
                , set->modelist[0], data->n_atoms, set->modelist[i], control
            );
            exit(EXIT_FAILURE);
        }
    }

// check for same masses in mode files (atom order/type)
    for(int i = 1; i < data->dimension; ++i){
        for(int j = 0; j < data->n_atoms; ++j){

            if( (mass[0][j] - mass[i][j])*(mass[0][j] - mass[i][j]) > set->threshold*set->threshold ){
                fprintf(stderr,
                        "\n (-) Error in reading mode files. Atom masses differ considerably"
                        "\n     between modes in entry number %d."
                        "\n             %s\t%s"
                        "\n     Value:  %le\t%le"
                        "\n     Aborting...\n\n"
                        , j, set->modelist[i-1], set->modelist[i], mass[i-1][j], mass[i][j]
                    );
                exit(EXIT_FAILURE);
            }
        }
    }

// if masses check succeeds free up memory:
//  so that only mass[0][entries] remains
    for(int i = 1; i < data->dimension; ++i){
        free( mass[i] ); mass[i] = NULL;
    }
// point data->atom_masses to the position of mass[0]
// now also mass, holding the pointers to mass[0] can be freed
    data->atom_masses = mass[0];
    free(mass);

// calculate the system's total mass
    data->tot_mass = 0.0;
    for(int i = 0; i < data->n_atoms; ++i){
        data->tot_mass += data->atom_masses[i];
    }


// Normalise mode files:
//------------------------------
    for(int i = 0; i < data->dimension; ++i){
        double norm = 0.0;
        for(int j = 0; j < data->n_atoms; ++j){
            norm += mode[i][j*3 + 0]*mode[i][j*3 + 0]; // x*x
            norm += mode[i][j*3 + 1]*mode[i][j*3 + 1]; // y*y
            norm += mode[i][j*3 + 2]*mode[i][j*3 + 2]; // z*z
        }

        norm = sqrt(norm);

    #ifdef debug_normalmode
        fprintf(stderr, "\t||%s|| = % lf\n", set->modelist[i], norm);
    #endif

        for(int j = 0; j < data->n_atoms; ++j){
            mode[i][j*3 + 0] /= norm;
            mode[i][j*3 + 1] /= norm;
            mode[i][j*3 + 2] /= norm;
        }
    }


#ifdef debug_normalmode
//{{{ print mode files
    fprintf(stderr, "Number of Atoms: %d\n", data->n_atoms);
    for(int i = 0; i < data->dimension; ++i){
        fprintf(stderr, "#%s:\n#", set->modelist[i]);
        fprintf(stderr, "\t    %s      ", "dx");
        fprintf(stderr, "\t    %s      ", "dy");
        fprintf(stderr, "\t    %s      ", "dz");
        fprintf(stderr, "\t    %s", "masses");
        fprintf(stderr, "\n");

        for(int j = 0; j < data->n_atoms; ++j){
            for(int k = 0; k < 3; ++k){
                fprintf(stderr, "\t% .8le", mode[i][j*3 + k]);
            }
            fprintf(stderr, "\t% .8le\n", data->atom_masses[j]);
        }
        fprintf(stderr, "\n");
    }
//}}}
#endif


//------------------------------------------------------------------------------
// Calculate Coriolis Coefficients Zeta    Calculate Coriolis Coefficients Zeta
//------------------------------------------------------------------------------
// for each direction of space {x,y,z} zeta contains n_modes x n_modes entries
// resulting in zeta[{x,y,z}][ n_modes*n_modes ]
    data->zeta = malloc( 3 * sizeof(double*) );
    if(data->zeta == NULL){ perror("data->zeta"); exit(errno); }
    for(int i = 0; i < 3; ++i){
        data->zeta[i] = calloc(data->dimension * data->dimension, sizeof(double));
        if(data->zeta[i] == NULL){ perror("data->zeta[i]"); exit(errno); }
    }

/* Calculate Coriolis Coefficients zeta
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

    Hence, only the upper triangle has to be calculated
//}}}*/
    for(int i = 0; i < data->dimension; ++i){
        for(int j = i+1; j < data->dimension; ++j){
            for(int k = 0; k < data->n_atoms; ++k){
                data->zeta[0][i*data->dimension + j] += mode[i][k*3 + 1]*mode[j][k*3 + 2] - mode[i][k*3 + 2]*mode[j][k*3 + 1];    // dy1*dz2 - dz1*dy2
                data->zeta[1][i*data->dimension + j] += mode[i][k*3 + 2]*mode[j][k*3 + 0] - mode[i][k*3 + 0]*mode[j][k*3 + 2];    // dz1*dx2 - dx1*dz2
                data->zeta[2][i*data->dimension + j] += mode[i][k*3 + 0]*mode[j][k*3 + 1] - mode[i][k*3 + 1]*mode[j][k*3 + 0];    // dx1*dy2 - dy1*dx2
            }

        }
    }

// modes are not required anymore, free up memory
    for(int i = 0; i < data->dimension; ++i){ free(mode[i]); mode[i] = NULL; }
    free(mode); mode = NULL;

// fill the lower triangle
    for(int i = 0; i < data->dimension; ++i){
        for(int j = i+1; j < data->dimension; ++j){

            data->zeta[0][j * data->dimension + i] = -data->zeta[0][i * data->dimension + j];
            data->zeta[1][j * data->dimension + i] = -data->zeta[1][i * data->dimension + j];
            data->zeta[2][j * data->dimension + i] = -data->zeta[2][i * data->dimension + j];

        }
    }
}
