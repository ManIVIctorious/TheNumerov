
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "data.h"
#include "settings.h"
#include "ConvertString.h"

void SetZetaAndDimension(settings* set, data* data){

    char ** zxarray = NULL;
    char ** zyarray = NULL;
    char ** zzarray = NULL;

    int nzx = splitstring_to_array(set->zeta_x, &zxarray, ":", 0);
    int nzy = splitstring_to_array(set->zeta_y, &zyarray, ":", 0);
    int nzz = splitstring_to_array(set->zeta_z, &zzarray, ":", 0);

// check if all Coriolis coefficients have the same number of fields
    if( nzx != nzy || nzx != nzz ){
        fprintf(stderr,
            "\n (-) Error: Different number of entries in the Coriolis coefficients."
            "\n     (n_zeta_x, n_zeta_y, n_zeta_z) = (%d, %d, %d)"
            "\n     Aborting...\n\n"
            , nzx, nzy, nzz
        );
        exit(EXIT_FAILURE);
    }

// since all coefficients are provided and have the same number of entries,
// calculate the dimension from the number of Coriolis coefficients.
    for(int D = 2, n_zeta = 0; n_zeta < 2*nzx; ++D){
        n_zeta = D*(D-1);
        data->dimension = D;
    }

    if( nzx*2 != data->dimension*(data->dimension-1) ){
        fprintf(stderr,
            "\n (-) Error: Invalid number of Coriolis coefficients (%d) provided."
            "\n     The number must be of the type <dimension>*(<dimension> - 1)."
            "\n     Aborting..."
            "\n\n", nzx
        );
        exit(EXIT_FAILURE);
    }

// for each direction of space {x,y,z} zeta contains n_modes x n_modes entries
// resulting in zeta[{x,y,z}][ n_modes*n_modes ]
    data->zeta = malloc( 3 * sizeof(double*) );
    if( !data->zeta ){ perror("data->zeta"); exit(errno); }

    data->zeta[0] = calloc( data->dimension*data->dimension, sizeof(double) );
    data->zeta[1] = calloc( data->dimension*data->dimension, sizeof(double) );
    data->zeta[2] = calloc( data->dimension*data->dimension, sizeof(double) );
    if( !(data->zeta[0]) ){ perror("data->zeta[0]"); exit(errno); }
    if( !(data->zeta[1]) ){ perror("data->zeta[1]"); exit(errno); }
    if( !(data->zeta[2]) ){ perror("data->zeta[2]"); exit(errno); }


    for(int i = 0, count = 0; i < data->dimension; ++i){
        for(int j = i+1; j < data->dimension; ++j, ++count){
        // upper triangle
            data->zeta[0][i*data->dimension + j] = convertstring_to_double(zxarray[count], "zxarray", NULL);
            data->zeta[1][i*data->dimension + j] = convertstring_to_double(zyarray[count], "zyarray", NULL);
            data->zeta[2][i*data->dimension + j] = convertstring_to_double(zzarray[count], "zzarray", NULL);
        // lower triangle
            data->zeta[0][j*data->dimension + i] = -data->zeta[0][i*data->dimension + j];
            data->zeta[1][j*data->dimension + i] = -data->zeta[1][i*data->dimension + j];
            data->zeta[2][j*data->dimension + i] = -data->zeta[2][i*data->dimension + j];
        }
    }

// free memory of z{x,y,z}array strings
    for(int i = 0; i < nzx; ++i){ free(zxarray[i]); }
    for(int i = 0; i < nzy; ++i){ free(zyarray[i]); }
    for(int i = 0; i < nzz; ++i){ free(zzarray[i]); }
    free(zxarray);
    free(zyarray);
    free(zzarray);
}
