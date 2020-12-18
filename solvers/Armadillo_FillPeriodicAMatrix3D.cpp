
#include <stdio.h>
#include <armadillo>

#include "settings.h"
#include "ArmadilloFillers.h"


arma::sp_mat FillPeriodicArmadillo_3D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil){

// Calculate the maximum number of non-zero entries in the A matrix
//  Should be <n_points - (n_stencil/2)*2> lines with <n_stencil> entries
    int n_points    = nq[0] * nq[1] * nq[2];
    int max_entries = (prefs->n_stencil*nq[0])
                     *(prefs->n_stencil*nq[1])
                     *(prefs->n_stencil*nq[2]);


// determine positions and values
    arma::umat locations(2, max_entries);
    arma::vec  values = arma::zeros(max_entries);
    unsigned int entry_index = 0;

    for(int i = 0; i < nq[0]; ++i){
    for(int j = 0; j < nq[1]; ++j){
    for(int k = 0; k < nq[2]; ++k){

        for(int xsh = 0; xsh < prefs->n_stencil; ++xsh){
        for(int ysh = 0; ysh < prefs->n_stencil; ++ysh){
        for(int zsh = 0; zsh < prefs->n_stencil; ++zsh){

            int xidx = ( (i + xsh) + (nq[0] - prefs->n_stencil/2) ) % nq[0];
            int yidx = ( (j + ysh) + (nq[1] - prefs->n_stencil/2) ) % nq[1];
            int zidx = ( (k + zsh) + (nq[2] - prefs->n_stencil/2) ) % nq[2];

        // auxiliary indices
            int row        = (  i  *nq[1] +  j   )*nq[2] +  k;
            int col        = ( xidx*nq[1] + yidx )*nq[2] + zidx;
            int stencilidx = (xsh*prefs->n_stencil + ysh)*prefs->n_stencil + zsh;

        // locations of stencil values (rows 0 and columns 1)
            locations(0, entry_index) = row;
            locations(1, entry_index) = col;
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //  The stencil values have to be divided by 2^(D-1)
            values(entry_index) = ekin_to_oue * 0.25*stencil[ stencilidx ];
        // increment number of entries
            ++entry_index;
        }
        }
        }

    }
    }
    }

// actual sparse matrix fill
    arma::sp_mat A(locations, values, n_points, n_points, true, true);

// add potential value
    for(int i = 0; i < n_points; ++i){ A(i,i) += v[i]; }

    printf("Matrix created, Potential added, %u entries\n", entry_index);
    return A;
}
