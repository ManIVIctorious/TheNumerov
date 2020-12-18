
#include <stdio.h>
#include <armadillo>

#include "settings.h"
#include "ArmadilloFillers.h"


arma::sp_mat FillPeriodicArmadillo_1D(int n_stencil,   int* nq, double* v, double ekin_to_oue, double* stencil){

// Calculate the maximum number of non-zero entries in the A matrix
//  Should be <n_points - (n_stencil/2)*2> lines with <n_stencil> entries
    int n_points    = nq[0];
    int max_entries = n_stencil*n_points;

// determine positions and values
    arma::umat locations(2, max_entries);
    arma::vec  values = arma::zeros(max_entries);
    unsigned int entry_index = 0;

    for(int i = 0; i < nq[0]; ++i){

        for(int xsh = 0; xsh < n_stencil; ++xsh){

            int xidx = ( (i + xsh) + (nq[0] - n_stencil/2) ) % nq[0];

        // locations of stencil values
            locations(0, entry_index) =  i;     // rows
            locations(1, entry_index) = xidx;   // columns
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            values(entry_index) = ekin_to_oue * stencil[xsh];
            ++entry_index;
        }
    }

// actual sparse matrix fill
    arma::sp_mat A(locations, values, n_points, n_points, true, true);

// add potential value
    for(int i = 0; i < n_points; ++i){ A(i,i) += v[i]; }

    printf("Matrix created, Potential added, %u entries\n", entry_index);
    return A;
}
