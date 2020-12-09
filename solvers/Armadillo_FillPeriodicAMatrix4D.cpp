
#include <stdio.h>
#include <armadillo>
#include "settings.h"

// Provided Prototypes
arma::sp_mat FillPeriodicArmadillo_4D(settings* prefs, int* nq, int n_points, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta);


// 3D fill
arma::sp_mat FillPeriodicArmadillo_4D(settings* prefs, int* nq, int n_points, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta){

// Calculate the maximum number of non-zero entries in the A matrix
//  Should be <n_points - (n_stencil/2)*2> lines with <n_stencil> entries
    int max_entries = (prefs->n_stencil*nq[0])
                     *(prefs->n_stencil*nq[1])
                     *(prefs->n_stencil*nq[2])
                     *(prefs->n_stencil*nq[3]);


// determine positions and values
    arma::umat locations(2, max_entries);
    arma::vec  values = arma::zeros(max_entries);
    unsigned int entry_index = 0;

    for(int i = 0; i < nq[0]; ++i){
    for(int j = 0; j < nq[1]; ++j){
    for(int k = 0; k < nq[2]; ++k){
    for(int l = 0; l < nq[3]; ++l){

        for(int ash = 0; ash < prefs->n_stencil; ++ash){
        for(int bsh = 0; bsh < prefs->n_stencil; ++bsh){
        for(int csh = 0; csh < prefs->n_stencil; ++csh){
        for(int dsh = 0; dsh < prefs->n_stencil; ++dsh){

            int aidx = ( (i + ash) + (nq[0] - prefs->n_stencil/2) ) % nq[0];
            int bidx = ( (j + bsh) + (nq[1] - prefs->n_stencil/2) ) % nq[1];
            int cidx = ( (k + csh) + (nq[2] - prefs->n_stencil/2) ) % nq[2];
            int didx = ( (l + dsh) + (nq[3] - prefs->n_stencil/2) ) % nq[3];

        // locations of stencil values (rows 0 and columns 1)
            locations (0, entry_index) = ( ( i  *nq[1] +  j)  *nq[2] +  k   )*nq[3] +  l;
            locations (1, entry_index) = ( (aidx*nq[1] + bidx)*nq[2] + cidx )*nq[3] + didx;
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            values(entry_index) = ekin_param * stencil[ ( (ash*prefs->n_stencil + bsh)*prefs->n_stencil + csh )*prefs->n_stencil + dsh ];
        // stencil entries have to be divided by 2^(D-1)
            values(entry_index) *= 0.125;
        // increment number of entries
            ++entry_index;
        }
        }
        }
        }

    }
    }
    }
    }

// actual sparse matrix fill
    arma::sp_mat A(locations, values, n_points, n_points, true, true);

// add potential value
    for(int i = 0; i < n_points; ++i){
        A(i,i) += v[i];
    }

    printf("Matrix created, Potential added, %u entries\n", entry_index);
    return A;
}
