
#include <stdio.h>
#include <armadillo>
#include "settings.h"

// Dependencies
extern "C" {
    void   init_watson_2d(settings* prefs);
    double exec_watson_2d(double*** mu, double** zeta, int* nq, double dq, double** q, int i, int j, int xsidx, int ysidx);
    void   free_watson_2d(void);
}

// Provided Prototypes
arma::sp_mat FillPeriodicArmadillo_2D(settings* prefs, int* nq, int n_points, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta);


// 2D fill
arma::sp_mat FillPeriodicArmadillo_2D(settings* prefs, int* nq, int n_points, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta){

// Calculate the maximum number of non-zero entries in the A matrix
//  Should be <n_points - (n_stencil/2)*2> lines with <n_stencil> entries
    int max_entries = (prefs->n_stencil*nq[0])
                     *(prefs->n_stencil*nq[1]);


// initialise for the calculation of the Watson Hamiltonian
    if( prefs->coriolis_file ){ init_watson_2d(prefs); }

// determine positions and values
    arma::umat locations(2, max_entries);
    arma::vec  values = arma::zeros(max_entries);
    unsigned int entry_index = 0;

    for(int i = 0; i < nq[0]; ++i){
    for(int j = 0; j < nq[1]; ++j){

        for(int xsh = 0; xsh < prefs->n_stencil; ++xsh){
        for(int ysh = 0; ysh < prefs->n_stencil; ++ysh){

            int xidx = ( (i + xsh) + (nq[0] - prefs->n_stencil/2) ) % nq[0];
            int yidx = ( (j + ysh) + (nq[1] - prefs->n_stencil/2) ) % nq[1];

        // locations of stencil values
            locations(0, entry_index) =  i  *nq[1] +  j;    // rows
            locations(1, entry_index) = xidx*nq[1] + yidx;  // columns
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            values(entry_index) = ekin_param * stencil[ xsh*prefs->n_stencil + ysh ];
        //  apply second term of Watson Hamiltonian
            if( prefs->coriolis_file ){
                values(entry_index) -= exec_watson_2d(mu, zeta, nq, dq, q, i, j, xsh, ysh);
            }
        //  The stencil values have to be divided by 2^(D-1)
            values(entry_index) *= 0.5;
        // increment number of entries
            ++entry_index;
        }
        }

    }
    }

// actual sparse matrix fill
    arma::sp_mat A(locations, values, n_points, n_points, true, true);

// free memory
    if( prefs->coriolis_file ){ free_watson_2d(); }

// add potential value
    for(int i = 0; i < n_points; ++i){ A(i,i) += v[i]; }

    printf("Matrix created, Potential added, %u entries\n", entry_index);
    return A;
}
