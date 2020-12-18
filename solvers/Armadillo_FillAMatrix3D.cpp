
#include <stdio.h>
#include <armadillo>

#include "settings.h"
#include "ArmadilloFillers.h"
extern "C" {
#include "Watson.h"
}

// Provided Prototypes
arma::sp_mat FillArmadillo_3D(settings* prefs, int* nq, int n_points, double* v, double ekin_to_oue, double* stencil, double** q, double dq, double*** mu, double** zeta);


arma::sp_mat FillArmadillo_3D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, double** q, double dq, double*** mu, double** zeta){

// Calculate the maximum number of non-zero entries in the A matrix
//  Should be <n_points - (n_stencil/2)*2> lines with <n_stencil> entries, the
//  first and last <n_stencil/2> lines are shortened via a triangular number sequence
    int n_points    = nq[0] * nq[1] * nq[2];
    int max_entries = (prefs->n_stencil*nq[0] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1))
                     *(prefs->n_stencil*nq[1] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1))
                     *(prefs->n_stencil*nq[2] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1));


// initialise for the calculation of the Watson Hamiltonian
    if( prefs->coriolis_file ){ init_watson(prefs); }

// determine positions and values
    arma::umat locations(2, max_entries);
    arma::vec  values = arma::zeros(max_entries);
    unsigned int entry_index = 0;

    for(int i = 0; i < nq[0]; ++i){
    for(int j = 0; j < nq[1]; ++j){
    for(int k = 0; k < nq[2]; ++k){

        for(int xsh = -(prefs->n_stencil/2); xsh <= (prefs->n_stencil/2); ++xsh){
        if( (i + xsh > -1) && (i + xsh < nq[0]) ){

        for(int ysh = -(prefs->n_stencil/2); ysh <= (prefs->n_stencil/2); ++ysh){
        if( (j + ysh > -1) && (j + ysh < nq[1]) ){

        for(int zsh = -(prefs->n_stencil/2); zsh <= (prefs->n_stencil/2); ++zsh){
        if( (k + zsh > -1) && (k + zsh < nq[2]) ){

        // auxiliary indices
            int s[3];
            s[0] = xsh + prefs->n_stencil/2;
            s[1] = ysh + prefs->n_stencil/2;
            s[2] = zsh + prefs->n_stencil/2;

            int row        = (   i    *nq[1] +   j     )*nq[2] +   k;
            int col        = ( (i+xsh)*nq[1] + (j+ysh) )*nq[2] + (k+zsh);
            int stencilidx = (s[0]*prefs->n_stencil + s[1])*prefs->n_stencil + s[2];

        // locations of stencil values (rows 0 and columns 1)
            locations(0, entry_index) = row;
            locations(1, entry_index) = col;
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //  The stencil values have to be divided by 2^(D-1)
            values(entry_index) = ekin_to_oue * 0.25*stencil[ stencilidx ];
//        //  apply second term of Watson Hamiltonian
//            if( prefs->coriolis_file ){
//                values(entry_index) -= exec_watson(mu, zeta, nq, dq, q, row, s);
//            }
        // increment number of entries
            ++entry_index;
        }}
        }}
        }}

    }
    }
    }

// actual sparse matrix fill
    arma::sp_mat A(locations, values, n_points, n_points, true, true);

// free memory
    if( prefs->coriolis_file ){ free_watson(); }

// add potential values to main diagonal
    for(int i = 0; i < n_points; ++i){ A(i,i) += v[i]; }

    printf("Matrix created, Potential added, %u entries\n", entry_index);
    return A;
}
