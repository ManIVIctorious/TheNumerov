
#include <stdio.h>
#include <armadillo>
#include "settings.h"

// Provided Prototypes
arma::sp_mat FillArmadillo_1D(int* nq, double* v, double ekin_param, double* stencil, int n_stencil);


// 1D fill
arma::sp_mat FillArmadillo_1D(int* nq, double* v, double ekin_param, double* stencil, int n_stencil){

// Calculate the matximum number of non-zero entries in the A matrix
//  Should be <n_points - (n_stencil/2)*2> lines with <n_stencil> entries, the
//  first and last <n_stencil/2> lines are shortened via a triangular number sequence
    int n_points = nq[0];
    int max_entries = n_stencil*n_points - (n_stencil/2)*(n_stencil/2 + 1);

// determine positions and values
    arma::umat locations(2, max_entries);
    arma::vec  values = arma::zeros(max_entries);
    unsigned int index = 0;

    for(int i = 0; i < nq[0]; ++i){

        for(int xsh = -(n_stencil/2); xsh < ((n_stencil/2) + 1); ++xsh){
        if( (i + xsh > -1) && (i + xsh < n_points) ){

        // locations of stencil values
            locations(0, index) =     i;       // rows
            locations(1, index) = ( i + xsh ); // columns
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            values(index) = ekin_param * stencil[xsh + n_stencil/2];
            ++index;
        }}

    }

// actual sparse matrix fill
    arma::sp_mat A(locations, values, nq[0], nq[0], true, true);

// add potential values to main diagonal
    for(int i = 0; i < nq[0]; ++i){
        A(i,i) += v[i];
    }

    printf("Matrix created, Potential added, %u entries\n", index);
    return A;
}
