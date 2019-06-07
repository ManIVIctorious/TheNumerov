
#include <stdio.h>
#include <armadillo>
#include "typedefinitions.h"

// Provided Prototypes
arma::sp_mat FillArmadillo_1D(int* nq, double* v, double ekin_param, double* stencil, int n_stencil);


// 1D fill
arma::sp_mat FillArmadillo_1D(int* nq, double* v, double ekin_param, double* stencil, int n_stencil){

    int i;
    int xsh;
    int max_entries = 0;
    double threshold = 1.0E-8;

// determine maximum number of entries (number of stencil entries inside the matrix greater than threshold)
    for(i = 0; i < nq[0]; ++i){

        for(xsh = -n_stencil/2; xsh <= n_stencil/2; ++xsh){
        if( ( (i+xsh) > -1 ) && ( (i+xsh) < nq[0] ) ){

            if( stencil[xsh + n_stencil/2] > threshold || stencil[xsh + n_stencil/2] < threshold ){
                max_entries++;
            }

        }}
    }

// determine positions and values
    arma::umat locations(2, max_entries);
    arma::vec  values = arma::zeros(max_entries);
    unsigned int index = 0;

    for(i = 0; i < nq[0]; ++i){

        for(xsh = -n_stencil/2; xsh <= n_stencil/2; ++xsh){
        if( ( (i+xsh) > -1 ) && ( (i+xsh) < nq[0] ) ){

            if( stencil[xsh + n_stencil/2] > threshold || stencil[xsh + n_stencil/2] < threshold ){
            // locations of stencil values around the main diagonal
                locations(0, index) =     i;       // rows
                locations(1, index) = ( i + xsh ); // columns
            //++++++++++++++++++++++++++++++++++++++++++++++++++
                values(index++) = stencil[xsh + n_stencil/2] * ekin_param;

            }
        }}
    }

// actual sparse matrix fill
    arma::sp_mat A(locations, values, nq[0], nq[0], true, true);

// add potential values to main diagonal
    for(i = 0; i < nq[0]; ++i){
        A(i,i) += v[i];
    }

    printf("Matrix created, Potential added, %u entries\n", index);
    return A;
}
