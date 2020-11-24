
#include <stdio.h>
#include <armadillo>
#include "settings.h"

// Provided Prototypes
arma::sp_mat FillArmadillo_3D(settings* prefs, int* nq, int n_points, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta);


// 3D fill
arma::sp_mat FillArmadillo_3D(settings* prefs, int* nq, int n_points, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta){

// Calculate the maximum number of non-zero entries in the A matrix
//  Should be <n_points - (n_stencil/2)*2> lines with <n_stencil> entries, the
//  first and last <n_stencil/2> lines are shortened via a triangular number sequence
    int max_entries = (prefs->n_stencil*nq[0] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1))
                     *(prefs->n_stencil*nq[1] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1))
                     *(prefs->n_stencil*nq[2] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1));


// determine positions and values
    arma::umat locations(2, max_entries);
    arma::vec  values = arma::zeros(max_entries);
    unsigned int index = 0;

    for(int i = 0; i < nq[0]; ++i){
        for(int j = 0; j < nq[1]; ++j){
            for(int k = 0; k < nq[2]; k++){

                for(int xsh = -(prefs->n_stencil/2); xsh <= (prefs->n_stencil/2); ++xsh){
                if( (i + xsh > -1) && (i + xsh < nq[0]) ){

                    for(int ysh = -(prefs->n_stencil/2); ysh <= (prefs->n_stencil/2); ++ysh){
                    if( (j + ysh > -1) && (j + ysh < nq[1]) ){

                        for(int zsh = -(prefs->n_stencil/2); zsh <= (prefs->n_stencil/2); ++zsh){
                        if( (k + zsh > -1) && (k + zsh < nq[2]) ){

                            int xsidx = xsh + prefs->n_stencil/2;
                            int ysidx = ysh + prefs->n_stencil/2;
                            int zsidx = zsh + prefs->n_stencil/2;

                        // locations of stencil values around the main diagonal
                            locations (0, index) = (    i     *nq[1] +    j      )*nq[2] +    k;
                            locations (1, index) = ( (i + xsh)*nq[1] + (j + ysh) )*nq[2] + (k + zsh);
                        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            values(index) = ekin_param * stencil[ (xsidx*prefs->n_stencil + ysidx*prefs->n_stencil)*prefs->n_stencil + zsidx ] / 4.0;
                            ++index;
                        }}
                    }}
                }}

            }
        }
    }

// actual sparse matrix fill
    arma::sp_mat A(locations, values, n_points, n_points, true, true);

// add potential value
    for(int i = 0; i < n_points; ++i){
        A(i,i) += v[i];
    }

    printf("Matrix created, Potential added, %u entries\n", index);
    return A;
}
