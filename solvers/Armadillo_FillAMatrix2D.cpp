
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
arma::sp_mat FillArmadillo_2D(settings* prefs, int* nq, int n_points, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta);


// 2D fill
arma::sp_mat FillArmadillo_2D(settings* prefs, int* nq, int n_points, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta){

// Calculate the maximum number of non-zero entries in the A matrix
//  Should be <n_points - (n_stencil/2)*2> lines with <n_stencil> entries, the
//  first and last <n_stencil/2> lines are shortened via a triangular number sequence
    int max_entries = (prefs->n_stencil*nq[0] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1))
                     *(prefs->n_stencil*nq[1] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1));


// initialise for the calculation of the Watson Hamiltonian
    if( prefs->coriolis_file ){ init_watson_2d(prefs); }

// determine positions and values
    arma::umat locations(2, max_entries);
    arma::vec  values = arma::zeros(max_entries);
    unsigned int entry_index = 0;

    for(int i = 0; i < nq[0]; ++i){
        for(int j = 0; j < nq[1]; ++j){

            for(int xsh = -(prefs->n_stencil/2); xsh < ((prefs->n_stencil/2) + 1); ++xsh){
            if( (i + xsh > -1) && (i + xsh < nq[0]) ){

                for(int ysh = -(prefs->n_stencil/2); ysh < ((prefs->n_stencil/2) + 1); ++ysh){
                if( (j + ysh > -1) && (j + ysh < nq[1]) ){

                    int xsidx = xsh + prefs->n_stencil/2;
                    int ysidx = ysh + prefs->n_stencil/2;

                // locations of stencil values
                    locations(0, entry_index) =    i     *nq[1] +    j;
                    locations(1, entry_index) = (i + xsh)*nq[1] + (j + ysh);
                //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    values(entry_index) = ekin_param * stencil[ xsidx*prefs->n_stencil + ysidx ];
                //  apply second term of Watson Hamiltonian
                    if( prefs->coriolis_file ){
                        values(entry_index) -= exec_watson_2d(mu, zeta, nq, dq, q, i, j, xsidx, ysidx);
                    }
                //  The stencil values have to be divided by a factor of two
                    values(entry_index) *= 0.5;

                    ++entry_index;
                }}
            }}

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
