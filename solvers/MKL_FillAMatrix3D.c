
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <mkl_solvers_ee.h>

#include "settings.h"
#include "MKLFillers.h"
#include "Watson.h"


int MKL_FillAMatrix3D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, double** q, double dq, double*** mu, double** zeta, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A){

// Calculate the maximum number of non-zero entries in the A matrix
//  Should be <n_points - (n_stencil/2)*2> lines with <n_stencil> entries, the
//  first and last <n_stencil/2> lines are shortened via a triangular number sequence
    int n_points    = nq[0] * nq[1] * nq[2];
    int max_entries = (prefs->n_stencil*nq[0] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1))
                     *(prefs->n_stencil*nq[1] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1))
                     *(prefs->n_stencil*nq[2] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1));

// The sparse matrix eigensolver of Intel MKL saves the positions of non empty entries in
//  two integer arrays (rows_A & cols_A), each non zero entry increments the counter by 1
//  rows_A contains the counter value of each first element in a line -> max entries n_points (+1 because Intel...)
//  cols_A contains the counter values of all non zero elements -> max entries n_stencil*n_points - borders
//  vals_A contains the values of all non zero entries -> max entries n_stencil*n_points - borders
    (*rows_A) = calloc(n_points + 1, sizeof(MKL_INT));
    (*cols_A) = calloc(max_entries,  sizeof(MKL_INT));
    (*vals_A) = calloc(max_entries,  sizeof(double));

    if( (*rows_A) == NULL ){ perror("3D MKL Fill rows_A"); exit(errno); }
    if( (*cols_A) == NULL ){ perror("3D MKL Fill cols_A"); exit(errno); }
    if( (*vals_A) == NULL ){ perror("3D MKL Fill vals_A"); exit(errno); }


// initialise for the calculation of the Watson Hamiltonian
    if( prefs->coriolis_file ){ init_watson(prefs, dq); }

// fill Numerov's A matrix
//  determine the non zero elements and store their positions in rows_A and cols_A
//  and their values in vals_A
    (*rows_A)[0] = 1;
    int entry_index = 0;

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

            int xidx = i + xsh;
            int yidx = j + ysh;
            int zidx = k + zsh;

            int row        = (i*nq[1] + j)*nq[2] + k;
            int stencilidx = (s[0]*prefs->n_stencil + s[1])*prefs->n_stencil + s[2];

        // set column index
            (*cols_A)[entry_index] = (xidx*nq[1] + yidx)*nq[2] + zidx + 1;

        // set matrix value, the stencil values have to be divided by 2^(D-1)
            (*vals_A)[entry_index] = ekin_to_oue * 0.25*stencil[ stencilidx ];
        //  apply second term of Watson Hamiltonian
            if( prefs->coriolis_file ){
                (*vals_A)[entry_index] += exec_watson_3d(mu, zeta, q, row, s);
            }

        // add potential to diagonal element
            if( (xsh == 0) && (ysh == 0) && (zsh == 0) ){
                (*vals_A)[entry_index] += v[ (i*nq[1] + j)*nq[2] + k ];
            }
            ++entry_index;
        }}
        }}
        }}
    // after inserting all entries in a row the total number of entries is inserted in the CSR format.
        (*rows_A)[((i*nq[1] + j)*nq[2] + k) + 1] = entry_index + 1;
    }
    }
    }

// free memory
    if( prefs->coriolis_file ){ free_watson(); }

    printf("Matrix created, Potential added, %d entries\n", entry_index);
    return 0;
}
