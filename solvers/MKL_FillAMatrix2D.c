
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <mkl_solvers_ee.h>

#include "settings.h"

// Dependencies
void   init_watson_2d(settings* prefs);
double exec_watson_2d(double*** mu, double** zeta, int* nq, double dq, double** q, int i, int j, int xsidx, int ysidx);
void   free_watson_2d(void);

// provided prototypes
int MKL_FillAMatrix2D(settings* prefs, int* nq, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);


// 2D fill
int MKL_FillAMatrix2D(settings* prefs, int* nq, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A){

// Calculate the maximum number of non-zero entries in the A matrix
//  Should be <n_points - (n_stencil/2)*2> lines with <n_stencil> entries, the
//  first and last <n_stencil/2> lines are shortened via a triangular number sequence
    int n_points = nq[0] * nq[1];
    int max_entries = (prefs->n_stencil*nq[0] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1))
                     *(prefs->n_stencil*nq[1] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1));

// The sparse matrix eigensolver of Intel MKL saves the positions of non empty entries in
//  two integer arrays (rows_A & cols_A), each non zero entry increments the counter by 1
//  rows_A contains the counter value of each first element in a line -> max entries n_points (+1 because Intel...)
//  cols_A contains the counter values of all non zero elements -> max entries n_stencil*n_points - borders
//  vals_A contains the values of all non zero entries -> max entries n_stencil*n_points - borders
    (*rows_A) = calloc(n_points + 1, sizeof(MKL_INT));
    (*cols_A) = calloc(max_entries,  sizeof(MKL_INT));
    (*vals_A) = calloc(max_entries,  sizeof(double));

    if( (*rows_A) == NULL ){ perror("2D MKL Fill rows_A"); exit(errno); }
    if( (*cols_A) == NULL ){ perror("2D MKL Fill cols_A"); exit(errno); }
    if( (*vals_A) == NULL ){ perror("2D MKL Fill vals_A"); exit(errno); }


// initialise for the calculation of the Watson Hamiltonian
    if( prefs->coriolis_file ){ init_watson_2d(prefs); }

// fill Numerov's A matrix
//  determine the non zero elements and store their positions in rows_A and cols_A
//  and their values in vals_A
    (*rows_A)[0] = 1;
    int entry_index = 0;

    for(int i = 0; i < nq[0]; ++i){
    for(int j = 0; j < nq[1]; ++j){

        for(int xsh = -(prefs->n_stencil/2); xsh < ((prefs->n_stencil/2) + 1); ++xsh){
        if( (i + xsh > -1) && (i + xsh < nq[0]) ){

            for(int ysh = -(prefs->n_stencil/2); ysh < ((prefs->n_stencil/2) + 1); ++ysh){
            if( (j + ysh > -1) && (j + ysh < nq[1]) ){

                int xsidx = xsh + prefs->n_stencil/2;    // stencil x index
                int ysidx = ysh + prefs->n_stencil/2;    // stencil y index

            // set column index
                (*cols_A)[entry_index] = (i + xsh)*nq[1] + (j + ysh) + 1;

            // set matrix value
                (*vals_A)[entry_index] = ekin_param * stencil[ xsidx*prefs->n_stencil + ysidx ];
            //  apply second term of Watson Hamiltonian
                if( prefs->coriolis_file ){
                    (*vals_A)[entry_index] -= exec_watson_2d(mu, zeta, nq, dq, q, i, j, xsidx, ysidx);
                }
            //  The stencil values have to be divided by 2^(D-1)
                (*vals_A)[entry_index] *= 0.5;

            // add potential to diagonal element
                if( (xsh == 0) && (ysh == 0) ){
                    (*vals_A)[entry_index] += v[i*nq[1] + j];
                }
                ++entry_index;
            }}
        }}
    // after inserting all entries in a row the total number of entries is inserted in the CSR format.
        (*rows_A)[i*nq[1] + (j + 1)] = entry_index + 1;
    }
    }

// free memory
    if( prefs->coriolis_file ){ free_watson_2d(); }

    printf("Matrix created, Potential added, %d entries\n", entry_index);
    return 0;
}
