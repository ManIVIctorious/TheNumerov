
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <mkl_solvers_ee.h>

#include "settings.h"
#include "MKLFillers.h"

// dependencies
void HeapSort(MKL_INT* array, double* values, int arraysize);


int MKL_FillPeriodicAMatrix2D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A){

// Calculate the maximum number of non-zero entries in the A matrix
//  Should be <n_points - (n_stencil/2)*2> lines with <n_stencil> entries, the
//  first and last <n_stencil/2> lines are shortened via a triangular number sequence
    int n_points    = nq[0] * nq[1];
    int midpoint    = prefs->n_stencil/2;
    int max_entries = (prefs->n_stencil*nq[0])
                     *(prefs->n_stencil*nq[1]);

// The sparse matrix eigensolver of Intel MKL saves the positions of non empty entries in
//  two integer arrays (rows_A & cols_A), each non zero entry increments the counter by 1
//  rows_A contains the counter value of each first element in a line -> max entries n_points (+1 because Intel...)
//  cols_A contains the counter values of all non zero elements -> max entries n_stencil*n_points
//  vals_A contains the values of all non zero entries -> max entries n_stencil*n_points
    (*rows_A) = calloc(n_points + 1, sizeof(MKL_INT));
    (*cols_A) = calloc(max_entries,  sizeof(MKL_INT));
    (*vals_A) = calloc(max_entries,  sizeof(double));

    if( (*rows_A) == NULL ){ perror("2D MKL Fill rows_A"); exit(errno); }
    if( (*cols_A) == NULL ){ perror("2D MKL Fill cols_A"); exit(errno); }
    if( (*vals_A) == NULL ){ perror("2D MKL Fill vals_A"); exit(errno); }


// fill Numerov's A matrix
//  determine the non zero elements and store their positions in rows_A and cols_A
//  and their values in vals_A
    (*rows_A)[0] = 1;
    int entry_index = 0;

    for(int i = 0; i < nq[0]; ++i){
    for(int j = 0; j < nq[1]; ++j){

        for(int xsh = 0; xsh < prefs->n_stencil; ++xsh){
        for(int ysh = 0; ysh < prefs->n_stencil; ++ysh){

            int xidx = ( (i + xsh) + (nq[0] - prefs->n_stencil/2) ) % nq[0];
            int yidx = ( (j + ysh) + (nq[1] - prefs->n_stencil/2) ) % nq[1];

            int stencilidx = xsh*prefs->n_stencil + ysh;

        // set column index
            (*cols_A)[entry_index] = xidx*nq[1] + yidx + 1;

        // set matrix value, the stencil values have to be divided by 2^(D-1)
            (*vals_A)[entry_index] = ekin_to_oue * 0.5*stencil[ stencilidx ];

        // add potential to diagonal element
            if( (xsh == midpoint) && (ysh == midpoint) ){
                (*vals_A)[entry_index] += v[i*nq[1] + j];
            }
            ++entry_index;
        }
        }
    // after inserting all entries in a row the total number of entries is inserted in the CSR format.
        (*rows_A)[(i*nq[1] + j) + 1] = entry_index + 1;
    }
    }

// The MKL CSR format requires the matrix entries to be in order, i.e. from left to right.
// Due to the modulo arithmetics this cannot be ascertained by the filling routine,
// hence values are sorted in retrospection
    for(int i = 0; i < n_points; ++i){
        int       width  = (*rows_A)[i+1] - (*rows_A)[i];
        MKL_INT * colidx = (*cols_A) + (*rows_A)[ i ] - 1;
        double  * values = (*vals_A) + (*rows_A)[ i ] - 1;

        HeapSort(colidx, values, width);
    }

    printf("Matrix created, Potential added, %d entries\n", entry_index);
    return 0;
}
