
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <mkl_solvers_ee.h>

#include "settings.h"

// provided prototypes
int FillMKL_1D(double* v, int* nq, double ekin_param, double* stencil, int n_stencil, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);


// 1D fill
int FillMKL_1D(double* v, int* nq, double ekin_param, double* stencil, int n_stencil, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A){

    int i, xsh;
    int element;
    int n_points = nq[0];
    int max_entries;
    int entry_index;

// calculate max_entries
    for(i = n_stencil, max_entries = 0; i > n_stencil/2; --i){
        max_entries += 2*i;
    }
    max_entries += n_stencil * (n_points - n_stencil + 1);

// The sparse matrix eigensolver of Intel MKL saves the positions of non empty entries in
//  two integer arrays (rows_A & cols_A), each non zero entry increments the counter by 1
//  rows_A contains the counter value of each first element in a line -> max entries n_points (+1 because Intel...)
//  cols_A contains the counter values of all non zero elements -> max entries n_stencil*n_points - borders
//  vals_A contains the values of all non zero entries -> max entries n_stencil*n_points - borders
    (*rows_A) = calloc(n_points + 1, sizeof(MKL_INT));
    (*cols_A) = calloc(max_entries,  sizeof(MKL_INT));
    (*vals_A) = calloc(max_entries,  sizeof(double));

    if( (*rows_A) == NULL ){ perror("1D MKL Fill rows_A"); exit(errno); }
    if( (*cols_A) == NULL ){ perror("1D MKL Fill cols_A"); exit(errno); }
    if( (*vals_A) == NULL ){ perror("1D MKL Fill vals_A"); exit(errno); }


// fill Numerov's A matrix
//  determine the non zero elements and store their positions in rows_A and cols_A
//  and their values in vals_A
    (*rows_A)[0] = 1;
    for(i = 0, entry_index = 0; i < n_points; ++i){
        for(xsh = -(n_stencil/2); xsh < ((n_stencil/2) + 1); ++xsh){

            element = i + xsh;

            if((element > -1) && (element < n_points)){

                (*cols_A)[entry_index] = element + 1;
                (*vals_A)[entry_index] = ekin_param * stencil[xsh + n_stencil/2];

            // add potential to diagonal element
                if(xsh == 0){
                    (*vals_A)[entry_index] += v[i];
                }
                ++entry_index;
            }
        }
    // after inserting all entries in a row the total number of entries is inserted in the CSR format.
        (*rows_A)[i+1] = entry_index + 1;
    }

    return 0;
}

