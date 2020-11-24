
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <mkl_solvers_ee.h>

// provided prototypes
int MKL_FillAMatrix1D(double* v, int* nq, double ekin_param, double* stencil, int n_stencil, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);


// 1D fill
int MKL_FillAMatrix1D(double* v, int* nq, double ekin_param, double* stencil, int n_stencil, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A){

// Calculate the maximum number of non-zero entries in the A matrix
//  Should be <n_points - (n_stencil/2)*2> lines with <n_stencil> entries, the
//  first and last <n_stencil/2> lines are shortened via a triangular number sequence
    int n_points = nq[0];
    int max_entries = n_stencil*n_points - (n_stencil/2)*(n_stencil/2 + 1);

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
    int entry_index = 0;

    for(int i = 0; i < n_points; ++i){

        for(int xsh = -(n_stencil/2); xsh < ((n_stencil/2) + 1); ++xsh){

            int element = i + xsh;

            if( (element > -1) && (element < n_points) ){

                (*vals_A)[entry_index] = ekin_param * stencil[xsh + n_stencil/2];
                (*cols_A)[entry_index] = element + 1;

            // add potential to diagonal element
                if( xsh == 0 ){
                    (*vals_A)[entry_index] += v[i];
                }
                ++entry_index;
            }
        }
    // after inserting all entries in a row the total number of entries is inserted in the CSR format.
        (*rows_A)[i+1] = entry_index + 1;
    }

    printf("Matrix created, Potential added, %d entries\n", entry_index);
    return 0;
}
