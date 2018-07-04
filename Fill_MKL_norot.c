
#include <stdio.h>
#include <stdlib.h>
#include <mkl_solvers_ee.h>

// provided prototypes
int FillMKL_2D_norot(double* v, int* nq, double ekin_param, double* stencil, int n_stencil, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);


// 2D fill without rotational terms
int FillMKL_2D_norot(double* v, int* nq, double ekin_param, double* stencil, int n_stencil, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A){

    int i, j;
    int n_entries = 0;
    int n_points = nq[0] * nq[1];
    int xsh, ysh;
    int element;

    // calculate max_entries
    int max_entries = 0;
    int sum_q1 = nq[0];
    int sum_q2 = nq[1];

    for(i = 1; i < (n_stencil/2 + 1); i++){
        sum_q1 += 2*(nq[0]-i);
        sum_q2 += 2*(nq[1]-i);
    }
    max_entries = sum_q1*sum_q2; // upper estimation for nnz entries in the matrix, but the easy way to code.

// The sparse matrix eigensolver of Intel MKL saves the positions of non empty entries in
//  two integer arrays (rows_A & cols_A), each non zero entry increments the counter by 1
//  rows_A contains the counter value of each first element in a line -> max entries n_points (+1 because Intel...)
//  cols_A contains the counter values of all non zero elements -> max entries n_stencil*n_points - borders
//  vals_A contains the values of all non zero entries -> max entries n_stencil*n_points - borders
    (*rows_A) = calloc(n_points + 1, sizeof(MKL_INT));
    (*cols_A) = calloc(max_entries,  sizeof(MKL_INT));
    (*vals_A) = calloc(max_entries,  sizeof(double));
    if((*rows_A) == NULL || (*cols_A) == NULL || (*vals_A) == NULL){
        fprintf(stderr,
            "\n (-) Error in memory allocation for row, column and/or value arrays"
            "\n     Aborting..."
            "\n\n"
        );
        exit(1);
    }

// fill Numerov's A matrix
//  determine the non zero elements and store their positions in rows_A and cols_A
//  and their values in vals_A
    for(i = 0; i < nq[0]; ++i){
        for(j = 0; j < nq[1]; ++j){
            for(xsh = -n_stencil/2; xsh < n_stencil/2 + 1; ++xsh){

                if( (i+xsh > -1) && (i+xsh < nq[0]) ){
                    for(ysh = -n_stencil/2; ysh < n_stencil/2 + 1; ++ysh){

                        if( (j+ysh > -1) && ( j+ysh < nq[1]) ){
                            element = (i + xsh)*nq[1] + j+ysh;
                            (*cols_A)[n_entries] = element+1; // wieso +1? weil intel!!

                        // stencil entries have to be divided by 2 to get the right result.
                        //  in three dimensions it should be a division by 4
                            (*vals_A)[n_entries] = ekin_param * stencil[(xsh+n_stencil/2)*n_stencil+ysh+n_stencil/2]/2;

                        // add potential to diagonal element
                            if(xsh == 0 && ysh ==0){
                                (*vals_A)[n_entries] += v[i*nq[1]+j];
                            }

                            n_entries++;
                        }
                    }
                }
            }
      // after inserting all entries in a row the total number of entries is inserted in the CSR format.
        (*rows_A)[i*nq[1] + j + 1] = n_entries + 1;
        }
    }
    (*rows_A)[0] = 1;

    return 0;
}
