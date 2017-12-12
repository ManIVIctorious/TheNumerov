
#include <stdio.h>
#include <stdlib.h>
#include <mkl_solvers_ee.h>

// provided prototypes
int EigensolverFEAST_MKL_2D(double *v, int n_points, int *nq, double *stencil, int n_stencil, double e_min, double e_max, double *E, double *X);


int EigensolverFEAST_MKL_2D(double *v, int n_points, int *nq, double *stencil, int n_stencil, double e_min, double e_max, double *E, double *X){

    int i, j;
    int n_entries = 0;
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
    MKL_INT * rows_A = calloc(n_points + 1, sizeof(MKL_INT));
    MKL_INT * cols_A = calloc(max_entries,  sizeof(MKL_INT));
    double  * vals_A = calloc(max_entries,  sizeof(double));
    if(rows_A == NULL || cols_A == NULL || vals_A == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for row, column and/or value arrays");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }

// fill Numerov's A matrix
//  determine the non zero elements and store their positions in rows_A and cols_A
//  and their values in vals_A
    for(i = 0; i < nq[0]; i++){
        for(j = 0; j < nq[1]; j++){
            for(xsh = -n_stencil/2; xsh < n_stencil/2 + 1; xsh++){

                if( (i+xsh > -1) && (i+xsh < nq[0]) ){
                    for(ysh = -n_stencil/2; ysh < n_stencil/2 + 1; ysh++){

                        if( (j+ysh > -1) && ( j+ysh < nq[1]) ){
                            element = (i + xsh)*nq[1] + j+ysh;
                            cols_A[n_entries] = element+1; // wieso +1? weil intel!!

                        // stencil entries have to be divided by 2 to get the right result.
                        //  in three dimensions it should be a division by 4
                            vals_A[n_entries] = stencil[(xsh+n_stencil/2)*n_stencil+ysh+n_stencil/2]/2;

                        // add potential to diagonal element
                            if(xsh == 0 && ysh ==0){
                                vals_A[n_entries] = vals_A[n_entries] + v[i*nq[1]+j];
                            }

                            n_entries ++;
                        }
                    }
                }
            }
      // after inserting all entries in a row the total number of entries is inserted in the CSR format.
        rows_A[i*nq[1] + j + 1] = n_entries + 1;
        }
    }
    rows_A[0] = 1;


// Start eigenstate calculation
    double * res = calloc(n_points, sizeof(double));    // Residual
    if(res == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation of residual array");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }

    char           UPLO     = 'F';
    const MKL_INT  N        = n_points;
    MKL_INT        fpm[128] = { 0 };        // Array to pass parameters to Intel MKL Extended Eigensolvers
    double         epsout   = 0.0;          // Relative error on the trace
    MKL_INT        loop     = 0;            // Number of refinement loop
    MKL_INT        M0       = n_points/2;   // Initial guess for subspace dimension to be used
    MKL_INT        n_out    = n_points/2;
    MKL_INT        info     = 0;            // Errors

    feastinit(fpm); // Initialize Eigensolver parameter array

    dfeast_scsrev(
        &UPLO,      // IN: UPLO = 'F', stores the full matrix
        &N,         // IN: Size of the problem
        vals_A,     // IN: CSR matrix A, values of non-zero elements
        rows_A,     // IN: CSR matrix A, index of the first non-zero element in row
        cols_A,     // IN: CSR matrix A, columns indices for each non-zero element
        fpm,        // IN/OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers
        &epsout,    // OUT: Relative error of on the trace
        &loop,      // OUT: Contains the number of refinement loop executed
        &e_min,     // IN: Lower bound of search interval
        &e_max,     // IN: Upper bound of search interval
        &M0,        // IN: The initial guess for subspace dimension to be used.
        E,          // OUT: The first M entries of Eigenvalues
        X,          // IN/OUT: The first M entries of Eigenvectors
        &n_out,     // OUT: The total number of eigenvalues found in the interval
        res,        // OUT: The first n_out components contain the relative residual vector
        &info       // OUT: Error code
    );

    // Error output
    if( info != 0 ){
        fprintf(stderr, "\n (-) Error in routine sfeast_scsrev");
        fprintf(stderr, "\n     Return code of ERROR: %d", (int)info);
        fprintf(stderr, "\n     Aborting...\n\n");
        exit((int)info);
    }

    free(res);      res     = NULL;
    free(rows_A);   rows_A  = NULL;
    free(cols_A);   cols_A  = NULL;
    free(vals_A);   vals_A  = NULL;

    return n_out;
}
