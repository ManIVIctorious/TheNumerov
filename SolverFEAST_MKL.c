
#include <stdio.h>
#include <stdlib.h>
#include <mkl_solvers_ee.h>
#include "typedefinitions.h"

// Dependencies (with rotation)
int FillMKL_2D(settings prefs, int* nq, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);
// Dependencies (without rotation)
int FillMKL_1D_norot(double* v, int* nq, double ekin_param, double* stencil, int n_stencil, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);
int FillMKL_2D_norot(double* v, int* nq, double ekin_param, double* stencil, int n_stencil, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);

// provided prototypes
int SolverFEAST_MKL(settings prefs, int* nq, double* v, double ekin_param, double* stencil, double* E, double* X, double** q, double dq, double*** mu, double** zeta);


int SolverFEAST_MKL(settings prefs, int* nq, double* v, double ekin_param, double* stencil, double* E, double* X, double** q, double dq, double*** mu, double** zeta){

    int i;
    int n_points;

// get number of points
    for(i = 0, n_points = 1; i < prefs.dimension; ++i){
        n_points *= nq[i];
    }

// declare arrays used by both, the fill and the evaluation function
    MKL_INT * rows_A = NULL;
    MKL_INT * cols_A = NULL;
    double  * vals_A = NULL;

// depending on the dimensionality use the appropriate filling routine
    switch(prefs.dimension){

        case 1:
            //if(prefs.coriolis_file != NULL){
            //// one dimensional filling routine with rotation enabled
            //    FillMKL_1D(prefs, nq, v, ekin_param, stencil, q, dq, mu, zeta, &rows_A, &cols_A, &vals_A);
            //}else{
            //// one dimensional filling routine without rotation enabled
                  FillMKL_1D_norot(v, nq, ekin_param, stencil, prefs.n_stencil, &rows_A, &cols_A, &vals_A);
            //}
            break;

        case 2:
            if(prefs.coriolis_file != NULL){
            // two dimensional filling routine with rotation enabled
                FillMKL_2D(prefs, nq, v, ekin_param, stencil, q, dq, mu, zeta, &rows_A, &cols_A, &vals_A);
            }else{
            // two dimensional filling routine without rotation enabled
                FillMKL_2D_norot(v, nq, ekin_param, stencil, prefs.n_stencil, &rows_A, &cols_A, &vals_A);
            }
            break;

        default:
            fprintf(stderr,
                "\n (-) Error: For the requested dimensionality of \"%d\""
                "\n     no adequate MKL FEAST fill routine exists."
                "\n     Please check your input. Aborting...\n\n"
                ,prefs.dimension
            );
            exit(-1);
    }

// Start eigenstate calculation
    double * res = calloc(n_points, sizeof(double));    // Residual
    if(res == NULL){
        fprintf(stderr,
            "\n (-) Error in memory allocation of residual array"
            "\n     Aborting..."
            "\n\n"
        );
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
        &UPLO,          // IN: UPLO = 'F', stores the full matrix
        &N,             // IN: Size of the problem
        vals_A,         // IN: CSR matrix A, values of non-zero elements
        rows_A,         // IN: CSR matrix A, index of the first non-zero element in row
        cols_A,         // IN: CSR matrix A, columns indices for each non-zero element
        fpm,            // IN/OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers
        &epsout,        // OUT: Relative error of on the trace
        &loop,          // OUT: Contains the number of refinement loop executed
        &prefs.e_min,   // IN: Lower bound of search interval
        &prefs.e_max,   // IN: Upper bound of search interval
        &M0,            // IN: The initial guess for subspace dimension to be used.
        E,              // OUT: The first M entries of Eigenvalues
        X,              // IN/OUT: The first M entries of Eigenvectors
        &n_out,         // OUT: The total number of eigenvalues found in the interval
        res,            // OUT: The first n_out components contain the relative residual vector
        &info           // OUT: Error code
    );

    // Error output
    if( info != 0 ){
        fprintf(stderr,
            "\n (-) Error in routine sfeast_scsrev"
            "\n     Return code of ERROR: %d"
            "\n     Aborting..."
            "\n\n"
            , (int)info
        );
        exit((int)info);
    }

// free memory of arrays which are not needed anymore
    free(res);      res     = NULL;
    free(rows_A);   rows_A  = NULL;
    free(cols_A);   cols_A  = NULL;
    free(vals_A);   vals_A  = NULL;


    return n_out;
}
