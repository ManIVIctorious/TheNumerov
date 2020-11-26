
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <mkl_solvers_ee.h>
#include "settings.h"

// Dependencies (with rotation)
int MKL_FillAMatrix1D(double* v, int* nq, double ekin_param, double* stencil, int n_stencil, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);
int MKL_FillAMatrix2D(settings* prefs, int* nq, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);
int MKL_FillAMatrix3D(settings* prefs, int* nq, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);

// provided prototypes
int SolverFEAST_MKL(settings* prefs, int* nq, double* v, double ekin_param, double* stencil, double** E, double** X, double** q, double dq, double*** mu, double** zeta);


int SolverFEAST_MKL(settings* prefs, int* nq, double* v, double ekin_param, double* stencil, double** E, double** X, double** q, double dq, double*** mu, double** zeta){

// get number of points
    int n_points = 1;
    for(int i = 0; i < prefs->dimension; ++i){
        n_points *= nq[i];
    }

// declare arrays used by both, the fill and the evaluation function
    MKL_INT * rows_A = NULL;
    MKL_INT * cols_A = NULL;
    double  * vals_A = NULL;

// depending on the dimensionality use the appropriate filling routine
    switch( prefs->dimension ){

        case 1:
        // one dimensional filling routine, the rotation terms are already set in main (Watson term)
            MKL_FillAMatrix1D(v, nq, ekin_param, stencil, prefs->n_stencil, &rows_A, &cols_A, &vals_A);
            break;

        case 2:
        // two dimensional filling routine, the second term of the Watson Hamiltonian is set in the filling
        //  routine if a Coriolis file is given, the third term is already set in main (Watson term)
            MKL_FillAMatrix2D(prefs, nq, v, ekin_param, stencil, q, dq, mu, zeta, &rows_A, &cols_A, &vals_A);
            break;

        case 3:
        // three dimensional filling routine:
        //  At the moment only the basic Hamiltonian is supported.
            if( prefs->coriolis_file ){
                fprintf(stderr, "At the moment only the basic Hamiltonian is implemented for this"
                                "\n3D problem. Nevertheless, the third term of the Watson Hamiltonian"
                                "\nis already set in main() => be prepared for some wrong results!\n\n"
                       );
            }
            MKL_FillAMatrix3D(prefs, nq, v, ekin_param, stencil, q, dq, mu, zeta, &rows_A, &cols_A, &vals_A);
            break;

        default:
            fprintf(stderr,
                "\n (-) Error: The requested MKL FEAST filling routine"
                "\n     for a dimensionality of \"%d\" is not implemented"
                "\n     Please check your input. Aborting...\n\n"
                ,prefs->dimension
            );
            exit(EXIT_FAILURE);
    }

// allocate memory for eigenvalues E and eigenvectors X
    (*E) = malloc(n_points          * sizeof(double));
    (*X) = malloc(n_points*n_points * sizeof(double));
    if( (*E) == NULL ){ perror("Eigenvalues" ); exit(errno); }
    if( (*X) == NULL ){ perror("Eigenvectors"); exit(errno); }

// Start eigenstate calculation
    double * res = calloc(n_points, sizeof(double));    // Residual
    if( res == NULL ){ perror("SolverFEAST_MKL \"residue\""); exit(errno); }

    char           UPLO     = 'F';
    const MKL_INT  N        = (MKL_INT)n_points;
    MKL_INT        fpm[128] = { 0 };                // Array to pass parameters to Intel MKL Extended Eigensolvers
    double         epsout   = 0.0;                  // Relative error on the trace
    MKL_INT        loop     = 0;                    // Number of refinement loop
    MKL_INT        M0       = (MKL_INT)n_points/2;  // Initial guess for subspace dimension to be used
    MKL_INT        n_out    = (MKL_INT)n_points/2;
    MKL_INT        info     = 0;                    // Errors

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
        &prefs->e_min,  // IN: Lower bound of search interval
        &prefs->e_max,  // IN: Upper bound of search interval
        &M0,            // IN: The initial guess for subspace dimension to be used.
        (*E),           // OUT: The first M entries of Eigenvalues
        (*X),           // IN/OUT: The first M entries of Eigenvectors
        &n_out,         // OUT: The total number of eigenvalues found in the interval
        res,            // OUT: The first n_out components contain the relative residual vector
        &info           // OUT: Error code
    );

// Error output
    if( info ){
        fprintf(stderr,
            "\n (-) Error, \"sfeast_scsrev\" returned code of error: %d"
            "\n     Aborting...\n\n"
            , (int)info
        );
        exit( (int)info );
    }

// reallocate to actual required memory amount
    (*E) = realloc( (*E), n_out          * sizeof(double));
    (*X) = realloc( (*X), n_out*n_points * sizeof(double));
    if( (*E) == NULL ){ perror("Eigenvalues" ); exit(errno); }
    if( (*X) == NULL ){ perror("Eigenvectors"); exit(errno); }


// Garbage Collection
    free(res);      res     = NULL;
    free(rows_A);   rows_A  = NULL;
    free(cols_A);   cols_A  = NULL;
    free(vals_A);   vals_A  = NULL;

    return n_out;
}
