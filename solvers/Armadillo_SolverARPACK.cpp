
#define ARMA_BLAS_LONG
#include <armadillo>
#include <stdio.h>

#include "settings.h"
#include "ArmadilloFillers.h"

// Provided Prototypes
extern "C" int SolverARPACK_Armadillo(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, double** E, double** X, double** q, double dq, double*** mu, double** zeta);


extern "C"{

    int SolverARPACK_Armadillo(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, double** E, double** X, double** q, double dq, double*** mu, double** zeta){

    // get number of points
        int n_points = 1;
        for(int i = 0; i < prefs->dimension; ++i){
            n_points *= nq[i];
        }

        arma::sp_mat A;

    // depending on the dimensionality use the appropriate filling routine
        switch( prefs->dimension ){

          case 1:
            if( prefs->periodic ){ A = FillPeriodicArmadillo_1D(prefs->n_stencil, nq, v, ekin_to_oue, stencil); }
            else{                  A =         FillArmadillo_1D(prefs->n_stencil, nq, v, ekin_to_oue, stencil); }
                break;

            case 2:
                if( prefs->periodic ){ A = FillPeriodicArmadillo_2D(prefs, nq, v, ekin_to_oue, stencil); }
                else{                  A =         FillArmadillo_2D(prefs, nq, v, ekin_to_oue, stencil, q, dq, mu, zeta); }
                break;

            case 3:
                if( prefs->periodic ){ A = FillPeriodicArmadillo_3D(prefs, nq, v, ekin_to_oue, stencil); }
                else{                  A =         FillArmadillo_3D(prefs, nq, v, ekin_to_oue, stencil, q, dq, mu, zeta); }
                break;

            case 4:
                if( prefs->periodic ){ A = FillPeriodicArmadillo_4D(prefs, nq, v, ekin_to_oue, stencil); }
                else{                  A =         FillArmadillo_4D(prefs, nq, v, ekin_to_oue, stencil, q, dq, mu, zeta); }
                break;

            default:
                fprintf(stderr,
                    "\n (-) Error: The requested ARPACK Armadillo filling routine"
                    "\n     for a dimensionality of \"%d\" is not implemented."
                    "\n     Please check your input. Aborting...\n\n"
                    ,prefs->dimension
                );
                exit(EXIT_FAILURE);
        }

    // start eigenstate calculation
    // the application of the Watson Hamiltonian breaks the matrix symmetry for
    // dimensions greater than two hence, a general solver implementation is required
        bool success = false;
        arma::cx_mat eigvec;
        arma::cx_vec eigval;

        success = eigs_gen(eigval, eigvec, A, prefs->n_out, "sm");
        if( !success ){
            fprintf(stderr,
                "\n (-) Error: Failed eigendecomposition."
                "\n     Aborting...\n\n"
            );
            exit(EXIT_FAILURE);
        }

    // allocate memory for eigenvalues E and eigenvectors X
        (*E) = (double*) malloc(eigval.n_elem                 * sizeof(double));
        (*X) = (double*) malloc(eigvec.n_rows * eigvec.n_cols * sizeof(double));
        if((*E) == NULL){ perror("Eigenvalues" ); exit(errno); }
        if((*X) == NULL){ perror("Eigenvectors"); exit(errno); }

    // store the real part of eigenvalues and eigenvectors in the E and X arrays
    // the imaginary part of these should be non-existent or at least negligible
        for(unsigned int i = 0; i < eigvec.n_cols; ++i){
            (*E)[i] = real(eigval[i]);

            for(unsigned int j = 0; j < eigvec.n_rows; ++j){
                (*X)[i*n_points + j] = real(eigvec(j,i));
            }
        }

        return eigval.n_elem;
    }
}
