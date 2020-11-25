
#define ARMA_BLAS_LONG
#include <armadillo>
#include <stdio.h>
#include "settings.h"

// Dependencies
arma::sp_mat FillArmadillo_1D(int* nq, double* v, double ekin_param, double* stencil, int n_stencil);
arma::sp_mat FillArmadillo_2D(settings* prefs, int* nq, int n_points, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta);
arma::sp_mat FillArmadillo_3D(settings* prefs, int* nq, int n_points, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta);

// Provided Prototypes
extern "C" int SolverARPACK_Armadillo(settings* prefs, int* nq, double* v, double ekin_param, double* stencil, double** E, double** X, double** q, double dq, double*** mu, double** zeta);


extern "C"{

    int SolverARPACK_Armadillo(settings* prefs, int* nq, double* v, double ekin_param, double* stencil, double** E, double** X, double** q, double dq, double*** mu, double** zeta){

    // get number of points
        int n_points = 1;
        for(int i = 0; i < prefs->dimension; ++i){
            n_points *= nq[i];
        }

        arma::sp_mat A;

    // depending on the dimensionality use the appropriate filling routine
        switch( prefs->dimension ){

            case 1:
            // one dimensional filling routine, the rotation terms are already set in main (Watson term)
                A = FillArmadillo_1D(nq, v, ekin_param, stencil, prefs->n_stencil);
                break;

            case 2:
            // two dimensional filling routine, the second term of the Watson Hamiltonian is set in the filling
            //  routine if a Coriolis file is given, the third term is already set in main (Watson term)
                A = FillArmadillo_2D(prefs, nq, n_points, v, ekin_param, stencil, q, dq, mu, zeta);
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
                A = FillArmadillo_3D(prefs, nq, n_points, v, ekin_param, stencil, q, dq, mu, zeta);
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
        bool success = false;
        arma::mat eigvec;
        arma::vec eigval;

        success = eigs_sym(eigval, eigvec, A, prefs->n_out, "sm");
        if( !success ){
            fprintf(stderr,
                "\n (-) Error: Failed eigen decomposition.\n"
                "\n     Aborting...\n\n"
            );
            exit(EXIT_FAILURE);
        }


    // allocate memory for eigenvalues E and eigenvectors X
        (*E) = (double*) malloc(eigval.n_elem                 * sizeof(double));
        (*X) = (double*) malloc(eigvec.n_rows * eigvec.n_cols * sizeof(double));
        if((*E) == NULL){ perror("Eigenvalues" ); exit(errno); }
        if((*X) == NULL){ perror("Eigenvectors"); exit(errno); }
    // fill eigenvalue E and eigenvector arrays X
        for(unsigned int i = 0; i < eigvec.n_cols; ++i){
            (*E)[i] = eigval[i];

            for(unsigned int j = 0; j < eigvec.n_rows; ++j){
                (*X)[i*n_points + j] = eigvec(j,i);
            }
        }

        return eigval.n_elem;
    }
}
