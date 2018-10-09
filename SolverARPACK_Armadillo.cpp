
#include <armadillo>
#include <stdio.h>
#include "typedefinitions.h"

// Dependencies
arma::sp_mat FillArmadillo_1D(int* nq, double* v, double ekin_param, double* stencil, int n_stencil);
arma::sp_mat FillArmadillo_2D(settings prefs, int* nq, int n_points, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta);

// Provided Prototypes
extern "C" int SolverARPACK_Armadillo(settings prefs, int* nq, double* v, double ekin_param, double* stencil, double* E, double* X, double** q, double dq, double*** mu, double** zeta);


extern "C"{

    int SolverARPACK_Armadillo(settings prefs, int* nq, double* v, double ekin_param, double* stencil, double* E, double* X, double** q, double dq, double*** mu, double** zeta){

        unsigned int i, j;
        int n_points;

    // get number of points
        for(i = 0, n_points = 1; i < (unsigned int)prefs.dimension; ++i){
            n_points *= nq[i];
        }

        arma::sp_mat A;

    // depending on the dimensionality use the appropriate filling routine
        switch(prefs.dimension){

            case 1:
            // one dimensional filling routine, the rotation terms are already set in main (Watson term)
                A = FillArmadillo_1D(nq, v, ekin_param, stencil, prefs.n_stencil);
                break;

            case 2:
            // two dimensional filling routine:
            //  The first and second terms of the Watson Hamiltonian are to be set in the filling routine
            //  when a Coriolis file is given, the third term is already set in main (Watson term)
                A = FillArmadillo_2D(prefs, nq, n_points, v, ekin_param, stencil, q, dq, mu, zeta);
                break;

            default:
                fprintf(stderr,
                    "\n (-) Error: The requested ARPACK Armadillo filling routine"
                    "\n     for a dimensionality of \"%d\" is not implemented."
                    "\n     Please check your input. Aborting...\n\n"
                    ,prefs.dimension
                );
                exit(EXIT_FAILURE);
        }

// start eigenstate calculation
    arma::mat eigvec;
    arma::vec eigval;

// In recent versions of armadillo the ARPACK is directly included
//  The eigendecomposition can either be invoked by the armadillo "high-level" eigs routine
//  or by using the underlying ARPACK implementation directly
//  To use the underlying ARPACK implementation some objects/classes have to be put into scope:
    using arma::newarp::SparseGenMatProd;
    using arma::newarp::SymEigsSolver;
    using arma::newarp::EigsSelect;

// Construct matrix operation object using the wrapper class SparseGenMatProd
    SparseGenMatProd<double> op(A);

// Construct eigen solver object, requesting the largest <prefs.n_out> eigenvalues
    SymEigsSolver< double, EigsSelect::SMALLEST_MAGN, SparseGenMatProd<double> > eigs(op, prefs.n_out, 4*prefs.n_out);

// Initialize and compute
    eigs.init();
    int n_conv = eigs.compute(10000, 1.0e-14);

    if(n_conv > 0){
        eigval = eigs.eigenvalues();
        eigvec = eigs.eigenvectors();
    }else{
        fprintf(stderr,
            "\n (-) Error: Failed eigen decomposition.\n"
            "\n     Aborting...\n\n"
        );
        exit(EXIT_FAILURE);
    }

    // fill eigenvalue E and eigenvector arrays X
        for(i = 0; i < eigvec.n_cols; ++i){
            E[i] = eigval[i];

            for(j = 0; j < eigvec.n_rows; ++j){
                X[i*n_points + j] = eigvec(j,i);
            }
        }

        return eigval.n_elem;
    }
}
