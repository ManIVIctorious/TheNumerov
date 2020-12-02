
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"

// Dependencies
#ifdef HAVE_MKL_INSTALLED
int        SolverFEAST_MKL(settings* prefs, int* nq, double* v, double ekin_param, double* stencil, double** E, double** X, double** q, double dq, double*** mu, double** zeta);
#endif
#ifdef HAVE_ARMA_INSTALLED
int SolverARPACK_Armadillo(settings* prefs, int* nq, double* v, double ekin_param, double* stencil, double** E, double** X, double** q, double dq, double*** mu, double** zeta);
#endif

// provided prototypes
int MetaEigensolver(settings* prefs, int* nq, double* v, double ekin_param, double* stencil, double* *E, double* *X, double** q, double dq, double*** mu, double** zeta);


int MetaEigensolver(settings* prefs, int* nq, double* v, double ekin_param, double* stencil, double* *E, double* *X, double** q, double dq, double*** mu, double** zeta){

    int n_out;

    switch( prefs->Eigensolver ){

#ifdef HAVE_MKL_INSTALLED
        case 1:
            n_out =        SolverFEAST_MKL(prefs, nq, v, ekin_param, stencil, E, X, q, dq, mu, zeta);
            break;
#endif

#ifdef HAVE_ARMA_INSTALLED
        case 2:
            n_out = SolverARPACK_Armadillo(prefs, nq, v, ekin_param, stencil, E, X, q, dq, mu, zeta);
            break;
#endif

        default:
            fprintf(stderr,
                "\n (-) Error: The required Eigensolver is not available."
                "\n     Please make sure it exists, also check your compile time flags."
                "\n     Aborting...\n\n"
            );
            exit(EXIT_FAILURE);
    }

    return n_out;
}
