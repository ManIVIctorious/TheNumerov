
#include <stdio.h>
#include <stdlib.h>
#include "typedefinitions.h"

// Dependencies
#ifdef HAVE_MKL_INSTALLED
int SolverFEAST_MKL(settings prefs, int *nq, double *v, double ekin_param, double *stencil, double *E, double *X, double **q, double dq, double ***mu, double ***zeta);
#endif

// Offered prototypes
int MetaEigensolver(settings prefs, int *nq, double *v, double ekin_param, double *stencil, double *E, double *X, double **q, double dq, double ***mu, double ***zeta);


int MetaEigensolver(settings prefs, int *nq, double *v, double ekin_param, double *stencil, double *E, double *X, double **q, double dq, double ***mu, double ***zeta){

    int n_out;

    switch(prefs.Eigensolver){

#ifdef HAVE_MKL_INSTALLED
        case 1:
            n_out = SolverFEAST_MKL(prefs, nq, v, ekin_param, stencil, E, X, q, dq, mu, zeta);
            break;
#endif

        default:
            fprintf(stderr,
                "\n (-) Error: The required Eigensolver is not available."
                "\n     Please make sure it exists, also check your compile time flags."
                "\n     Aborting...\n\n"
            );
            exit(-1);
    }

    return n_out;
}
