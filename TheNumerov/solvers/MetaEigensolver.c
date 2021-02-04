
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"

// Dependencies
int        SolverFEAST_MKL(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, double** E, double** X, double** q, double dq, double*** mu, double** zeta);
int SolverARPACK_Armadillo(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, double** E, double** X, double** q, double dq, double*** mu, double** zeta);

// provided prototypes
int MetaEigensolver(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, double* *E, double* *X, double** q, double dq, double*** mu, double** zeta);


int MetaEigensolver(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, double* *E, double* *X, double** q, double dq, double*** mu, double** zeta){

    int n_out;

    switch( prefs->Eigensolver ){
    // Eigensolver bit-mask
    // 1    Intel MKL FEAST
    // 2    Armadillo ARPACK
    // 4    Not implemented yet

        case 1:
            n_out =        SolverFEAST_MKL(prefs, nq, v, ekin_to_oue, stencil, E, X, q, dq, mu, zeta);
            break;

        case 2:
            n_out = SolverARPACK_Armadillo(prefs, nq, v, ekin_to_oue, stencil, E, X, q, dq, mu, zeta);
            break;

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
