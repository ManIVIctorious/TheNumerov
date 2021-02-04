
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "settings.h"

// provided prototypes
void PrintEKin(FILE* fd, settings *prefs, int n_out, int n_points, int* nq, double dq, double* X, double* stencil, double ekin_to_oue);

// Dependencies
double Integrate(int dimension, int* nq, double dx, double* integrand);


void PrintEKin(FILE* fd, settings *prefs, int n_out, int n_points, int* nq, double dq, double* X, double* stencil, double ekin_to_oue){

// allocate memory for integrand
    double * integrand = malloc(n_points * sizeof(double));
    if(integrand == NULL){ perror("Integrand"); exit(errno); }

// print header and index line
    fprintf(fd, "#\n#\n# Kinetic Energy:\n#\n#");
    for(int i = 0; i < n_out; ++i){
        fprintf(fd, "       %7d", i);
    }

// integration block, calculate <X|mu|X> for all pairs of X[i] and X[j]
//  => i, j are the eigenvector indices, k loops over all entries
    for(int i = 0; i < n_out; ++i){
        fprintf(fd, "\n# %3d", i);

        for(int j = 0; j < (i+1); ++j){

            for(int k = 0; k < nq[0]; k++){
                for(int l = 0; l < nq[1]; l++){
                    integrand[k*nq[1] + l] = 0;
//---------------------------------------------------------------------------------------------------------
                    for(int xsh = -prefs->n_stencil/2; xsh < (prefs->n_stencil/2 + 1); xsh++){

                        if( (k+xsh > -1) && (k+xsh < nq[0]) ){
                            for(int ysh = -prefs->n_stencil/2; ysh < (prefs->n_stencil/2 + 1); ysh++){

                                if( (l+ysh > -1) && (l+ysh < nq[1]) ){
                                    int element = (k + xsh)*nq[1] + l + ysh;

                                // integrand has to be divided by d^2,
                                //  but the division is already set in the "ekin_to_oue" parameter
                                    integrand[k*nq[1] + l] += X[element + i*n_points] * ekin_to_oue * stencil[(xsh + prefs->n_stencil/2)*prefs->n_stencil + ysh + prefs->n_stencil/2]/2;
                                }
                            }
                        }
                    }
//---------------------------------------------------------------------------------------------------------
                    integrand[k*nq[1] + l] *= X[(k*nq[1] + l) + j*n_points];
                }
            }
            fprintf(fd, "  % 12.5e", Integrate(2, nq, dq, integrand));
        }
    }
    fprintf(fd, "\n");

    free(integrand); integrand = NULL;
}
