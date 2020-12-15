
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

// provided prototypes
void PrintOrthonormality(FILE* fd, int dimensionality, int n_out, int n_points, int* nq, double dq, double* X);

// Dependencies
double Integrate(int dimension, int* nq, double dx, double* integrand);


// check for ortho-normality of evaluated wave functions:
void PrintOrthonormality(FILE* fd, int dimensionality, int n_out, int n_points, int* nq, double dq, double* X){

// allocate memory for integrand
    double * integrand = malloc(n_points * sizeof(double));
    if(integrand == NULL){ perror("Integrand"); exit(errno); }

// print header and index line
    fprintf(fd, "\n# Orthonormality:\n#\n#");
    for(int i = 0; i < n_out; ++i){
        fprintf(fd, "       %7d", i);
    }

// calculate int Psi_i*Psi_j dτ (i.e. <X[i]|X[j]>)
    for(int i = 0; i < n_out; ++i){
        fprintf(fd, "\n# %3d", i);

    // first entry of i-th eigenvector
        int istart = i*n_points;

        for(int j = 0; j <= i; ++j){

        // first entry of j-th eigenvector
            int jstart = j*n_points;

        // prepare and perform integration
            for(int k = 0; k < n_points; ++k){
                integrand[k] = X[istart + k]*X[jstart + k];
            }
            fprintf(fd, "  % 12.5e", Integrate(dimensionality, nq, dq, integrand));
        }
    }
    fprintf(fd, "\n#\n#");

    free(integrand); integrand = NULL;
}
