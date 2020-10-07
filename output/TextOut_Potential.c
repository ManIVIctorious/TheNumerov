
#include <stdio.h>
 
// Dependencies
double Integrate(int dimension, int* nq, double dx, double* integrand);

// provided prototypes
int TextOut_Potential(FILE* fd, int dimensionality, int n_out, int n_points, int* nq, double* integrand, double dq, double* X, double* v);


// Potential energy output
//  calculate int Psi_i*V*Psi_j dτ (i.e. <X[i]|V|X[j]>)
int TextOut_Potential(FILE* fd, int dimensionality, int n_out, int n_points, int* nq, double* integrand, double dq, double* X, double* v){

    int i, j, k;

// print header and index line
    fprintf(fd, "\n# Potential Energy:\n#\n#");
    for(i = 0; i < n_out; ++i){
        fprintf(fd, "       %7d", i);
    }
    for(i = 0; i < n_out; ++i){
        fprintf(fd, "\n# %3d", i);

        for(j = 0; j < (i+1); ++j){
            for(k = 0; k < n_points; ++k){
            // generate integrand
                integrand[k] = X[k + i*n_points]*X[k + j*n_points] * v[k];
            }
            fprintf(fd, "  % 12.5e", Integrate(dimensionality, nq, dq, integrand));
        }
    }
    fprintf(fd, "\n#\n#");

    return 0;
}
