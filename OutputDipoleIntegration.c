
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedefinitions.h"
#include "constants.h"

// Dependencies
double Integrate(int dimension, int* nq, double dx, double* integrand);

// provided prototypes
int OutputDipoleIntegration(settings prefs, int* nq, int n_out, double dq, double** dip, double* E, double* X, FILE* fd);

int OutputDipoleIntegration(settings prefs, int* nq, int n_out, double dq, double** dip, double* E, double* X, FILE* fd){

    int i, j, k, m;
    int n_points, element;

    double    integral;
    double *  integrand = NULL; // array to be integrated, freed
    double ** ts_dip    = NULL; // transition dipole moment integral values, freed

// calculate number of points
    for(i = 0, n_points = 1; i < prefs.dimension; ++i){
        n_points *= nq[i];
    }

// allocate memory for arrays
    integrand = malloc(n_points * sizeof(double));
    if(integrand == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation of dipole integrand");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }

// To represent the lower triangle a number of <n_out * (n_out-1) / 2> data points is needed,
//  but due to limitations introduced by integer division the term was increased to n_out * n_out / 2
    ts_dip = malloc(3 * sizeof(double*));
    for(m = 0; m < 3; ++m){
        ts_dip[m] = malloc(n_out*n_out/2 * sizeof(double));
    }
    if(ts_dip == NULL || ts_dip[0] == NULL || ts_dip[1] == NULL || ts_dip[2] == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for transition dipole moment array");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }



// Evaluate the {x,y,z}-components of the IR-intensity, by calculation of
//  <X[i]|\mu_m|X[j]>, for m ∈ {x,y,z}  (i.e. int Psi_i * mu_m * Psi_j dxdy)
// Here m loops over {x,y,z}
    for(m = 0; m < 3; ++m){

    // print header and index line
        fprintf(fd, "\n# Dipole - %c-component:\n#\n#", "xyz"[m]);
        for(i = 0; i < n_out; ++i){
            fprintf(fd, "       %7d", i);
        }


    // integration block, calculate <X|mu|X> for all pairs of X[i] and X[j]
    //  => i, j represent the respective eigenvectors
    //     k  loops over the entries within
        for(i = 0, element = 0; i < n_out; ++i){
            fprintf(fd, "\n# %3d", i);
            for(j = 0; j <= i; ++j){

            // generate integrand and perform integration
                for(k = 0; k < n_points; ++k){
                    integrand[k] = X[k + i*n_points]*X[k + j*n_points] * dip[m][k];
                }
                integral = Integrate(prefs.dimension, nq, dq, integrand);

            // output integral
                fprintf(fd, "  % 12.5e", integral);

            // ts_dip is needed for the calculation of the oscillator strength
            //  since the oscillator strength is only defined between two different
            //  states the diagonal elements have to be excluded
                if(i != j){
                    ts_dip[m][element++] = integral;
                }

            }
        }
        fprintf(fd, "\n#\n#");

    }
// free integrand
    free(integrand); integrand = NULL;



// Evaluate the oscillator strength according to the following equation
//
//  f_ij  =   (4*pi * elmass) / (3 * hbar * elcharge^2)         // 1 / (A^2 s m^2)
//          * || sum_k { <Psi_i|mu_k|Psi_j> } ||^2 * DipToAsm^2 // * A^2 s^2 m^2
//          * (E_j - E_i) / (planck * avogadro / 1000)          // * 1 / s
//
// [f_ij] = 1

// print header and index line
    fprintf(fd, "\n# Oscillator strength:\n#\n#");
    for(i = 0; i < (n_out - 1); ++i){
        fprintf(fd,"       %7d", i);
    }

// calculate and immediately print oscillator strength f_ij
    for(i = 1, element = 0; i < n_out; ++i){
        fprintf(fd, "\n# %3d", i);

        for(j = 0; j < i; ++j){

            fprintf(fd, "  % 12.5e",

                  (4*M_PI * elmass) / (3 * planck/(2*M_PI) * elementarycharge*elementarycharge)
                * prefs.DipToAsm*prefs.DipToAsm
                * (
                      ts_dip[0][element]*ts_dip[0][element]
                    + ts_dip[1][element]*ts_dip[1][element]
                    + ts_dip[2][element]*ts_dip[2][element]
                  )
                * ((E[i] - E[j]) / prefs.ekin_factor) / (planck * avogadro / 1000)

            );

            ++element;
        }
    }
    fprintf(fd, "\n#\n#");

// free ts_dip 2D array
    for(m = 0; m < 3; ++m){
        free(ts_dip[m]); ts_dip[m] = NULL;
    }
    free(ts_dip); ts_dip = NULL;

    return 0;
}
