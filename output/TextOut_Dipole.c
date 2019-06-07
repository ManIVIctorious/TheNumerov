
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

#include "typedefinitions.h"
#include "constants.h"
 
// Dependencies
double Integrate(int dimension, int* nq, double dx, double* integrand);

// provided prototypes
int TextOut_Dipole(FILE* fd, settings prefs, int n_out, int n_points, int* nq, double* integrand, double dq, double* E, double* X, double** dip);


// perform dipole integration to get intensities and oscillator strengths
int TextOut_Dipole(FILE* fd, settings prefs, int n_out, int n_points, int* nq, double* integrand, double dq, double* E, double* X, double** dip){

    int i, j, k, m, shift;
    double ** ts_dip = NULL;

// allocate memory for the transition moment integrals in x, y and z direction
    ts_dip = malloc(3 * sizeof(double*));
    if(ts_dip == NULL){ perror("Transition dipole moment (ts_dip)"); exit(errno); }

    for(m = 0; m < 3; ++m){
        ts_dip[m] = malloc( (n_out*(n_out-1) + n_out) * sizeof(double));
    }
    if(ts_dip[0] == NULL){ perror("Transition dipole moment (ts_dip[0])"); exit(errno); }
    if(ts_dip[1] == NULL){ perror("Transition dipole moment (ts_dip[1])"); exit(errno); }
    if(ts_dip[2] == NULL){ perror("Transition dipole moment (ts_dip[2])"); exit(errno); }


// Evaluate the {x,y,z}-components of the IR-intensity, by calculation of
//  <X[i]|\mu_m|X[j]>, for m ∈ {x,y,z}  (i.e. int Psi_i * mu_m * Psi_j dτ)
    for(m = 0; m < 3; ++m){

    // print header and index line
        fprintf(fd, "\n# Dipole - %c-component:\n#\n#", "xyz"[m]);
        for(i = 0; i < n_out; ++i){
            fprintf(fd, "       %7d", i);
        }

    // integration block, calculate <X|mu|X> for all pairs of X[i] and X[j]
    //  => i, j are the eigenvector indices, k loops over all entries
        for(i = 0, shift = 0; i < n_out; ++i){
            fprintf(fd, "\n# %3d", i);

            for(j = 0; j <= i; ++j){
                for(k = 0; k < n_points; ++k){
                // generate integrand
                    integrand[k] = X[i*n_points + k] * X[j*n_points + k] * dip[m][k];
                }
                ts_dip[m][shift+j] = Integrate(prefs.dimension, nq, dq, integrand);
                fprintf(fd, "  % 12.5e", ts_dip[m][shift+j]);
            }
            shift += (i+1);
        }
        fprintf(fd, "\n#\n#");
    }


// Evaluate oscillator strength f_ij ( [f_ij] = 1 ):
//
//  f_ij  =   (4*pi * elmass) / (3 * hbar * elcharge^2)             // 1 / (A^2 s m^2)
//          * || sum_k { <Psi_i|mu_k|Psi_j> } ||^2 * DipToAsm^2     // * A^2 s^2 m^2
//          * (E_j - E_i) / (planck * avogadro / 1000)              // * 1 / s
//
const double f_osc_prefactor = (8000.0/3.0 * M_PI*M_PI * elmass) * prefs.DipToAsm*prefs.DipToAsm
                             / (planck*planck * avogadro * elementarycharge*elementarycharge);

// print header and index line
    fprintf(fd, "\n# Oscillator strength:\n#\n#");
    for(i = 0; i < (n_out - 1); ++i){
        fprintf(fd,"       %7d", i);
    }

// calculate and immediately print oscillator strength f_ij
    for(i = 1, shift = 1; i < n_out; ++i){
        fprintf(fd, "\n# %3d", i);

        for(j = 0; j < i; ++j){

            fprintf(fd, "  % 12.5e",

                f_osc_prefactor * ((E[i] - E[j]) / prefs.ekin_factor) * (
                      ts_dip[0][shift+j]*ts_dip[0][shift+j] // <psi|mu_x|psi>
                    + ts_dip[1][shift+j]*ts_dip[1][shift+j] // <psi|mu_y|psi>
                    + ts_dip[2][shift+j]*ts_dip[2][shift+j] // <psi|mu_z|psi>
                )

            );
        }
        shift += (i+1);
    }
    fprintf(fd, "\n#\n#");


// free ts_dip 2D array
    for(m = 0; m < 3; ++m){
        free(ts_dip[m]); ts_dip[m] = NULL;
    }
    free(ts_dip); ts_dip = NULL;

    return 0;
}
