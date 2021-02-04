
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

#include "settings.h"
#include "constants.h"

// provided prototypes
void PrintDipole(FILE* fd, settings *prefs, int n_out, int n_points, int* nq, double dq, double* E, double* X, double** dip);

// Dependencies
double Integrate(int dimension, int* nq, double dx, double* integrand);


// perform dipole integration to get intensities and oscillator strengths
void PrintDipole(FILE* fd, settings *prefs, int n_out, int n_points, int* nq, double dq, double* E, double* X, double** dip){

// allocate memory for the transition moment integrals in x, y and z direction
    double ** ts = malloc(3 * sizeof(double*));
    if(ts == NULL){ perror("Transition dipole moment (ts)"); exit(errno); }

    for(int m = 0; m < 3; ++m){ ts[m] = malloc( (n_out*(n_out+1))/2 * sizeof(double)); }
    if(ts[0] == NULL){ perror("Transition dipole moment (ts[0])"); exit(errno); }
    if(ts[1] == NULL){ perror("Transition dipole moment (ts[1])"); exit(errno); }
    if(ts[2] == NULL){ perror("Transition dipole moment (ts[2])"); exit(errno); }

// allocate memory for integrand
    double * integrand = malloc(n_points * sizeof(double));
    if(integrand == NULL){ perror("Integrand"); exit(errno); }



// Evaluate the {x,y,z}-components of the IR-intensity
//--------------------------------------------------------------------------------
    for(int m = 0; m < 3; ++m){

    // print header and index line
        fprintf(fd, "#\n#\n# Dipole - %c-component in input units:\n#\n#", "xyz"[m]);
        for(int i = 0; i < n_out; ++i){
            fprintf(fd, "       %7d", i);
        }

    // calculate int Psi_i * dip_m * Psi_j dτ (i.e. <X[i]|mu[m]|X[j]>)
    // entry counter for ts
        int count = 0;
        for(int i = 0; i < n_out; ++i){
            fprintf(fd, "\n# %3d", i);

        // first entry of i-th eigenvector
            int istart = i*n_points;

            for(int j = 0; j <= i; ++j){

            // first entry of j-th eigenvector
                int jstart = j*n_points;

            // prepare and perform integration
                for(int k = 0; k < n_points; ++k){
                    integrand[k] = X[istart + k]*X[jstart + k] * dip[m][k];
                }
                ts[m][count] = Integrate(prefs->dimension, nq, dq, integrand);
                fprintf(fd, "  % 12.5e", ts[m][count]);
                ++count;
            }
        }
        fprintf(fd, "\n");
    }

    free(integrand); integrand = NULL;



// Evaluate and print oscillator strengths f_ij ( [f_ij] = 1 ):
//--------------------------------------------------------------------------------
/*
 *  f_ij = (4*pi * elmass) / (3 * hbar * elcharge^2)            // 1 / (A^2 s m^2)
 *       * || sum_m { <Psi_i|mu_m|Psi_j> } ||^2 * dip_to_Asm^2  // * A^2 s^2 m^2
 *       * (E_j - E_i) / (planck * avogadro / 1000)             // * 1 / s
 */
    const double f_perDipsq = (8.0 * M_PI*M_PI * elmass) / (3.0 * planck * elcharge*elcharge);
    const double E_to_freq  = 1000.0 / (avogadro * planck) / prefs->kJpermol_to_oue;
    const double prefactor  = f_perDipsq * E_to_freq * prefs->dip_to_Asm*prefs->dip_to_Asm;

// print header and index line
    fprintf(fd, "#\n#\n# Oscillator strength:\n#\n#");
    for(int i = 0; i < (n_out - 1); ++i){
        fprintf(fd,"       %7d", i);
    }

// calculate and immediately print oscillator strength f_ij
// entry counter for ts
    int count = 0;
    for(int i = 1; i < n_out; ++i){
        fprintf(fd, "\n# %3d", i);

    // jump over entries of the main diagonal
        count += 1;

        for(int j = 0; j < i; ++j){

            fprintf(fd, "  % 12.5e",

                prefactor * (E[i] - E[j]) * (
                      ts[0][count]*ts[0][count] // <psi|mu_x|psi>
                    + ts[1][count]*ts[1][count] // <psi|mu_y|psi>
                    + ts[2][count]*ts[2][count] // <psi|mu_z|psi>
                )

            );
            ++count;
        }
    }
    fprintf(fd, "\n");


// free ts 2D array
    for(int m = 0; m < 3; ++m){ free(ts[m]); ts[m] = NULL; }
    free(ts); ts = NULL;
}
