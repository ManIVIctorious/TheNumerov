
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

#include "typedefinitions.h"
#include "constants.h"

// Dependencies
double Integrate(int dimension, int* nq, double dx, double* integrand);

// provided prototypes
int    TextOut_Frequencies(FILE* fd, settings prefs, int n_out, double* E);
int TextOut_Orthonormality(FILE* fd, settings prefs, int n_out, int n_points, int* nq, double* integrand, double dq, double* X);
int      TextOut_Potential(FILE* fd, settings prefs, int n_out, int n_points, int* nq, double* integrand, double dq, double* X, double* v);
int       TextOut_EKinetic(FILE* fd, settings prefs, int n_out, int n_points, int* nq, double* integrand, double dq, double* X, double* stencil, double ekin_param);
int         TextOut_Dipole(FILE* fd, settings prefs, int n_out, int n_points, int* nq, double* integrand, double dq, double* E, double* X, double** dip);
int   TextOut_Eigenvectors(FILE* fd, settings prefs, int n_out, int n_points, int* nq, double** q, double* v, double*** mu, double* X, double** dip);


// Output eigenvalues and calculate Frequencies
int TextOut_Frequencies(FILE* fd, settings prefs, int n_out, double* E){

    int i, j;
    const double kJpermolToWavenumber = 10.0 / (avogadro*planck*lightspeed); // cm^-1 / (kJ/mol)

// output eigenvalues
    fprintf(fd, "# Eigenvalues:");
    for(i = 0; i < n_out; ++i){
        fprintf(fd, " %24.16lf", E[i]);
    }

// and output frequencies
    fprintf(fd, "\n#\n# Frequencies:\n#\n#");
    for(i = 0; i < (n_out - 1); ++i){
        fprintf(fd, "       %7d", i);
    }
    for(i = 1; i < n_out; i++){
        fprintf(fd, "\n# %3d", i);

        for(j = 0; j < i; j++){
            fprintf(fd, "  % 12.5e", (E[i] - E[j]) * kJpermolToWavenumber / prefs.ekin_factor);
        }
    }
    fprintf(fd, "\n#\n#");

    return 0;
}


// check for ortho-normality of evaluated wave functions:
//  calculate int Psi_i*Psi_j dτ (i.e. <X[i]|X[j]>)
int TextOut_Orthonormality(FILE* fd, settings prefs, int n_out, int n_points, int* nq, double* integrand, double dq, double* X){

    int i, j, k;

// print header and index line
    fprintf(fd, "\n# Orthonormality:\n#\n#");
    for(i = 0; i < n_out; ++i){
        fprintf(fd, "       %7d", i);
    }
    for(i = 0; i < n_out; ++i){
        fprintf(fd, "\n# %3d", i);

        for(j = 0; j < (i+1); ++j){
            for(k = 0; k < n_points; ++k){
            // generate integrand
                integrand[k] = X[k + i*n_points]*X[k + j*n_points];
            }
            fprintf(fd, "  % 12.5e", Integrate(prefs.dimension, nq, dq, integrand));
        }
    }
    fprintf(fd, "\n#\n#");

    return 0;
}


// Potential energy output
//  calculate int Psi_i*V*Psi_j dτ (i.e. <X[i]|V|X[j]>)
int TextOut_Potential(FILE* fd, settings prefs, int n_out, int n_points, int* nq, double* integrand, double dq, double* X, double* v){

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
            fprintf(fd, "  % 12.5e", Integrate(prefs.dimension, nq, dq, integrand));
        }
    }
    fprintf(fd, "\n#\n#");

    return 0;
}


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


int TextOut_Eigenvectors(FILE* fd, settings prefs, int n_out, int n_points, int* nq, double** q, double* v, double*** mu, double* X, double** dip){

    int i, j, k;

// print header
    fprintf(fd, "\n# Potential and Eigenfunctions: %d data-points", n_points);

// output size per dimension
    fprintf(fd, "\n  N");
    for(i = 0; i < prefs.dimension; ++i){
        fprintf(fd, " %4d", nq[i]);
    }
    fprintf(fd, "\n");

// output key
    fprintf(fd, "\n#");
    for(i = 0; i < prefs.dimension; ++i){
        fprintf(fd, "\t          q[%d]          ", i);
    }
    fprintf(fd, "\t           v(q)         ");

    if(prefs.dipole){
        fprintf(fd, "\t          dip_x         ");
        fprintf(fd, "\t          dip_y         ");
        fprintf(fd, "\t          dip_z         ");
    }

    if(prefs.coriolis_file_set){
        fprintf(fd, "\tv(q) - sum_i(mu[i][i])/8");
    }

    for(i = 0; i < n_out; ++i){
        fprintf(fd, "\t        Psi[%d]          ", i);
    }
    fprintf(fd, "\n");


// output data
    for(i = 0; i < n_points; ++i){
    // add a newline every time a index jumps (from max to min)
        for(j = (prefs.dimension - 1), k = 1; j >= 1; --j){
            k *= nq[j];
            if(i%k == 0){
                fprintf(fd, "\n");
            }
        }

    // output coordinates q[j][i] and potential v[i]
        for(j = 0; j < prefs.dimension; ++j){
            fprintf(fd, "\t% 24.16lf", q[j][i]);
        }
        if(prefs.coriolis_file_set){
            fprintf(fd, "\t% 24.16lf", v[i] + ((mu[0][0][i] + mu[1][1][i] + mu[2][2][i])/8.0 * (prefs.mu_factor * prefs.ekin_factor)));
        }else{
            fprintf(fd, "\t% 24.16lf", v[i]);
        }

    // output dipole moment components
        if(prefs.dipole != 0){
            fprintf(fd, "\t% 24.16lf", dip[0][i]);
            fprintf(fd, "\t% 24.16lf", dip[1][i]);
            fprintf(fd, "\t% 24.16lf", dip[2][i]);
        }


    // output potential after addition of Watson potential term
        if(prefs.coriolis_file_set){
            fprintf(fd, "\t% 24.16lf", v[i]);
        }

    // output wave functions
        for(j = 0; j < n_out; ++j){
            fprintf(fd, "\t% 24.16lf", X[i + j*n_points]);
        }

        fprintf(fd, "\n");
    }

    return 0;
}


int TextOut_EKinetic(FILE* fd, settings prefs, int n_out, int n_points, int* nq, double* integrand, double dq, double* X, double* stencil, double ekin_param){

    int i, j, k, l;
    int xsh, ysh;
    int element;

// print header and index line
    fprintf(fd, "\n# Kinetic Energy:\n#\n#");
    for(i = 0; i < n_out; ++i){
        fprintf(fd, "       %7d", i);
    }

// integration block, calculate <X|mu|X> for all pairs of X[i] and X[j]
//  => i, j are the eigenvector indices, k loops over all entries
    for(i = 0; i < n_out; ++i){
        fprintf(fd, "\n# %3d", i);

        for(j = 0; j < (i+1); ++j){

            for(k = 0; k < nq[0]; k++){
                for(l = 0; l < nq[1]; l++){
                    integrand[k*nq[1] + l] = 0;
//---------------------------------------------------------------------------------------------------------
                    for(xsh = -prefs.n_stencil/2; xsh < (prefs.n_stencil/2 + 1); xsh++){

                        if( (k+xsh > -1) && (k+xsh < nq[0]) ){
                            for(ysh = -prefs.n_stencil/2; ysh < (prefs.n_stencil/2 + 1); ysh++){

                                if( (l+ysh > -1) && (l+ysh < nq[1]) ){
                                    element = (k + xsh)*nq[1] + l + ysh;

                                // integrand has to be divided by d^2,
                                //  but the division is already set in the "ekin_param" parameter
                                    integrand[k*nq[1] + l] += X[element + i*n_points] * ekin_param * stencil[(xsh + prefs.n_stencil/2)*prefs.n_stencil + ysh + prefs.n_stencil/2]/2;
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
    fprintf(fd, "\n#\n#");

    return 0;
}
