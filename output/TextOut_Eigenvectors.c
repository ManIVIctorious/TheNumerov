
#include <stdio.h>
#include "settings.h"

// Dependencies
double Integrate(int dimension, int* nq, double dx, double* integrand);

// provided prototypes
int TextOut_Eigenvectors(FILE* fd, settings* prefs, int n_out, int n_points, int* nq, double** q, double* v, double*** mu, double* X, double** dip);


int TextOut_Eigenvectors(FILE* fd, settings* prefs, int n_out, int n_points, int* nq, double** q, double* v, double*** mu, double* X, double** dip){

    int i, j, k;

// print header
    fprintf(fd, "\n# Potential and Eigenfunctions: %d data-points", n_points);

// output size per dimension
    fprintf(fd, "\n  N");
    for(i = 0; i < prefs->dimension; ++i){
        fprintf(fd, " %4d", nq[i]);
    }
    fprintf(fd, "\n");

// output key
    fprintf(fd, "\n#");
    for(i = 0; i < prefs->dimension; ++i){
        fprintf(fd, "\t          q[%d]          ", i);
    }
    fprintf(fd, "\t           v(q)         ");

    if(prefs->dipole){
        fprintf(fd, "\t          dip_x         ");
        fprintf(fd, "\t          dip_y         ");
        fprintf(fd, "\t          dip_z         ");
    }

    if( prefs->coriolis_file ){
        fprintf(fd, "\tv(q) - sum_i(mu[i][i])/8");
    }

    for(i = 0; i < n_out; ++i){
        fprintf(fd, "\t        Psi[%d]          ", i);
    }
    fprintf(fd, "\n");


// output data
    for(i = 0; i < n_points; ++i){
    // add a newline every time a index jumps (from max to min)
        for(j = (prefs->dimension - 1), k = 1; j >= 1; --j){
            k *= nq[j];
            if(i%k == 0){
                fprintf(fd, "\n");
            }
        }

    // output coordinates q[j][i] and potential v[i]
        for(j = 0; j < prefs->dimension; ++j){
            fprintf(fd, "\t% 24.16lf", q[j][i]);
        }
        if( prefs->coriolis_file ){
            fprintf(fd, "\t% 24.16lf", v[i] + ((mu[0][0][i] + mu[1][1][i] + mu[2][2][i])/8.0 * (prefs->mu_factor * prefs->ekin_factor)));
        }else{
            fprintf(fd, "\t% 24.16lf", v[i]);
        }

    // output dipole moment components
        if(prefs->dipole != 0){
            fprintf(fd, "\t% 24.16lf", dip[0][i]);
            fprintf(fd, "\t% 24.16lf", dip[1][i]);
            fprintf(fd, "\t% 24.16lf", dip[2][i]);
        }


    // output potential after addition of Watson potential term
        if( prefs->coriolis_file ){
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
