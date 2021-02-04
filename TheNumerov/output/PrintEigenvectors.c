
#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "settings.h"

// provided prototypes
void PrintEigenvectors(FILE* fd, settings* prefs, int n_out, int n_points, int* nq, double** q, double* v, double*** mu, double* X, double** dip);


void PrintEigenvectors(FILE* fd, settings* prefs, int n_out, int n_points, int* nq, double** q, double* v, double*** mu, double* X, double** dip){

// print header
    fprintf(fd, "#\n#\n# Potential and Eigenfunctions: %d data-points", n_points);

// output size per dimension
    fprintf(fd, "\n  N");
    for(int i = 0; i < prefs->dimension; ++i){
        fprintf(fd, " %4d", nq[i]);
    }
    fprintf(fd, "\n");


// output key
    fprintf(fd, "\n#");
    for(int i = 0; i < prefs->dimension; ++i){
        fprintf(fd, "\t   q[%d]              ", i);
    }
    fprintf(fd, "\t   v(q)              ");

  // given dipole moment
    if( prefs->dipole ){
        fprintf(fd, "\t   dip_x             ");
        fprintf(fd, "\t   dip_y             ");
        fprintf(fd, "\t   dip_z             ");
    }

  // if a coriolis file is given print the modified potential
    if( prefs->coriolis_file ){
        fprintf(fd, "\t   v - sum_i(mu_ii)/8");
    }

  // output of wave functions
    for(int i = 0; i < n_out; ++i){
        fprintf(fd, "\t   Psi[%d]           ", i);
    }
    fprintf(fd, "\n");


// output data
    for(int i = 0; i < n_points; ++i){
    // add a newline every time a index jumps (from max to min)
        int k = 1;
        for(int j = (prefs->dimension - 1); j >= 1; --j){
            k *= nq[j];
            if(i%k == 0){
                fprintf(fd, "\n");
            }
        }

    // output coordinates q[j][i] and potential v[i]
        for(int j = 0; j < prefs->dimension; ++j){
            fprintf(fd, "\t% .14le", q[j][i]);
        }

    // if Coriolis file is set first print the potential without/undo the modifications
        if( prefs->coriolis_file ){
            double prefactor = 1.25E19 * hbar*hbar*avogadro*avogadro
                                       * prefs->InvInertia_to_molpergAasq
                                       * prefs->kJpermol_to_oue;
            fprintf(fd, "\t% .14le", v[i] + prefactor * (mu[0][0][i] + mu[1][1][i] + mu[2][2][i]));
        }else{
            fprintf(fd, "\t% .14le", v[i]);
        }

    // output dipole moment components
        if( prefs->dipole ){
            fprintf(fd, "\t% .14le", dip[0][i]);
            fprintf(fd, "\t% .14le", dip[1][i]);
            fprintf(fd, "\t% .14le", dip[2][i]);
        }


    // output potential after addition of (negative) Watson potential term
        if( prefs->coriolis_file ){
            fprintf(fd, "\t% .14le", v[i]);
        }

    // output wave functions
        for(int j = 0; j < n_out; ++j){
            fprintf(fd, "\t% .14le", X[i + j*n_points]);
        }

        fprintf(fd, "\n");
    }
}
