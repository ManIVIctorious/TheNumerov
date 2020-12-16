
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#include "constants.h"
#include "settings.h"

// Dependencies
int FirstDerivative (int n_stencil, double*  first_derivative);
int SecondDerivative(int n_stencil, double* second_derivative);

// Provided prototypes
void   init_watson_2d(settings* prefs);
double exec_watson_2d(double*** mu, double** zeta, int* nq, double dq, double** q, int i, int j, int xsidx, int ysidx);
void   free_watson_2d(void);

// The rotational terms of the molecular Hamiltionian require a first-
//  and second-derivative of the inverse moment of inertia tensor mu
//  Their respective finite difference stencils are stored in the
//  double arrays fst_deriv and sec_deriv.
//  Additionally, the auxiliary array "donothing" ensures the one-dimensional
//  derivatives are only applied to one row/column/etc.
static double * donothing = NULL;
static double * fst_deriv = NULL;
static double * sec_deriv = NULL;

static double conversion_factor = 0.0;

void init_watson_2d(settings* prefs){
//{{{
// allocate memory for derivative stencils
    donothing = calloc(prefs->n_stencil,  sizeof(double));
    fst_deriv = malloc(prefs->n_stencil * sizeof(double));
    sec_deriv = malloc(prefs->n_stencil * sizeof(double));

    if(donothing == NULL){ perror("2D MKL \"donothing\""); exit(errno); }
    if(fst_deriv == NULL){ perror("2D MKL \"fst_deriv\""); exit(errno); }
    if(sec_deriv == NULL){ perror("2D MKL \"sec_deriv\""); exit(errno); }

// set central point of "donothing" to 1.0
    donothing[prefs->n_stencil/2] = 1.0;

// fill first and second derivative stencils
    int errone = FirstDerivative (prefs->n_stencil, fst_deriv);
    int errtwo = SecondDerivative(prefs->n_stencil, sec_deriv);
    if( errone || errtwo ){
        fprintf(stderr,
            "\n (-) Error initialising derivative stencil parameters."
            "\n     Aborting..."
            "\n\n"
        );
        exit(EXIT_FAILURE);
    }

// Pre-calculate conversion factor from Watson-term to output unit of energy
/*--------------------------------------------------------------------------
 *  [zeta]                              1
 *  [mu]                                g/mol/angstrom^2 (per default)
 *  [prefs->InvInertia_to_molpergAasq]  kJ/mol / [mu]
 *  [prefs->kJpermol_to_oue]            (output unit of energy) / (kJ/mol)
 */
    conversion_factor = 1.0E20 * hbar*hbar*avogadro*avogadro //   kJ/mol . g.angstrom^2/mol
                      * prefs->InvInertia_to_molpergAasq     // * mol/(g.angstrom^2) / [mu]
                      * prefs->kJpermol_to_oue;              // * oue / (kJ/mol)
}//}}}

double exec_watson_2d(double*** mu, double** zeta, int* nq, double dq, double** q, int i, int j, int xsidx, int ysidx){
//{{{

// Memory layout of mu and zeta:
// mu[3][3][index]      is a 3x3 matrix containing the [x;y;z] * [x y z] terms of the
//                      reciprocal moment of inertia tensor for each configuration
// zeta[3][D*(D-1)/2]   contains all D*(D-1)/2 combinations of normal modes in {x,y,z}

// calculate pre-factor
    double prefactor = 0.0;
    for(int n = 0; n < 3; ++n){
        for(int m = 0; m < 3; ++m){
            prefactor -= zeta[n][0]*zeta[m][0] * mu[n][m][i*nq[1] + j];
        }
    }

    prefactor *= conversion_factor;

    double watson2d = prefactor * (
                q[0][i*nq[1] + j] *                     fst_deriv[xsidx] * donothing[ysidx] / dq
              + q[1][i*nq[1] + j] *                     donothing[xsidx] * fst_deriv[ysidx] / dq
              - q[0][i*nq[1] + j] * q[0][i*nq[1] + j] * sec_deriv[xsidx] * donothing[ysidx] / dq / dq
              - q[1][i*nq[1] + j] * q[1][i*nq[1] + j] * donothing[xsidx] * sec_deriv[ysidx] / dq / dq
        + 2.0 * q[0][i*nq[1] + j] * q[1][i*nq[1] + j] * fst_deriv[xsidx] * fst_deriv[ysidx] / dq / dq
    );

    return watson2d;
}//}}}

void free_watson_2d(void){
//{{{

    free(donothing); donothing = NULL;
    free(fst_deriv); fst_deriv = NULL;
    free(sec_deriv); sec_deriv = NULL;

}//}}}
