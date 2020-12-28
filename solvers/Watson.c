
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#include "constants.h"
#include "settings.h"

// Provided prototypes
#include "Watson.h"

// Dependencies
int  FirstDerivative(int n_stencil, double*  first_derivative);
int SecondDerivative(int n_stencil, double* second_derivative);


/* Triplet of init_ exec_ free_
 *
 * init_watson:
 *   allocate memory and set values for:
 *     * fst_deriv: array containing the first  derivative stencil
 *     * sec_deriv: array containing the second derivative stencil
 *     * donothing: auxiliary array only containing a central point, to make the
 *                  application of one-dimensional derivatives more consistent.
 *
 *   pre-calculate conv_factor:
 *     This factor represents the conversion from the inverse moment of inertia
 *     to the output unit of energy utilising Planck's constant.
 *
 * free_watson:
 *   free allocated arrays
 *
 * exec_watson:
 *   calculate first term of Watson molecular Hamiltonian H
 *   The second term is already calculated in the main() function, while the
 *   third and fourth represent the standard, non-rotational Hamiltonian.
 *
 *      H =   hbar^2/2  * sum_ab pi_a * mu_ab * pi_b
 *          - hbar^2/8  * sum_a  mu_aa
 *          - hbar^2/2m * sum_r  d^2/dq_r^2
 *          + V(q)
 *
 *      pi_a = -i sum_rs zeta_rs^a * q_r * d/dq_s
 *
 *      i ... imaginary unit
 *      a,b in {x,y,z}
 *      r,s in {modes}
 *
 */
static double * donothing = NULL;
static double * fst_deriv = NULL;
static double * sec_deriv = NULL;
static double conv_factor = 0.0;


void init_watson(settings* prefs, double dq){
// allocate memory for derivative stencils
    donothing = calloc(prefs->n_stencil,  sizeof(double));
    fst_deriv = malloc(prefs->n_stencil * sizeof(double));
    sec_deriv = malloc(prefs->n_stencil * sizeof(double));

    if(donothing == NULL){ perror("\"donothing\""); exit(errno); }
    if(fst_deriv == NULL){ perror("\"fst_deriv\""); exit(errno); }
    if(sec_deriv == NULL){ perror("\"sec_deriv\""); exit(errno); }

// set central point of "donothing" to 1.0
    donothing[prefs->n_stencil/2] = 1.0;

// fill first and second derivative stencils
    int errone =  FirstDerivative(prefs->n_stencil, fst_deriv);
    int errtwo = SecondDerivative(prefs->n_stencil, sec_deriv);
    if( errone || errtwo ){
        fprintf(stderr,
            "\n (-) Error initialising derivative stencil parameters."
            "\n     Aborting..."
            "\n\n"
        );
        exit(EXIT_FAILURE);
    }

// the filled derivative stencils do not contain the necessary division by dq
    for(int i = 0; i < prefs->n_stencil; ++i){
        fst_deriv[i] /=   dq;
        sec_deriv[i] /= (dq*dq);
    }

// Conversion factor from inverse moment of inertia to output unit of energy
// [conv_factor] = (output unit of energy) / [mu]
// Actually, it is the pre-factor of hbar^2/2 * (sum_ab pi_a * mu_ab * pi_b)
// the part in the parenthesis is calculated in the exec_watson() function
    conv_factor = 0.5 * 1.0E20 * hbar*hbar*avogadro*avogadro
                * prefs->InvInertia_to_molpergAasq
                * prefs->kJpermol_to_oue;
}


/* Memory layout of mu and zeta (z):
 * mu[3][3][index]  is a 3x3 matrix containing the [x;y;z] * [x y z] terms of the
 *                  reciprocal moment of inertia tensor for each configuration
 * z[3][D*(D-1)/2]  contains the Coriolis factors of all D*(D-1)/2 mode combinations
 */
double exec_watson(double*** mu, double** z, double** q, int index, int* shift){
//{{{

// calculate pre-factor
    double prefactor = 0.0;
    for(int a = 0; a < 3; ++a){
        for(int b = 0; b < 3; ++b){
            prefactor += mu[a][b][index] * z[a][0]*z[b][0];
        }
    }

    prefactor *= conv_factor;

    double watson = prefactor * (
                  q[0][index] *               fst_deriv[shift[0]] * donothing[shift[1]]
                - q[0][index] * q[0][index] * donothing[shift[0]] * sec_deriv[shift[1]]
          + 2.0 * q[0][index] * q[1][index] * fst_deriv[shift[0]] * fst_deriv[shift[1]]
                - q[1][index] * q[1][index] * sec_deriv[shift[0]] * donothing[shift[1]]
                + q[1][index] *               donothing[shift[0]] * fst_deriv[shift[1]]
    );

    return watson;
}//}}}


double exec_watson_3d(double*** mu, double** z, double** q, int index, int* shift){
//{{{

    double watson3d = 0.0;

    for(int a = 0; a < 3; ++a){
        for(int b = 0; b < 3; ++b){

            watson3d += mu[a][b][index] * (
            // fst total
                  z[a][0] * z[b][0] * q[0][index] * fst_deriv[shift[0]] * donothing[shift[1]] * donothing[shift[2]]
                + z[a][0] * z[b][0] * q[1][index] * donothing[shift[0]] * fst_deriv[shift[1]] * donothing[shift[2]]
                + z[a][0] * z[b][1] * q[1][index] * donothing[shift[0]] * donothing[shift[1]] * fst_deriv[shift[2]]
                - z[a][0] * z[b][2] * q[0][index] * donothing[shift[0]] * donothing[shift[1]] * fst_deriv[shift[2]]

                + z[a][1] * z[b][0] * q[2][index] * donothing[shift[0]] * fst_deriv[shift[1]] * donothing[shift[2]]
                + z[a][1] * z[b][1] * q[0][index] * fst_deriv[shift[0]] * donothing[shift[1]] * donothing[shift[2]]
                + z[a][1] * z[b][1] * q[2][index] * donothing[shift[0]] * donothing[shift[1]] * fst_deriv[shift[2]]
                + z[a][1] * z[b][2] * q[0][index] * donothing[shift[0]] * fst_deriv[shift[1]] * donothing[shift[2]]

                - z[a][3] * z[b][0] * q[2][index] * fst_deriv[shift[0]] * donothing[shift[1]] * donothing[shift[2]]
                + z[a][3] * z[b][1] * q[1][index] * fst_deriv[shift[0]] * donothing[shift[1]] * donothing[shift[2]]
                + z[a][3] * z[b][2] * q[1][index] * donothing[shift[0]] * fst_deriv[shift[1]] * donothing[shift[2]]
                + z[a][3] * z[b][2] * q[2][index] * donothing[shift[2]] * donothing[shift[1]] * fst_deriv[shift[2]]

            // second total
                - z[a][0] * z[b][0] * q[1][index] * q[1][index] * sec_deriv[shift[0]] * donothing[shift[1]] * donothing[shift[2]]
                - z[a][0] * z[b][0] * q[0][index] * q[0][index] * donothing[shift[0]] * sec_deriv[shift[1]] * donothing[shift[2]]
                - z[a][0] * z[b][1] * q[1][index] * q[2][index] * sec_deriv[shift[0]] * donothing[shift[1]] * donothing[shift[2]]
                + z[a][0] * z[b][2] * q[0][index] * q[2][index] * donothing[shift[0]] * sec_deriv[shift[1]] * donothing[shift[2]]

                - z[a][1] * z[b][0] * q[1][index] * q[2][index] * sec_deriv[shift[0]] * donothing[shift[1]] * donothing[shift[2]]
                - z[a][1] * z[b][1] * q[2][index] * q[2][index] * sec_deriv[shift[0]] * donothing[shift[1]] * donothing[shift[2]]
                - z[a][1] * z[b][1] * q[0][index] * q[0][index] * donothing[shift[0]] * donothing[shift[1]] * sec_deriv[shift[2]]
                - z[a][1] * z[b][2] * q[0][index] * q[1][index] * donothing[shift[0]] * donothing[shift[1]] * sec_deriv[shift[2]]

                + z[a][2] * z[b][0] * q[0][index] * q[2][index] * donothing[shift[0]] * sec_deriv[shift[1]] * donothing[shift[2]]
                - z[a][2] * z[b][1] * q[0][index] * q[1][index] * donothing[shift[0]] * donothing[shift[1]] * sec_deriv[shift[2]]
                - z[a][2] * z[b][2] * q[2][index] * q[2][index] * donothing[shift[0]] * sec_deriv[shift[1]] * donothing[shift[2]]
                - z[a][2] * z[b][2] * q[1][index] * q[1][index] * donothing[shift[0]] * donothing[shift[1]] * sec_deriv[shift[2]]

            // cross singled
                + z[a][0] * z[b][1] * q[0][index] * q[2][index] * fst_deriv[shift[0]] * fst_deriv[shift[1]] * donothing[shift[2]]
                + z[a][0] * z[b][1] * q[0][index] * q[1][index] * fst_deriv[shift[0]] * donothing[shift[1]] * fst_deriv[shift[2]]
                - z[a][0] * z[b][1] * q[0][index] * q[0][index] * donothing[shift[0]] * fst_deriv[shift[1]] * fst_deriv[shift[2]]
                - z[a][0] * z[b][2] * q[1][index] * q[2][index] * fst_deriv[shift[0]] * fst_deriv[shift[1]] * donothing[shift[2]]
                + z[a][0] * z[b][2] * q[1][index] * q[1][index] * fst_deriv[shift[0]] * donothing[shift[1]] * fst_deriv[shift[2]]
                - z[a][0] * z[b][2] * q[0][index] * q[1][index] * donothing[shift[0]] * fst_deriv[shift[1]] * fst_deriv[shift[2]]

                + z[a][1] * z[b][0] * q[0][index] * q[2][index] * fst_deriv[shift[0]] * fst_deriv[shift[1]] * donothing[shift[2]]
                + z[a][1] * z[b][0] * q[0][index] * q[1][index] * fst_deriv[shift[0]] * donothing[shift[1]] * fst_deriv[shift[2]]
                - z[a][1] * z[b][0] * q[0][index] * q[0][index] * donothing[shift[0]] * fst_deriv[shift[1]] * fst_deriv[shift[2]]
                - z[a][1] * z[b][2] * q[2][index] * q[2][index] * fst_deriv[shift[0]] * fst_deriv[shift[1]] * donothing[shift[2]]
                + z[a][1] * z[b][2] * q[1][index] * q[2][index] * fst_deriv[shift[0]] * donothing[shift[1]] * fst_deriv[shift[2]]
                + z[a][1] * z[b][2] * q[0][index] * q[2][index] * donothing[shift[0]] * fst_deriv[shift[1]] * fst_deriv[shift[2]]

                - z[a][2] * z[b][0] * q[1][index] * q[2][index] * fst_deriv[shift[0]] * fst_deriv[shift[1]] * donothing[shift[2]]
                + z[a][2] * z[b][0] * q[1][index] * q[1][index] * fst_deriv[shift[0]] * donothing[shift[1]] * fst_deriv[shift[2]]
                - z[a][2] * z[b][0] * q[0][index] * q[1][index] * donothing[shift[0]] * fst_deriv[shift[1]] * fst_deriv[shift[2]]
                - z[a][2] * z[b][1] * q[2][index] * q[2][index] * fst_deriv[shift[0]] * fst_deriv[shift[1]] * donothing[shift[2]]
                + z[a][2] * z[b][1] * q[1][index] * q[2][index] * fst_deriv[shift[0]] * donothing[shift[1]] * fst_deriv[shift[2]]
                + z[a][2] * z[b][1] * q[0][index] * q[2][index] * donothing[shift[0]] * fst_deriv[shift[1]] * fst_deriv[shift[2]]

            // cross doubled
                + 2 * z[a][0] * z[b][0] * q[0][index] * q[1][index] * fst_deriv[shift[0]] * fst_deriv[shift[1]] * donothing[shift[2]]
                + 2 * z[a][1] * z[b][1] * q[0][index] * q[2][index] * fst_deriv[shift[0]] * donothing[shift[1]] * fst_deriv[shift[2]]
                + 2 * z[a][2] * z[b][2] * q[1][index] * q[2][index] * donothing[shift[0]] * fst_deriv[shift[1]] * fst_deriv[shift[2]]
            );

        }
    }
    watson3d *= conv_factor;

    return watson3d;
}//}}}


void free_watson(void){

    free(donothing); donothing = NULL;
    free(fst_deriv); fst_deriv = NULL;
    free(sec_deriv); sec_deriv = NULL;

}
