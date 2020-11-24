
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <mkl_solvers_ee.h>

#include "settings.h"

// Dependencies
int FirstDerivative (int n_stencil, double*  first_derivative);
int SecondDerivative(int n_stencil, double* second_derivative);

// provided prototypes
int MKL_FillAMatrix2D(settings* prefs, int* nq, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);


// 2D fill
int MKL_FillAMatrix2D(settings* prefs, int* nq, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A){

// Calculate the maximum number of non-zero entries in the A matrix
//  Should be <n_points - (n_stencil/2)*2> lines with <n_stencil> entries, the
//  first and last <n_stencil/2> lines are shortened via a triangular number sequence
    int n_points = nq[0] * nq[1];
    int max_entries = (prefs->n_stencil*nq[0] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1))
                     *(prefs->n_stencil*nq[1] - (prefs->n_stencil/2)*(prefs->n_stencil/2 + 1));

// The sparse matrix eigensolver of Intel MKL saves the positions of non empty entries in
//  two integer arrays (rows_A & cols_A), each non zero entry increments the counter by 1
//  rows_A contains the counter value of each first element in a line -> max entries n_points (+1 because Intel...)
//  cols_A contains the counter values of all non zero elements -> max entries n_stencil*n_points - borders
//  vals_A contains the values of all non zero entries -> max entries n_stencil*n_points - borders
    (*rows_A) = calloc(n_points + 1, sizeof(MKL_INT));
    (*cols_A) = calloc(max_entries,  sizeof(MKL_INT));
    (*vals_A) = calloc(max_entries,  sizeof(double));

    if( (*rows_A) == NULL ){ perror("2D MKL Fill rows_A"); exit(errno); }
    if( (*cols_A) == NULL ){ perror("2D MKL Fill cols_A"); exit(errno); }
    if( (*vals_A) == NULL ){ perror("2D MKL Fill vals_A"); exit(errno); }


// variables required for the calculation of the Watson Hamiltonian rotational terms
    double * nothing   = NULL;  // "nothing" ensures the differentiation is only applied to one row/column/etc.
    double * fst_deriv = NULL;  // first derivative stencil
    double * sec_deriv = NULL;  // second derivative stencil

    if( prefs->coriolis_file ){
    // allocate memory for derivative stencils
        nothing   = calloc(prefs->n_stencil,  sizeof(double));
        fst_deriv = malloc(prefs->n_stencil * sizeof(double));
        sec_deriv = malloc(prefs->n_stencil * sizeof(double));

        if(nothing   == NULL){ perror("2D MKL \"nothing\""  ); exit(errno); }
        if(fst_deriv == NULL){ perror("2D MKL \"fst_deriv\""); exit(errno); }
        if(sec_deriv == NULL){ perror("2D MKL \"sec_deriv\""); exit(errno); }

    //  set central point of "nothing" to 1.0
        nothing[prefs->n_stencil/2] = 1.0;

    //  fill first and second derivative stencils
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
    }


// fill Numerov's A matrix
//  determine the non zero elements and store their positions in rows_A and cols_A
//  and their values in vals_A
    (*rows_A)[0] = 1;
    int entry_index = 0;

    for(int i = 0; i < nq[0]; ++i){
        for(int j = 0; j < nq[1]; ++j){

            for(int xsh = -(prefs->n_stencil/2); xsh < ((prefs->n_stencil/2) + 1); ++xsh){
            if( (i + xsh > -1) && (i + xsh < nq[0]) ){

                for(int ysh = -(prefs->n_stencil/2); ysh < ((prefs->n_stencil/2) + 1); ++ysh){
                if( (j + ysh > -1) && (j + ysh < nq[1]) ){

                    int xsidx = xsh + prefs->n_stencil/2;    // stencil x index
                    int ysidx = ysh + prefs->n_stencil/2;    // stencil y index

                // stencil entries have to be divided by 2 to get the right result
                    (*vals_A)[entry_index] = ekin_param * stencil[ xsidx*prefs->n_stencil + ysidx ] / 2.0;
                    (*cols_A)[entry_index] = (i + xsh)*nq[1] + (j + ysh) + 1;

                // enable second term of Watson Hamiltonian when Coriolis file is set {{{
                    if( prefs->coriolis_file ){
                    // calculate pre-factor
                        double prefactor = 0.0;
                        for(int n = 0; n < 3; ++n){
                            for(int m = 0; m < 3; ++m){
                                prefactor -= zeta[n][0]*zeta[m][0] * mu[n][m][i*nq[1] + j];
                            }
                        }
                    //  zeta is normalized      =>          non-dimensional
                    //  mu                  is given in     g/mol/angstrom^2
                    //  prefs->mu_factor     is given in     kJ/mol / [mu]
                    //  prefs->ekin_factor   is given in     (output unit of energy) / (kJ/mol)
                        prefactor *= ((prefs->mu_factor * prefs->ekin_factor)/2.0);

                        (*vals_A)[entry_index] -= prefactor * (
                                    q[0][i*nq[1] + j] *                     fst_deriv[xsidx] *   nothing[ysidx] / dq
                                  + q[1][i*nq[1] + j] *                       nothing[xsidx] * fst_deriv[ysidx] / dq
                                  - q[0][i*nq[1] + j] * q[0][i*nq[1] + j] * sec_deriv[xsidx] *   nothing[ysidx] / dq / dq
                                  - q[1][i*nq[1] + j] * q[1][i*nq[1] + j] *   nothing[xsidx] * sec_deriv[ysidx] / dq / dq
                            + 2.0 * q[0][i*nq[1] + j] * q[1][i*nq[1] + j] * fst_deriv[xsidx] * fst_deriv[ysidx] / dq / dq
                        );
                    }
                //}}}

                // add potential to diagonal element
                    if( (xsh == 0) && (ysh == 0) ){
                        (*vals_A)[entry_index] += v[i*nq[1] + j];
                    }
                    ++entry_index;
                }}
            }}
        // after inserting all entries in a row the total number of entries is inserted in the CSR format.
            (*rows_A)[i*nq[1] + (j + 1)] = entry_index + 1;
        }
    }

// free memory
    if( prefs->coriolis_file ){
        free(nothing);   nothing   = NULL;
        free(fst_deriv); fst_deriv = NULL;
        free(sec_deriv); sec_deriv = NULL;
    }

    printf("Matrix created, Potential added, %d entries\n", entry_index);
    return 0;
}
