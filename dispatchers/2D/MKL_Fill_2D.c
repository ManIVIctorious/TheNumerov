
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <mkl_solvers_ee.h>

#include "settings.h"

// Dependencies
int FirstDerivative (int n_stencil, double*  first_derivative);
int SecondDerivative(int n_stencil, double* second_derivative);

// provided prototypes
int FillMKL_2D(settings prefs, int* nq, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);


// 2D fill
int FillMKL_2D(settings prefs, int* nq, double* v, double ekin_param, double* stencil, double** q, double dq, double*** mu, double** zeta, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A){

    int i, j;
    int n_entries = 0;
    int n_points = nq[0] * nq[1];
    int xsh, ysh;
    int element;

    // calculate max_entries
    int max_entries = 0;
    int sum_q1 = nq[0];
    int sum_q2 = nq[1];

    for(i = 1; i < (prefs.n_stencil/2 + 1); ++i){
        sum_q1 += 2*(nq[0]-i);
        sum_q2 += 2*(nq[1]-i);
    }
    max_entries = sum_q1*sum_q2; // upper estimation for nnz entries in the matrix, but the easy way to code.

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


// variables needed for the calculation of the Watson Hamiltonian rotational terms
    int m, n;
    double   prefactor = 0.0;

//  derivative stencils for the calculation of rotational terms
    double * nothing   = NULL;  // "nothing" ensures the differentiation is only applied to one row/column/etc.
    double * fst_deriv = NULL;  // first derivative stencil
    double * sec_deriv = NULL;  // second derivative stencil


    if(prefs.coriolis_file_set){
    // allocate memory for derivative stencils
        nothing   = calloc(prefs.n_stencil,  sizeof(double));
        fst_deriv = malloc(prefs.n_stencil * sizeof(double));
        sec_deriv = malloc(prefs.n_stencil * sizeof(double));

        if(nothing   == NULL){ perror("2D MKL \"nothing\""  ); exit(errno); }
        if(fst_deriv == NULL){ perror("2D MKL \"fst_deriv\""); exit(errno); }
        if(sec_deriv == NULL){ perror("2D MKL \"sec_deriv\""); exit(errno); }


    //  set central point of "nothing" to 1.0
        nothing[prefs.n_stencil/2] = 1.0;

    //  fill first and second derivative stencils
        m = FirstDerivative (prefs.n_stencil, fst_deriv);
        n = SecondDerivative(prefs.n_stencil, sec_deriv);
        if(m != 0 || n != 0){
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
//  as well as their values in vals_A
    for(i = 0; i < nq[0]; ++i){
        for(j = 0; j < nq[1]; ++j){
            for(xsh = -prefs.n_stencil/2; xsh < prefs.n_stencil/2 + 1; ++xsh){

                if( (i+xsh > -1) && (i+xsh < nq[0]) ){
                    for(ysh = -prefs.n_stencil/2; ysh < prefs.n_stencil/2 + 1; ++ysh){

                        if( (j+ysh > -1) && (j+ysh < nq[1]) ){
                            element = (i + xsh)*nq[1] + j+ysh;
                            (*cols_A)[n_entries] = element+1; // wieso +1? weil intel!!


                        // stencil entries have to be divided by 2 to get the right result.
                        //  in three dimensions it should be a division by 4
                            (*vals_A)[n_entries] = ekin_param * stencil[(xsh + prefs.n_stencil/2)*prefs.n_stencil + ysh + prefs.n_stencil/2]/2.0;

                        // enable second term of Watson Hamiltonian when Coriolis file is set
                        //####################################################################################################
                            if(prefs.coriolis_file_set){
                            // calculate pre-factor
                                for(n = 0, prefactor = 0.0; n < 3; ++n){
                                    for(m = 0; m < 3; ++m){
                                        prefactor -= zeta[n][0]*zeta[m][0] * mu[n][m][i*nq[1] + j];
                                    }
                                }
                            //  zeta is normalized      =>          non-dimensional
                            //  mu                  is given in     g/mol/angstrom^2
                            //  prefs.mu_factor     is given in     kJ/mol / [mu]
                            //  prefs.ekin_factor   is given in     (output unit of energy) / (kJ/mol)
                                prefactor *= ((prefs.mu_factor * prefs.ekin_factor)/2.0);

                                (*vals_A)[n_entries] -= prefactor * (
                                            q[0][i*nq[1]+j] *                   fst_deriv[(prefs.n_stencil/2) + xsh] *   nothing[(prefs.n_stencil/2) + ysh] / dq
                                          + q[1][i*nq[1]+j] *                   fst_deriv[(prefs.n_stencil/2) + ysh] *   nothing[(prefs.n_stencil/2) + xsh] / dq
                                          - q[0][i*nq[1]+j] * q[0][i*nq[1]+j] * sec_deriv[(prefs.n_stencil/2) + ysh] *   nothing[(prefs.n_stencil/2) + xsh] / dq / dq
                                          - q[1][i*nq[1]+j] * q[1][i*nq[1]+j] * sec_deriv[(prefs.n_stencil/2) + xsh] *   nothing[(prefs.n_stencil/2) + ysh] / dq / dq
                                    + 2.0 * q[0][i*nq[1]+j] * q[1][i*nq[1]+j] * fst_deriv[(prefs.n_stencil/2) + xsh] * fst_deriv[(prefs.n_stencil/2) + ysh] / dq / dq
                                );
                            }
                        //####################################################################################################


                        // add potential to diagonal element
                            if(xsh == 0 && ysh == 0){
                                (*vals_A)[n_entries] += v[i*nq[1] + j];
                            }

                            n_entries++;
                        }
                    }
                }
            }
      // after inserting all entries in a row the total number of entries is inserted in the CSR format.
        (*rows_A)[i*nq[1] + j + 1] = n_entries + 1;
        }
    }
    (*rows_A)[0] = 1;

// free memory
    if(prefs.coriolis_file_set){
        free(nothing);   nothing   = NULL;
        free(fst_deriv); fst_deriv = NULL;
        free(sec_deriv); sec_deriv = NULL;
    }

    return 0;
}
