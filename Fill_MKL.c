
#include <stdio.h>
#include <stdlib.h>
#include <mkl_solvers_ee.h>
#include "typedefinitions.h"

// Offered prototypes
int FillMKL_2D(settings prefs, int *nq, double *v, double ekin_param, double *stencil, double **q, double dq, double ***mu, double ***zeta, MKL_INT **rows_A, MKL_INT **cols_A, double **vals_A);


int FillMKL_2D(settings prefs, int *nq, double *v, double ekin_param, double *stencil, double **q, double dq, double ***mu, double ***zeta, MKL_INT **rows_A, MKL_INT **cols_A, double **vals_A){

double onestencil[52]={0.0, 0.0, 0.0, 0.0, 0.0, -1.0/2.0, 0.0, 1.0/2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0/60.0, 3.0/20.0, -3.0/4.0, 0.0, 3.0/4.0, -3.0/20.0, 1.0/60.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/280.0, -4.0/105.0, 1.0/5.0, -4.0/5.0, 0.0, 4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0, 0.0, 0.0};
double second_der[78]={0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, -1.0/12.0, 4.0/3.0, -5.0/2.0, 4.0/3.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 1.0/90.0, -3.0/20.0, 3.0/2.0, -49.0/18.0, 3.0/2.0, -3.0/20.0, 1.0/90.0, 0.0, 0.0, 0.0,0.0, 0.0, -1.0/560.0, 8.0/315.0, -1.0/5.0,  8.0/5.0,-205.0/72.0, 8.0/5.0, -1.0/5.0, 8.0/315.0, -1.0/560.0, 0.0, 0.0,0.0, 1.0/3150.0, -5.0/1008.0, 5.0/126.0, -5.0/21.0, 5.0/3.0, -5296.0/1800.0, 5.0/3.0, -5.0/21.0, 5.0/126.0, -5.0/1008.0, 1.0/3150.0, 0.0,-1.0/16632.0, 2.0/1925.0, -1.0/112.0, 10.0/189.0, -15.0/56.0, 12.0/7.0, -5369.0/1800.0, 12.0/7.0, -15.0/56.0,10.0/189.0, -1.0/112.0, 2.0/1925.0, -1.0/16632.0};

int m,n;
double vorfaktor=0.0;
double nothing[13]={0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0};// für die einzelnen Ableitungen.

int sec_st;
if(prefs.n_stencil<10)
{sec_st=(prefs.n_stencil-1)/2-1;}
else
{sec_st=3;}

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
    if((*rows_A) == NULL || (*cols_A) == NULL || (*vals_A) == NULL){
        fprintf(stderr,
            "\n (-) Error in memory allocation for row, column and/or value arrays"
            "\n     Aborting..."
            "\n\n"
        );
        exit(1);
    }

// fill Numerov's A matrix
//  determine the non zero elements and store their positions in rows_A and cols_A
//  as well as their values in vals_A
    for(i = 0; i < nq[0]; i++)
    {
        for(j = 0; j < nq[1]; j++)
       {

            for(xsh = -prefs.n_stencil/2; xsh < prefs.n_stencil/2 + 1; xsh++)
        {
                if( (i+xsh > -1) && (i+xsh < nq[0]) )
            {

                    for(ysh = -prefs.n_stencil/2; ysh < prefs.n_stencil/2 + 1; ysh++)
            {
                        if( (j+ysh > -1) && ( j+ysh < nq[1]) )
            {

                            element = (i + xsh)*nq[1] + j+ysh;
                            (*cols_A)[n_entries] = element+1; // wieso +1? weil intel!!


                        // stencil entries have to be divided by 2 to get the right result.
                        //  in three dimensions it should be a division by 4
                            vorfaktor=0.0;

                              for (n=0;n<3;n++)
                                 {
                                 for(m=0;m<3;m++)
                                 {

                                  vorfaktor=vorfaktor -   zeta[n][0][i*nq[1]+j]*zeta[m][0][i*nq[1]+j]*mu[n][m][i*nq[1]+j];

                                 }// for m
                                }// for n


                              (*vals_A)[n_entries] = ekin_param * stencil[(xsh+prefs.n_stencil/2)*prefs.n_stencil+ysh+prefs.n_stencil/2]/2.0;

     (*vals_A)[n_entries] -= ((prefs.mu_factor * prefs.ekin_factor)/2.0) * vorfaktor * (q[0][i*nq[1]+j] *                   onestencil[sec_st*13+6+xsh] * nothing[6+ysh]              / dq
                                                                                   + q[1][i*nq[1]+j] *                   onestencil[sec_st*13+6+ysh] * nothing[6+xsh]              / dq
                                                                                   - q[0][i*nq[1]+j] * q[0][i*nq[1]+j] * second_der[sec_st*13+6+ysh] * nothing[6+xsh]              / dq / dq
                                                                                   - q[1][i*nq[1]+j] * q[1][i*nq[1]+j] * second_der[sec_st*13+6+xsh] * nothing[6+ysh]              / dq / dq
                                                                                +2.0*q[0][i*nq[1]+j] * q[1][i*nq[1]+j] * onestencil[sec_st*13+6+xsh] * onestencil[sec_st*13+6+ysh] / dq / dq);


                                if(xsh == 0 && ysh ==0)                        // add potential to diagonal element
                                {
                                 (*vals_A)[n_entries] += v[i*nq[1]+j] - ((mu[0][0][i*nq[1]+j] + mu[1][1][i*nq[1]+j] + mu[2][2][i*nq[1]+j]) * (prefs.mu_factor * prefs.ekin_factor) / 8.0);// *dq*dq*mass*1.0/4.0; // watson pot

                                 }// end if xsh=ysh=0


                            n_entries ++;
                        }// end if 0<=ysh<nq[1]
                    }// end for ysh
                } // end if 0<=xsh<nq[0]
            }// end for xsh
      // after inserting all entries in a row the total number of entries is inserted in the CSR format.
        (*rows_A)[i*nq[1]+j+1] = n_entries+1;
        } // end for j
    }//end for i
    (*rows_A)[0] = 1;


    return 0;
}
