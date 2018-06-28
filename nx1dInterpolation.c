
#include <stdio.h>
#include <stdlib.h>

// provided prototypes
int nx1dInterpolation(double** v, int* nq_in, double dq, int dimension, int n_spline);

int nx1dInterpolation(double** v, int* nq_in, double dq, int dimension, int n_spline){

/*  One dimensional cubic spline interpolation:
//{{{
    The 1D cubic spline interpolation is based on the calculation of a third
    order (cubic) polynomial between each subset of neighbouring grid data points.

        f(h(i))   =   a  +   b h(i)  +   c h(i)^2  + d h(i)^3
        f'(h(i))  =   b  + 2 c h(i)  + 3 d h(i)^2
        f''(h(i)) = 2 c  + 6 d h(i)

    Where h is the grid spacing (h_n = x_n+1 - x_n). Using the first and second
    derivatives and the natural boundary conditions f''(h_0) = 0 and f''(h_n) = 0
    allows for the derivation of the following recursion formulas:

        a(i+1) = a(i) +   b(i) h(i) +   c(i) h(i)^2 + d(i) h(i)^3
        b(i+1) = b(i) + 2 c(i) h(i) + 3 d(i) h(i)^2
        c(i+1) = c(i) + 3 d(i) h(i)

    Here a directly corresponds to the original data points, giving 3n equations for
    the calculation of 3n unknown parameters. Rearrangement leads to

        3 ( (a(i+1) - a(i)) / h(i) - (a(i) - a(i-1)) / h(i-1)
        =
        h(i-1) c(i-1) + 2 (h(i) + h(i-1)) c(i) + h(i) c(i+1)

    On an equispaced grid (as required by the Numerov method) this simplifies to

        3/h^2 ( a(i-1) - 2 a(i) + a(i+1) ) = c(i-1) + 4 c(i) + c(i+1)

    In matrix notation (Ax = b) this can be represented as tridiagonal matrix equation

        | 4  1            | |  c_1  |          |             0              |
        | 1  4  1         | |  c_2  |          | a(0)   - 2 a(1)   + a(2)   |
        |    .  .  .      | |   .   |  = 3/h/h |             .              |
        |      .  .  .    | |   .   |          |             .              |
        |         1  4  1 | | c_n-2 |          | a(n-3) - 2 a(n-2) + a(n-1) |
        |            1  4 | | c_n-1 |          | a(n-2) - 2 a(n-1) + a(n)   |

    which can be solved using the Thomas (tridiagonal matrix) algorithm (vide infra).

    With this information the b and d coefficients can be calculated

        d(i) = 1/(3 h) * (c(i+1) - c(i))
        b(i) = 1/h     * (a(i+1) - a(i)) - h/3 * (2 c(i) + c(i+1))

//}}}*/

/*  Thomas algorithm for solving tridiagonal matrix problems:
//{{{
    The Thomas algorithm is a very efficient way to solve equations of the type

        a_i x_i-1 + b_i x_i + c_i x_i+1 = d_i

    with the two constraints a_1 = 0 and c_n = 0

    In matrix notation (Ax = b) this corresponds to

       | b_1    c_1                        | | x_1 |     | d_1 |
       | a_2    b_2    c_2                 | | x_2 |     | d_2 |
       |         .      .      .           | |  .  |  =  |  .  |
       |                .      .      .    | |  .  |     |  .  |
       |                      a_n    b_n   | | x_n |     | d_n |

    The algorithm starts with a forward sweep to eliminate the
    a coefficients by modifying the coefficients c (to c') and
    d (to d')

        c'_{i} =          c_{i}           / (b_{i} - a_{i} c'_{i-1})
        d'_{i} = (d_{i} - a_{i} d'_{i-1}) / (b_{i} - a_{i} c'_{i-1})

    Note, since a_1 = 0:
        c'_1 = c_1 / b_1
        d'_1 = d_1 / b_1

    The solution is then obtained by back substitution

        x_{i} = d'_{i} - c'_{i} x_{i+1}

    Since c_n = 0:
        x_{n} = d'_{n}

//}}}*/

/*
   In this particular case
       ∀ i ∊ (1,n] : a_i = 1
       ∀ i ∊ [1,n] : b_i = 4
       ∀ i ∊ [1,n) : c_i = 1

   allowing for further simplifications
     i = 1
       c'_{1} = 1/4
       d'_{1} = 0

     i ∊ [2,n)
       c'_{i} = 1 / (4 - c'_{i-1})
       d'_{i} = (d_{i} - d'_{i-1}) * c'_{i}

     i = n
       c'_{n} = 0
       d'_{n} = (d_{n} - d'_{n-1}) / (4 - c'_{n-1})
//*/

//----------------------------------------------------------------------
//  Declaration   Declaration   Declaration   Declaration   Declaration
//----------------------------------------------------------------------
// For better readability the parameters of the Thomas algorithm will
//  be prefixed with "matrix_" and the interpolation parameters will
//  be prefixed with "inter_"

    int i, j, k, l, m;
    int maxdim    = 0;
    int n_points  = 0;
    int nn_points = 0;
    double newdq  = 0;

// nq array containing number of points per dimension
//  the new nq array is needed to prevent overwriting of the original
//  value in the superordinate function stack
    int * nq = NULL;            // freed

// tridiagonal matrix algorithm coefficients
    double * matrix_c = NULL;   // freed
    double * matrix_d = NULL;   // freed
// interpolation coefficients
    double * inter_c  = NULL;   // freed
    double   inter_b;
    double   inter_d;
// interpolated data array
    double * yy = NULL;         // freed

// auxiliary pointers to allow freeing and swapping
    double * aux_v  = (*v);     // freed
    double * aux_yy = NULL;     // at the end points to new *v => don't free!!!

// variables for dimensional loop
//  to break down n-dimensional arrays to a set of 1D arrays
    int n_iteration;      // Number of 1D runs over ND array
    int jump;             // Index difference between the nth and (n+1)th entry for a given dimension
    int d_jump;           // The number of jumps because of superordinate dimensions
    int offset = 0;       // Offset added after each <d_jump> position, corresponds to previous jump
    int inter_offset = 0; // Offset on new interpolated array


//----------------------------------------------------------------------
//  Initialisation   Initialisation   Initialisation   Initialisation
//----------------------------------------------------------------------
// determine the size of the dimension with the most entries
//  and calculate n_points
    for(i = 0, maxdim = 0, n_points = 1, nn_points = 1; i < dimension; ++i){
        if(nq_in[i] > maxdim){
            maxdim = nq_in[i];
        }
        n_points *= nq_in[i];
        nn_points *= ((nq_in[i] - 1) * (n_spline + 1) + 1);
    }

// allocate memory for matrix_c, matrix_d and inter_c arrays
//  as well as the new array containing interpolated data yy
    matrix_c = malloc((maxdim-1) * sizeof(double)); // freed
    matrix_d = malloc((maxdim-1) * sizeof(double)); // freed
    inter_c  = malloc((maxdim-1) * sizeof(double)); // freed
    yy       = malloc(nn_points  * sizeof(double)); // freed
    aux_yy   = malloc(nn_points  * sizeof(double)); // points to new *v (old is freed by freeing aux_v)

//  allocate memory for nq and fill it with nq_in values
    nq = malloc(dimension * sizeof(int));
    for(i = 0; i < dimension; ++i){
        nq[i] = nq_in[i];
    }

// define new delta q (newdq) which accounts for interpolated points
    newdq = dq / ((double)n_spline + 1.0);

//----------------------------------------------------------------------------------------------------
//  Dimensional loop    Dimensional loop    Dimensional loop    Dimensional loop    Dimensional loop
//----------------------------------------------------------------------------------------------------
/* Break down n-dimensional arrays to a set of 1D arrays
//{{{
    Example: 3x2x4 3D dataset (n_points = 24)

        q0  q1  q2    index
        -----------------   The example cuboid is built upon 8 (2*4) one dimensional arrays,
        0   0   0   |   0   each containing 3 entries. The number of iterations needed to cover
        0   0   1   |   1   each element of this n-dimensional by a 1D array is therefore 8
        0   0   2   |   2
        0   0   3   |   3       n_iteration = n_points / nq[0] = 24/3 = 8

        0   1   0   |   4   The index difference between the nth and (n+1)th entry
        0   1   1   |   5   for the first dimension is 8.
        0   1   2   |   6
        0   1   3   |   7       jump = n_points / nq[0] = 24/3 = 8

        1   0   0   |   8   The first dimension can be mapped by a simple iteration over the first
        1   0   1   |   9   <n_iteration> points and an iteration over <nq[0]> times the jump.
        1   0   2   |  10
        1   0   3   |  11       0,8,16; 1,9,17; 2,10,18; etc.

        1   1   0   |  12   But by doing the same for the second dimension one gets
        1   1   1   |  13
        1   1   2   |  14       0, 4;  1, 5;  2, 6;  3, 7;  which are correct and
        1   1   3   |  15       4, 8;  5, 9;  6,10;  7,11;  etc. which are utter garbage...
                                8,12;  9,13; 10,14; 11,15;  would be correct
        2   0   0   |  16
        2   0   1   |  17   There are as many additional jump positions (n_iteration/d_jump),
        2   0   2   |  18   requiring an offset of the length of the previous jump,
        2   0   3   |  19   as  there are superordinate entries.

        2   1   0   |  20       q0: d_jump = 1               offset = n * 0
        2   1   1   |  21       q1: d_jump = nq[0]           offset = n * jump[0]
        2   1   2   |  22       q2: d_jump = nq[0] * nq[1]   offset = n * jump[1]   n ∊ [0,d_jump[
        2   1   3   |  23       etc.

    The following code segment outputs the indices of above 3D array as a set of 26 1D arrays
//*/
/*
    dimension=3;
    nq[0] = 3;
    nq[1] = 2;
    nq[2] = 4;
    n_points = 2*3*4;

    for(k = 0, jump = n_points, d_jump = 1; k < dimension; ++k){

        n_iteration = n_points/nq[k];
        jump /= nq[k];

        for(l = 0; l < d_jump; ++l){
            for(m = 0; m < n_iteration/d_jump; ++m){

                for(n = 0; n < nq[k]; ++n){
                    printf("\t%d", m + n*jump + l*offset);
                }
                printf("\n");
            }
        }
        printf("\n");

        offset = jump;
        d_jump *= nq[k];
    }
    return 0;
//}}}*/

    for(k = 0, jump = n_points, d_jump = 1; k < dimension; ++k){

        n_iteration = n_points/nq[k];
        jump /= nq[k];

        for(l = 0; l < d_jump; ++l){
            for(m = 0; m < n_iteration/d_jump; ++m){

            //----------------------------------------------------------------------
            //   Forward sweep   Forward sweep   Forward sweep   Forward sweep
            //----------------------------------------------------------------------
            // fill matrix_d array with tridiagonal matrix d elements
                matrix_d[0] = 0.0;
                for(i = 1; i < nq[k]-1; ++i){
                    matrix_d[i] = 3/dq/dq * (   (*v)[(i-1)*jump + m + l*offset]
                                             -2*(*v)[ i   *jump + m + l*offset]
                                             +  (*v)[(i+1)*jump + m + l*offset]);
                }

            // calculate tridiagonal matrix elements c' and d' and save them to
            //  matrix_c and matrix_d, respectively (overwrite matrix_d)
                matrix_c[0] = 1.0/4.0;
                for(i = 1; i < nq[k]-2; ++i){
                    matrix_c[i] = 1.0 / (4.0 - matrix_c[i-1]);
                    matrix_d[i] = (matrix_d[i] - matrix_d[i-1]) * matrix_c[i];
                }


            //----------------------------------------------------------------------
            // Back substitution to determine inter_c coefficients and interpolation
            //----------------------------------------------------------------------
            // solve the last entry of the matrix problem
                inter_c[i] = (matrix_d[i] - matrix_d[i-1]) / (4.0 - matrix_c[i-1]);

            // calculate inter_b and inter_d values from inter_c
                inter_d = - inter_c[i] / (3.0 * dq);
                inter_b = ((*v)[(i+1)*jump + m + l*offset] - (*v)[i*jump + m + l*offset]) / dq - (2.0 / 3.0 * dq * inter_c[i]);

            // interpolation procedure for the last n_spline+1 points
                for(j = n_spline+1; j >= 0; --j){
                    yy[(i * (n_spline+1) + j)*jump + m + l*inter_offset]
                                            = (*v)[i*jump + m + l*offset]
                                            + inter_b    * (newdq*j)
                                            + inter_c[i] * (newdq*j) * (newdq*j)
                                            + inter_d    * (newdq*j) * (newdq*j) * (newdq*j);
                }

                for(i = nq[k]-3; i >= 0; --i){
                // solve the rest of the matrix problem from back to front
                    inter_c[i] = matrix_d[i] - matrix_c[i] * inter_c[i+1];

                // calculate inter_b and inter_d values from inter_c
                    inter_d = (inter_c[i+1] - inter_c[i]) / (3.0 * dq);
                    inter_b = ((*v)[(i+1)*jump + m + l*offset] - (*v)[i*jump + m + l*offset]) / dq - dq / 3.0 * (2.0 * inter_c[i] + inter_c[i+1]);


                // interpolation procedure for the remaining points (all but the last set)
                    for(j = n_spline; j >= 0; --j){
                        yy[(i * (n_spline+1) + j)*jump + m + l*inter_offset]
                                                = (*v)[i*jump + m + l*offset]
                                                + inter_b    * (newdq*j)
                                                + inter_c[i] * (newdq*j) * (newdq*j)
                                                + inter_d    * (newdq*j) * (newdq*j) * (newdq*j);
                    }
                }
            }
        }

    // update *v to point to yy, swap yy and aux_yy
    //  and initialize new yy to zero.
        (*v) = yy;

        yy = aux_yy;
        for(i = 0; i < nn_points; ++i){
            yy[i] = 0;
        }
        aux_yy = (*v);

    // update size of dimension and n_points
        n_points /= nq[k];
        nq[k]     = (nq[k]-1) * (n_spline+1) + 1;
        n_points *= nq[k];

        offset = jump;
        inter_offset = (jump-1) * (n_spline+1) + 1;

        d_jump *= nq[k];
    }


// free memory
    free(matrix_c); matrix_c = NULL;
    free(matrix_d); matrix_d = NULL;
    free(inter_c);  inter_c  = NULL;
    free(aux_v);    aux_v    = NULL;
    free(yy);       yy       = NULL;
    free(nq);       nq       = NULL;

    return n_points;

}
