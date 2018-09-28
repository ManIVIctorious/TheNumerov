
/*  Finite difference stencils:
//{{{
    
    For the calculation of the finite difference stencils of e.g. the fourth derivative
    one starts with the expansion to a Taylor series

    Taylor series:
    
        f(b) = ∑_{n=0}^∞ f⁽ⁿ⁾(a) / n! · (b-a)^n = f(a) + f'(a)·(b-a) + f''(a)/2 · (b-a)^2 + f'''(a)/6 · (b-a)^3 +···

    Example:

        ∂^4/∂x^4 ≈ A·f(x-2h) + B·f(x-h) + C·f(x) + D·f(x+h) + E·f(x+2h)

    Next insert the values for a and b where a is the starting point i.e. a = x
    and b is the function argument, e.g.: A: b = x-2h; B: b = x-h; C: b = x; etc.

        f(x-2h) ≈ f(x) + f'(x)·(-2h) + f''(x)·( 4h²)/2 + f'''(x)·(-8h³)/6 + f⁽⁴⁾(x)·(16h⁴)/24 +···
        f(x-1h) ≈ f(x) + f'(x)·(-1h) + f''(x)·(  h²)/2 + f'''(x)·(-1h³)/6 + f⁽⁴⁾(x)·( 1h⁴)/24 +···
        f(x)    = f(x) + f'(x)·( 0h) + f''(x)·( 0h²)/2 + f'''(x)·( 0h³)/6 + f⁽⁴⁾(x)·( 0h⁴)/24 +···
        f(x+1h) ≈ f(x) + f'(x)·( 1h) + f''(x)·(  h²)/2 + f'''(x)·( 1h³)/6 + f⁽⁴⁾(x)·( 1h⁴)/24 +···
        f(x+2h) ≈ f(x) + f'(x)·( 2h) + f''(x)·( 4h²)/2 + f'''(x)·( 8h³)/6 + f⁽⁴⁾(x)·(16h⁴)/24 +···

    writing down above equations in matrix form yields

        | 1   -2h    4h²/2   -8h³/6   16h⁴/24 |   | f⁽⁰⁾(x) |     | f(x-2h) |
        | 1   -1h     h²/2   -1h³/6    1h⁴/24 |   | f⁽¹⁾(w) |     | f(x-1h) |
        | 1    0     0        0        0      | · | f⁽²⁾(x) |  =  | f(x)    |
        | 1    1h     h²/2    1h³/6    1h⁴/24 |   | f⁽³⁾(x) |     | f(x+1h) |
        | 1    2h    4h²/2    8h³/6   16h⁴/24 |   | f⁽⁴⁾(x) |     | f(x+2h) |

    which is equivalent to

        | 1   -2    4/2   -8/6   16/24 |   | 1  · f⁽⁰⁾(x) |     | f(x-2h) |            | 1  · f⁽⁰⁾(x) |     | f(x-2h) |
        | 1   -1    1/2   -1/6    1/24 |   | h  · f⁽¹⁾(w) |     | f(x-1h) |            | h  · f⁽¹⁾(w) |     | f(x-1h) |
        | 1    0     0      0      0   | · | h² · f⁽²⁾(x) |  =  | f(x)    |   or   A · | h² · f⁽²⁾(x) |  =  | f(x)    |
        | 1    1    1/2    1/6    1/24 |   | h³ · f⁽³⁾(x) |     | f(x+1h) |            | h³ · f⁽³⁾(x) |     | f(x+1h) |
        | 1    2    4/2    8/6   16/24 |   | h⁴ · f⁽⁴⁾(x) |     | f(x+2h) |            | h⁴ · f⁽⁴⁾(x) |     | f(x+2h) |

    Note that the entries of the A matrix all follow the formula (b-a)^n / n!
    with n incrementing from left (n = 0) to right (in this particular example n = 4)

    This equation system can be solved to

        | 1  · f⁽⁰⁾(x) |     |   0       0       1       0       0   |   | f(x-2h) |              | f(x-2h) |
        | h  · f⁽¹⁾(w) |     |  1/12   -8/12     0     -8/12    1/12 |   | f(x-1h) |              | f(x-1h) |
        | h² · f⁽²⁾(x) |  =  | -1/12   16/12  -30/12   16/12   -1/12 | · | f(x)    |  =  inv(A) · | f(x)    |
        | h³ · f⁽³⁾(x) |     | -1/ 2    2/ 2     0     -2/ 2    1/ 2 |   | f(x+1h) |              | f(x+1h) |
        | h⁴ · f⁽⁴⁾(x) |     |   1      -4       6      -4       1   |   | f(x+2h) |              | f(x+2h) |

    Each line (i ∈ [0,n-1]) of the matrix inv(A) now contains the finite difference stencil of the i-th derivative
    beginning with the identity (zeroth derivative) to the n-1th derivative.
    In the particular example of the matrix contains the identity, first, second, third and fourth derivative stencils.
    For its application each stencil has to be divided by h^<derivative>
    
//}}}*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permute.h>
#include <gsl/gsl_linalg.h>


// Provided prototypes
int FiniteDifferenceStencil(double* stencil, unsigned int size, unsigned int derivative);
static unsigned int factorial(unsigned int n);


int FiniteDifferenceStencil(double* stencil, unsigned int size, unsigned int derivative){

    if(derivative > size){
        fprintf(stderr,
            "\n (-) Error: The calculation of finite difference stencils"
            "\n     requires at least as many points in stencil size (%d)"
            "\n     as the requested derivative order (%d)."
            "\n     Aborting...\n\n", size, derivative
        );
        exit(EXIT_FAILURE);
    }

    int signum;
    unsigned int i, j;
    gsl_matrix      * A     = NULL;
    gsl_matrix      * invA  = NULL;
    gsl_permutation * p     = NULL;

// gsl has its own error handling
    A = gsl_matrix_alloc(size, size);
    p = gsl_permutation_alloc(size);

// set values of first and second column
    for(i = 0; i < size; ++i){
        gsl_matrix_set(A, i, 0, 1.0);
        gsl_matrix_set(A, i, 1, stencil[i]);
    }

// fill remaining columns
    for(i = 2; i < size; ++i){
        for(j = 0; j < size; ++j){
            gsl_matrix_set(A, j, i, pow(stencil[j], i) / (double)factorial(i));
        }
    }

// calculate LU decomposition of matrix A. The LU matrix is stored in matrix A
//  using the upper triangle and main diagonal as U matrix
//  and the lower triangle as L matrix
    gsl_linalg_LU_decomp(A, p, &signum);

// allocate memory for inverse of matrix A
    invA = gsl_matrix_alloc(size,size);

// compute the inverse of matrix A from its LU decomposition (LU, p)
//  The inverse is computed by solving the system A x = b for each column of
//  the identity matrix. It is preferable to avoid direct use of the inverse
//  whenever possible, as the linear solver functions can obtain the same
//  result more efficiently and reliably
    gsl_linalg_LU_invert(A, p, invA);

// free memory of A matrix and permutation p
    gsl_matrix_free(A);
    gsl_permutation_free(p);

// store the requested stencil values in stencil
    for(i = 0; i < size; ++i){
        stencil[i] = gsl_matrix_get(invA, derivative, i);
    }

// free memory of invA
    gsl_matrix_free(invA);

    return 0;
}


static unsigned int factorial(unsigned int n){

    unsigned int factorial = 1;

    if(n == 0){ return  1; }

    while(n > 0){
        factorial *= n;
        --n;
    }

    return factorial;
}
