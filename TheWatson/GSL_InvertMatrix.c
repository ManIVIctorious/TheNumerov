
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

void InvertMatrix(gsl_matrix* Matrix, gsl_matrix* InvMatrix, int dimension);

void InvertMatrix(gsl_matrix* Matrix, gsl_matrix* InvMatrix, int dimension){

    gsl_matrix      *AUX = gsl_matrix_calloc(dimension, dimension);
    gsl_permutation   *p = gsl_permutation_alloc(dimension);

    // copy matrix to auxiliary matrix
    for(int i = 0; i < dimension; ++i){
        for(int j = 0; j < dimension; ++j){
            gsl_matrix_set(AUX, i, j, gsl_matrix_get(Matrix, i, j));
        }
    }

    int signum;
    gsl_linalg_LU_decomp(AUX, p, &signum);
    gsl_linalg_LU_invert(AUX, p, InvMatrix);

    gsl_matrix_free(AUX);
    gsl_permutation_free(p);
}

