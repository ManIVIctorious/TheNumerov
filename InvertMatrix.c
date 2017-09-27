
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

int InvertMatrix(gsl_matrix *Matrix, gsl_matrix *InvMatrix, int dimension);

int InvertMatrix(gsl_matrix *Matrix, gsl_matrix *InvMatrix, int dimension){

    int i, j, signum;
    gsl_matrix      *AUX = gsl_matrix_calloc(dimension, dimension);
    gsl_permutation   *p = gsl_permutation_alloc(dimension);

    // error handling
    if(AUX == NULL || p == NULL){
        fprintf(stderr, "ERROR: Memory for auxiliary matrix \"AUX\" and/or gsl_permutation \"p\"\n");
        fprintf(stderr, "       in InvertMatrix couldn't get allocated.\n");
        exit(-1);
    }

    // copy matrix to auxiliary matrix
    for(i=0; i<dimension; ++i){
        for(j=0; j<dimension; ++j){
            gsl_matrix_set(AUX, i, j, gsl_matrix_get(Matrix, i, j));
        }
    }

    gsl_linalg_LU_decomp(AUX, p, &signum);
    gsl_linalg_LU_invert(AUX, p, InvMatrix);

    gsl_matrix_free(AUX);
    gsl_permutation_free(p);

    return 0;
}

