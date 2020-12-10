#include <mkl_solvers_ee.h>

// The MKL CSR format requires the matrix entries to be in order, i.e. from left to right.
// In cases where this is not done directly by the filling routine the following, simple
// matrix line sorter can be utilised.
void mkl_sort_matrix_rows(int n_rows, MKL_INT* row_start_index, MKL_INT* column_index, double* value){

// iterate over each row
    for(int i = 0; i < n_rows; ++i){

        int start = row_start_index[ i ] - 1;
        int end   = row_start_index[i+1] - 1;

    // sort column indices and values
        int swapped;
        do{

            swapped = 0;
            for(int j = start; j < end-1; ++j){

                if( column_index[j] > column_index[j+1] ){
                // swap column indices
                    long long int tmp_column = column_index[j];
                    column_index[ j ] = column_index[j+1];
                    column_index[j+1] = tmp_column;
                // swap values
                    double tmp_value = value[j];
                    value[ j ] = value[j+1];
                    value[j+1] = tmp_value;
                // indicate swap happened
                    swapped++;
                }
            }
        }while( swapped );

    }
}
