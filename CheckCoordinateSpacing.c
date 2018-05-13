
#include <stdio.h>
#include <stdlib.h>

// provided prototypes
double CheckCoordinateSpacing(double** q, int* nq, double threshold, int dimension);

double CheckCoordinateSpacing(double** q, int* nq, double threshold, int dimension){

    int i, j, k, l;
    int n_points;
    double dq;
    double dq_new = 0;          // new dq to check dq against (in coordinate check)
    int n_blocks, block_size;   // integers for coordinate spacing
    int index1, index2;         // indices for better readability

    for(i = 0, n_points = 1; i < dimension; ++i){
        n_points *= nq[i];
    }

/* Breaking down a ND array to 4 loops
//{{{
    Example: 3x2x4 3D dataset (n_points = 24)

        q0  q1  q2    index
        -------------------     Each dimensional position index shows up n_points/nq[i] times,
        0   0   0   |   0       giving the same number of starting points <n_starts>.
        0   0   1   |   1
        0   0   2   |   2           n_starts = n_points / nq[i]
        0   0   3   |   3
                                Additionally for all dimensions (except the last one) each
        0   1   0   |   4       of these start indices occur subsequently multiple times.
        0   1   1   |   5       The number of subsequent occurrences of the same position
        0   1   2   |   6       index is referred to as <block_size>
        0   1   3   |   7
                                    block_size = n_starts / n_blocks
        1   0   0   |   8                      = n_points / nq[i] / n_blocks
        1   0   1   |   9
        1   0   2   |  10       Where the number of blocks <n_blocks> is 1 for the first and
        1   0   3   |  11       n_points for the last dimension. The representation will be
                                directed from the last dimension to the first (D, D-1, ..., 1)
        1   1   0   |  12       and therefore <n_blocks> breaks down to the initialisation
        1   1   1   |  13
        1   1   2   |  14           n_blocks = n_points
        1   1   3   |  15
                                and the new assignment for each dimension i
        2   0   0   |  16
        2   0   1   |  17           n_blocks = n_blocks / nq[i]
        2   0   2   |  18
        2   0   3   |  19       With the last two variables <block_size> and <n_blocks> as well
                                as the number of entries for all dimensions <nq[i]> and the
        2   1   0   |  20       resulting total number of points <n_points> it should be possible
        2   1   1   |  21       to represent all wished array indices.
        2   1   2   |  22       Below is an example code to represent the running index, and all
        2   1   3   |  23       dimensions q.

//  Example code
//{{{
int main(int argc, char **argv){

   int i,j,k,l;
   int n_blocks;
   int block_size;
   int index;

// prerequisites:
int dimension = 3;
int n_points;
int * nq = malloc(dimension * sizeof(int));
nq[0] = 3;
nq[1] = 2;
nq[2] = 4;

for(i = 0, n_points = 1; i < dimension; ++i){
    n_points *= nq[i];
}

    for(i = (dimension-1), n_blocks = n_points; i >= 0; --i){

        n_blocks   /= nq[i];
        block_size  = n_points / nq[i] / n_blocks;

        printf("\nDimension = %2d\n\n", i);

        for(j = 0; j < n_blocks; ++j){
            for(k = 0; k < nq[i]; ++k){
                for(l = 0; l < block_size; ++l){

                    index = l + k*block_size + j*block_size*nq[i];

                    printf("\t%2d", index);
                    printf("\t%2d", k);
                    printf("\n");
                }
            }
        }
        printf("\n");
    }
    return 0;

}
//}}}
//}}}*/

// get initial value to compare with:
    dq = q[dimension-1][1] - q[dimension-1][0];

// first check spacing within each dimension:
//  q[i][index2] - q[i][index1] must be same as dq

    for(i = (dimension-1), n_blocks = n_points; i >= 0; --i){

        n_blocks   /= nq[i];
        block_size  = n_points/nq[i]/n_blocks;

        for(j = 0; j < n_blocks; ++j){
            for(k = 0; k < nq[i]-1; ++k){
                for(l = 0; l < block_size; ++l){

                    index1 = l +     k*block_size + j*block_size*nq[i];
                    index2 = l + (k+1)*block_size + j*block_size*nq[i];

                    dq_new = q[i][index2] - q[i][index1];

                    if( (dq - dq_new) * (dq - dq_new) > threshold*threshold ){
                        fprintf(stderr,
                            "\n (-) Error in input file."
                            "\n     Coordinate spacing (check no. 1) not equivalent:"
                            "\n     \ttarget = % .12lf"
                            "\n     \tvalue  = % .12lf"
                            "\n     Aborting - please check your input..."
                            "\n\n"
                            ,dq, dq_new
                        );
                        exit(-1);
                    }

                }
            }
        }
    }

// now check if all start points within a dimension have been the same
//  q[i][index2] - q[i][index1] must be zero
    for(i = (dimension-1), n_blocks = n_points; i >= 0; --i){

        n_blocks   /= nq[i];
        block_size  = n_points/nq[i]/n_blocks;

        for(j = 0; j < (n_blocks-1); ++j){
            for(k = 0; k < block_size; ++k){

                index1 = k +   j  *nq[i]*block_size;
                index2 = k + (j+1)*nq[i]*block_size;

                dq_new = q[i][index2] - q[i][index1];

                if( dq_new * dq_new > threshold*threshold ){
                    fprintf(stderr,
                        "\n (-) Error in input file."
                        "\n     Coordinate spacing (check no. 2) not equivalent:"
                        "\n     \ttarget = % .12lf"
                        "\n     \tvalue  = % .12lf"
                        "\n     Aborting - please check your input..."
                        "\n\n"
                        ,dq, dq_new
                    );
                    exit(-1);
                }

            }
        }
    }

// above loop does not iterate over q[0] since in this case n_blocks = 1
//  therefore check q[0] start points explicitly
    block_size = n_points/nq[0];
    for(i = 0; i < (block_size-1); ++i){

        dq_new = q[0][i+1] - q[0][i];

        if( dq_new * dq_new > threshold*threshold ){
            fprintf(stderr,
                "\n (-) Error in input file."
                "\n     Coordinate spacing (check no. 3) not equivalent:"
                "\n     \ttarget = % .12lf"
                "\n     \tvalue  = % .12lf"
                "\n     Aborting - please check your input..."
                "\n\n"
                ,dq, dq_new
            );
            exit(-1);
        }

    }

    return dq;
}
