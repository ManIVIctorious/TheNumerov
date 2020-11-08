
#include <stdio.h>
#include <stdlib.h>

/* Breaking down an nD array
//{{{
  Example: 3x2x4 3D dataset (n_points = 24)

  q0  q1  q2    index
  -------------------   For each dimension the individual position indices are arranged
  0   0   0   |   0     in a number of blocks, where each block represents a frame of
  0   0   1   |   1     iso-space of the previous dimension.
  0   0   2   |   2     The first dimension only consists of one block, while for
  0   0   3   |   3     subsequent dimensions it is represented by the product of the
                        number of unique points of the previous dimensions.
  0   1   0   |   4     The representation will be directed from the last dimension to
  0   1   1   |   5     the first (D, D-1, ..., 1) and therefore <n_blocks> breaks down
  0   1   2   |   6     to the initialisation and the new assignment for each dimension i
  0   1   3   |   7
                            n_blocks    = n_points
  1   0   0   |   8         n_blocks[i] = n_blocks / gridsize[i]
  1   0   1   |   9
  1   0   2   |  10     With <gridsize[i]> being the number of unique entries per
  1   0   3   |  11     dimension.

  1   1   0   |  12     The size of these blocks can be calculated via
  1   1   1   |  13
  1   1   2   |  14         blocksize    = n_points / n_blocks
  1   1   3   |  15                      = gridsize[last] * ... * gridsize[i]
                            blocksize[0] = n_points
  2   0   0   |  16
  2   0   1   |  17     For all dimensions (except the last one) each of the unique
  2   0   2   |  18     indices occurs subsequently multiple times. Hence, the distance
  2   0   3   |  19     between two subsequent unique entries has to be determined.

  2   1   0   |  20         distance       = blocksize / gridsize[i]
  2   1   1   |  21                        = gridsize[last] * ... * gridsize[i+1]
  2   1   2   |  22         distance[last] = 1
  2   1   3   |  23

  In above example:
      n_points  = 24        // total number of points 3 x 2 x 4
      gridsize  =  3, 2, 4  // number of unique points per dimension
      n_blocks  =  1, 3, 6  // number of blocks to be looped over
      blocksize = 24, 8, 4  // size of each individual block
      distance  =  8, 4, 1  // distance between consecutive points
//}}}*/

// provided prototypes
double CheckCoordinateSpacing(double** q, int* gridsize, double threshold, int dimension);

double CheckCoordinateSpacing(double** q, int* gridsize, double threshold, int dimension){

// precalculate threshold squared
    double threshold_sq = threshold*threshold;

// get total number of points
    int n_points = 1;
    for(int i = 0; i < dimension; ++i){ n_points *= gridsize[i]; }

// get initial spacing value
// calcluated between the first two entries of the last coordinate column
    double initial_dq = q[dimension-1][1] - q[dimension-1][0];


// First check the spacing within each dimension:
//  q[i][index2] - q[i][index1] must be same as dq
    int n_blocks = n_points;
    for(int i = (dimension-1); i >= 0; --i){

        n_blocks     /= gridsize[i];
        int blocksize = n_points  / n_blocks;
        int distance  = blocksize / gridsize[i];

    // iterate over each individual block
        for(int j = 0; j < n_blocks; ++j){
            for(int k = 0; k < (blocksize-distance); ++k){

                int index2 = j*blocksize + k + distance;
                int index1 = j*blocksize + k;

                double dq = q[i][index2] - q[i][index1];

                if( (dq - initial_dq)*(dq - initial_dq) > threshold_sq ){
                    fprintf(stderr,
                        "\n (-) Error: Coordinate spacing not equivalent:"
                        "\n       target = % .12lf"
                        "\n       value  = % .12lf"
                        "\n     Aborting - please check your input...\n\n"
                        , initial_dq, dq
                    );
                    exit(EXIT_FAILURE);
                }

            }
        }
    }


// Then check if all start points within a dimension have been the same
//  q[i][index2] - q[i][index1] must be zero
    n_blocks = n_points;
    for(int i = (dimension-1); i >= 0; --i){

        n_blocks     /= gridsize[i];
        int blocksize = n_points / n_blocks;
        int distance  = blocksize / gridsize[i];

        for(int j = 1; j < n_blocks; ++j){
            for(int k = 0; k < distance; ++k){

                double dq = q[i][j*blocksize+k] - q[i][k];

                if( (dq*dq) > threshold_sq ){
                    fprintf(stderr,
                        "\n (-) Error: Initial block-coordinates not equivalent:"
                        "\n       target = % .12lf"
                        "\n       value  = % .12lf"
                        "\n     Aborting - please check your input...\n\n"
                        , q[i][0], q[i][j*blocksize]
                    );
                    exit(EXIT_FAILURE);
                }

            }
        }
    }

// Finally, check start points of primary dimension
// above loop does not iterate over q[0] due to n_blocks being 1
    for(int i = 1; i < (n_points / gridsize[0]); ++i){

        double dq = q[0][i] - q[0][0];

        if( (dq*dq) > threshold_sq ){
            fprintf(stderr,
                "\n (-) Error: Initial block-coordinates not equivalent:"
                "\n       target = % .12lf"
                "\n       value  = % .12lf"
                "\n     Aborting - please check your input...\n\n"
                , q[0][0], q[0][i]
            );
            exit(EXIT_FAILURE);
        }

    }


    return initial_dq;
}
