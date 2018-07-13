
#include <stdio.h>
#include <stdlib.h>

// provided prototypes
int FirstDerivative (int n_stencil, double*  first_derivative);
int SecondDerivative(int n_stencil, double* second_derivative);


int FirstDerivative(int n_stencil, double* first_derivative){

// fill first derivative stencil
    switch(n_stencil){

        case 3:
            first_derivative[0] = -1.0/2.0;

            first_derivative[1] =    0.0;

            first_derivative[2] = -first_derivative[0];
            return 0;

        case 5:
            first_derivative[0] =  1.0/12.0;
            first_derivative[1] = -8.0/12.0;

            first_derivative[2] =    0.0;

            first_derivative[3] = -first_derivative[1];
            first_derivative[4] = -first_derivative[0];
            return 0;

        case 7:
            first_derivative[0] =  -1.0/60.0;
            first_derivative[1] =   9.0/60.0;
            first_derivative[2] = -45.0/60.0;

            first_derivative[3] =     0.0;

            first_derivative[4] = -first_derivative[2];
            first_derivative[5] = -first_derivative[1];
            first_derivative[6] = -first_derivative[0];
            return 0;

        case 9:
            first_derivative[0] =    3.0/840.0;
            first_derivative[1] =  -32.0/840.0;
            first_derivative[2] =  168.0/840.0;
            first_derivative[3] = -672.0/840.0;

            first_derivative[4] =      0.0;

            first_derivative[5] = -first_derivative[3];
            first_derivative[6] = -first_derivative[2];
            first_derivative[7] = -first_derivative[1];
            first_derivative[8] = -first_derivative[0];
            return 0;

        case 11:
            first_derivative[ 0] =    -2.0/2520.0;
            first_derivative[ 1] =    25.0/2520.0;
            first_derivative[ 2] =  -150.0/2520.0;
            first_derivative[ 3] =   600.0/2520.0;
            first_derivative[ 4] = -2100.0/2520.0;

            first_derivative[ 5] =       0.0;

            first_derivative[ 6] = -first_derivative[4];
            first_derivative[ 7] = -first_derivative[3];
            first_derivative[ 8] = -first_derivative[2];
            first_derivative[ 9] = -first_derivative[1];
            first_derivative[10] = -first_derivative[0];
            return 0;

        case 13:
            first_derivative[ 0] =      5.0/27720.0;
            first_derivative[ 1] =    -72.0/27720.0;
            first_derivative[ 2] =    495.0/27720.0;
            first_derivative[ 3] =  -2200.0/27720.0;
            first_derivative[ 4] =   7425.0/27720.0;
            first_derivative[ 5] = -23760.0/27720.0;

            first_derivative[ 6] =        0.0;

            first_derivative[ 7] = -first_derivative[5];
            first_derivative[ 8] = -first_derivative[4];
            first_derivative[ 9] = -first_derivative[3];
            first_derivative[10] = -first_derivative[2];
            first_derivative[11] = -first_derivative[1];
            first_derivative[12] = -first_derivative[0];
            return 0;

        default:
        // if stencil is not implemented print an error message and exit program
            fprintf(stderr,
                "\n (-) Error: First derivative %d-point stencil (1D) not implemented"
                "\n     Aborting..."
                "\n\n"
                , n_stencil
            );
            return (-1);
    }
}


int SecondDerivative(int n_stencil, double* second_derivative){

// fill second derivative stencil
    switch(n_stencil){

        case 3:
            second_derivative[0] =  1.0;

            second_derivative[1] = -2.0;

            second_derivative[2] = second_derivative[0];
            return 0;

        case 5:
            second_derivative[0] =  -1.0/12.0;
            second_derivative[1] =  16.0/12.0;

            second_derivative[2] = -30.0/12.0;

            second_derivative[3] = second_derivative[1];
            second_derivative[4] = second_derivative[0];
            return 0;

        case 7:
            second_derivative[0] =    2.0/180.0;
            second_derivative[1] =  -27.0/180.0;
            second_derivative[2] =  270.0/180.0;

            second_derivative[3] = -490.0/180.0;

            second_derivative[4] = second_derivative[2];
            second_derivative[5] = second_derivative[1];
            second_derivative[6] = second_derivative[0];
            return 0;

        case 9:
            second_derivative[0] =     -9.0/5040.0;
            second_derivative[1] =    128.0/5040.0;
            second_derivative[2] =  -1008.0/5040.0;
            second_derivative[3] =   8064.0/5040.0;

            second_derivative[4] = -14350.0/5040.0;

            second_derivative[5] = second_derivative[3];
            second_derivative[6] = second_derivative[2];
            second_derivative[7] = second_derivative[1];
            second_derivative[8] = second_derivative[0];
            return 0;

        case 11:
            second_derivative[ 0] =      8.0/25200.0;
            second_derivative[ 1] =   -125.0/25200.0;
            second_derivative[ 2] =   1000.0/25200.0;
            second_derivative[ 3] =  -6000.0/25200.0;
            second_derivative[ 4] =  42000.0/25200.0;

            second_derivative[ 5] = -73766.0/25200.0;

            second_derivative[ 6] = second_derivative[4];
            second_derivative[ 7] = second_derivative[3];
            second_derivative[ 8] = second_derivative[2];
            second_derivative[ 9] = second_derivative[1];
            second_derivative[10] = second_derivative[0];
            return 0;

        case 13:
            second_derivative[ 0] =      -50.0/831600.0;
            second_derivative[ 1] =      864.0/831600.0;
            second_derivative[ 2] =    -7425.0/831600.0;
            second_derivative[ 3] =    44000.0/831600.0;
            second_derivative[ 4] =  -222750.0/831600.0;
            second_derivative[ 5] =  1425600.0/831600.0;

            second_derivative[ 6] = -2480478.0/831600.0;

            second_derivative[ 7] = second_derivative[5];
            second_derivative[ 8] = second_derivative[4];
            second_derivative[ 9] = second_derivative[3];
            second_derivative[10] = second_derivative[2];
            second_derivative[11] = second_derivative[1];
            second_derivative[12] = second_derivative[0];
            return 0;

        default:
        // if stencil is not implemented print an error message and exit program
            fprintf(stderr,
                "\n (-) Error: Second derivative %d-point stencil (1D) not implemented"
                "\n     Aborting..."
                "\n\n"
                , n_stencil
            );
            return (-1);
    }
}
