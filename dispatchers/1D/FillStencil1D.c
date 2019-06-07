
#include <stdio.h>
#include <stdlib.h>

// provided prototypes
int FillStencil1D(double* stencil, int n_stencil);
 
int FillStencil1D(double* stencil, int n_stencil){

// fill the stencil array
    switch(n_stencil){
        case 3:
            stencil[0] =  1.0;
            stencil[1] = -2.0;
            stencil[2] =  1.0;
            return 0;

        case 5:
            stencil[0] =   -1.0/12.0;
            stencil[1] =   16.0/12.0;
            stencil[2] =  -30.0/12.0;
            stencil[3] =   16.0/12.0;
            stencil[4] =   -1.0/12.0;
            return 0;

        case 7:
            stencil[0] =    2.0/180.0;
            stencil[1] =  -27.0/180.0;
            stencil[2] =  270.0/180.0;
            stencil[3] = -490.0/180.0;
            stencil[4] =  270.0/180.0;
            stencil[5] =  -27.0/180.0;
            stencil[6] =    2.0/180.0;
            return 0;

        case 9:
            stencil[0] =     -9.0/5040.0;
            stencil[1] =    128.0/5040.0;
            stencil[2] =  -1008.0/5040.0;
            stencil[3] =   8064.0/5040.0;
            stencil[4] = -14350.0/5040.0;
            stencil[5] =   8064.0/5040.0;
            stencil[6] =  -1008.0/5040.0;
            stencil[7] =    128.0/5040.0;
            stencil[8] =     -9.0/5040.0;
            return 0;

        case 11:
            stencil[0]  =      8.0/25200.0;
            stencil[1]  =   -125.0/25200.0;
            stencil[2]  =   1000.0/25200.0;
            stencil[3]  =  -6000.0/25200.0;
            stencil[4]  =  42000.0/25200.0;
            stencil[5]  = -73766.0/25200.0;
            stencil[6]  =  42000.0/25200.0;
            stencil[7]  =  -6000.0/25200.0;
            stencil[8]  =   1000.0/25200.0;
            stencil[9]  =   -125.0/25200.0;
            stencil[10] =      8.0/25200.0;
            return 0;

        case 13:
            stencil[0]  =      -50.0/831600.0;
            stencil[1]  =      864.0/831600.0;
            stencil[2]  =    -7425.0/831600.0;
            stencil[3]  =    44000.0/831600.0;
            stencil[4]  =  -222750.0/831600.0;
            stencil[5]  =  1425600.0/831600.0;
            stencil[6]  = -2480478.0/831600.0;
            stencil[7]  =  1425600.0/831600.0;
            stencil[8]  =  -222750.0/831600.0;
            stencil[9]  =    44000.0/831600.0;
            stencil[10] =    -7425.0/831600.0;
            stencil[11] =      864.0/831600.0;
            stencil[12] =      -50.0/831600.0;
            return 0;

        default:
        // if stencil is not implemented print an error message and exit program
            fprintf(stderr,
                "\n (-) Error no data for %d-point stencil available."
                "\n     Aborting - please check your input...\n\n"
                , n_stencil
            );
            exit(EXIT_FAILURE);
    }
}
