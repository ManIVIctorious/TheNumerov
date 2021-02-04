
#include <stdio.h>
#include <stdlib.h>

// Dependencies
int FillStencil1D(double* stencil, int n_stencil);
int FillStencil2D(double* stencil, int n_stencil);
int FillStencil3D(double* stencil, int n_stencil);
int FillStencil4D(double* stencil, int n_stencil);


// provided prototypes
void MetaGetStencil(double* stencil, int n_stencil, int dimension);

void MetaGetStencil(double* stencil, int n_stencil, int dimension){

// Assignment of stencil functions
    switch(dimension){

        case 1:
            FillStencil1D(stencil, n_stencil);
            break;

        case 2:
            FillStencil2D(stencil, n_stencil);
            break;

        case 3:
            FillStencil3D(stencil, n_stencil);
            break;

        case 4:
            FillStencil4D(stencil, n_stencil);
            break;

        default:
            fprintf(stderr,
                "\n (-) Error: Currently there is no stencil of dimension %d implemented"
                "\n     Aborting...\n\n"
                , dimension
            );
            exit(EXIT_FAILURE);
            break;
    }

}
