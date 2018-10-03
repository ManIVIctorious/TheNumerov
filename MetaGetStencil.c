
#include <stdio.h>
#include <stdlib.h>

// Dependencies
int FillStencil1D(double* stencil, int n_stencil);
int FillStencil2D(double* stencil, int n_stencil);


// provided prototypes
int MetaGetStencil(double* stencil, int n_stencil, int dimension);

int MetaGetStencil(double* stencil, int n_stencil, int dimension){

    int control = -1;

// Assignment of stencil functions
    switch(dimension){

        case 1:
            control = FillStencil1D(stencil, n_stencil);
            break;

        case 2:
            control = FillStencil2D(stencil, n_stencil);
            break;

        default:
            fprintf(stderr,
                "\n (-) Error currently there is no stencil of dimension %d implemented"
                "\n     Aborting...\n\n"
                , dimension
            );
            exit(EXIT_FAILURE);
            break;
    }

    return control;
}
