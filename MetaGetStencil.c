
#include <stdio.h>
#include <stdlib.h>

// Dependencies
int FillStencil1D(double* stencil, int n_stencil);
int FillStencil2D(double* stencil, int n_stencil);


// Offered prototypes
int MetaGetStencil(double* stencil, int n_stencil, int dimension);

int MetaGetStencil(double* stencil, int n_stencil, int dimension){

    int control = -1;

// dimension has to be greater than zero
    if(dimension <= 0){
        fprintf(stderr,
            "\n (-) Error dimension has to be greater than zero"
            "\n     Aborting..."
            "\n\n"
        );
        exit(-1);
    }

// stencil has to have an odd number of entries (at least in current implementation)
    if(n_stencil%2 == 0){
        fprintf(stderr,
            "\n (-) Stencil size is given as even (%d), but must be an odd number."
            "\n     Aborting - please check your input..."
            "\n\n"
            , n_stencil
        );
        exit(-1);
    }


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
                "\n     Aborting..."
                "\n\n"
                , dimension
            );
            exit(-1);
            break;
    }

    return control;
}
