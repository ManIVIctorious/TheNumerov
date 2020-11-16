
#include <stdio.h>
#include <stdlib.h>

// dependencies
double Integrate1D(int n, double dx, double* integrand);
double Integrate2D(int nq1, int nq2, double dx, double* integrand);
double Integrate3D(int nq1, int nq2, int nq3, double dx, double* integrand);
double Integrate4D(int nq1, int nq2, int nq3, int nq4, double dx, double* integrand);


// provided prototypes
double Integrate(int dimension, int* nq, double dx, double* integrand);

double Integrate(int dimension, int* nq, double dx, double* integrand){

    double integral = 0.0;

    switch(dimension){

        case 1:
            integral = Integrate1D(nq[0], dx, integrand);
            break;

        case 2:
            integral = Integrate2D(nq[0], nq[1], dx, integrand);
            break;

        case 3:
            integral = Integrate3D(nq[0], nq[1], nq[2], dx, integrand);
            break;

        case 4:
            integral = Integrate4D(nq[0], nq[1], nq[2], nq[3], dx, integrand);
            break;

        default:
            fprintf(stderr,
                "\n (-) Error: The requested integration routine is not implemented"
                "\n     at the moment, Aborting...\n\n"
            );
            exit(EXIT_FAILURE);
    }

    return integral;
}
