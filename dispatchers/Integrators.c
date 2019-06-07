
#include <stdio.h>
#include <stdlib.h>

// dependencies
double integrate_1d(int n, double dx, double* integrand);
double integrate_2d(int nq1, int nq2, double dx, double* integrand);
double integrate_3d(int nq1, int nq2, int nq3, double dx, double* integrand);
double integrate_4d(int nq1, int nq2, int nq3, int nq4, double dx, double* integrand);


// provided prototypes
double Integrate(int dimension, int* nq, double dx, double* integrand);

double Integrate(int dimension, int* nq, double dx, double* integrand){

    double integral = 0.0;

    switch(dimension){

        case 1:
            integral = integrate_1d(nq[0], dx, integrand);
            break;

        case 2:
            integral = integrate_2d(nq[0], nq[1], dx, integrand);
            break;

        case 3:
            integral = integrate_3d(nq[0], nq[1], nq[2], dx, integrand);
            break;

        case 4:
            integral = integrate_4d(nq[0], nq[1], nq[2], nq[3], dx, integrand);
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
