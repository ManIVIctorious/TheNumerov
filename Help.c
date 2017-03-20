
#include <stdio.h>

// Offered prototypes
int Help(char *filename);


int Help(char *filename){

    printf("\n%s\t [OPTIONS] -i INPUT-FILENAME -o OUTPUT-FILENAME", filename);
    printf("\n\t");
    printf("\n\t"); printf("-h, --help              Show this help dialogue");

// Numerov
    printf("\n");
    printf("\n\t"); printf("-m, --mass              Set reduced mass of involved modes in g/mol");
    printf("\n\t"); printf("-k, --fkin              Conversion factor from kJ/mol to desired output energy unit");
    printf("\n\t"); printf("                        e.g. 1.0/4.184 (kcal/kJ) gives output in kcal/mol (default)");
    printf("\n\t"); printf("-v, --fpot              Conversion factor input -> output unit of energy");
    printf("\n\t"); printf("-n, --n-stencil         Set stencil size");

// Matrix solver
    printf("\n");
    printf("\n\t"); printf("-l, --lower-bound       Set lower energy bound of calculated eigenstates");
    printf("\n\t"); printf("-u, --upper-bound       Set upper energy bound of calculated eigenstates");
    printf("\n\t"); printf("-s, --spline            Set number of interpolated points");

// I/O
    printf("\n");
    printf("\n\t"); printf("-t, --dq-threshold      Maximal value of abs(q[i]-q[i+1])");
    printf("\n\t"); printf("-a, --analyze           Output additional information containing");
    printf("\n\t"); printf("                            Orthonormality, Potential, kinetic energy and coupling");
    printf("\n\t"); printf("-d, --dipole            Expect 6 columns of input (q1, q2, potential, dipolemoment{x,y,z})");
    printf("\n\t"); printf("-P, --pipe              Read input from pipe instead of input file");
    printf("\n\t"); printf("-i, --input-file        Set name of input file");
    printf("\n\t"); printf("-o, --output-file       Set name of output file");

    printf("\n\n");
    return 0;
}
