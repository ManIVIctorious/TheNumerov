
#include <stdio.h>

// Offered prototypes
int Help(char *filename);


int Help(char *filename){

    printf("\n%s\t [OPTIONS] -i INPUT-FILENAME [-o OUTPUT-FILENAME]", filename);

// Help flag
    printf("\n");
    printf("\n\t"); printf("-h, --help");
    printf("\n\t"); printf("    Show this help dialogue");

// Numerov
    printf("\n");
    printf("\n\t"); printf("-m, --mass");
    printf("\n\t"); printf("    Set reduced mass of involved modes in g/mol");

    printf("\n");
    printf("\n\t"); printf("-k, --fkin");
    printf("\n\t"); printf("    Conversion factor from kJ/mol to desired output energy unit");
    printf("\n\t"); printf("    e.g. 1.0/4.184 (kcal/kJ) gives output in kcal/mol (default)");
    printf("\n\t"); printf("    assuming the mass input is given in g/mol and the coordinate");
    printf("\n\t"); printf("    input is given in angstrom");

    printf("\n");
    printf("\n\t"); printf("-v, --fpot");
    printf("\n\t"); printf("    Conversion factor from the input unit of energy to the desired output");
    printf("\n\t"); printf("    unit of energy, e.g. 2625.49962 (kJ/mol)/Hartree.");
    printf("\n\t"); printf("    Be aware that for output energy units different from kcal/mol the");
    printf("\n\t"); printf("    kinetic energy factor [-k, -fkin] has to be set too.");
    printf("\n\t"); printf("    E.g. above example requires the option \"-k 1.0\".");

    printf("\n");
    printf("\n\t"); printf("-M, --fmu");
    printf("\n\t"); printf("    Conversion factor for the \"effective reciprocal inertia tensor\".");
    printf("\n\t"); printf("    The conversion is performed from the input dimension to kJ/mol,");
    printf("\n\t"); printf("    further conversion to the desired output unit of energy is handled");
    printf("\n\t"); printf("    by the kinetic energy factor [-k, --fkin].");

    printf("\n");
    printf("\n\t"); printf("-n, --n-stencil");
    printf("\n\t"); printf("    Set the one dimensional stencil size, e.g. \"-n 11\" gives a");
    printf("\n\t"); printf("    11x11 stencil in the two dimensional Numerov.");

#ifdef HAVE_OPT_SPLINE
    printf("\n");
    printf("\n\t"); printf("-s, --spline");
    printf("\n\t"); printf("    Enable spline interpolation of potential points and set the Number");
    printf("\n\t"); printf("    of points interpolated between each point-pair.");
    printf("\n\t"); printf("    A spline of 1 doubles the number of data points for each dimension,");
    printf("\n\t"); printf("    n_points times (spline_points + 1)**dimension.");
#endif

// I/O
    printf("\n");
    printf("\n"); printf("    Input/Output options");

    printf("\n");
    printf("\n\t"); printf("-i, --input-file");
    printf("\n\t"); printf("    Path to input file. The program expects the following structure of input:");
    printf("\n\t"); printf("    <dimension> columns containing the coordinates and one column containing");
    printf("\n\t"); printf("    the potential energy values. (e.g. 2D expects 3 columns)");

    printf("\n");
    printf("\n\t"); printf("-d, --dipole");
    printf("\n\t"); printf("    Expect 3 additional columns containing the dipole moment in {x,y,z} direction");
    printf("\n\t"); printf("    after the potential energy (e.g. 2D: q1, q2, potential, dipolemoment{x,y,z})");

    printf("\n");
    printf("\n\t"); printf("-P, --pipe");
    printf("\n\t"); printf("    Read input from standard input instead of an input file");

    printf("\n");
    printf("\n\t"); printf("-t, --dq-threshold");
    printf("\n\t"); printf("    The multi dimensional Numerov procedure requires a equi spaced grid");
    printf("\n\t"); printf("    This option sets the maximal variation of the grid spacing.");
    printf("\n\t"); printf("    Default is \"-t 1.0E-12\"");

    printf("\n");
    printf("\n\t"); printf("-a, --analyze");
    printf("\n\t"); printf("    Perform additional calculations, giving insight to");
    printf("\n\t"); printf("    Orthonormality, Potential, kinetic energy and coupling");

    printf("\n");
    printf("\n\t"); printf("-c, --coriolis-input");
    printf("\n\t"); printf("    Path to an optional file containing Coriolis coefficients and the \"effective");
    printf("\n\t"); printf("    reciprocal inertia tensor\" for each grid point to allow for the calculation of");
    printf("\n\t"); printf("    rotational terms.");

    printf("\n");
    printf("\n\t"); printf("-o, --output-file");
    printf("\n\t"); printf("    Set the path to the output file.");

// Matrix solver
#ifdef HAVE_MKL_INSTALLED
    printf("\n");
    printf("\n"); printf("    Options for the Intel Math Kernel Library FEAST eigensolver:");

    printf("\n");
    printf("\n\t"); printf("-l, --lower-bound       Set lower energy bound of calculated eigenstates");
    printf("\n\t"); printf("-u, --upper-bound       Set upper energy bound of calculated eigenstates");
#endif

    printf("\n\n");
    return 0;
}
