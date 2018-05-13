
#include <stdio.h>
#include "typedefinitions.h"

// provided prototypes
int Help(char* filename, settings defaults);


int Help(char* filename, settings defaults){

    printf("\n%s\t [OPTIONS] -i INPUT-FILENAME [-o OUTPUT-FILENAME]", filename);

// Help flag
    printf("\n"
           "\n\t-h, --help"
           "\n\t    Show this help dialogue"
    );

// Numerov
    printf("\n"
           "\n\t-m, --mass"
           "\n\t    Set reduced mass of involved modes in g/mol"
           "\n\t    For dimensions > 1 this should always be set to 1.0"
           "\n\t    Default:\t% lf", defaults.mass
    );

    printf("\n"
           "\n\t-k, --fkin"
           "\n\t    Conversion factor from kJ/mol to desired output energy unit"
           "\n\t    e.g. 1.0/4.184 (kcal/kJ) gives output in kcal/mol (default)"
           "\n\t    assuming the mass input is given in g/mol and the coordinate"
           "\n\t    input is given in angstrom"
           "\n\t    Default:\t% lf", defaults.ekin_factor
    );

    printf("\n"
           "\n\t-v, --fpot"
           "\n\t    Conversion factor from the input unit of energy to the desired output"
           "\n\t    unit of energy, e.g. 2625.49962 (kJ/mol)/Hartree."
           "\n\t    Be aware that for output energy units different from kcal/mol the"
           "\n\t    kinetic energy factor [-k, -fkin] has to be set too."
           "\n\t    E.g. above example requires the option \"-k 1.0\"."
           "\n\t    Default:\t% lf", defaults.epot_factor
    );

    printf("\n"
           "\n\t-M, --fmu"
           "\n\t    Conversion factor for the \"effective reciprocal inertia tensor\"."
           "\n\t    The conversion is performed from the input dimension to kJ/mol,"
           "\n\t    further conversion to the desired output unit of energy is handled"
           "\n\t    by the kinetic energy factor [-k, --fkin]."
           "\n\t    Default:\t% lf", defaults.mu_factor
    );

    printf("\n"
           "\n\t-n, --n-stencil"
           "\n\t    Set the one dimensional stencil size, e.g. \"-n 11\" gives a"
           "\n\t    11x11 stencil in the two dimensional Numerov."
           "\n\t    Default:\t%d", defaults.n_stencil
    );

    printf("\n"
           "\n\t-D, --dimension"
           "\n\t    Set the dimensionality of the applied Numerov procedure"
           "\n\t    Default:\t%d", defaults.dimension
    );

    printf("\n"
           "\n\t-s, --spline"
           "\n\t    Enable spline interpolation of potential points and set the Number"
           "\n\t    of points interpolated between each point-pair."
           "\n\t    A spline of 1 effectively doubles the number of data points for each"
           "\n\t    dimension, n_points times (spline_points + 1)**dimension."
           "\n\t    Default:\t%d", defaults.n_spline
    );

// I/O
    printf("\n");
    printf("\n"); printf("    Input/Output options");

    printf("\n"
           "\n\t-i, --input-file"
           "\n\t    Path to input file. The program expects the following structure of input:"
           "\n\t    <dimension> columns containing the coordinates and one column containing"
           "\n\t    the potential energy values. (e.g. 2D expects 3 columns)"
    );

    printf("\n"
           "\n\t-d, --dipole"
           "\n\t    Expect 3 additional columns containing the dipole moment in {x,y,z} direction"
           "\n\t    after the potential energy (e.g. 2D: q1, q2, potential, dipolemoment{x,y,z})"
           "\n\t    Default:\t%s", defaults.dipole ? "true" : "false"
    );

    printf("\n"
           "\n\t-P, --pipe"
           "\n\t    Read input from standard input instead of an input file"
           "\n\t    Default:\t%s", "false"
    );

    printf("\n"
           "\n\t-t, --dq-threshold"
           "\n\t    The multi dimensional Numerov procedure requires a equi spaced grid"
           "\n\t    This option sets the maximal variation of the grid spacing."
           "\n\t    Default:\t% le", defaults.threshold
    );

    printf("\n"
           "\n\t-T, --no-spacing-check"
           "\n\t    Disable the coordinate spacing check completely"
           "\n\t    This option does not disable the check if the Coriolis files"
           "\n\t    contain the same coordinates as the standard input file,"
           "\n\t    so threshold (-t) may still be needed"
           "\n\t    Default:\t%s", defaults.check_spacing ? "false" : "true"
    );

    printf("\n"
           "\n\t-a, --analyze"
           "\n\t    Perform additional calculations, giving insight to"
           "\n\t    Orthonormality, Potential, kinetic energy and coupling"
           "\n\t    Default:\t%s", defaults.analyze ? "true" : "false"
    );

    printf("\n"
           "\n\t-c, --coriolis-input"
           "\n\t    Path to an optional file containing Coriolis coefficients and the \"effective"
           "\n\t    reciprocal inertia tensor\" for each grid point to allow for the calculation of"
           "\n\t    rotational terms."
    );

    printf("\n"
           "\n\t-o, --output-file"
           "\n\t    Set the path to the output file."
    );

// Matrix solver
    printf("\n\n");
    printf("    Matrix Eigen-solvers:");

    printf("\n"
           "\n\tThis program includes multiple implementations for the calculation"
           "\n\tof the matrix eigenstates. These implementations require different"
           "\n\ttypes of information and therefore different flags apply for their"
           "\n\tsettings. The following lines show these individual flags for the"
           "\n\trespective calculation methods"
    );

#ifdef HAVE_MKL_INSTALLED
    printf("\n\n");
    printf("\tOptions for the Intel Math Kernel Library FEAST eigen solver:");

    printf("\n"
           "\n\t    --mkl               Use the Intel Math Kernel Library FEAST eigen solver"
           "\n\t    -l, --lower-bound   Set lower energy bound of calculated eigenstates"
           "\n\t    -u, --upper-bound   Set upper energy bound of calculated eigenstates"
    );
#endif

#ifdef HAVE_ARMA_INSTALLED
    printf("\n\n");
    printf("\tOptions for the Armadillo ARPACK eigensolver:");

    printf("\n"
           "\n\t    --armadillo         Use the Armadillo ARPACK eigen solver"
           "\n\t    -N, --nout          Set number of eigenstates to be calculated"
    );
#endif

    printf("\n\n");
    return 0;
}
