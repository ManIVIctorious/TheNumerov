
#define _GNU_SOURCE
#include <stdio.h>
#include <errno.h>
#include "settings.h"
#include "gitversion.h"

// provided prototypes
void usage(void);

// dependencies
settings SetDefaultSettings();

void usage(void){

    settings defaults = SetDefaultSettings();

    printf("\n  %s\t[OPTIONS] -i INPUT-FILENAME [-o OUTPUT-FILENAME]", program_invocation_short_name);

// Description
    printf("\n"
           "\n    Numerical solver for the time independent Schrödinger equation based on the"
           "\n    well known Numerov procedure for the calculation of second order differential"
           "\n    without any terms of first order."
           "\n    This implementation is specifically designed for the use in computational chemistry"
           "\n    enhancing the original procedure by the inclusion of higher dimensions as well as"
           "\n    the consideration of the molecular Hamiltonian first described by Watson."
           "\n\n    Version: \"%s\"\n", gitversion
    );

// Help flag
    printf("\n"
           "\n\t-h, --help"
           "\n\t    Show this help dialogue"
    );

// Control file
    printf("\n"
           "\n\t-C, --control-file"
           "\n\t    Define a control file, containing keyword-value pairs of the form"
           "\n\t    <keyword>=<value>;. Used command line flags will always overwrite"
           "\n\t    the settings presented in the control file."
           "\n\t    This allows to use the control file as basis to pass local defaults,"
           "\n\t    but preserves the flexibility introduced with the command line flags."
    );

// Numerov
    printf("\n"
           "\n\t-m, --masses"
           "\n\t    Set reduced mass of involved modes in g/mol"
           "\n\t    This setting is only required if the input is not given in mass weighted"
           "\n\t    coordinates. For dimensions greater than 1 the reduced masses are passed"
           "\n\t    as a colon separated array of reduced masses."
           "\n\t    The default values can be restored by the string \"default\"."
           "\n\t    Keyword:\tReduced_Masses"
           "\n\t    Default:\t1.0:1.0"
    );

    printf("\n"
           "\n\t-k, --fkin"
           "\n\t    Conversion factor from kJ/mol to desired output energy unit"
           "\n\t    e.g. 1.0/4.184 (kcal/kJ) gives output in kcal/mol (default)"
           "\n\t    assuming the mass input is given in g/mol and the coordinate"
           "\n\t    input is given in angstrom"
           "\n\t    Keyword:\tkJpermolToOUE"
           "\n\t    Default:\t% le", defaults.kJpermol_to_oue
    );

    printf("\n"
           "\n\t-v, --fpot"
           "\n\t    Conversion factor from the input unit of energy to the desired output"
           "\n\t    unit of energy, e.g. 2625.49962 (kJ/mol)/Hartree."
           "\n\t    Be aware that for output energy units different from kcal/mol the"
           "\n\t    kinetic energy factor [-k, -fkin] has to be set too."
           "\n\t    E.g. above example requires the option \"-k 1.0\"."
           "\n\t    Keyword:\tEPotToOUE"
           "\n\t    Default:\t% le", defaults.epot_to_oue
    );

    printf("\n"
           "\n\t-f, --fdipole"
           "\n\t    Conversion factor of the \"dipole moment\"."
           "\n\t    The conversion is performed from the input dimension to Asm"
           "\n\t    (Ampere Second Meter), further conversions are handled by the"
           "\n\t    kinetic energy factor [-k, --fkin]."
           "\n\t    Keyword:\tDipToAsm"
           "\n\t    Default:\t% le", defaults.dip_to_Asm
    );

    printf("\n"
           "\n\t-M, --fmu"
           "\n\t    Conversion factor for the \"effective reciprocal inertia tensor\"."
           "\n\t    The conversion is performed from the input dimension to mol/(g.Å^2)."
           "\n\t    Keyword:\tIMOITomolpergAasq"
           "\n\t    Default:\t% le", defaults.InvInertia_to_molpergAasq
    );

    printf("\n"
           "\n\t-w, --watson-threshold"
           "\n\t    Threshold for the entries of the \"effective reciprocal moment of"
           "\n\t    inertia tensor\". When set to a positive number all entries outside"
           "\n\t    the region [-threshold,threshold] are set to the appropriate"
           "\n\t    boundary. Setting it to a negative number disables the threshold."
           "\n\t    Keyword:\tIMOI_Threshold"
           "\n\t    Default:\t% le", defaults.InvInertiaThreshold
    );

    printf("\n"
           "\n\t-n, --n-stencil"
           "\n\t    Set the one dimensional stencil size, e.g. \"-n 11\" gives a"
           "\n\t    11x11 stencil in the two dimensional Numerov."
           "\n\t    Keyword:\tStencil_Size"
           "\n\t    Default:\t%d", defaults.n_stencil
    );

    printf("\n"
           "\n\t-D, --dimension"
           "\n\t    Set the dimensionality of the applied Numerov procedure"
           "\n\t    Keyword:\tDimensionality"
           "\n\t    Default:\t%d", defaults.dimension
    );

    printf("\n"
           "\n\t-s, --spline"
           "\n\t    Enable spline interpolation of potential points and set the Number"
           "\n\t    of points interpolated between each point-pair."
           "\n\t    A spline of 1 effectively doubles the number of data points for each"
           "\n\t    dimension, n_points times (spline_points + 1)**dimension."
           "\n\t    Keyword:\tInterpolation_points"
           "\n\t    Default:\t%d", defaults.n_spline
    );

// I/O
    printf("\n");
    printf("\n"); printf("    Input/Output options");

    printf("\n"
           "\n\t-i, --input-file"
           "\n\t    Path to input file. The program expects the following structure of input:"
           "\n\t    <dimension> columns containing the coordinates in Ångstrom and one column"
           "\n\t    containing the potential energy values. (e.g. 2D expects 3 columns)"
           "\n\t    Keyword:\tInput_File"
    );

    printf("\n"
           "\n\t-d, --dipole"
           "\n\t    Expect 3 additional columns containing the dipole moment in {x,y,z} direction"
           "\n\t    after the potential energy (e.g. 2D: q1, q2, potential, dipolemoment{x,y,z})"
           "\n\t    Keyword:\tDipole"
           "\n\t    Default:\t%s", defaults.dipole ? "true" : "false"
    );

    printf("\n"
           "\n\t-e, --ext-dipole-file"
           "\n\t    Instead of reading the dipole moments from the primary input file, they can"
           "\n\t    also be passed on in an explicit external file consisting of <dimension> columns"
           "\n\t    containing the coordinates and additional three columns containing the dipole"
           "\n\t    moment in x, y and z direction (e.g. 2D: q1, q2, dip_x, dip_y, dip_z)."
           "\n\t    Keyword:\tExternal_Dipole_File"
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
           "\n\t    Keyword:\tSpacing_Threshold"
           "\n\t    Default:\t% le", defaults.threshold
    );

    printf("\n"
           "\n\t-T, --no-spacing-check"
           "\n\t    Disable the check for equi-distant coordinate spacing."
           "\n\t    This option does not disable the coordinate comparison between"
           "\n\t    different input files. In this case setting the threshold"
           "\n\t    (-t argument) to a sufficiently high value could be required."
           "\n\t    Keyword:\tCheck_Spacing"
           "\n\t    Default:\t%s", defaults.check_spacing ? "false" : "true"
    );

    printf("\n"
           "\n\t-p, --periodic"
           "\n\t    For periodic problems (e.g. rotations) use a periodic implementation"
           "\n\t    of the matrix filling routine."
           "\n\t    Keyword:\tPeriodic"
           "\n\t    Default:\t%s", defaults.periodic ? "true" : "false"
    );

    printf("\n"
           "\n\t-F, --frequencies"
           "\n\t    Calculate the energy differences between all eigenvalues in wavenumbers"
           "\n\t    Keyword:\tFrequencies"
           "\n\t    Default:\t%s", defaults.frequencies ? "true" : "false"
    );

    printf("\n"
           "\n\t-a, --analyse"
           "\n\t    Perform additional calculations, giving insight to"
           "\n\t    Orthonormality, Potential, kinetic energy and coupling"
           "\n\t    Keyword:\tAnalyse"
           "\n\t    Default:\t%s", defaults.analyse ? "true" : "false"
    );

    printf("\n"
           "\n\t-c, --coriolis-input"
           "\n\t    Path to an optional file containing Coriolis coefficients and the \"effective"
           "\n\t    reciprocal inertia tensor\" for each grid point to allow for the calculation of"
           "\n\t    rotational terms."
           "\n\t    Keyword:\tCoriolis_File"
    );

    printf("\n"
           "\n\t-o, --output-file"
           "\n\t    Set the path to the output file."
           "\n\t    Keyword:\tOutput_File"
           "\n\t    Default:\t%s", defaults.output_file
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
    printf("\tOptions for the Intel Math Kernel Library FEAST eigensolver:");

    printf("\n"
           "\n\t    Flags               Keywords          Description"
           "\n\t    --mkl               Intel_MKL_FEAST   Use the Intel Math Kernel Library FEAST eigensolver"
           "\n\t    -l, --lower-bound   Lower_Bound       Set lower energy bound of calculated eigenstates"
           "\n\t    -u, --upper-bound   Upper_Bound       Set upper energy bound of calculated eigenstates"
    );
#endif

#ifdef HAVE_ARMA_INSTALLED
    printf("\n\n");
    printf("\tOptions for the Armadillo ARPACK eigensolver:");

    printf("\n"
           "\n\t    Flags               Keywords          Description"
           "\n\t    --armadillo         Armadillo_ARPACK  Use the Armadillo ARPACK eigensolver"
           "\n\t    -N, --n-out         N_Eigenstates     Set number of eigenstates to be calculated"
    );
#endif

    printf("\n\n");
}
