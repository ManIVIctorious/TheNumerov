
#define _GNU_SOURCE
#include <stdio.h>
#include <errno.h>
#include "settings.h"
#include "gitversion.h"

// dependencies
settings SetDefaultSettings(void);

// provided prototypes
void usage(void);


void usage(void){

// first get default settings to print them alongside the individual options
    settings defaults = SetDefaultSettings();

// Usage
    printf("\n  %s  [OPTIONS] -m <Mode file A> ... -m <Mode file N> -i <Input file>", program_invocation_short_name);


// Description
    printf("\n"
           "\n  Calculate the system's Coriolis-coefficients \"zeta\" from the mode files"
           "\n  and the reciprocal moment of inertia tensor for each configuration."
           "\n"
           "\n  The program takes a number of mode files (all prepended with -m|-M) as"
           "\n  input (this also sets the problem's dimension) and calculates the zeta"
           "\n  values. Due to the identities"
           "\n                            zeta_ii =    0"
           "\n                            zeta_ij = -zeta_ji"
           "\n  only the upper triangle (without the main diagonal) of the zeta values"
           "\n  is stored in the output file."
           "\n  Furthermore, the program takes an input file containing <dimension> columns"
           "\n  representing the deviation from the minimum geometry and an additional"
           "\n  column with the file-path to a configuration file. This file will then"
           "\n  be processed yielding the upper triangle (including the main diagonal)"
           "\n  of the reciprocal moment of inertia tensor."
           "\n\n  Version: \"%s\"\n", gitversion
    );


// Help flag
    printf("\n"
           "\n\t-h, --help"
           "\n\t    Show this help dialogue"
    );

// Threshold
    printf("\n"
           "\n\t-t, --threshold"
           "\n\t    Whenever two floating point numbers are compared, their equality is"
           "\n\t    measured as their absolute difference being smaller than the absolute"
           "\n\t    value of a given threshold. This parameter is set by this flag."
           "\n\t    Default:\t%le"
           , defaults.threshold
    );

// Mode file
    printf("\n"
           "\n\t-m|M, --mode-file"
           "\n\t    Set a file containing the vector of a normal mode."
           "\n\t    This file should have the following 4 column structure:"
           "\n\t    The first three columns represent the deviation of each atom relative"
           "\n\t    to the minimum geometry, while the 4th column contains the respective"
           "\n\t    particle mass. This flag can be used multiple times (e.g. two times"
           "\n\t    for two, or three times for three dimensional problems)."
    );

// Input file
    printf("\n"
           "\n\t-i, --input-file"
           "\n\t    The reciprocal effective moment of inertia has to be determined for"
           "\n\t    every single grid point. Hence, the system geometry of each grid point"
           "\n\t    as well as the deviation from minimum, leading to this geometry has to"
           "\n\t    be known."
           "\n\t    Instead of setting these values on the command line, a input file is"
           "\n\t    used containing <dimension> + 1 columns. The first <dimension> columns"
           "\n\t    represent the deviation from the minimum geometry alongside the normal"
           "\n\t    modes A to N and the last column is a path to the respective input"
           "\n\t    files (e.g. *.com files) containing the geometry."
    );

// Pipe
    printf("\n"
           "\n\t-P, --pipe-input"
           "\n\t    Expect the input file being provided from stdin."
    );

// Output file
    printf("\n"
           "\n\t-o, --output-file"
           "\n\t    Name of the output file. The output file is structured in two blocks:"
           "\n\t    First three lines beginning with Zeta_{x,y,z} and followed by"
           "\n\t    (<dimension>*<dimension> - <dimension>) / 2 data columns."
           "\n\t    Theses columns represent the Coriolis coefficients between two modes"
           "\n\t    (upper triangle without main diagonal)."
           "\n\t    The second block are <n_grid_points> lines with <dimension> + 6 columns"
           "\n\t    The first <dimension> columns represent the deviations as given in"
           "\n\t    input file and the latter 6 columns contain the upper triangle (with"
           "\n\t    main diagonal) of the reciprocal effective moment of inertia tensor."
    );

    printf("\n\n");
}
