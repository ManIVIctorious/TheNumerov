
#define _GNU_SOURCE
#include <stdio.h>
#include <errno.h>
#include "settings.h"

// dependencies
settings SetDefaultSettings();

// provided prototypes
void Help();


void Help(){

    settings defaults = SetDefaultSettings();

    printf("\n  %s  [OPTIONS] -m <Mode file A> ... -m <Mode file N> -i <Input file>", program_invocation_short_name);

// Description


// Help flag
    printf("\n"
           "\n\t-h, --help"
           "\n\t    Show this help dialogue"
    );

// Threshold
    printf("\n"
           "\n\t-t, --threshold"
           "\n\t    Whenever two floating point numbers are compared, their"
           "\n\t    equality is measured as their absolute difference being"
           "\n\t    smaller than the absolute value of a given threshold."
           "\n\t    This parameter is set by this flag."
           "\n\t    Default:\t%le"
           , defaults.threshold
    );

// Mode file
    printf("\n"
           "\n\t-m|M, --mode-file"
           "\n\t    Set a file containing the vector of a normal mode."
           "\n\t    This file should have the following 4 column structure:"
           "\n\t    The first three columns represent the deviation of each"
           "\n\t    atom relative to the minimum geometry, while the 4th"
           "\n\t    column contains the respective particle mass."
           "\n\t    This flag can be used multiple times (e.g. two times for"
           "\n\t    two, or three times for three dimensional problems)."
    );

// Input file
    printf("\n"
           "\n\t-i, --input-file"
           "\n\t    The reciprocal effective moment of inertia has to be determined"
           "\n\t    for every single grid point. Therefore, the system geometry of"
           "\n\t    each grid point as well as the deviation from minimum, leading"
           "\n\t    to this geometry has to be known. Instead of setting these"
           "\n\t    values on the command line, a input file is used containing"
           "\n\t    <dimension> + 1 columns. The first columns represent the deviation"
           "\n\t    from the minimum geometry alongside the normal modes A to N"
           "\n\t    and the last column is a path to the respective input files (e.g."
           "\n\t    .com files) containing the geometry."
    );

// Pipe
    printf("\n"
           "\n\t-P, --pipe-input"
           "\n\t    Instead of providing an input file the exact same structure can"
           "\n\t    be created on the command line and piped to the program."
           "\n\t    To indicate that the input should be read from stdin instead of"
           "\n\t    a dedicated input file this flag can be used."
    );

// Output file
    printf("\n"
           "\n\t-o, --output-file"
           "\n\t    Name of the output file. The output file is structured in two"
           "\n\t    blocks: first three lines beginning with Zeta_{x,y,z} and"
           "\n\t    ( <dimension>*<dimension> - <dimension> ) / 2 data columns."
           "\n\t    Theses columns represent the Coriolis coefficients between two"
           "\n\t    modes. Only the upper triangle without main diagonal is printed"
           "\n\t    since zeta(AA) = 0 and zeta(AB) = -zeta(BA)."
           "\n\t    The second block are <n_grid_points> lines with <dimension> + 6"
           "\n\t    columns. The first <dimension> columns contain the deviations as"
           "\n\t    given in input file and the latter 6 columns contain the upper"
           "\n\t    triangle (with main diagonal) of the reciprocal effective moment"
           "\n\t    of inertia tensor. This tensor is symmetric therefore, only the"
           "\n\t    upper triangle is printed."
    );

    printf("\n\n");
}
