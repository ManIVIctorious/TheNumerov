
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
           "\n\t-m, --mode-file"
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

// Masses file
    printf("\n"
           "\n\t-M, --masses-file"
           "\n\t    Instead of expecting a fourth column in the mode files representing the"
           "\n\t    particle mass just read the first three columns and get the masses from"
           "\n\t    the provided masses file."
    );

// Coriolis coefficients from command line
    printf("\n"
           "\n\t-X, --zeta-x"
           "\n\t-Y, --zeta-y"
           "\n\t-Z, --zeta-z"
           "\n\t    Instead of calculating the Coriolis coefficients on each run from the"
           "\n\t    mode files, they can also be provided as a colon separated list, where"
           "\n\t    the fields are the upper triangle of the coefficient matrices without"
           "\n\t    the main diagonal, i.e."
           "\n\t            a b c        j k l        s t u  --zeta_x b:c:f"
           "\n\t      z_x = d e f  z_y = m n o  z_z = v w x  --zeta_y k:l:o"
           "\n\t            g h i        p q r        y z 1  --zeta_z t:u:x"
    );

// Pipe
    printf("\n"
           "\n\t-P, --pipe-input"
           "\n\t    Expect the structure of the main input file to be provided from stdin."
           "\n\t    In combination with masses file, coriolis coefficients and append mode" "\n\t    this can be used to process each geometry file individually, e.g."
           "\n"
           "\n\t      for i in $(seq $start_1 $dq_1 $stop_1); do"
           "\n\t      for j in $(seq $start_2 $dq_2 $stop_2); do"
           "\n\t      for k in $(seq $start_3 $dq_3 $stop_3); do"
           "\n"
           "\n\t        file=\"scan_3d_${i}_${j}_${k}.com\""
           "\n\t        echo \"$i $j $k $file\" | TheWatson -M masses\\"
           "\n\t                                          -X b:c:f -Y k:l:o -Z t:u:x\\"
           "\n\t                                          -a Watson.out -P --"
           "\n"
           "\n\t      done"
           "\n\t      done"
           "\n\t      done"

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

    printf("\n"
           "\n\t-a, --append-file"
           "\n\t    Same as -o, but if the set file already exists append to it instead of"
           "\n\t    re-writing it."
           "\n\t    The output routine still checks if the given file exists and if not"
           "\n\t    writes headers/keys and Coriolis coefficients, so these are only once"
           "\n\t    present in the generated output file."
    );

    printf("\n\n");
}
