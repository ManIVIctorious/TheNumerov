#ifndef _SETTINGS_H
#define _SETTINGS_H

typedef struct settings {

// files
// the main input file consists of dimension + 1 columns:
//  the first <dimension> columns contain the deviations from the minimum
//  geometry and the last column represents the path to a file containing
//  the system geometry (e.g. a .com file) at the respective deviation.
    char ** modelist;           // array of mode-files
    char *  masses_file;        // optional file to read atomic masses
    char    masses_from_modes;  // wheter masses should be read from mode files
    char *  input_coordinates;  // main input file

// output
    char *  output_file;        // standard output file
    char *  output_fmode;       // fopen file-mode, e.g. "w" or "a"
    FILE *  fdout;              // output filestream

// threshold for floating point number comparison
    double  threshold;

// Coriolis coefficient colon separated strings
// instead of calculating the zeta values they can also be provided
// on the command line as colon separated lists.
    char * zeta_x;
    char * zeta_y;
    char * zeta_z;

} settings;

#endif
