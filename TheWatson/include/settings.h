#ifndef _SETTINGS_H
#define _SETTINGS_H

typedef struct settings {

    double threshold;  // threshold for number comparison

    char ** modelist;   // list of input mode files
    int dimension;      // number of included modes

// files
    // input file consists of dimension + 1 columns:
    //  the first <dimension> columns contain the deviations from
    //  the minimum geometry and the last column represents the path
    //  to a file containing the system geometry at the respective deviation.
    char * input_coordinates;
    char   input_coordinates_set;

    // standard output file
    char * output_file;
    FILE * fdout;


// data:
    int       n_atoms;
    double    tot_mass;
    double *  atom_masses;
    double ** zeta;

} settings;



#endif
