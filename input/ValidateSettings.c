
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "ConvertString.h"
#include "settings.h"

void ValidateSettings(settings *set){

//Reduced masses:
// The reduced masses are given as a colon separated string array with the number of
//  entries being the dimensionality. To ensure the dimensionality is set beforehand
//  the string is saved temporarily and the conversion to a double array is
//  performed after the initial input

// allocate memory array and initialize all entries to 1.0 g/mol
    set->masses = malloc(set->dimension * sizeof(double));
    if(set->masses == NULL){ perror("Reduced masses"); exit(errno); }

    for(int i = 0; i < set->dimension; ++i){ set->masses[i] = 1.0; }

// get reduced masses from masses_string
    if( set->masses_string_set ){

        char * stringp = set->masses_string;
        char * token   = NULL;

        for(int i = 0; i < set->dimension; ++i){
            token = strsep(&stringp, ":");
            set->masses[i] = convertstring_to_double(token, "Reduced mass", NULL);
        }
    }


// Eigensolver:
//  Valid settings for the used Eigensolver are
//      0         unknown
//      1   aka   Intel_MKL_FEAST
//      2   aka   Armadillo_ARPACK
    if( set->Eigensolver > 2 ){
        set->Eigensolver = 0;
    }


// External dipole file:
//  For the calculation of intensities / oscillator strengths the system's dipole moment in
//  {x,y,z}-direction has to be provided.  Either these dipole moments are directly provided
//  in the main data input file or via an additional dipole moments file.
//  Hence, disable reading dipole from primary input file when external dipole is used
    if( set->ext_dip_file ){
        set->dipole = 0;
    }

}
