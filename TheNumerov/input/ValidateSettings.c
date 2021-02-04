
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "ConvertString.h"
#include "settings.h"

void ValidateSettings(settings *set){

/*Eigensolver:
 * create bit-mask of available eigensolvers
 * 0:   none    0 0 0 0 0 0 0 0
 * 1:   MKL     0 0 0 0 0 0 0 1
 * 2:   ARMA    0 0 0 0 0 0 1 0
 * 4:   GSL     0 0 0 0 0 1 0 0 (not implemented)
 */
    unsigned int available_eigensolvers = 0;
#ifdef HAVE_MKL_INSTALLED
    available_eigensolvers += 1;
#endif
#ifdef HAVE_ARMA_INSTALLED
    available_eigensolvers += 2;
#endif

  // check if set->Eigensolver is a power of two
  //  if x is a power of two, (x-1) must be its ones-complement
    if( !(set->Eigensolver) || (set->Eigensolver & (set->Eigensolver-1)) ){
        fprintf(stderr,
            "\n (-) Error: Eigensolver (number %u) is invalid."
            "\n     Its number must be a power of two, greater than zero."
            "\n     Aborting...\n\n"
            , set->Eigensolver
        );
        exit(EXIT_FAILURE);
    }

  // check if the set eigensolver is actually available (compare with bit-mask)
    if( !(set->Eigensolver & available_eigensolvers) ){

        char * solvername = "unknown";
        if     ( set->Eigensolver == 1 ){ solvername = "Intel_MKL_FEAST";  }
        else if( set->Eigensolver == 2 ){ solvername = "Armadillo_ARPACK"; }

        fprintf(stderr,
            "\n (-) Error: The set eigensolver \"%s\" (number %d) is not available."
            "\n     Please make sure to compile the program with the appropriate"
            "\n     makro definitions, e.g."
            "\n         -D HAVE_MKL_INSTALLED"
            "\n         -D HAVE_ARMA_INSTALLED"
            "\n     Aborting...\n\n"
            , solvername, set->Eigensolver
        );
        exit(EXIT_FAILURE);
    }


// Rotation and periodicity are mutually exclusive
    if( set->periodic && set->coriolis_file ){
        fprintf(stderr,
            "\n (-) Error: Rotation and periodicity are mutually exclusive."
            "\n     Either disable the periodicity setting or unset the"
            "\n     Coriolis file."
            "\n     Aborting...\n\n"
        );
        exit(EXIT_FAILURE);
    }

// Bit-mask of general solvers
    unsigned int general_solvers = 0;
    general_solvers += 2;   // Armadillo ARPACK gen solver

// for dimensions higher than two the matrix symmetry is broken via the terms of
// the Watson Hamiltonian resulting in the requirenment of a general solver.
    if( (set->dimension > 2) && set->coriolis_file && !(general_solvers & set->Eigensolver) ){
        fprintf(stderr,
            "\n (-) Error: The application of the molecular Hamiltonian including"
            "\n     rotational terms in dimensions greater than two requires the"
            "\n     utilisation of a general matrix solving algorithm."
            "\n     Please choose an appropriate eigensolver."
            "\n     Aborting...\n\n"
        );
        exit(EXIT_FAILURE);
    }



//Reduced masses:
// The reduced masses are given as a colon separated string array with the number of
//  entries being the dimensionality. To ensure the dimensionality is set beforehand
//  the string is saved temporarily and the conversion to a double array is
//  performed after the initial input

// allocate memory array and initialise all entries to 1.0 g/mol
    set->masses = malloc(set->dimension * sizeof(double));
    if(set->masses == NULL){ perror("Reduced masses"); exit(errno); }

    for(int i = 0; i < set->dimension; ++i){ set->masses[i] = 1.0; }

// if masses string is set and not set to "default" get reduced masses
    if( set->masses_string && (strcasecmp(set->masses_string, "default") != 0) ){

        char * stringp = set->masses_string;
        char * token   = NULL;

        for(int i = 0; i < set->dimension; ++i){
            token = strsep(&stringp, ":");
            if( !token ){
                fprintf(stderr,
                    "\n (-) Error: Invalid masses string detected."
                    "\n     Only %d of %d expected entries found."
                    "\n     Aborting...\n\n"
                    , i+1, set->dimension
                );
                exit(EXIT_FAILURE);
            }
            set->masses[i] = convertstring_to_double(token, "Reduced mass", NULL);
        }
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
