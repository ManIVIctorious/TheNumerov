
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

// Required prototypes
char * PreprocessBuffer(char* inputfile, int linenumber, char* buffer, const char* comment);

// Provided prototypes
int InputNormalMode(char* inputfile, double* *mode, char read_masses, double* *mass);

int InputNormalMode(char* inputfile, double* *mode, char read_masses, double* *mass){

// if read_masses is true expect four columns, else expect three
    int expected_column_count;
    if( read_masses ){ expected_column_count = 4; }
    else{              expected_column_count = 3; }

// open input file read only
    FILE * fd = fopen(inputfile, "r");
    if( fd == NULL ){ perror(inputfile); exit(errno); }

// define comment character
    const char * comment   = "#%\n";

// allocate memory of size _MaxLineLength_ for buffer
    char * pos    = NULL;
    char * buffer = malloc((_MaxLineLength_) * sizeof(char));
    if( buffer == NULL ){ perror("InputNormalMode buffer"); exit(errno); }


// start file parsing
    int entry_rows = 0;
    int linenumber = 0;
    while( fgets(buffer, _MaxLineLength_, fd) ){

    // pre-process buffer
    // strip buffer off comments, leading white spaces and do some error handling
        linenumber++;
        pos = PreprocessBuffer(inputfile, linenumber, buffer, comment);
        if(pos == NULL){ continue; }

// pos now points to the first non-white-space character of buffer
//-----------------------------------------------------------------------------------

    // memory allocation for mode coordinates
    // data format: mode[x0,y0,z0,x1,y1,z1,...]
        (*mode) = realloc( (*mode), 3*(entry_rows + 1) * sizeof(double) );
        if( (*mode) == NULL ){ perror("Mode"); exit(errno); }

        int control = sscanf(pos, "%lf  %lf  %lf",
                        &(*mode)[3*entry_rows + 0],
                        &(*mode)[3*entry_rows + 1],
                        &(*mode)[3*entry_rows + 2]
                      );

    // if mass should be read reallocate masses array and read its values
    // data format: mass[0,1,...]
        if( read_masses ){
            (*mass) = realloc( (*mass), (entry_rows + 1) * sizeof(double) );
            if( (*mass) == NULL ){ perror("Mass"); exit(errno); }

            control += sscanf(pos, "%*f  %*f  %*f  %lf", &(*mass)[entry_rows]);
        }

        if( control != expected_column_count ){
            fprintf(stderr,
                "\n (-) Error reading data from input-file \"%s\"."
                "\n     Too few entries in input line number %d (%d of expected %d)"
                "\n     A proper mode line contains {x,y,z}-displacements and, if no"
                "\n     masses file is given, a column representing the atomic_masses."
                "\n\n", inputfile, linenumber, control, expected_column_count
            );
            exit(EXIT_FAILURE);
        }

        ++entry_rows;
    }
    fclose(fd);   fd     = NULL;
    free(buffer); buffer = NULL;

    return entry_rows;
}
