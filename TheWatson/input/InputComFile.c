
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
int InputComFile(char* inputfile, double* x, double* y, double* z, int max_lines);

int InputComFile(char* inputfile, double* x, double* y, double* z, int max_lines){

// open input file read only
    FILE * fd = fopen(inputfile, "r");
    if( fd == NULL ){ perror(inputfile); exit(errno); }

// define comment characters
    const char * comment = "#%\n";

// allocate memory of size _MaxLineLength_ for buffer
    char * pos     = NULL;
    char * buffer  = malloc((_MaxLineLength_) * sizeof(char));
    if( buffer == NULL ){ perror("InputComFile buffer"); exit(errno); }


// start parsing file
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

    // throw a warning if file contains more non-comment lines than number of atoms
        if(entry_rows == max_lines){
            fprintf(stderr,
                "Warning: More than n_atoms lines exist in file \"%s\", please check your input\n"
                , inputfile
            );
        }

    // get {x,y,z}-data from line
    // only read lines with exactly four entries
        int control = sscanf(pos, "%*s  %lf  %lf  %lf", &(x[entry_rows]), &(y[entry_rows]), &(z[entry_rows]));
        if( control != 3 ){ continue; }

        ++entry_rows;
    }
    fclose(fd); fd = NULL;

    return entry_rows;
}
