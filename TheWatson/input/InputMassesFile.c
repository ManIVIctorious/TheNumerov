
#define _GNU_SOURCE
#define _MaxLineLength_ 512

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

// Required prototypes
char * PreprocessBuffer(char* inputfile, int linenumber, char* buffer, const char* comment);

// Provided prototypes
int InputMassesFile(char* inputfile, double* *m);

int InputMassesFile(char* inputfile, double* *m){

// open input file read only
    FILE * fd = fopen(inputfile, "r");
    if( fd == NULL ){ perror(inputfile); exit(errno); }

// define comment characters
    const char * comment = "#%\n";

// allocate memory of size _MaxLineLength_ for buffer
    char * pos     = NULL;
    char * buffer  = malloc((_MaxLineLength_) * sizeof(char));
    if( buffer == NULL ){ perror("InputMassesFile buffer"); exit(errno); }


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

    // get individual atom masses from line
    // the program expects exactly one entry per line,
    // containing the atomic mass of the entry_row's particle
        (*m) = realloc((*m), entry_rows+1);
        int control = sscanf(pos, "%lf", &(*m)[entry_rows]);
        if( control != 1 ){ continue; }

        ++entry_rows;
    }
    fclose(fd); fd = NULL;

    return entry_rows;
}
