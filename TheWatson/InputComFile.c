
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

// Offered prototypes
int InputComFile(char* inputfile, double* x, double* y, double* z, int max_lines);

int InputComFile(char* inputfile, double* x, double* y, double* z, int max_lines){

    FILE * fd     = NULL;
    char * token  = NULL;
    char * pos    = NULL;
    char * buffer = NULL;
    const char * comment = "#%\n";

    int i, control;
    int linenumber;
    int entry_rows;

// open input file read only
    fd = fopen(inputfile, "r");
    if( fd == NULL ){ perror(inputfile); exit(errno); }

// allocate memory of size _MaxLineLength_ for buffer
    buffer = malloc((_MaxLineLength_) * sizeof(char));
    if( buffer == NULL ){ perror("InputComFile buffer"); exit(errno); }

// start parsing file
    entry_rows = 0;
    linenumber = 0;
    while( fgets(buffer, _MaxLineLength_, fd) != NULL ){

        linenumber++;

    // check for existence of newline character. If not found the line is not
    //  fully inside of the buffer and therefore, exceeding line length.
        for(i = 0, control = 0; i < (int)strlen(buffer); ++i){
            if(buffer[i] == '\n'){
                control = 1;
            }
        }
        if(control == 0){
            fprintf(stderr,
                "\n (-) Error in input file \"%s\", line \"%d\" is too long."
                "\n     Aborting..."
                "\n\n", inputfile, linenumber
            );
            exit(EXIT_FAILURE);
        }

    // strip buffer from comments
        token = buffer;
        pos = strsep(&token, comment);
        if(pos == NULL) continue;

    // remove leading white spaces and empty lines
        while( isspace(*pos) && *pos != '\0' ){ pos++; }
        if(strlen(pos) == 0) continue;

// buffer now contains a full (non empty) line of the input file, stripped of
//  comments and *pos points to the first, non white space character of buffer
//-----------------------------------------------------------------------------------

    if(entry_rows == max_lines){
        fprintf(stderr,
                "\n (-) Error: Reached maximum number of lines in file \"%s\"\n"
                "\n     but more lines do exist, please check your input."
                "\n     Aborting...\n\n"
                , inputfile
            );
    }

        control = sscanf(buffer, "%*s  %lf  %lf  %lf", &(x[entry_rows]), &(y[entry_rows]), &(z[entry_rows]));
        if( control != 3){ continue; }

        ++entry_rows;
    }
    fclose(fd); fd = NULL;

    return entry_rows;
}
