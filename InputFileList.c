
#define _GNU_SOURCE
#define _MaxLineLength_ 2048
#define _PATH_MAX_ 4096

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

// provided prototypes
int InputFileList(char* inputfile, int dimension, double** *Q, char** *coordslist);

int InputFileList(char* inputfile, int dimension, double** *Q, char** *coordslist){

    FILE * fd     = NULL;
    char * token  = NULL;
    char * pos    = NULL;
    char * buffer = NULL;
    const char * comment = "#%\n";
    const char * delim   = " \t";

    int i, control;
    int linenumber;
    int entry_rows;

// open input file read only
    fd = fopen(inputfile, "r");
    if( fd == NULL ){ perror(inputfile); exit(errno); }

// allocate memory of size _MaxLineLength_ for buffer
    buffer = malloc((_MaxLineLength_) * sizeof(char));
    if( buffer == NULL ){ perror("InputFileList buffer"); exit(errno); }

// start file parsing
    entry_rows = 0;
    linenumber = 0;
    while(fgets(buffer, _MaxLineLength_, fd) != NULL){

        linenumber++;

    // check for existence of newline character. If not found the line is not
    //  fully inside of the buffer and therefore, exceeding maximum line length
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

        token = pos;

        (*coordslist) = realloc( *(coordslist), (entry_rows + 1) * sizeof(char*) );
        if( (*coordslist) == NULL ){ perror("coordslist"); exit(errno); }

        (*coordslist)[entry_rows] = malloc( _PATH_MAX_ * sizeof(char) );
        if( (*coordslist)[entry_rows] == NULL ){ perror("coordslist[]"); exit(errno); }

        for(i = 0; i < dimension; ++i){
            (*Q)[i] = realloc( (*Q)[i], (entry_rows + 1) * sizeof(double) );
            if( (*Q)[i] == NULL ){ perror("Q[i]"); exit(errno); }
        }


    // read deviation from minimum geometry Q from input file
        i = 0;
        while(i < dimension){
            pos = strsep( &token, delim );

        // throw an error if no data found
            if(pos == NULL){
                fprintf(stderr,
                        "\n (-) Error reading data from input file \"%s\"."
                        "\n     Too few lines in input line number %d (only found %d of the expected %d columns)"
                        "\n     Aborting...\n\n"
                        , inputfile, linenumber, i, dimension+1
                    );
                exit(EXIT_FAILURE);
            }

        // ignore adjacent delimiting characters
            if(*pos == '\0'){ continue; }

            (*Q)[i++][entry_rows] = atof(pos);
        }


    // read coordinate file names from input file
        while( i < (dimension + 1) ){
            pos = strsep( &token, delim );

        // throw an error if no data found
            if(pos == NULL){
                fprintf(stderr,
                        "\n (-) Error reading data from input file \"%s\"."
                        "\n     Too few lines in input line number %d (only found %d of the expected %d columns)"
                        "\n     Aborting...\n\n"
                        , inputfile, linenumber, i, dimension+1
                    );
                exit(EXIT_FAILURE);
            }

        // ignore adjacent delimiting characters
            if(*pos == '\0'){ continue; }

            strncpy( (*coordslist)[entry_rows], pos, _PATH_MAX_ );
            ++i;
        }

        ++entry_rows;

    }
    fclose(fd); fd = NULL;

    return entry_rows;
}
