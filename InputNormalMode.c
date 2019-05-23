
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

// Offered prototypes
int InputNormalMode(char* inputfile, double* *mode, double* *mass);

int InputNormalMode(char* inputfile, double* *mode, double* *mass){

    FILE * fd     = NULL;
    char * token  = NULL;
    char * pos    = NULL;
    char * buffer = NULL;
    const char * comment   = "#%\n";

    int i, control;
    int linenumber;
    int entry_rows;

// open input file read only
    fd = fopen(inputfile, "r");
    if( fd == NULL ){ perror(inputfile); exit(errno); }

// allocate memory of size _MaxLineLength_ for buffer
    buffer = malloc((_MaxLineLength_) * sizeof(char));
    if( buffer == NULL ){ perror("InputNormalMode buffer"); exit(errno); }

// start file parsing
    entry_rows = 0;
    linenumber = 0;
    while( fgets(buffer, _MaxLineLength_, fd) != NULL ){

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

    // memory allocation for input lines
        (*mode) = realloc( (*mode), 3*(entry_rows + 1) * sizeof(double) );
        if( (*mode) == NULL ){ perror("Mode"); exit(errno); }
        (*mass) = realloc( (*mass),   (entry_rows + 1) * sizeof(double) );
        if( (*mass) == NULL ){ perror("Mass"); exit(errno); }

        control = sscanf(buffer, "%lf  %lf  %lf  %lf",
                            &(*mode)[3*entry_rows    ],
                            &(*mode)[3*entry_rows + 1],
                            &(*mode)[3*entry_rows + 2],
                            &(*mass)[  entry_rows    ]
                        );

        if(control != 4){
            fprintf(stderr,
                    "\n (-) Error reading data from input-file \"%s\"."
                    "\n     Too few entries in input line number %d (only found %d of the expected %d columns)"
                    , inputfile, linenumber, control, 4
                );
            exit(EXIT_FAILURE);
        }

        ++entry_rows;
    }
    fclose(fd); fd = NULL;

    return entry_rows;
}
