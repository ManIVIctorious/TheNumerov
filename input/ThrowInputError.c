
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

// Provided prototypes
void ThrowInputError(char* inputfile, int linenumber, char* format,...);

void ThrowInputError(char* inputfile, int linenumber, char* format,...){

// use a variable argument list to print errors:
//  has to be initialised by va_start() and closed by va_end() macros
//  afterwards, each call to ap gives the subsequent argument
    va_list ap;

    fprintf(stderr,
        "\n (-) Error in input file \"%s\", line %d."
        , inputfile, linenumber
    );

    if( (format != NULL) && (strlen(format) > 0) ){
        fprintf(stderr, "\n     ");
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
    }

    fprintf(stderr, "\n     Aborting...\n\n");

    exit(EXIT_FAILURE);
}
