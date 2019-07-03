
#define _GNU_SOURCE
#define _MaxLineLength_ 2048
#define _PATH_MAX_ 4096

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "settings.h"

// dependencies
void EffectiveReciprocalMomentofInertia(settings prefs, double* q, char* coordsfile);

// provided prototypes
int ProcessFileList(settings prefs);

int ProcessFileList(settings prefs){

    FILE * fd     = NULL;
    char * token  = NULL;
    char * pos    = NULL;
    char * buffer = NULL;
    const char * comment = "#%\n";
    const char * delim   = " \t";

    int i, control;
    int linenumber;
    int entry_rows;

// data in FileList-file (to be further processed)
//--------------------------------------------------
    char coordsfile[_PATH_MAX_];
    double * q = malloc(prefs.dimension * sizeof(double));
    if(q == NULL){ perror("q in ProcessFileList"); exit(errno); }

// Begin with file input
//--------------------------------------------------
// open input file read only
    fd = fopen(prefs.input_coordinates, "r");
    if( fd == NULL ){ perror(prefs.input_coordinates); exit(errno); }

// allocate memory of size _MaxLineLength_ for buffer
    buffer = malloc((_MaxLineLength_) * sizeof(char));
    if( buffer == NULL ){ perror("ProcessFileList buffer"); exit(errno); }

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
                "\n\n", prefs.input_coordinates, linenumber
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


    // read deviation from minimum geometry Q from input file
        i = 0;
        while(i < prefs.dimension){
            pos = strsep( &token, delim );

        // throw an error if no data found
            if(pos == NULL){
                fprintf(stderr,
                        "\n (-) Error reading data from input file \"%s\"."
                        "\n     Too few lines in input line number %d (only found %d of the expected %d columns)"
                        "\n     Aborting...\n\n"
                        , prefs.input_coordinates, linenumber, i, prefs.dimension+1
                    );
                exit(EXIT_FAILURE);
            }

        // ignore adjacent delimiting characters
            if(*pos == '\0'){ continue; }

            q[i++] = atof(pos);
        }


    // read coordinate file names from input file
        while( i < (prefs.dimension + 1) ){
            pos = strsep( &token, delim );

        // throw an error if no data found
            if(pos == NULL){
                fprintf(stderr,
                        "\n (-) Error reading data from input file \"%s\"."
                        "\n     Too few lines in input line number %d (only found %d of the expected %d columns)"
                        "\n     Aborting...\n\n"
                        , prefs.input_coordinates, linenumber, i, prefs.dimension+1
                    );
                exit(EXIT_FAILURE);
            }

        // ignore adjacent delimiting characters
            if(*pos == '\0'){ continue; }

            strncpy( coordsfile, pos, _PATH_MAX_ );
            ++i;
        }

    // After each read data line start directly with processing
    //  Determine Effective Reciprocal Moment of Inertia Tensor mu
    //------------------------------------------------------------
    #ifdef debug_coords
    //{{{
    // output coordinates for control
        fprintf(stderr, "\n%d\t%s\n", entry_rows+1, coordsfile);
        fprintf(stderr,
                "\nNumber of Atoms:\t%d"
                "\nInput coordinates:\n"
                "\t         x         "
                "\t         y         "
                "\t         z         "
                "\t    atomic mass    "
                "\n"
                , prefs.n_atoms
            );
    //}}}
    #endif

        EffectiveReciprocalMomentofInertia(prefs, q, coordsfile);

        ++entry_rows;
    }
    fclose(fd); fd = NULL;

    return entry_rows;
}
