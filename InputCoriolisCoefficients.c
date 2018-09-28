
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>

// provided prototypes
int InputCoriolisCoefficients(char* inputfile, double** *q, double** zeta, double*** *mu, int dimension);

int InputCoriolisCoefficients(char* inputfile, double** *q, double** zeta, double*** *mu, int dimension){

    FILE * fd      = NULL;
    char * buffer  = NULL;
    char * stringp = NULL;
    char * pos     = NULL;
    const char * comment   = "#%\n";
    const char * delimiter = " \t";

    int i, m, n, control;
    int linenumber;
    int entry_rows;

// open input file
    fd = fopen(inputfile, "r");
    if(fd == NULL){ perror(inputfile); exit(errno); }

// allocate memory for buffer
    buffer = malloc((_MaxLineLength_) * sizeof(char));
    if(buffer == NULL){ perror("Input Coriolis file buffer"); exit(errno); }

// start file parsing
    entry_rows = 0;
    linenumber = 0;
    while( fgets(buffer, _MaxLineLength_, fd) != NULL ){

        linenumber++;

    // check for existence of newline character, if not found the line
    //  exceeds maximum line length
        for(i = 0, control = 0; i < (int)strlen(buffer); ++i){
            if(buffer[i] == '\n'){
                control = 1;
            }
        }
        if(control == 0){
            fprintf(stderr,
                "\n (-) Error in Coriolis input file \"%s\", line \"%d\" is too long."
                "\n     Aborting..."
                "\n\n", inputfile, linenumber
            );
            exit(EXIT_FAILURE);
        }

    // strip buffer from comments
        stringp = buffer;
        pos = strsep(&stringp, comment);
        if(pos == NULL){ continue; }

    // remove leading white spaces and skip empty lines
        while( isspace(*pos) && *pos != '\0' ){ pos++; }
        if(strlen(pos) == 0) continue;

// buffer now contains a (non empty) line of the input file, stripped from comments
//  and *pos points to the first non white space character of buffer
//-----------------------------------------------------------------------------------

    // point stringp to first non white space character in buffer
        stringp = pos;

    // Coordinates:
    //  reallocate memory
        for(i = 0; i < dimension; ++i){
            (*q)[i] = realloc((*q)[i], (entry_rows + 1)*sizeof(double));
            if( (*q)[i] == NULL){ perror("Coriolis input function q[i]"); exit(errno); }
        }

    //  read data from input file and store it in coordinates
        i = 0;
        while(i < dimension){
            pos = strsep(&stringp, delimiter);

        // throw an error if no data found
            if(pos == NULL){
                fprintf(stderr,
                    "\n (-) Error reading data from input-file \"%s\"."
                    "\n     Too few entries in input line number %d (only found %d of expected %d columns)"
                    "\n     Aborting, please check your input...\n\n"
                    , inputfile, linenumber, i, (dimension + 3 * ((dimension * (dimension - 1))/2) + 6)
                );
                exit(EXIT_FAILURE);
            }

        // ignore adjacent delimiting characters
            if(*pos == '\0'){ continue; }

            (*q)[i][entry_rows] = atof(pos);
            ++i;
        }
        n = i;


    // Zeta values:
    //  read data from input file and store it in zeta
        for(m = 0; m < 3; ++m){

            i = 0;
            while(i < (dimension*(dimension - 1)) / 2){

                pos = strsep(&stringp, delimiter);

            // throw an error if no data found
                if(pos == NULL){
                    fprintf(stderr,
                        "\n (-) Error reading data from input-file \"%s\"."
                        "\n     Too few entries in input line number %d (only found %d of expected %d columns)"
                        "\n     Aborting, please check your input...\n\n"
                        , inputfile, linenumber, n, (dimension + 3 * ((dimension * (dimension - 1))/2) + 6)
                    );
                    exit(EXIT_FAILURE);
                }

            // ignore adjacent delimiting characters
                if(*pos == '\0'){ continue; }

                if(zeta[m][i] == 0.0){
                    zeta[m][i] = atof(pos);

                }else if( zeta[m][i] != atof(pos) ){
                    fprintf(stderr,
                        "\n (-) Error reading data from input-file \"%s\"."
                        "\n     Since the zeta parameters are only dependent on the mode files"
                        "\n     they have to be constant in regard of the deviation from minimum"
                        "\n     Line\told zeta\tnew zeta:"
                        "\n      %d \t% .12le \t% .12le"
                        "\n     Aborting - please check your input..."
                        "\n\n", inputfile, linenumber, zeta[m][i], atof(pos)
                    );
                    exit(EXIT_FAILURE);
                }
                ++i;
                ++n;
            }
        }


    // "Effective reciprocal inertia tensor" mu:
    //  reallocate memory
        for(i = 0; i < 3; ++i){
            for(m = 0; m < 3; ++m){
                (*mu)[i][m] = realloc((*mu)[i][m], (entry_rows + 1) * sizeof(double));
                if( (*mu)[i][m] == NULL){ perror("Coriolis input function mu[m][n]"); exit(errno); }
            }
        }

    //  read data from input file and store it in mu
        for(m = 0; m < 3; ++m){

            i = 0;
            while( (i+m) < 3 ){

                pos = strsep(&stringp, delimiter);

            // throw an error if no data found
                if(pos == NULL){
                    fprintf(stderr,
                        "\n (-) Error reading data from input-file \"%s\"."
                        "\n     Too few entries in input line number %d (only found %d of expected %d columns)"
                        "\n     Aborting, please check your input...\n\n"
                        , inputfile, linenumber, n, (dimension + 3 * ((dimension * (dimension - 1))/2) + 6)
                    );
                    exit(EXIT_FAILURE);
                }

            // ignore adjacent delimiting characters
                if(*pos == '\0'){ continue; }

                (*mu)[m][i+m][entry_rows] = atof(pos);
                (*mu)[i+m][m][entry_rows] = (*mu)[m][i+m][entry_rows];
                ++i;
                ++n;
            }
        }

    // increment number of rows by 1
        ++entry_rows;
    }
    fclose(fd); fd = NULL;

    return entry_rows;
}
