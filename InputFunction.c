
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>

// provided prototypes
int InputFunction(char* inputfile, double** *q, int* nq, double* *v, double** *mu, int dimension, char potential_true, char dipole_true);


int InputFunction(char* inputfile, double** *q, int* nq, double* *v, double** *mu, int dimension, char potential_true, char dipole_true){

    int expected_column_count = dimension;
    if(potential_true){ expected_column_count += 1; }
    if(dipole_true){    expected_column_count += 3; }

    FILE * fd     = NULL;
    char * token  = NULL;
    char * pos    = NULL;
    char * buffer = NULL;
    const char * comment   = "#%\n";
    const char * delimiter = " \t";

    int i, control;
    int linenumber;
    int entry_rows;

// open input file read only
    fd = fopen(inputfile, "r");
    if(fd == NULL){ perror(inputfile); exit(errno); }

// allocate memory of size _MaxLineLength_ for buffer
    buffer = malloc((_MaxLineLength_) * sizeof(char));
    if(buffer == NULL){ perror("Input function buffer"); exit(errno); }

// start parsing of file
    entry_rows = 0;
    linenumber = 0;
    while( fgets(buffer, _MaxLineLength_, fd) != NULL ){

        linenumber++;

    // check for existence of newline character. If not found the line is not
    //  fully inside of the buffer and, therefore exceeding maximum line length
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
        while(isspace(*pos) != 0 && *pos != '\0') { pos++; }
        if(strlen(pos) == 0) continue;

// buffer now contains a full (non empty) line of the input file, stripped of
//  comments and *pos points to the first, non white space character of buffer
//--------------------------------------------------------------------------------

    // point token to buffer
        token = pos;

    // if the first non-white-space char of a line is either n or N:
    //  get the number of coordinate entries
        if( *pos == 'N' || *pos == 'n' ){

        // skip reading entries when nq == NULL
            if( nq != NULL ){

            // tokenize buffer, get first entry which ain't a white space (the N)
                pos = strsep(&token, delimiter);

            // there must be <dimension> entries after the N flag:
                i = 0;
                while(i < dimension){
                    pos = strsep(&token, delimiter);

                // print an error if less than <dimension> entries are found
                    if(pos == NULL){
                        fprintf(stderr,
                            "\n (-) Error reading data from input-file \"%s\"."
                            "\n     The N line doesn't contain %d entries."
                            "\n     Aborting - please check your input..."
                            "\n\n"
                            , inputfile, dimension
                        );
                        exit(EXIT_FAILURE);
                    }

                // ignore adjacent delimiting characters
                    if(*pos == '\0'){ continue; }

                // store data in integer array
                    nq[i] = atoi(pos);
                    ++i;
                }
            }
            continue;
        }


    // Start reading data:
    // coordinates:
    //  first reallocate memory for coordinates q
        for(i = 0; i < dimension; ++i){
            (*q)[i] = realloc( (*q)[i], (entry_rows + 1) * sizeof(double) );
            if( (*q)[i] == NULL){ perror("Input function q[i]"); exit(errno); }
        }

    // read data from input file and store it in coordinates
        i = 0;
        while(i < dimension){
            pos = strsep(&token, delimiter);

        // throw an error if no data found
            if(pos == NULL){
                fprintf(stderr,
                    "\n (-) Error reading data from input-file \"%s\"."
                    "\n     Too few entries in input line number %d (only found %d of the expected %d columns)"
                    "\n     Aborting - please check your input...\n\n"
                    , inputfile, linenumber, i, expected_column_count
                );
                exit(EXIT_FAILURE);
            }

        // ignore adjacent delimiting characters
            if(*pos == '\0'){ continue; }

            (*q)[i][entry_rows] = atof(pos);
            ++i;
        }


    // if potential is to be read
        if(potential_true){
        //  reallocate memory for potential v
            (*v) = realloc( (*v), (entry_rows + 1) * sizeof(double) );
            if( (*v) == NULL){ perror("Input function v"); exit(errno); }

        // get token and convert to double
            i = 0;
            while(i < 1){
                pos = strsep(&token, delimiter);

            // throw an error if no data found
                if(pos == NULL){
                    fprintf(stderr,
                        "\n (-) Error reading data from input-file \"%s\"."
                        "\n     Too few entries in input line number %d (only found %d of the expected %d columns)"
                        "\n     Aborting - please check your input...\n\n"
                        , inputfile, linenumber, i+dimension, expected_column_count
                    );
                    exit(EXIT_FAILURE);
                }

            // ignore adjacent delimiting characters
                if(*pos == '\0'){ continue; }

                (*v)[entry_rows] = atof(pos);
                ++i;
            }
        }


    // if dipole flag equals 1 read dipole data, else increment number of entry rows and continue
        if(dipole_true){

        // dipole:
        //  reallocate memory for mu[0] to mu[2]
            for(i = 0; i < 3; ++i){
                (*mu)[i] = realloc( (*mu)[i], (entry_rows + 1) * sizeof(double) );
                if( (*mu)[i] == NULL ){ perror("Input function mu[i]"); exit(errno); }
            }

        // get tokens and convert to double
            i = 0;
            while(i < 3){
                pos = strsep(&token, delimiter);

            // throw an error if no data found
                if(pos == NULL){
                    fprintf(stderr,
                        "\n (-) Error reading data from input-file \"%s\"."
                        "\n     Too few entries in input line number %d (only found %d of the expected %d columns)"
                        "\n     Aborting - please check your input...\n\n"
                        , inputfile, linenumber, i+dimension+1, expected_column_count
                    );
                    exit(EXIT_FAILURE);
                }

            // ignore adjacent delimiting characters
                if(*pos == '\0'){ continue; }

                (*mu)[i][entry_rows] = atof(pos);
                ++i;
            }
        }

    // increment number of entry rows and move on to next line
        ++entry_rows;
    }
    free(buffer); buffer = NULL;
    fclose(fd);   fd     = NULL;

    return entry_rows;
}
