
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>

#include "ConvertString.h"

// Required prototypes
char * PreprocessBuffer(char* inputfile, int linenumber, char* buffer, const char* comment);
void ThrowInputError(char* inputfile, int linenumber, char* format,...);

// Provided prototypes
int InputDataFile(char* inputfile, double** *q, int* nq, double* *v, double** *mu, int dimension, char read_potential, char read_dipole);


int InputDataFile(char* inputfile, double** *q, int* nq, double* *v, double** *mu, int dimension, char read_potential, char read_dipole){

// open input file read only
    FILE * fd = fopen(inputfile, "r");
    if(fd == NULL){ perror("Data inputfile"); exit(errno); }

// define delimiters
    char * comment = "#%\n";  // strip buffer off anything after these chars
    char * delimit = " \t";   // column delimiters

// allocate memory of size _MaxLineLength_ for buffer
    char * pos     = NULL;
    char * stringp = NULL;
    char * buffer  = malloc((_MaxLineLength_) * sizeof(char));
    if(buffer == NULL){ perror("Input function buffer"); exit(errno); }

// calculated expected column count
    int expected_column_count = dimension;
    if(read_potential){ expected_column_count += 1; }
    if(read_dipole){    expected_column_count += 3; }


// start parsing of file
    int entryrows  = 0;
    int linenumber = 0;
    while( fgets(buffer, _MaxLineLength_, fd) != NULL ){

    // pre-process buffer
    // strip buffer off comments, leading white spaces and do some error handling
        linenumber++;
        pos = PreprocessBuffer(inputfile, linenumber, buffer, comment);
        if(pos == NULL){ continue; }

        stringp = pos;

    // Parse the N line (but only if nq != NULL)
    //--------------------------------------------------------------------------
        if( nq && strncasecmp(pos, "n", 1) == 0 ){

        // tokenise buffer, get first entry which ain't a white space (the N)
            do{
                pos = strsep(&stringp, delimit);
            }while( (pos != NULL) && (*pos == '\0') );

        // case insensitive check of the first item
            if( strcasecmp(pos, "N") != 0 ){
                ThrowInputError(inputfile, linenumber,
                    "Invalid string detected (\"%s\"), should be \"N\".", pos
                );
            }

        // there must be <dimension> entries after the N flag
            for(int i = 0; i < dimension; ){

                do{
                    pos = strsep(&stringp, delimit);
                }while ( (pos != NULL) && (*pos == '\0') );

            // print an error if less than <dimension> entries are found
                if( pos == NULL ){
                    ThrowInputError(inputfile, linenumber,
                        "The N line doesn't contain %d entries."
                        , dimension
                    );
                }

            // store data in integer array
                nq[i++] = (int)convertstring_to_long(pos, "nq", NULL);

            }

            continue;
        }


    // Start reading data:
    //--------------------------------------------------------------------------
    // Coordinates q
    //----------------
    //  reallocate memory for coordinates
        for(int i = 0; i < dimension; ++i){
            (*q)[i] = realloc( (*q)[i], (entryrows + 1) * sizeof(double) );
            if( (*q)[i] == NULL){ perror("Input data function q[i]"); exit(errno); }
        }

    // read data from input file and store it in coordinates
        for(int i = 0; i < dimension; ){

        // get next token
            do{
                pos = strsep(&stringp, delimit);
            }while( (pos != NULL) && (*pos == '\0') );

        // throw an error if no data found
            if( pos == NULL ){
                ThrowInputError(inputfile, linenumber,
                    "\n     Too few entries in input line number %d"
                    "(only found %d of the expected %d columns)"
                    , i, expected_column_count
                );
            }

        // save value to coordinate array
            (*q)[i++][entryrows] = convertstring_to_double(pos, "Coordinates q", NULL);
        }


    // Potential v
    //--------------
        if( read_potential ){
        //  reallocate memory for potential v
            (*v) = realloc( (*v), (entryrows + 1) * sizeof(double) );
            if( (*v) == NULL){ perror("Input data function v"); exit(errno); }

        // get next token
            do{
                pos = strsep(&stringp, delimit);
            }while( (pos != NULL) && (*pos == '\0') );

        // throw an error if no data found
            if( pos == NULL ){
                ThrowInputError(inputfile, linenumber,
                    "\n     Too few entries in input line number %d"
                    "(only found %d of the expected %d columns)"
                    , dimension, expected_column_count
                );
            }

        // save value to potential array
            (*v)[entryrows] = convertstring_to_double(pos, "Potential v", NULL);
        }


    // Dipole moment mu
    //-------------------
        if( read_dipole ){
        //  reallocate memory for dipole moment mu[3] array
            for(int i = 0; i < 3; ++i){
                (*mu)[i] = realloc( (*mu)[i], (entryrows + 1) * sizeof(double) );
                if( (*mu)[i] == NULL ){ perror("Input function mu[i]"); exit(errno); }
            }

        // read data from input file and store it in dipole moment array
            for(int i = 0; i < 3; ){

            // get next token
                do{
                    pos = strsep(&stringp, delimit);
                }while( (pos != NULL) && (*pos == '\0') );

            // throw an error if no data found
                if( pos == NULL ){
                    ThrowInputError(inputfile, linenumber,
                        "\n     Too few entries in input line number %d"
                        "(only found %d of the expected %d columns)"
                        , i, expected_column_count
                    );
                }

            // save value to coordinate array
                (*mu)[i++][entryrows] = convertstring_to_double(pos, "Dipole moment mu", NULL);
            }
        }

    // increment number of entry rows and move on to next line
        ++entryrows;
    }
    free(buffer); buffer = NULL;
    fclose(fd);   fd     = NULL;

    return entryrows;
}
