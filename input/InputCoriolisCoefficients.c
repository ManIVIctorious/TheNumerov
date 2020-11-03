
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
int InputCoriolisCoefficients(char* inputfile, double** *q, double** zeta, double*** *mu, int dimension);

int InputCoriolisCoefficients(char* inputfile, double** *q, double** zeta, double*** *mu, int dimension){

// open input file read only
    FILE * fd = fopen(inputfile, "r");
    if(fd == NULL){ perror("Coriolis input file"); exit(errno); }

// define delimiters
    char * comment = "#%\n";  // strip buffer off anything after these chars
    char * delimit = " \t";   // column delimiters

// allocate memory of size _MaxLineLength_ for buffer
    char * pos     = NULL;
    char * stringp = NULL;
    char * buffer  = malloc((_MaxLineLength_) * sizeof(char));
    if(buffer == NULL){ perror("Input Coriolis file buffer"); exit(errno); }

// expected column counts
    int zeta_column_count = ( dimension * (dimension-1) ) / 2;

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


    // Parse Zeta_{x,y,z} values:
    //--------------------------------------------------------------------------
    //  lines containing zeta values start with Zeta_{x,y,z}: followed by
    //  (dim*(dim-1))/2 columns containing the respective Coriolis coefficients

    // check if line starts with "Zeta_"
        if( strncasecmp(stringp, "Zeta_", 5) == 0 ){

        // get first token of buffer
            do{
                pos = strsep(&stringp, delimit);
            }while( (pos != NULL) && (*pos == '\0') );

        // define list of allowed first tokens
            const char * zeta_line_inits[3] = { "Zeta_x:" , "Zeta_y:" , "Zeta_z:" };

            for(int i = 0; i < 3; ++i){
                int recognised = 0;

            // compare first token to entries of token list
                if( strcasecmp(pos, zeta_line_inits[i]) == 0 ){
                    recognised = 1;

                // get zeta values
                    for(int j = 0; j < zeta_column_count; ){

                        do{
                            pos = strsep(&stringp, delimit);
                        }while( (pos != NULL) && (*pos == '\0') );

                    // print an error if less than (dim*(dim-1)/2) entries are found
                        if( pos == NULL ){
                            ThrowInputError(inputfile, linenumber,
                                "Too few entries in Zeta_%c: line"
                                "(only found %d of expected %d columns)"
                                , "xyz"[i], j+1, zeta_column_count
                            );
                        }

                    // store data in double array
                        zeta[i][j] = convertstring_to_double(pos, "zeta", NULL);
                    }
                }

                if( !recognised ){
                    ThrowInputError(inputfile, linenumber, "Unrecognised keyword %s" , pos);
                }
            }
        }


    // Effective Reciprocal Moment of Inertia Tensor lines:
    //--------------------------------------------------------------------------
    //  these lines start with <dim> coordinates on the potential energy surface
    //  followed by the upper triangle (with main diagonal) of the symmetric tensor
    //  <dim1>...<dimN> <mu_xx> <mu_xy> <mu_xz> <mu_yy> <mu_yz> <mu_zz>

    // Coordinates q
    //----------------
    //  reallocate memory for coordinates
        for(int i = 0; i < dimension; ++i){
            (*q)[i] = realloc( (*q)[i], (entryrows + 1) * sizeof(double) );
            if( (*q)[i] == NULL ){ perror("Coriolis input function q[i]"); exit(errno); }
        }

    //  read data from input file and store it in coordinates
        for(int i = 0; i < dimension; ){

        // get next token
            do{
                pos = strsep(&stringp, delimit);
            }while( (pos != NULL) && (*pos == '\0') );

        // throw an error if no data found
            if( pos == NULL ){
                ThrowInputError(inputfile, linenumber,
                    "Too few entries in input line"
                    "(only found %d of expected %d columns)"
                    , i, (dimension+6)
                );
            }

        // save value to coordinate array
            (*q)[i++][entryrows] = convertstring_to_double(pos, "Coordinates q", NULL);
        }


    // "Effective reciprocal inertia tensor" mu:
    //--------------------------------------------
    //  reallocate memory, only the upper triangle of the 3x3 tensor is actually allocated memory,
    //  the lower tringle just points to the values of the upper one
        for(int i = 0; i < 3; ++i){
            for(int j = i; j < 3; ++j){
                (*mu)[i][j] = realloc((*mu)[i][j], (entryrows + 1) * sizeof(double));
                if( (*mu)[i][j] == NULL){ perror("Coriolis input function mu[m][n]"); exit(errno); }
            }
        }

    //  read data from input file and store it in mu
        int count = 0;
        for(int i = 0; i < 3; ++i){
            for(int j = i; j < 3; ++j){
            // get next token
                do{
                    pos = strsep(&stringp, delimit);
                }while( (pos != NULL) && (*pos == '\0') );

            // throw an error if no data found
                if( pos == NULL ){
                    ThrowInputError(inputfile, linenumber,
                        "Too few entries in input line"
                        "(only found %d of expected %d columns)"
                        , dimension+count, (dimension+6)
                    );
                }

            // save value to reciprocal moment of inertia tensor
                (*mu)[i][j][entryrows] = convertstring_to_double(pos, "Coordinates q", NULL);

            // increment count
                count++;
            }
        }


    // increment number of rows by 1
        ++entryrows;
    }
    fclose(fd); fd = NULL;

// point lower triangle of reciprocal moment of inertia 3x3 tensor to upper triangle
    (*mu)[1][0] = (*mu)[0][1];
    (*mu)[2][0] = (*mu)[0][2];
    (*mu)[2][1] = (*mu)[1][2];

    return entryrows;
}
