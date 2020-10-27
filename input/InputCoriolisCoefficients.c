
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


    // zeta values:
    //  lines containing zeta values start with Zeta_{x,y,z} followed by
    //  (dim*(dim-1))/2 columns containing the respective Coriolis coefficients
        for(i = 0; i < 3; ++i){
            control = 0;

        //  check if line starts with Zeta_{x,y,z}:
            const char * zeta_line_inits[3] = {"Zeta_x:", "Zeta_y:", "Zeta_z:" };
            if( strncasecmp(stringp, zeta_line_inits[i], 7) == 0 ){
                control = 1;

            // set stringp to second field
                strsep(&stringp, delimiter);

                m = 0;
                while( m < (dimension * (dimension - 1))/2 ){
                // set pos to next field
                    pos = strsep(&stringp, delimiter);

                // throw an error if no data found
                    if(pos == NULL){
                        fprintf(stderr,
                            "\n (-) Error reading data from input-file \"%s\"."
                            "\n     Too few entries in input line number %d (only found %d of expected %d columns)"
                            "\n     Aborting, please check your input...\n\n"
                            , inputfile, linenumber, m, (dimension * (dimension - 1))/2
                        );
                        exit(EXIT_FAILURE);
                    }

                // ignore adjacent delimiters
                    if(*pos == '\0'){ continue; }

                    zeta[i][m++] = atof(pos);
                }
                break;
            }
        }
        if(control == 1){ continue; }


    // Effective Reciprocal Moment of Inertia Tensor lines:
    //  these lines start with <dimension> coordinates on the potential energy surface
    //  followed by the upper triangle (with main diagonal) of the symmetric tensor
    //  <dim1>...<dimN> <mu_xx> <mu_xy> <mu_xz> <mu_yy> <mu_yz> <mu_zz>

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
                    , inputfile, linenumber, i, (dimension + 6)
                );
                exit(EXIT_FAILURE);
            }

        // ignore adjacent delimiting characters
            if(*pos == '\0'){ continue; }

            (*q)[i++][entry_rows] = atof(pos);
        }
        n = i;


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
                        , inputfile, linenumber, n, (dimension + 6)
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
