
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// provided prototypes
int InputCoriolisCoefficients(char* inputfile, double** *q, double*** *zeta, double*** *mu, int dimension);

int InputCoriolisCoefficients(char* inputfile, double** *q, double*** *zeta, double*** *mu, int dimension){

    int linenumber;
    int rows, comment_flag;
    int m, n;
    unsigned int i;
    char * comment = "#%\n";
    char * line    = NULL;
    char * token = NULL;
    char   buffer[_MaxLineLength_] = "";
    FILE * fd;

    fd = fopen(inputfile, "r");
    if(fd == NULL){
        fprintf(stderr,
            "\n(-) ERROR opening input-file: \"%s\""
            "\n    Exiting..."
            "\n\n"
            , inputfile
        );
        exit(1);
    }

    rows = 0;
    linenumber = 0;
    while(fgets(buffer, sizeof(buffer), fd) != NULL){

        ++linenumber;

    // check if the first character in buffer is a comment char,
    //  if yes jump to next line
        comment_flag = 0;
        for(i = 0; i < strlen(comment); ++i){
            if(buffer[0] == comment[i]){
                comment_flag = 1;
                break;
            }
        }
        if(comment_flag == 1) continue;

    // copy "buffer" with stripped comments to new buffer "line"
        line = strtok(buffer, comment);
        if(line == NULL) continue;

    // remove leading white spaces and tabulators
    //  and skip empty lines
        while(isspace(*line)) line++;
        if(strlen(line) == 0) continue;

    // At this point the requested input line is stripped of
    //  comments and blank lines. From here on the parsing starts:
    //    printf("%s\n", line);
//-----------------------------------------------------------------------------------

    // Start reading data:
    //  Break down string array "line" to <dimension + 3> string tokens "token"

    // Get first token, convert it to double and store it in q[0]
        token = strtok(line, " \t");
        if(token != NULL){
        // reallocate memory for q array
            for(n = 0; n < dimension; ++n){
                (*q)[n] = realloc((*q)[n], (rows + 1)*sizeof(double));
                if( (*q)[n] == NULL){
                    fprintf(stderr,
                        "\n(-) ERROR in reallocation of %s"
                        "\n    Aborting..."
                        "\n\n"
                        , "coordinate"
                    );
                    exit(2);
                }
            }
            (*q)[0][rows] = atof(token);
        }

    // store the other <dimension - 1> coordinate entries
        for(m = 1; m < dimension; ++m){
            token = strtok(NULL, " \t");
            if(token == NULL){
                fprintf(stderr,
                    "\n(-) ERROR reading data from input-file \"%s\"."
                    "\n    Too few entries in input line number %d"
                    "\n    Aborting - please check your input..."
                    "\n\n"
                    , inputfile, linenumber
                );
                exit(1);
            }
            (*q)[m][rows] = atof(token);
        }

    // store potential values:
    //  increase array size:
        for(m = 0; m < 3; ++m){
            for(n = 0; n < ((dimension*(dimension - 1))/2); ++n){
                (*zeta)[m][n] = realloc((*zeta)[m][n], (rows + 1) * sizeof(double));

                if( (*zeta)[m][n] == NULL){
                    fprintf(stderr,
                        "\n(-) ERROR in reallocation of %s"
                        "\n    Aborting..."
                        "\n\n"
                        , "Coriolis coefficient"
                    );
                    exit(2);
                }
            }
        }
    // get token and convert to double
        for(m = 0; m < 3; ++m){
            for(n = 0; n < ((dimension*(dimension - 1))/2); ++n){
                token = strtok(NULL, " \t");
                if(token == NULL){
                    fprintf(stderr,
                        "\n(-) ERROR reading data from input-file \"%s\"."
                        "\n    Too few entries in input line number %d"
                        "\n    Aborting - please check your input..."
                        "\n\n", inputfile, linenumber
                    );
                    exit(1);
                }
                (*zeta)[m][n][rows] = atof(token);
            }
        }

    // store the "effective reciprocal inertia tensor" mu
    //  increase array size:
        for(m = 0; m < 3; ++m){
            for(n = 0; n < 3; ++n){
                (*mu)[m][n] = realloc((*mu)[m][n], (rows + 1) * sizeof(double));

                if( (*mu)[m][n] == NULL){
                    fprintf(stderr,
                        "\n(-) ERROR in reallocation of %s"
                        "\n    Aborting..."
                        "\n\n"
                        , "effective reciprocal inertia tensor"
                    );
                    exit(2);
                }
            }
        }
    // get token and convert to double
        for(m = 0; m < 3; ++m){
            for(n = m; n < 3; ++n){
                token = strtok(NULL, " \t");
                if(token == NULL){
                    fprintf(stderr,
                        "\n(-) ERROR reading data from input-file \"%s\"."
                        "\n    Too few entries in input line number %d"
                        "\n    Aborting - please check your input..."
                        "\n\n", inputfile, linenumber
                    );
                    exit(1);
                }
                (*mu)[m][n][rows] = atof(token);
                (*mu)[n][m][rows] = atof(token);
            }
        }

    // increment number of rows by 1
        ++rows;

    }
    fclose(fd); fd = NULL;

    return rows;
}
