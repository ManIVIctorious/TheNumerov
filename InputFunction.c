
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// Offered prototypes
int InputFunction(char *inputfile, double ***q, int *nq, double **V, int dimension);

int InputFunction(char *inputfile, double ***q, int *nq, double **V, int dimension){

    int rows, comment_flag;
    int n;
    unsigned int i;
    char * comment = "#%\n";
    char * line    = NULL;
    char * token = NULL;
    char   buffer[_MaxLineLength_] = "";
    FILE * fd;

    fd = fopen(inputfile, "r");
    if(fd == NULL){
        fprintf(stderr, "\n(-) ERROR opening input-file: \"%s\"", inputfile);
        fprintf(stderr, "\n    Exiting...\n\n");
        exit(1);
    }

    rows = 0;
    while(fgets(buffer, sizeof(buffer), fd) != NULL){

    // check if the first character in buffer is a comment char,
    //  if yes jump to next line
        comment_flag = 0;
        for(i=0; i<strlen(comment); ++i){
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

    // If the first non-white-space char of a line is either n or N get number of coordinate entries
        if(strcasestr(line, "N") != NULL){
        // Tokenize line, get first entry which ain't a white space
            token = strtok(line, " \t");

        // read the first <dimension> entries after the N flag
            for(n = 0; n < dimension; ++n){
                token = strtok(NULL, " \t");
            // if there are less than <dimension> entries print an error
                if(token == NULL){
                    fprintf(stderr, "\n(-) ERROR reading data from input-file \"%s\".", inputfile);
                    fprintf(stderr, "\n    The N line doesn't contain %d entries.", dimension);
                    fprintf(stderr, "\n    Aborting - please check your input...\n\n");
                    exit(1);
                }
                nq[n] = atoi(token);
            }
            continue;
        }


    // Start reading data:
    //  Break down string array "line" to <dimension + 3> string tokens "token"

    // Get first token, convert it to double and store it in q[0]
        token = strtok(line, " \t");
        if(token != NULL){
        // reallocate memory for q array
            for(n = 0; n < dimension; ++n){
                (*q)[n] = realloc((*q)[n], (rows + 1)*sizeof(double));
                if( (*q)[n] == NULL){
                    fprintf(stderr, "\n(-) ERROR in reallocation of %s", "coordinate");
                    fprintf(stderr, "\n    Aborting...\n\n");
                    exit(2);
                }
            }
            (*q)[0][rows] = atof(token);
        }

    // store the other <dimension - 1> coordinate entries
        for(n = 1; n < dimension; ++n){
            token = strtok(NULL, " \t");
            if(token == NULL){
                fprintf(stderr, "\n(-) ERROR reading data from input-file \"%s\".", inputfile);
                fprintf(stderr, "\n    Too few entries in input line number %d", rows);
                fprintf(stderr, "\n    Aborting - please check your input...\n\n");
                exit(1);
            }
            (*q)[n][rows] = atof(token);
        }

    // store potential values:
    //  increase array size:
        (*V) = realloc((*V), (rows + 1)*sizeof(double));
        if( (*V) == NULL){
            fprintf(stderr, "\n(-) ERROR in reallocation of %s", "potential");
            fprintf(stderr, "\n    Aborting...\n\n");
            exit(2);
        }
    // get token and convert to double
        token = strtok(NULL, " \t");
        if(token == NULL){
            fprintf(stderr, "\n(-) ERROR reading data from input-file \"%s\".", inputfile);
            fprintf(stderr, "\n    Too few entries in input line number %d", rows);
            fprintf(stderr, "\n    Aborting - please check your input...\n\n");
            exit(1);
        }
        (*V)[rows] = atof(token);

    // increment number of rows by 1
        ++rows;

    }
    fclose(fd); fd = NULL;

    return rows;
}
