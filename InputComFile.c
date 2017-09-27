
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// Offered prototypes
int InputComFile(char *inputfile, double **x, double **y, double **z);

int InputComFile(char *inputfile, double **x, double **y, double **z){

    int i, rows, comment_flag, control;
    char * comment = "#%\n";
    char * line    = NULL;
    char   buffer[_MaxLineLength_] = "";
    FILE *fd;

    double aux1 = 0;
    double aux2 = 0;
    double aux3 = 0;

    fd = fopen(inputfile, "r");
    if(fd == NULL){
        fprintf(stderr, "\n(-) ERROR opening input-file: \"%s\"", inputfile);
        fprintf(stderr, "\n    Exiting...\n\n");
        return(-1);
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
    //printf("%s\n", line);
//-----------------------------------------------------------------------------------

        control = sscanf(line, "%*s  %lf  %lf  %lf", &aux1, &aux2, &aux3);
        if(control == 3){
            (*x)  = realloc((*x),  (rows + 1) * sizeof(double));
            (*y)  = realloc((*y),  (rows + 1) * sizeof(double));
            (*z)  = realloc((*z),  (rows + 1) * sizeof(double));

            (*x)[rows] = aux1;
            (*y)[rows] = aux2;
            (*z)[rows] = aux3;

            ++rows;
        }
    }
    fclose(fd); fd = NULL;

    return rows;
}
