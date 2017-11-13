
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// Offered prototypes
int InputMasses(char *inputfile, double **m);

int InputMasses(char *inputfile, double **m){

    int rows, comment_flag, control;
    unsigned int i;
    char * comment = "#%\n";
    char * line    = NULL;
    char   buffer[_MaxLineLength_] = "";
    FILE *fd;

    double aux;

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

        control = sscanf(line, "%lf", &aux);
        if(control == 1){
            (*m)  = realloc((*m),  (rows + 1) * sizeof(double));
            (*m)[rows] = aux;

            ++rows;
        }
    }
    fclose(fd); fd = NULL;

    return rows;
}
