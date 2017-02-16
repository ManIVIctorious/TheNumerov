
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int InputFunction(char *inputfile, double **q1, double **q2, double **V, int *nq1, int *nq2){

    char *comment = "#%\n";
    int i, rows, comment_flag, control;
    char *line, buffer[_MaxLineLength_];
    FILE *fd;

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
        //printf("%s\n", line);
//-----------------------------------------------------------------------------------

        if(strcasestr(line, "N") != NULL){
            control = sscanf(line, "%*s  %d  %d", &(*nq1), &(*nq2));
            if(control != 2){
                fprintf(stderr, "\n(-) ERROR reading data from input-file \"%s\".", inputfile);
                fprintf(stderr, "\n    Aborting - please check your input...\n\n");
                exit(1);
            }
            continue;
        }

        (*q1)  = realloc((*q1),  (rows + 1) * sizeof(double));
        (*q2)  = realloc((*q2),  (rows + 1) * sizeof(double));
        (*V)   = realloc((*V),   (rows + 1) * sizeof(double));

        control = sscanf(line, "%lf  %lf  %lf", &(*q1)[rows], &(*q2)[rows], &(*V)[rows]);
        if(control != 3){
          fprintf(stderr, "\n(-) ERROR reading data from input-file \"%s\".", inputfile);
          fprintf(stderr, "\n    Aborting - please check your input...\n\n");
          exit(1);
        }

        ++rows;
    }
    fclose(fd);

    return rows;
}
