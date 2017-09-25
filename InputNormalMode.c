
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// Offered prototypes
int InputNormalMode(char *inputfile, int start, double **modedisplacement, double **mass);

int InputNormalMode(char *inputfile, int start, double **modedisplacement, double **mass){

    int i, rows, comment_flag, control;
    char * comment = "#%\n";
    char * line    = NULL;
    char   buffer[_MaxLineLength_] = "";
    FILE *fd;

    fd = fopen(inputfile, "r");
    if(fd == NULL){
        fprintf(stderr, "\n(-) ERROR opening input-file: \"%s\"", inputfile);
        fprintf(stderr, "\n    Exiting...\n\n");
        return(-1);
    }

    rows = start;
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

        (*modedisplacement) = realloc((*modedisplacement), 3*(rows + 1) * sizeof(double));
        (*mass)             = realloc((*mass),               (rows + 1) * sizeof(double));

        control = sscanf(line, "%lf  %lf  %lf  %lf",
                                &(*modedisplacement)[3*rows],
                                &(*modedisplacement)[3*rows + 1],
                                &(*modedisplacement)[3*rows + 2],
                                &(*mass)[rows]
                        );
        if(control != 4){
          fprintf(stderr, "\n(-) ERROR reading data from input-file \"%s\".", inputfile);
          fprintf(stderr, "\n    Aborting - please check your input...\n\n");
          return(-1);
        }

        ++rows;
    }
    fclose(fd); fd = NULL;

    return rows;
}
