
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// Offered prototypes
int InputCoriolisCoefficients(char *inputfile, double ****zeta, int dimension);

int InputCoriolisCoefficients(char *inputfile, double ****zeta, int dimension){

    int rows, comment_flag, control;
    unsigned int i;
    char * comment = "#%\n";
    char * line    = NULL;
    char   buffer[_MaxLineLength_] = "";
    FILE *fd;

    int  index1 = 0;
    int  index2 = 0;
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

        control = sscanf(line, "%d  %d  %lf  %lf  %lf", &index1, &index2, &aux1, &aux2, &aux3);
        if(control == 5){

            if(index1 > (dimension-1) || index2 > (dimension-1)){
                fprintf(stderr, "\n (-) Error in input file \"%s\"", inputfile);
                fprintf(stderr, "\n     Indices in file (%d|%d) exceed maximal value (%d)", index1, index2, dimension-1);
                fprintf(stderr, "\n     Aborting...\n\n");
                exit (3);
            }


            (*zeta)[index1][index2][0] = aux1;
            (*zeta)[index1][index2][1] = aux2;
            (*zeta)[index1][index2][2] = aux3;

            ++rows;

            if(rows*3 > dimension*dimension*3){
                fprintf(stderr, "\n (-) Error in input file \"%s\"", inputfile);
                fprintf(stderr, "\n     Number of entries exceeds expectations (%d)", dimension*dimension*3);
                fprintf(stderr, "\n     Please check your input (e.g. do dimensionalities match?)");
                fprintf(stderr, "\n     Aborting...\n\n");
                exit (3);
            }
        }
    }
    fclose(fd); fd = NULL;

    return rows;
}
