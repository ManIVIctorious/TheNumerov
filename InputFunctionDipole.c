
#define _GNU_SOURCE
#define _MaxLineLength_ 2048

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// Offered prototypes
int InputFunctionDipole(char *inputfile, double **q1, double **q2, double **V, double **mux, double **muy, double **muz, int *nq1, int *nq2);


int InputFunctionDipole(char *inputfile, double **q1, double **q2, double **V, double **mux, double **muy, double **muz, int *nq1, int *nq2){

    int i, rows, comment_flag, control;
    char * comment = "#%\n";
    char * line    = NULL;
    char   buffer[_MaxLineLength_] = "";
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
        (*mux) = realloc((*mux), (rows + 1) * sizeof(double));
        (*muy) = realloc((*muy), (rows + 1) * sizeof(double));
        (*muz) = realloc((*muz), (rows + 1) * sizeof(double));

        control = sscanf(line, "%lf  %lf  %lf  %lf  %lf  %lf", &(*q1)[rows], &(*q2)[rows], &(*V)[rows], &(*mux)[rows], &(*muy)[rows], &(*muz)[rows]);
        if(control != 6){
          fprintf(stderr, "\n(-) ERROR reading data from input-file \"%s\".", inputfile);
          fprintf(stderr, "\n    Aborting - please check your input...\n\n");
          exit(1);
        }

        ++rows;
    }
    fclose(fd); fd = NULL;

    return rows;
}


//int main(int argc, char **argv){
//
//    double *q1, *q2, *v, *mux, *muy, *muz;
//    int n_pot;
//    int i;
//    char *inputfile = "testfile";
//    int nq1, nq2;
//    double v_min = 1.0E100;
//    
//    q1  = malloc(sizeof(double));
//    q2  = malloc(sizeof(double));
//    v   = malloc(sizeof(double));
//    mux = malloc(sizeof(double));
//    muy = malloc(sizeof(double));
//    muz = malloc(sizeof(double));
//
//    n_pot = InputFunctionDipole(inputfile, &q1, &q2, &v, &mux, &muy, &muz, &nq1, &nq2, &v_min);
//
//    printf("Number of Lines:\t%d\n", n_pot);
//    printf("%d\t%d\n", nq1, nq2);
//    printf("%lf\n", v_min);
//    for(i = 0; i < 1; ++i){
//        printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", q1[i], q2[i], v[i], mux[i], muy[i], muz[i]);
//    }
//
//    return 0;
//}
