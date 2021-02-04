
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>

#include "controlfile.h"

// Provided Prototypes
void free_keywordlistvalues(struct keywords* kwlist);
static void rmSurroundingWhitespaces(char* string);

// Dependencies
char * PreprocessBuffer(char* inputfile, int linenumber, char* buffer, const char* comment);

void free_keywordlistvalues(struct keywords* kwlist){

    for(int i = 0; kwlist[i].keyword; ++i){

        if(kwlist[i].value != NULL){
            for(int j = 0; j < kwlist[i].set; ++j){
                free(kwlist[i].value[j]);
                kwlist[i].value[j] = NULL;
            }
            free(kwlist[i].value);
            kwlist[i].value = NULL;
        }
    }
}


int TokeniseControlFile(char* filename, struct keywords* kwlist, int ExitOnError){

    const char * comment = "%#\n";

// open input file and create buffer
    FILE * fd = fopen(filename, "r");
    if( !fd ){ perror(filename); exit(errno); }

    char * buffer = malloc((_MaxLineLength_) * sizeof(char));
    if( !buffer ){ perror("Control file parser buffer"); exit(errno); }

// start reading input file
    int linenumber     = 0;
    int n_keywords_set = 0;
    while( fgets(buffer, _MaxLineLength_, fd) ){

    // pre-process buffer
    // strip buffer of comments, leading white spaces and do some error handling
        linenumber++;
        char * pos = PreprocessBuffer(filename, linenumber, buffer, comment);
        if( pos == NULL ){ continue; }

// buffer now contains a full (non empty) line of the control file, stripped of
//  comments and *pos points to the first, non white space character of buffer
//--------------------------------------------------------------------------------

    // split at every occurrence of semi-colon
        char * stringp = pos;
        while( (pos = strsep(&stringp, ";")) ){

        // ignore leading whitespaces and check string length
            while( isspace(*pos) && (*pos != '\0') ){ pos++; }
            if( strlen(pos) == 0 ){ continue; }

        // split keyword = value pair at equality sign
            char * value = pos;
            char * kword = strsep(&value, "=");
            if( value == NULL ){
                fprintf(stderr, "Error: Keyval pair \"%s\" (line %d) lacking equality sign\n", kword, linenumber);
                if( ExitOnError ){ exit(EXIT_FAILURE); }
                continue;
            }

        // strip keyword and value from leading and trailing whitespaces
            rmSurroundingWhitespaces(kword);
            if( kword == NULL ){
                fprintf(stderr, "Error: Empty keyword on line %d\n", linenumber);
                if( ExitOnError ){ exit(EXIT_FAILURE); }
                continue;
            }

            rmSurroundingWhitespaces(value);
            if( value == NULL ){
                fprintf(stderr, "Error: Empty value on line %d\n", linenumber);
                if( ExitOnError ){ exit(EXIT_FAILURE); }
                continue;
            }
            if( value[0] == '=' ){
                fprintf(stderr, "Error: Keyword \"%s\" on line %d, value starts with equality sign\n", kword, linenumber);
                if( ExitOnError ){ exit(EXIT_FAILURE); }
                continue;
            }

        // iterate over list of keywords and check if found word is keyword

            int keyword_recognised = 0;
            for(int i = 0; kwlist[i].keyword; ++i){

                if( strcasecmp(kwlist[i].keyword, kword) == 0 ){

                // mark keyword as recognised
                    ++keyword_recognised;

                // increment set number
                    ++kwlist[i].set;

                // allocate memory for value buffer
                    kwlist[i].value = realloc(kwlist[i].value, kwlist[i].set * sizeof(char*));
                    if(kwlist[i].value == NULL){ perror("kwlist->value"); exit(errno); }

                    kwlist[i].value[kwlist[i].set - 1] = malloc( (_MaxValLength_) * sizeof(char));
                    if(kwlist[i].value[kwlist[i].set - 1] == NULL){ perror("kwlist->value[i]"); exit(errno); }

                // copy value to kwlist array and ensure '\0' termination
                    strncpy( (kwlist[i].value[kwlist[i].set - 1]), value, (_MaxValLength_));
                    if( kwlist[i].value[kwlist[i].set - 1][_MaxValLength_ - 1] != '\0'){
                        fprintf(stderr, "Error: in keyword \"%s\". Value exceeds maximum length of %d characters\n"
                                , kword, (_MaxValLength_-1)
                        );
                        exit(EXIT_FAILURE);
                    }

                // increment number of successfully set keywords and break loop
                    ++n_keywords_set;
                    break;
                }
            }

        // print a warning if keyword was not recognised
            if( !keyword_recognised ){
                fprintf(stderr, "Ignoring unrecognised keyword \"%s\" on line %d\n", kword, linenumber);
                if( ExitOnError ){ exit(EXIT_FAILURE); }
            }

        }
    }

// close buffers and free memory
    fclose(fd);   fd     = NULL;
    free(buffer); buffer = NULL;

    return n_keywords_set;
}


//#include <ctype.h>
//#include <string.h>
static void rmSurroundingWhitespaces(char* string){
// remove leading whitespaces:
//  store index of first non-whitespace character and shift all characters
//  overwrite last character with '\0' (increment happens at end of loop body)
    int index, i;
    for(index = 0; isspace(string[index]); ++index);
    for(  i   = 0; string[index+i]; ++i){ string[i] = string[index+i]; }
    string[i] = '\0';

// overwrite trailing whitespaces with '\0'
    for(int end = i-1; isspace(string[end]) && (end >= 0); --end){ string[end] = '\0'; }

// if string is empty set it to NULL
    if( strlen(string) == 0 ){ string = NULL; }
}
