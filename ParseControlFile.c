
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "typedefinitions.h"
#include "controlfile.h"

// Provided prototypes
int ControlFileParser(char* filename, struct keywords* keywordlist, int verbosity_flag);


int ControlFileParser(char* filename, struct keywords* keywordlist, int verbosity_flag){

    long unsigned int i, j;
    FILE * fd     = NULL;
    char * buffer = NULL;
    const char * comment = "%#\n";
    char *token, *pos1, *pos2, *pos3;

    int recognized_keyword;
    int n_keywords_set;
    int linenumber;
    int control;

    fd = fopen(filename, "r");
    if(fd == NULL){
        fprintf(stderr,
            "\n (-) Error opening input-file: \"%s\""
            "\n     Aborting..."
            "\n\n"
            , filename
        );
        exit(EXIT_FAILURE);
    }

    buffer = malloc((_MaxLineLength_) * sizeof(char));
    if(buffer == NULL){
        fprintf(stderr,
            "\n (-) Error in memory allocation of buffer in control file input"
            "\n     Aborting..."
            "\n\n"
        );
        exit(EXIT_FAILURE);
    }


    linenumber = 0;
    n_keywords_set = 0;
    while( (fgets(buffer, _MaxLineLength_, fd)) != NULL ){

        linenumber++;

    // check for existence of newline character,
    //  if not found the line is not fully inside of the buffer
    //  and therefore exceeding maximum line length
        for(i = 0, control = 0; i < strlen(buffer); ++i){
            if(buffer[i] == '\n'){
                control = 1;
            }
        }
        if(control == 0){
            fprintf(stderr,
                "\n (-) Error in control file \"%s\", line \"%d\" is too long."
                "\n     Aborting..."
                "\n\n", filename, linenumber
            );
            exit(1);
        }

    // strip buffer from comments
        token  = buffer;
        pos1 = strsep(&token, comment);
        if(pos1 == NULL) continue;

    // remove leading white spaces and empty lines
        while(isspace(*pos1) != 0 && *pos1 != '\0') { pos1++; }
        if(strlen(pos1) < 1) continue;

// buffer now contains a full (non empty) line of the control file, stripped of
//  comments and *pos1 points to the first, non white space character of buffer
//--------------------------------------------------------------------------------


    // split at every occurrence of semi colon
    //  &token --> token --> char array;
    //                       ↑
    //                      pos1
        token = pos1;

        while( (pos1 = strsep(&token, ";")) != NULL ){

        // remove leading white spaces
            while(isspace(*pos1) != 0) { pos1++; }

            if(strlen(pos1) > 0){

            // search for an equality sign after word:
                pos3 = pos2 = pos1;
                while( *pos3 != '=' && strlen(pos3) > 0 ){
                // first set pos2 to start of word and iterate to its end
                    pos2 = pos1;
                    while( isspace(*pos2) == 0 && *pos2 != '=' && *pos2 != '\0') { pos2++; }

                // then move pos3 forwards over eventual white spaces
                    pos3 = pos2;
                    while( isspace(*pos3) != 0 && *pos3 != '=' && *pos3 != '\0') { pos3++; }

                // and check if pos3 now points to an equality sign
                //  if not ignore the word and check if there is a keyword
                //  afterwards (set pos1 to new word start and reiterate loop)
                    if( *pos3 != '=' ){ pos1 = pos3; }
                }

                if( *pos3 == '='){
                // set char after end of word to \0 (delimit *pos1)
                    *pos2 = '\0';

                    i = 0;
                    recognized_keyword = 0;
                // iterate over list of keywords and check if found word is keyword
                    while( (keywordlist [i] .keyword) != NULL ){

                    // first check (case insensitive) name of keyword
                        if(strncasecmp((keywordlist[i].keyword), pos1, strlen(keywordlist[i].keyword)) == 0){

                        // additionally verify word length and keyword length are the same
                            if( strlen(keywordlist [i] .keyword) == strlen(pos1) ){

                            //  mark word as recognized (set recognized flag to 1)
                                recognized_keyword = 1;

                            // check if keyword has already been set, if yes throw a warning
                                if(keywordlist[i].set > 0){
                                    fprintf(stderr,
                                        " (-) Warning: Control file \"%s\", line \"%d\": Keyword \"%s\" has already been set, overwriting content\n"
                                        , filename, linenumber, pos1
                                    );
                                }

                            // remove leading and trailing white spaces from value:
                                if(strlen(pos3) > 0){
                                // trailing: point pos2 to end of pos3, iterate backwards
                                //  until non white space is found (since pos3 points to an
                                //  equality sign a non white space character must exist)
                                //  and set the next char to \0
                                    pos2 = pos3 + strlen(pos3) - 1;
                                    while( isspace(*pos2) != 0 && pos2 > pos3 ) { pos2--; }
                                    ++pos2;
                                    *pos2 = '\0';

                                // leading: move pos3 to start of value
                                    ++pos3;
                                    while( isspace(*pos3) != 0 && *pos3 != '\0' ) { pos3++; }

                                }

                            // make sure value (pointed to by pos3) is not empty
                                if( strlen(pos3) < 1 ){
                                    fprintf(stderr,
                                        "\n (-) Error in control file \"%s\", line \"%d\""
                                        "\n     Keyword \"%s\" does not contain any information."
                                        "\n     Please check your input."
                                        "\n     Aborting..."
                                        "\n\n", filename, linenumber, keywordlist[i].keyword
                                    );
                                    exit(1);
                                }


                            // check if <value> contains any white spaces and throw a warning if true
                                for(j = 0, control = 0; j < strlen(pos3); ++j){
                                    if(isspace(pos3[j]) != 0){
                                        control = 1;
                                    }
                                }
                                if(control == 1){
                                    fprintf(stderr,
                                        " (-) Warning: Control file \"%s\", line \"%d\": Value of keyword \"%s\" contains white spaces\n"
                                        , filename, linenumber, keywordlist[i].keyword
                                    );
                                }

                            // pos3 now contains the value in raw character array form
                            //  copy value to keywordlist array and ensure \0 termination
                                strncpy( (keywordlist [i] .value), pos3, (_MaxLineLength_ - 1));
                                keywordlist[i].value[_MaxLineLength_ - 1] = '\0';


                            // if verbosity flag is set: output found keywords and values
                                if(verbosity_flag != 0){
                                    fprintf(stderr, "\t%s\t\t= \"%s\"\n", pos1, pos3);
                                }


                            // since keyword is found:
                            //  mark keyword as set (increment keywordlist[i].set)
                                ++keywordlist[i].set;
                            //  increment number of found keywords and
                                ++n_keywords_set;
                            //  stop iteration over keywords
                                break;
                            }
                        }
                        ++i;
                    }

                // check if keyword was recognized, else throw a warning
                    if(recognized_keyword == 0){
                        fprintf(stderr,
                            " (-) Warning: Control file \"%s\", line \"%d\": Ignoring unrecognized keyword \"%s\"\n"
                            , filename, linenumber, pos1
                        );
                    }
                }
            }
        }
    }

// free memory
    free(buffer); buffer = NULL;

    return n_keywords_set;
}
