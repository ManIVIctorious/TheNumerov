
#include <ctype.h>
#include <string.h>

void ThrowInputError(char* inputfile, int linenumber, char* format,...);

// Provided prototypes
char * PreprocessBuffer(char* inputfile, int linenumber, char* buffer, const char* comment);

char * PreprocessBuffer(char* inputfile, int linenumber, char* buffer, const char* comment){

// check if line is complete (search for newline char)
    int control = 0;
    for(long unsigned i = 0; i < strlen(buffer); ++i){
        if(buffer[i] == '\n') { control = 1; }
    }
    if(control == 0){
        ThrowInputError(inputfile, linenumber,
            "Line too long (no newline char '\\n' found)."
        );
    }

// strip buffer from comments
    char * stringp = buffer;
    char * pos     = strsep(&stringp, comment);
    if(pos == NULL){ return pos; }

// remove leading whitespaces
    while( isspace( *pos ) && (*pos != '\0') ){ pos++; }
    if(strlen(pos) < 1){ pos = NULL; }

    return pos;
}
