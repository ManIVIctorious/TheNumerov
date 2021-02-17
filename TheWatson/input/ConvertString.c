
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <limits.h>
#include <math.h>

#include "ConvertString.h"

/*  Convert string to bool, long integer or double
 *
 *    All three functions get two char pointers and an integer pointer
 *    as arguments (optarg, varname and control)
 *
 *    optarg:   the string to be converted to a number
 *
 *    varname:  optional variable name, if not NULL print it in front
 *              of possible errors for easier identification.
 *
 *    control:  if control is NULL let the function handle the error,
 *              else set its value to errvar (either errno or):
 *                  1     No conversion
 *                  2     Only part of string was converted to number
 *
 *    Depending on control the function either returns or exits (NULL)
 *    the program on first encountered error therefore, if control is
 *    set varname has no effect.
 */

long int convertstring_to_bool(char* optarg, char* varname, int *control){

// set value to 1 if optarg equals to "true" (case insensitive)
// set it to 0 if optarg equals to "false" and else assume it to
// be set as a number. In this case do a integer comparison and
// return the appropriate value.
    if      ( strcasecmp( optarg, "true"  ) == 0 ){ return 1; }
    else if ( strcasecmp( optarg, "yes"   ) == 0 ){ return 1; }
    else if ( strcasecmp( optarg, "false" ) == 0 ){ return 0; }
    else if ( strcasecmp( optarg, "no"    ) == 0 ){ return 0; }
    else{ return convertstring_to_long(optarg, varname, control); }

}


long int convertstring_to_long(char* optarg, char* varname, int *control){

    char * endptr  = NULL;
// set error variables
    errno = 0;
    int exitonerror;
    if( control ){ exitonerror = 0; *control = 0; }
    else         { exitonerror = 1; }

    long int value = strtol(optarg, &endptr, 10);

/* check for various possible errors */

// general error (also catches out of bound errors)
    if( errno != 0 ){
        if( !exitonerror ){ *control = errno; return value; }

        if( varname ){ perror(varname);  }
        else         { perror("strtol"); }
        exit(errno);
    }

// check if conversion did actually happen
    if( endptr == optarg ){
        if( !exitonerror ){ *control = 1; return value; }

        if( varname ){ fprintf(stderr, "%s:  ", varname); }
        fprintf(stderr, "No digits were found\n");
        exit(EXIT_FAILURE);
    }

// check if whole string was converted, i.e. until '\0' termination
    if( *endptr != '\0' ){
        if( !exitonerror ){ *control = 2; return value; }

        if( varname ){ fprintf(stderr, "%s:  ", varname); }
        fprintf(stderr, "Further characters after number\n");
        exit(EXIT_FAILURE);
    }

    return value;
}


double convertstring_to_double(char* optarg, char* varname, int *control){

    char * endptr = NULL;
// set error variables
    errno = 0;
    int exitonerror;
    if( control ){ exitonerror = 0; *control = 0; }
    else         { exitonerror = 1; }

    double value = strtod(optarg, &endptr);

/* check for various possible errors */

// general error (also catches out of bound errors)
    if( errno != 0 ){
        if( !exitonerror ){ *control = errno; return value; }

        if( varname ){ perror(varname);  }
        else         { perror("strtod"); }
        exit(errno);
    }

// check if conversion did actually happen
    if( endptr == optarg ){
        if( !exitonerror ){ *control = 1; return value; }

        if( varname ){ fprintf(stderr, "%s:  ", varname); }
        fprintf(stderr, "No digits were found\n");
        exit(EXIT_FAILURE);
    }

// check if whole string was converted, i.e. until NULL termination
    if( *endptr != '\0' ){
        if( !exitonerror ){ *control = 2; return value; }

        if( varname ){ fprintf(stderr, "%s:  ", varname); }
        fprintf(stderr, "Further characters after number\n");
        exit(EXIT_FAILURE);
    }

    return value;
}


int splitstring_to_array(char* array, char** *stringlist, const char* delimit, size_t preserve_array){

    char * pos = array;
    char * stringp = pos;

// if preserve_array is set copy the content of array into an auxiliary array
    if( preserve_array ){
        pos = malloc( preserve_array * sizeof(char) );
        if( !pos ){ perror("auxiliary array in split_string_array"); exit(errno); }
        strncpy(pos, array, (preserve_array-1) );
    }

    int entries = 0;

    while( 1 ){
    // get token
        do{
            pos = strsep(&stringp, delimit);
        }while ( (pos != NULL) && (*pos == '\0') );

        if( pos == NULL ){ break; }

    // get size of element and increase counter
        size_t elementlength = strlen(pos) + 1;
        entries++;

    // add additional field to stringlist
        (*stringlist) = realloc( (*stringlist), entries*sizeof(char*) );
        if( !(*stringlist) ){ perror("stringlist"); exit(errno); }

        (*stringlist)[entries-1] = malloc( elementlength*sizeof(char) );
        if( (*stringlist)[entries-1] == NULL ){ perror("stringlist[i]"); exit(errno); }

    // copy entry into stringlist element
        strncpy( (*stringlist)[entries-1], pos, (elementlength-1) );
    }

    if( preserve_array ){ free(pos); }

    return entries;
}
