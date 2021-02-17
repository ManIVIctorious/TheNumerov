
#define _GNU_SOURCE
#define _MaxLineLength_ 2048
#define _PATH_MAX_ 4096

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#include "data.h"
#include "settings.h"

// dependencies
char * PreprocessBuffer(char* inputfile, int linenumber, char* buffer, const char* comment);
void ThrowInputError(char* inputfile, int linenumber, char* format,...);
double convertstring_to_double(char* optarg, char* varname, int *control);
void EffectiveReciprocalMomentofInertia(settings *set, data* data, double* q, char* coordsfile);

// provided prototypes
int ProcessFileList(settings *set, data* data);

int ProcessFileList(settings *set, data* data){

// open input file read only
    FILE * fd = fopen(set->input_coordinates, "r");
    if( fd == NULL ){ perror(set->input_coordinates); exit(errno); }

// define comment characters
    const char * comment = "#%\n";
    const char * delimit = " \t";

// allocate memory of size _MaxLineLength_ for buffer
    char * pos     = NULL;
    char * stringp = NULL;
    char * buffer  = malloc((_MaxLineLength_) * sizeof(char));
    if( buffer == NULL ){ perror("ProcessFileList buffer"); exit(errno); }


// start parsing file
    int entry_rows = 0;
    int linenumber = 0;
    while( fgets(buffer, _MaxLineLength_, fd) ){

    // pre-process buffer
    // strip buffer off comments, leading white spaces and do some error handling
        linenumber++;
        pos = PreprocessBuffer(set->input_coordinates, linenumber, buffer, comment);
        if(pos == NULL){ continue; }

// pos now points to the first non-white-space character of buffer
//-----------------------------------------------------------------------------------

        stringp = pos;

    // allocate memory for data in FileList-line. This data will be further
    // processed in the EffectiveReciprocalMomentofInertia() function and
    // will be freed there.
        double * q = malloc(data->dimension * sizeof(double));
        char * coordsfile = malloc( (_PATH_MAX_) * sizeof(char) );

        if(q == NULL){ perror("q in ProcessFileList"); exit(errno); }
        if(coordsfile == NULL){ perror("coordsfile");  exit(errno); }


    // read deviation from minimum geometry Q from input file
        for(int i = 0; i < data->dimension; ){

        // get next token
            do{
                pos = strsep(&stringp, delimit);
            }while( (pos != NULL) && (*pos == '\0') );

        // throw an error if no data found
            if( pos == NULL ){
                ThrowInputError(set->input_coordinates, linenumber,
                    "\n     Too few entries in input line "
                    "(only found %d of the expected %d columns)"
                    , i, data->dimension + 1
                );
            }

        // save value to coordinate array
            q[i++] = convertstring_to_double(pos, "Coordinates q", NULL);
        }


    // get last token (coordinate filepath)
        do{
            pos = strsep(&stringp, delimit);
        }while( (pos != NULL) && (*pos == '\0') );

    // throw an error if column is not available
        if( pos == NULL ){
            ThrowInputError(set->input_coordinates, linenumber,
                "\n     Too few entries in input line "
                "(only found %d of the expected %d columns)"
                , data->dimension, data->dimension + 1
            );
        }

    // copy file path to coordsfile array
        strncpy( coordsfile, pos, _PATH_MAX_ );


    // After each read data line start directly with processing
    //  Determine Effective Reciprocal Moment of Inertia Tensor mu
    //------------------------------------------------------------
    #ifdef debug_coords
    //{{{
    // output coordinates for control
        fprintf(stderr, "\n%d\t%s\n", entry_rows+1, coordsfile);
        fprintf(stderr,
                "\nNumber of Atoms:\t%d"
                "\nInput coordinates:\n"
                "\t         x         "
                "\t         y         "
                "\t         z         "
                "\t    atomic mass    "
                "\n"
                , set->n_atoms
            );
    //}}}
    #endif

        EffectiveReciprocalMomentofInertia(set, data, q, coordsfile);

        ++entry_rows;
    }
    fclose(fd); fd = NULL;

    return entry_rows;
}
