
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>

#include "data.h"
#include "settings.h"
#include "gitversion.h"

// dependencies
settings SetDefaultSettings(void);
void GetSettingsGetopt(settings* prefs, data* data, int argc, char** argv);

void SetZetaAndDimension(settings* set, data* data);
int InputMassesFile(char* inputfile, double* *m);
int ProcessModeFiles(settings *prefs, data *data);
int ProcessFileList(settings *prefs, data *data);

// output
void PrintCoriolisCoefficients(settings* set, data* data);


int main(int argc, char **argv){

    struct data data;
    data.dimension = 0;
    settings prefs = SetDefaultSettings();
    GetSettingsGetopt(&prefs, &data, argc, argv);

    if( prefs.zeta_x ){
        SetZetaAndDimension(&prefs, &data);
        if( !prefs.masses_file ){
            fprintf(stderr,
                "\n (-) Error: When the Coriolis coefficients are directly"
                "\n     provided via the command line the atomic masses"
                "\n     have to be made available this way as well."
                "\n     Aborting...\n\n"
            );
            exit(EXIT_FAILURE);
        }
    }

// return early if dimension is too small
    if(data.dimension < 2){
        fprintf(stderr,
                "\n (-) Error: At least two modes have to be considered"
                "\n     Please check your input, aborting...\n\n"
            );
        exit(EXIT_FAILURE);
    }

// if a masses file is set get the atom_masses from there
    if( prefs.masses_file ){
        data.atom_masses = NULL;
        data.n_atoms = InputMassesFile( prefs.masses_file, &data.atom_masses );
    }

//--------------------------------------------------------------------------------------
//  Coriolis Coefficients Zeta  Coriolis Coefficients Zeta  Coriolis Coefficients Zeta
//--------------------------------------------------------------------------------------
// input mode files:
//--------------------------------------------------
//  Read mode files, get n_atoms, atom_masses and total_mass then
//  normalise the mode files and calculate the Coriolis coefficients zeta
    if( !prefs.zeta_x ){
        ProcessModeFiles(&prefs, &data);
    }

// calculate the system's total mass
    data.tot_mass = 0.0;
    for(int i = 0; i < data.n_atoms; ++i){
        data.tot_mass += data.atom_masses[i];
    }

//--------------------------------------------------------------------------------------
//       Output preparations        Output preparations        Output preparations
//--------------------------------------------------------------------------------------
// write Coriolis coefficients and header/key only if:
//  1. output file-mode is write, or
//  2. the file does not exist yet
// This way, even in append mode, the header should only be written once
    int print_header = 0;
    if(    (strcmp(prefs.output_fmode, "w") == 0)
        || (access(prefs.output_file, F_OK) != 0)
    ){
        print_header = 1;
    }

// if output file is set open it, else print to default (stdout)
    if( prefs.output_file ){ prefs.fdout = fopen(prefs.output_file, prefs.output_fmode); }
    if( prefs.fdout == NULL ){ perror(prefs.output_file); exit(errno); }


// output Coriolis coefficients and return earyl
// if no geometry/coordinate files were provided
//--------------------------------------------------
    if( print_header ){
    // Print version information to header
        fprintf(prefs.fdout, "#Version: \"%s\"\n", gitversion);
        PrintCoriolisCoefficients(&prefs, &data);
    }

    if( !prefs.input_coordinates ){
        fclose(prefs.fdout); prefs.fdout = NULL;
        for(int i = 0; i < 3; ++i){ free(data.zeta[i]); data.zeta[i] = NULL; } free(data.zeta);
        data.zeta = NULL;

        return EXIT_SUCCESS;
    }


//--------------------------------------------------------------------------------------
//   Moment of inertia    Moment of inertia    Moment of inertia    Moment of inertia
//--------------------------------------------------------------------------------------
// print header/key
//--------------------------------------------------
    if( print_header ){
        fprintf(prefs.fdout, "#");
        for(int i = 0; i < data.dimension; ++i){
            fprintf(prefs.fdout, "\t q[%2d]              ", i);
        }
        for(int i = 0; i < 3; ++i){
        for(int j = i; j < 3; ++j){
            fprintf(prefs.fdout, "\t mu_%c%c             ", "xyz"[i], "xyz"[j]);
        }
        }
        fprintf(prefs.fdout, "\n");
    }

/* read and process input coordinates list:
//--------------------------------------------------
//{{{
    This file is of the following structure:
        <dimension> columns representing the displacement from the minimum geometry
            1       column  containing the path to a file containing the system's geometry

    For each file in the list the Effective Reciprocal Moment of Inertia Tensor
    is evaluated and printed to the output file.

//}}}*/
    ProcessFileList(&prefs, &data);

// close files and free unused memory
    fclose(prefs.fdout);     prefs.fdout = NULL;
    for(int i = 0; i < 3; ++i){ free(data.zeta[i]); data.zeta[i] = NULL; }
    free(data.zeta);        data.zeta        = NULL;
    free(data.atom_masses); data.atom_masses = NULL;

    return EXIT_SUCCESS;
}
