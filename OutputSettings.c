
#include <stdio.h>
#include "typedefinitions.h"

// provided prototypes
int OutputSettings(FILE *fd, settings prefs);


int OutputSettings(FILE *fd, settings prefs){

    char * Eigensolver = NULL;

    switch(prefs.Eigensolver){

        case 1:
            Eigensolver = "Intel_MKL_FEAST";
            break;

        case 2:
            Eigensolver = "ARMADILLO_ARPACK";
            break;

        default:
            Eigensolver = "unknown";
    }
    

    fprintf(fd,
        "#------------------------------------------------------------------------------\n"
        "#  Settings   Settings   Settings   Settings   Settings   Settings   Settings\n"
        "#------------------------------------------------------------------------------\n"
    );

fprintf(fd,
"## General settings:\n");
// integer values
    fprintf(fd, "#\tDimensionality = % d;\n", prefs.dimension);
    fprintf(fd, "#\tStencil_Size   = % d;\n", prefs.n_stencil);
    if(prefs.n_spline == 0){
        fprintf(fd, "#\tInterpolation  = %s;\n", "false");
    }else{
        fprintf(fd, "#\tN_Inter_Points = % d;\n", prefs.n_spline);
    }
// double values
    fprintf(fd, "#\n");
    fprintf(fd, "#\tReduced_Mass      = % 12.6lf; # in g/mol\n", prefs.mass);
    fprintf(fd, "#\tKin_E_Factor      = % 12.6lf; # in x per kJ/mol\n", prefs.ekin_factor);
    fprintf(fd, "#\tPot_E_Factor      = % 12.6lf; # input -> output unit of energy\n", prefs.epot_factor);
    fprintf(fd, "#\tIMOI_Mu_Factor    = % 12.6lf; # input (mu) -> kJ/mol\n", prefs.mu_factor);
    fprintf(fd, "#\tSpacing_Threshold = % 12.6lf;\n", prefs.threshold);
    fprintf(fd, "#\n#\n");


fprintf(fd,
"## Flags:\n");
    fprintf(fd, "#\tAnalyze        = %s;\n", prefs.analyze        ? "true" : "false");
    fprintf(fd, "#\tDipole         = %s;\n", prefs.dipole         ? "true" : "false");
    fprintf(fd, "#\tCheck_Spacing  = %s;\n", prefs.check_spacing  ? "true" : "false");
    fprintf(fd, "#\n#\n");


fprintf(fd,
"## Eigensolver specific settings:\n");
    fprintf(fd, "#\tEigensolver    = %s;\n", Eigensolver);
    if(prefs.Eigensolver == 1){
        fprintf(fd, "#\tLower_Bound    = % 12.6lf;\n", prefs.e_min);
        fprintf(fd, "#\tUpper_Bound    = % 12.6lf;\n", prefs.e_max);
    }
    if(prefs.Eigensolver == 2){
        fprintf(fd, "#\tN_Eigenstates  = %d;\n", prefs.n_out);
    }
    fprintf(fd, "#\n#\n");


fprintf(fd,
"## Files:\n");
    fprintf(fd, "#\tInput_File     = %s;\n", prefs.input_file);
    if(prefs.coriolis_file != NULL){
        fprintf(fd, "#\tCoriolis_File  = %s;\n", prefs.coriolis_file);
    }
    fprintf(fd, "#\tOutput_File    = %s;\n", prefs.output_file);
    fprintf(fd, "#\n#\n");


    fprintf(fd,
        "#------------------------------------------------------------------------------\n"
        "#   End Settings   End Settings   End Settings   End Settings   End Settings\n"
        "#------------------------------------------------------------------------------\n"
        "#\n#\n"
    );

    return 0;
}
