
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <pwd.h>
#include <time.h>
#include "typedefinitions.h"
#include "gitversion.h"

// provided prototypes
int OutputSettings(FILE *fd, settings prefs);


int OutputSettings(FILE *fd, settings prefs){

    char * Eigensolver = NULL;
    time_t current_time = time(NULL);
    char * Hostname = malloc(16 * sizeof(char));


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
        "#-----------------------------------------------------------------------------------------\n"
        "#  Settings   Settings   Settings   Settings   Settings   Settings   Settings   Settings\n"
        "#-----------------------------------------------------------------------------------------\n"
    );

// General Information
fprintf(fd,
"## General Information:\n");
    if(strlen(gitversion) > 0){
        fprintf(fd, "##\tGit Revision:\t%s\n", gitversion);
    }

    if(0 == gethostname(Hostname, 16)){
        struct passwd *p = getpwuid(getuid());
        fprintf(fd, "##\tSystem:      \t%s@%s (UID:%d)\n", p->pw_name, Hostname, getuid());
    }

    fprintf(fd, "##\tUnix Epoch:  \t%ld\n", current_time);
    fprintf(fd, "##\tDate & Time: \t%s#\n", ctime(&current_time));
    fprintf(fd, "#\n#\n");


fprintf(fd,
"## General settings:\n");
// integer values
    fprintf(fd, "#\tDimensionality       = % d;\n", prefs.dimension);
    fprintf(fd, "#\tStencil_Size         = % d;\n", prefs.n_stencil);
    if(prefs.n_spline == 0){
    fprintf(fd, "#\tInterpolation_Points = %s;\n", "none");
    }else{
        fprintf(fd, "#\tInterpolation_Points = % d;\n", prefs.n_spline);
    }
// double values
    fprintf(fd, "#\n");
    fprintf(fd, "#\tReduced_Mass         = % 12.6lf; # in g/mol\n", prefs.mass);
    fprintf(fd, "#\tKin_E_Factor         = % 12.6lf; # in x per kJ/mol\n", prefs.ekin_factor);
    fprintf(fd, "#\tPot_E_Factor         = % 12.6lf; # input (v)  -> output unit of energy\n", prefs.epot_factor);
    fprintf(fd, "#\tIMOI_Factor          = % 12.6lf; # input (mu) -> kJ/mol\n", prefs.mu_factor);
    fprintf(fd, "#\tSpacing_Threshold    = % 12.5le;\n", prefs.threshold);
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
        "#-----------------------------------------------------------------------------------------\n"
        "# End Settings   End Settings   End Settings   End Settings   End Settings   End Settings\n"
        "#-----------------------------------------------------------------------------------------\n"
        "#\n#\n"
    );

    return 0;
}
