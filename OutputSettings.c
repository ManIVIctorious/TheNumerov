
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pwd.h>
#include <time.h>
#include "typedefinitions.h"
#include "gitversion.h"

// provided prototypes
int OutputSettings(FILE* fd, settings prefs);


int OutputSettings(FILE* fd, settings prefs){

    int i;
    char * Eigensolver = NULL;
    time_t current_time = time(NULL);
    char * Hostname = malloc(16 * sizeof(char));


    switch(prefs.Eigensolver){

        case 1:
            Eigensolver = "Intel_MKL_FEAST";
            break;

        case 2:
            Eigensolver = "Armadillo_ARPACK";
            break;

        default:
            Eigensolver = "unknown";
    }


fprintf(fd,
    "##-----------------------------------------------------------------------------------------\n"
    "##  Settings   Settings   Settings   Settings   Settings   Settings   Settings   Settings\n"
    "##-----------------------------------------------------------------------------------------\n"
);

// General Information
fprintf(fd,
"## General Information:");
  if(strlen(gitversion) > 0){
    fprintf(fd, "\n##\tGit Revision:\t%s", gitversion);
  }

  if(gethostname(Hostname, 16) == 0){
    struct passwd *p = getpwuid(getuid());
    fprintf(fd, "\n##\tSystem:      \t%s@%s (UID:%d)", p->pw_name, Hostname, getuid());
  }

    fprintf(fd, "\n##\tUnix Epoch:  \t%ld", current_time);
    fprintf(fd, "\n##\tDate & Time: \t%s#", ctime(&current_time));


fprintf(fd,
"\n#\n#\n"
"## General settings:");
// integer values
    fprintf(fd, "\n#\tDimensionality       = % d;", prefs.dimension);
    fprintf(fd, "\n#\tStencil_Size         = % d;", prefs.n_stencil);
  if(prefs.n_spline == 0){
    fprintf(fd, "\n#\tInterpolation_Points = %s;", "none");
  }else{
    fprintf(fd, "\n#\tInterpolation_Points = % d;", prefs.n_spline);
  }
// double values
    fprintf(fd, "\n#");

    fprintf(fd, "\n#\tKin_E_Factor         = % 12.8lf;\t# in x per kJ/mol", prefs.ekin_factor);
    fprintf(fd, "\n#\tPot_E_Factor         = % 12.8lf;\t# input (v)  -> output unit of energy", prefs.epot_factor);
    fprintf(fd, "\n#\tIMOI_Factor          = % 12.8le;\t# input (mu) -> kJ/mol", prefs.mu_factor);
    fprintf(fd, "\n#\tSpacing_Threshold    = % 12.8le;", prefs.threshold);
    fprintf(fd, "\n#");

    fprintf(fd, "\n#\t# Reduced mass of each dimension,");
    fprintf(fd, "\n#\t# given as colon separated string, with entries in g/mol");
    fprintf(fd, "\n#\tReduced_Masses = %.12lf", prefs.masses[0]);
    for(i = 1; i < prefs.dimension; ++i){
        fprintf(fd, ":%.12lf", prefs.masses[i]);
    }
    fprintf(fd, ";");


fprintf(fd,
"\n#\n#\n"
"## Flags:");
    fprintf(fd, "\n#\tAnalyze        = %s;", prefs.analyze        ? "true" : "false");
    fprintf(fd, "\n#\tDipole         = %s;", prefs.dipole         ? "true" : "false");
    fprintf(fd, "\n#\tCheck_Spacing  = %s;", prefs.check_spacing  ? "true" : "false");


fprintf(fd,
"\n#\n#\n"
"## Eigensolver specific settings:");
    fprintf(fd, "\n#\tEigensolver    = %s;", Eigensolver);
  if(prefs.Eigensolver == 1){
    fprintf(fd, "\n#\tLower_Bound    = % 12.8lf;", prefs.e_min);
    fprintf(fd, "\n#\tUpper_Bound    = % 12.8lf;", prefs.e_max);
  }
  if(prefs.Eigensolver == 2){
    fprintf(fd, "\n#\tN_Eigenstates  = %d;", prefs.n_out);
  }


fprintf(fd,
"\n#\n#\n"
"## Files:");
    fprintf(fd, "\n#\tInput_File           = %s;", prefs.input_file);
  if(prefs.ext_dip_file_set){
    fprintf(fd, "\n#\tExternal_Dipole_File = %s;", prefs.ext_dip_file);
  }
  if(prefs.coriolis_file_set){
    fprintf(fd, "\n#\tCoriolis_File        = %s;", prefs.coriolis_file);
  }
    fprintf(fd, "\n#\tOutput_File          = %s;", prefs.output_file);


fprintf(fd,
    "\n#\n#\n"
    "##-----------------------------------------------------------------------------------------\n"
    "## End Settings   End Settings   End Settings   End Settings   End Settings   End Settings\n"
    "##-----------------------------------------------------------------------------------------\n"
    "#\n#\n"
);

    return 0;
}
