
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <pwd.h>
#include <time.h>

#include "settings.h"
#include "gitversion.h"

// provided prototypes
void PrintSettings(settings* prefs, FILE *fd);


void PrintSettings(settings* prefs, FILE *fd){

    time_t current_time = time(NULL);
    char * Hostname = malloc(16 * sizeof(char));

    if(fd == NULL){
    // open file for writing
        fd = fopen(prefs->output_file, "w");
        if(fd == NULL){ perror(prefs->output_file); exit(errno); }
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
    fprintf(fd, "\n#\tDimensionality       = % d;", prefs->dimension);
    fprintf(fd, "\n#\tStencil_Size         = % d;", prefs->n_stencil);
  if(prefs->n_spline == 0){
    fprintf(fd, "\n#\tInterpolation_Points = %s;", "none");
  }else{
    fprintf(fd, "\n#\tInterpolation_Points = % d;", prefs->n_spline);
  }
// double values
    fprintf(fd, "\n#");

    fprintf(fd, "\n#\tEPotToOUE            = % 12.8lf;\t# input (v)  -> output unit of energy (oue)", prefs->epot_to_oue);
    fprintf(fd, "\n#\tkJpermolToOUE        = % 12.8lf;\t# in x per kJ/mol", prefs->kJpermol_to_oue);
    fprintf(fd, "\n#\tIMOITomolpergAasq    = % 12.8le;\t# input (mu) -> mol/(g.Å^2)", prefs->InvInertia_to_molpergAasq);
if( prefs->InvInertiaThreshold > 0.0 ){
    fprintf(fd, "\n#\tIMOI_Threshold       = % 12.8le;\t# in [input mu]", prefs->InvInertiaThreshold);
}else{
    fprintf(fd, "\n#\tIMOI_Threshold       = none;\t# in [input mu]");
}
    fprintf(fd, "\n#\tDipToAsm             = % 12.8le;\t# input (dipole moment) -> A.s.m", prefs->dip_to_Asm);
    fprintf(fd, "\n#\tSpacing_Threshold    = % 12.8le;", prefs->threshold);
    fprintf(fd, "\n#");

    fprintf(fd, "\n#\t# Reduced mass of each dimension,");
    fprintf(fd, "\n#\t# given as colon separated string, with entries in g/mol");
    fprintf(fd, "\n#\tReduced_Masses = %.12lf", prefs->masses[0]);
    for(int i = 1; i < prefs->dimension; ++i){
        fprintf(fd, ":%.12lf", prefs->masses[i]);
    }
    fprintf(fd, ";");


// Boolean values
fprintf(fd,
"\n#\n#\n"
"## Flags:");
    fprintf(fd, "\n#\tFrequencies    = %s;", prefs->frequencies    ? "true" : "false");
    fprintf(fd, "\n#\tAnalyse        = %s;", prefs->analyse        ? "true" : "false");
    fprintf(fd, "\n#\tCheck_Spacing  = %s;", prefs->check_spacing  ? "true" : "false");
    fprintf(fd, "\n#\tDipole         = %s;", prefs->dipole         ? "true" : "false");
    fprintf(fd, "\n#\tPeriodic       = %s;", prefs->periodic       ? "true" : "false");


// Eigensolver
// 1    Intel MKL FEAST
// 2    Armadillo ARPACK
// 4    Not implemented yet
    char * solvername = "unknown";
    if     ( prefs->Eigensolver == 1 ){ solvername = "Intel_MKL_FEAST";  }
    else if( prefs->Eigensolver == 2 ){ solvername = "Armadillo_ARPACK"; }

fprintf(fd,
"\n#\n#\n"
"## Eigensolver specific settings:");
    fprintf(fd, "\n#\tEigensolver    = %s;", solvername);
// Intel MKL FEAST specific settings
  if(prefs->Eigensolver == 1){
    fprintf(fd, "\n#\tLower_Bound    = % 12.8lf;", prefs->e_min);
    fprintf(fd, "\n#\tUpper_Bound    = % 12.8lf;", prefs->e_max);
  }
// Armadillo ARPACK specific settings
  else if(prefs->Eigensolver == 2){
    fprintf(fd, "\n#\tN_Eigenstates  = %d;", prefs->n_out);
  }


// Files
fprintf(fd,
"\n#\n#\n"
"## Files:");
    fprintf(fd, "\n#\tInput_File           = %s;", prefs->input_file);
  if( prefs->ext_dip_file ){
    fprintf(fd, "\n#\tExternal_Dipole_File = %s;", prefs->ext_dip_file);
  }
  if( prefs->coriolis_file ){
    fprintf(fd, "\n#\tCoriolis_File        = %s;", prefs->coriolis_file);
  }
    fprintf(fd, "\n#\tOutput_File          = %s;", prefs->output_file);


fprintf(fd,
    "\n#\n#\n"
    "##-----------------------------------------------------------------------------------------\n"
    "## End Settings   End Settings   End Settings   End Settings   End Settings   End Settings\n"
    "##-----------------------------------------------------------------------------------------\n"
);

// close file
    fclose(fd); fd = NULL;
}
