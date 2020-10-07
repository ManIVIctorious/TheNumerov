
#include <stdio.h>

#include "constants.h"

// provided prototypes
int TextOut_Frequencies(FILE* fd, double ekin_factor, int n_out, double* E);
 
// Output eigenvalues and calculate Frequencies
int TextOut_Frequencies(FILE* fd, double ekin_factor, int n_out, double* E){

    int i, j;
    const double kJpermolToWavenumber = 10.0 / (avogadro*planck*lightspeed); // cm^-1 / (kJ/mol)

// output eigenvalues
    fprintf(fd, "# Eigenvalues:");
    for(i = 0; i < n_out; ++i){
        fprintf(fd, " %24.16lf", E[i]);
    }

// and output frequencies
    fprintf(fd, "\n#\n# Frequencies:\n#\n#");
    for(i = 0; i < (n_out - 1); ++i){
        fprintf(fd, "       %7d", i);
    }
    for(i = 1; i < n_out; i++){
        fprintf(fd, "\n# %3d", i);

        for(j = 0; j < i; j++){
            fprintf(fd, "  % 12.5e", (E[i] - E[j]) * kJpermolToWavenumber / ekin_factor);
        }
    }
    fprintf(fd, "\n#\n#");

    return 0;
}
