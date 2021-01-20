
#include <stdio.h>

#include "constants.h"

// provided prototypes
void PrintFrequencies(FILE* fd, double kJpermol_to_oue, int n_out, double* E);

// Output eigenvalues and calculate Frequencies
void PrintFrequencies(FILE* fd, double kJpermol_to_oue, int n_out, double* E){

// calculate conversion factor from output unit of energy (oue) to wavenumbers
    const double kJpermol_to_wavenumber = 10.0 / (avogadro*planck*lightspeed); // cm^-1 / (kJ/mol)
    double oue_to_wavenumber = kJpermol_to_wavenumber / kJpermol_to_oue;       // cm^-1 / oue

// output eigenvalues in output unit of energy (oue)
    fprintf(fd, "#\n#\n# Eigenvalues:");
    for(int i = 0; i < n_out; ++i){
        fprintf(fd, " %24.16lf", E[i]);
    }

// and output frequencies in wavenumbers
    fprintf(fd, "\n#\n#\n# Frequencies in 1/cm:\n#\n#");
    for(int i = 0; i < (n_out - 1); ++i){
        fprintf(fd, "       %7d", i);
    }
    for(int i = 1; i < n_out; i++){
        fprintf(fd, "\n# %3d", i);

        for(int j = 0; j < i; j++){
            fprintf(fd, "  % 12.5e", (E[i] - E[j]) * oue_to_wavenumber);
        }
    }
    fprintf(fd, "\n");
}
