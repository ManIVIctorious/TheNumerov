
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "settings.h"

settings SetDefaultSettings(void);

settings SetDefaultSettings(void){

    settings defaults = (const struct settings) {
        .dimension = 2,           // dimension of the problem
        .n_stencil = 9,           // 1D stencil size
        .n_spline  = 0,           // Number of interpolation points

        .masses_string = NULL,    // String containing reduced masses
        .masses = NULL,           // Actual reduced masses array

        .kJpermol_to_oue = 1.0/4.184,       // kJ/mol -> output unit of energy
        .epot_to_oue = 1.0,                 // E_pot  -> output unit of energy
        .InvInertia_to_molpergAasq = 1.0,   // [reciprocal moment of inertia] -> mol/(g.Å^2)
        .InvInertiaThreshold = 0.0,         // if not zero set all mu above and below to +/-Threshold
        .dip_to_Asm  = 3.3356E-30,          // Conversion from input dipole moment to A.s.m
        .threshold   = 1.0E-10,             // Spacing threshold for double comparison

    // Flags
        .frequencies   = 1,
        .analyse       = 0,
        .dipole        = 0,
        .check_spacing = 1,
        .periodic      = 0,

    // Eigensolver specific values
        .Eigensolver = 2,     // Which eigensolver to use (1 FEAST, 2 ARPACK)
        .n_out = 8,           // Number of output eigenstates (ARPACK)
        .e_min =   0.0,       // Minimal energy of eigenstates (FEAST)
        .e_max = 400.0,       // Maximal energy of eigenstates (FEAST)

    // Files
    //  input file
        .input_file        = NULL,
    //  external dipole file
        .ext_dip_file      = NULL,
    //  coriolis file
        .coriolis_file     = NULL,
    //  output file
        .output_file       = "/dev/stdout",
    };

    return defaults;
}
