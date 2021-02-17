
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"

// provided prototypes
settings SetDefaultSettings();


settings SetDefaultSettings(){

    settings defaults = (const struct settings) {

        .modelist          = NULL,
        .masses_file       = NULL,
        .masses_from_modes = 1,
        .input_coordinates = NULL,

        .output_file  = NULL,
        .output_fmode = "w",
        .fdout        = stdout,

        .threshold    = 1.0E-10,

        .zeta_x = NULL,
        .zeta_y = NULL,
        .zeta_z = NULL

    };

    return defaults;
}
