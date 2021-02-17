
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"

// provided prototypes
settings SetDefaultSettings();


settings SetDefaultSettings(){

    settings defaults = (const struct settings) {

        .threshold    = 1.0E-10,
        .modelist     = NULL,
        .input_coordinates = NULL,
        .output_file  = NULL,
        .output_fmode = "w",
        .fdout        = stdout,

    };

    return defaults;
}
