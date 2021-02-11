
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"

// provided prototypes
settings SetDefaultSettings();


settings SetDefaultSettings(){

    settings defaults = (const struct settings) {

        .threshold   = 1.0E-10,
        .modelist    = NULL,
        .dimension   = 0,
        .input_coordinates = NULL,
        .output_file = NULL,
        .fdout       = stdout,

    };

    return defaults;
}
