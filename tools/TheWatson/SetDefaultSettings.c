
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"

// provided prototypes
settings SetDefaultSettings();


settings SetDefaultSettings(){

    settings defaults = (const struct settings) {

        .dimension = 0,
        .threshold = 10E-10,

        .modelist = NULL,

        .output_file = NULL,
        .output_file_set = 0,
        .fdout = stdout,

    };

    return defaults;
}
