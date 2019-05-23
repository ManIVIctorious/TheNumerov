
#include <stdlib.h>
#include "settings.h"

// provided prototypes
settings SetDefaultSettings();


settings SetDefaultSettings(){

    settings defaults = (const struct settings) {

        .dimension = 0,
        .threshold = 10E-10,

        .modelist = NULL,

        .output_file = "/dev/stdout",
        .output_file_set = 0,

    };

    return defaults;
}
