
#include <stdio.h>
#include <stdlib.h>

#include "data.h"
#include "settings.h"

// Provided prototypes
void PrintCoriolisCoefficients(settings* set, data* data);

void PrintCoriolisCoefficients(settings* set, data* data){

// print header identifying the mode combinations
    if( set->modelist ){
        fprintf(set->fdout, "#Modes:");
        for(int i = 0; i < data->dimension; ++i){
            for(int j = i+1; j < data->dimension; ++j){
                fprintf(set->fdout, "\t %s|%s", set->modelist[i], set->modelist[j]);
            }
        }
        fprintf(set->fdout, "\n");
    }

// print zeta values
    for(int k = 0; k < 3; ++k){
    // identifier
        fprintf(set->fdout, "Zeta_%c:", "xyz"[k]);

    // upper triangle of values
        for(int i = 0; i < data->dimension; ++i){
            for(int j = i+1; j < data->dimension; ++j){
                fprintf(set->fdout, "\t% .12le", data->zeta[k][i*data->dimension + j]);
            }
        }
        fprintf(set->fdout, "\n");
    }
}
