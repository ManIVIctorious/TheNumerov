
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"

// Provided prototypes
void PrintCoriolisCoefficients(settings* set);

void PrintCoriolisCoefficients(settings* set){

// print header identifying the mode combinations
    fprintf(set->fdout, "#Modes:");
    for(int i = 0; i < set->dimension; ++i){
        for(int j = i+1; j < set->dimension; ++j){
            fprintf(set->fdout, "\t %s|%s", set->modelist[i], set->modelist[j]);
        }
    }
    fprintf(set->fdout, "\n");

// print zeta values
    for(int k = 0; k < 3; ++k){
    // identifier
        fprintf(set->fdout, "Zeta_%c:", "xyz"[k]);

    // upper triangle of values
        for(int i = 0; i < set->dimension; ++i){
            for(int j = i+1; j < set->dimension; ++j){
                fprintf(set->fdout, "\t% .12le", set->zeta[k][i*set->dimension + j]);
            }
        }
        fprintf(set->fdout, "\n");
    }
}
