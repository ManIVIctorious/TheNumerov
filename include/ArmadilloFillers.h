#ifndef _ARMADILLO_FILLERS_H
#define _ARMADILLO_FILLERS_H

#include "settings.h"

// periodic
arma::sp_mat FillPeriodicArmadillo_1D(int n_stencil,   int* nq, double* v, double ekin_to_oue, double* stencil);
arma::sp_mat FillPeriodicArmadillo_2D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil);
arma::sp_mat FillPeriodicArmadillo_3D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil);
arma::sp_mat FillPeriodicArmadillo_4D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil);

#endif
