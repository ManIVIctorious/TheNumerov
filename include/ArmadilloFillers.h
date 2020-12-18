#ifndef _ARMADILLO_FILLERS_H
#define _ARMADILLO_FILLERS_H

#include "settings.h"

// non-periodic
arma::sp_mat FillArmadillo_1D(int n_stencil,   int* nq, double* v, double ekin_to_oue, double* stencil);
arma::sp_mat FillArmadillo_2D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, double** q, double dq, double*** mu, double** zeta);
arma::sp_mat FillArmadillo_3D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, double** q, double dq, double*** mu, double** zeta);
arma::sp_mat FillArmadillo_4D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, double** q, double dq, double*** mu, double** zeta);

// periodic
arma::sp_mat FillPeriodicArmadillo_1D(int n_stencil,   int* nq, double* v, double ekin_to_oue, double* stencil);
arma::sp_mat FillPeriodicArmadillo_2D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil);
arma::sp_mat FillPeriodicArmadillo_3D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil);
arma::sp_mat FillPeriodicArmadillo_4D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil);

#endif
