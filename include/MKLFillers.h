#ifndef _MKL_FILLERS_H
#define _MKL_FILLERS_H

#include "settings.h"

// non-periodic
int MKL_FillAMatrix1D(int n_stencil,   int* nq, double* v, double ekin_to_oue, double* stencil, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);
int MKL_FillAMatrix2D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, double** q, double dq, double*** mu, double** zeta, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);
int MKL_FillAMatrix3D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, double** q, double dq, double*** mu, double** zeta, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);

// periodic
int MKL_FillPeriodicAMatrix1D(int n_stencil,   int* nq, double* v, double ekin_to_oue, double* stencil, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);
int MKL_FillPeriodicAMatrix2D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);
int MKL_FillPeriodicAMatrix3D(settings* prefs, int* nq, double* v, double ekin_to_oue, double* stencil, MKL_INT* *rows_A, MKL_INT* *cols_A, double* *vals_A);

#endif
