#ifndef _TYPEDEFINITIONS_H
#define _TYPEDEFINITIONS_H

#define _MaxSettingsStringLength_ 1024

typedef struct settings {
  int dimension;      // dimension of the problem
  int n_stencil;      // 1D stencil size
  int n_spline;       // Number of interpolation points

  char   * masses_string;   // String containing reduced masses
  int      masses_string_set;
  double * masses;          // Actual reduced masses array

  double ekin_factor; // Kinetic     energy factor
  double epot_factor; // Potential   energy factor
  double mu_factor;   // Watson term energy factor
  double DipToAsm;    // Conversion from input dipole moment to A.s.m
  double threshold;   // Spacing threshold for double comparison

// Flags
  char analyze;
  char dipole;
  char check_spacing;

// Eigensolver specific values
  int Eigensolver;  // Which eigensolver to use (1 FEAST, 2 ARPACK)
  int    n_out;     // Number of output eigenstates (ARPACK)
  double e_min;     // Minimal energy of eigenstates (FEAST)
  double e_max;     // Maximal energy of eigenstates (FEAST)

// Files
  char * input_file;
  char   input_file_set;

  char * ext_dip_file;
  char   ext_dip_file_set;

  char * coriolis_file;
  char   coriolis_file_set;

  char * output_file;
  char   output_file_set;

} settings;

#endif