#ifndef _TYPEDEFINITIONS_H
#define _TYPEDEFINITIONS_H

#define _MaxSettingsStringLength_ 1024

typedef struct settings {
  int dimension;        // dimension of the problem
  int n_stencil;        // stencil size per dimension
  int n_spline;         // Number of interpolation points

  char   * masses_string;   // String containing reduced masses
  double * masses;          // Actual reduced masses array

  double ekin_factor;   // Kinetic energy factor
  double epot_to_oue;   // Conversion factor of potential to output unit of energy
  double mu_factor;     // Watson term energy factor
  double DipToAsm;      // Conversion from input dipole moment to A.s.m
  double threshold;     // Spacing threshold for double comparison

// Flags
  char frequencies;     // output frequencies
  char analyse;         // output additional information
  char dipole;          // read dipole moment and output transition moment integral and oscillator strength
  char check_spacing;   // check coordinate spacing
  char periodic;        // use matrix fillers for periodic boundary conditions

// Eigensolver specific values
  int Eigensolver;      // Which eigensolver to use (1 FEAST, 2 ARPACK)
  int    n_out;         // Number of output eigenstates (ARPACK)
  double e_min;         // Minimal energy of eigenstates (FEAST)
  double e_max;         // Maximal energy of eigenstates (FEAST)

// Files
  char * input_file;    // data input file

  // optional files
  char * ext_dip_file;  // coordinates and dipole moment in external file
  char * coriolis_file; // Coriolis coefficients and inverse moment of inertia tensor for each data point
  char * output_file;   // output file

} settings;

#endif
