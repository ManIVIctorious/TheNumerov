#ifndef _SETTINGS_H
#define _SETTINGS_H

#define _MaxSettingsStringLength_ 1024

typedef struct settings {
  int dimension;        // dimension of the problem
  int n_stencil;        // stencil size per dimension
  int n_spline;         // Number of interpolation points

  char   * masses_string;   // String containing reduced masses
  double * masses;          // Actual reduced masses array

  double kJpermol_to_oue;           // kJ/mol -> output unit of energy
  double epot_to_oue;               // Conversion factor of potential to output unit of energy
  double InvInertia_to_molpergAasq; // reciprocal moment of inertia -> mol/(g.Å^2)
  double InvInertiaThreshold;       // if not zero set all mu above and below to +/-Threshold
  double dip_to_Asm;                // Conversion from input dipole moment to A.s.m
  double threshold;                 // Spacing threshold for double comparison

// Flags
  char frequencies;     // output frequencies
  char analyse;         // output additional information
  char dipole;          // read dipole moment and output transition moment integral and oscillator strength
  char check_spacing;   // check coordinate spacing
  char periodic;        // use matrix fillers for periodic boundary conditions

// Eigensolver specific values
// Eigensolvers are encoded in a bit-mask, i.e.
//  1   Intel MKL FEAST
//  2   Armadillo ARPACK
//  4   Not yet implemented
  unsigned int Eigensolver; // Which eigensolver to use
  int    n_out;             // Number of output eigenstates (ARPACK)
  double e_min;             // Minimal energy of eigenstates (FEAST)
  double e_max;             // Maximal energy of eigenstates (FEAST)

// Files
  char * input_file;    // data input file

  // optional files
  char * ext_dip_file;  // coordinates and dipole moment in external file
  char * coriolis_file; // Coriolis coefficients and inverse moment of inertia tensor for each data point
  char * output_file;   // output file

} settings;

#endif
