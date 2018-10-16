
#ifndef DEFAULT_SETTINGS
#define DEFAULT_SETTINGS

settings defaults = (const struct settings) {
  .dimension = 2,           // dimension of the problem
  .n_stencil = 9,           // 1D stencil size
  .n_spline  = 0,           // Number of interpolation points

  .masses_string = NULL,    // String containing reduced masses
  .masses_string_set = 0,
  .masses = NULL,           // Actual reduced masses array

  .ekin_factor = 1.0/4.184, // Kinetic     energy factor
  .epot_factor = 1.0,       // Potential   energy factor
  .mu_factor   = 1.0E20 * AVOGADRO*AVOGADRO * PLANCK*PLANCK/(4.0*M_PI*M_PI),   // Watson term energy factor
  .DipToAsm    = 3.3356E-30,// Conversion from input dipole moment to A.s.m
  .threshold   = 1.0E-10,   // Spacing threshold for double comparison

// Flags
  .analyze = 0,
  .dipole  = 0,
  .check_spacing = 1,

// Eigensolver specific values
  .Eigensolver = 2,     // Which eigensolver to use (1 FEAST, 2 ARPACK)
  .n_out = 8,           // Number of output eigenstates (ARPACK)
  .e_min =   0.0,       // Minimal energy of eigenstates (FEAST)
  .e_max = 400.0,       // Maximal energy of eigenstates (FEAST)

// Files
  .input_file        = NULL,
  .input_file_set    = 0,

  .coriolis_file     = NULL,
  .coriolis_file_set = 0,

  .output_file       = "/dev/stdout",
  .output_file_set   = 0,
};

#endif
