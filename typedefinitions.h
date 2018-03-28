
typedef struct settings {
  int dimension;      // dimension of the problem
  int n_stencil;      // 1D stencil size
  int n_spline;       // Number of interpolation points

  double mass;        // Reduced mass

  double ekin_factor; // Kinetic     energy factor
  double epot_factor; // Potential   energy factor
  double mu_factor;   // Watson term energy factor
  double threshold;   // Spacing threshold for double comparison

// Flags
  int analyze;
  int dipole;
  int check_spacing;

// Eigensolver specific values
  int Eigensolver;  // Which eigensolver to use (1 FEAST, 2 ARPACK)
  int    n_out;     // Number of output eigenstates (ARPACK)
  double e_min;     // Minimal energy of eigenstates (FEAST)
  double e_max;     // Maximal energy of eigenstates (FEAST)

// Files
  char * input_file;
  char * coriolis_file;
  char * output_file;
} settings;