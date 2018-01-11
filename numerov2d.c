#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

// input/output functions
int InputFunction(char *inputfile, double ***q, int *nq, double **V, int dimension);
int InputFunctionDipole(char *inputfile, double ***q, int *nq, double **V, double ***mu, int dimension);
int InputCoriolisCoefficients(char *inputfile, double ***q, double ****zeta, double ****mu, int dimension);
int Help(char *filename);

// meta functions
int MetaGetStencil(double *stencil, int n_stencil, int dimension);

// other
double integrate_1d(int n, double dx, double integrand[]);
double integrate_2d(int nx, int ny, double dx, double integrand[]);

// functions requiring compile time flags
#ifdef HAVE_GSL_INSTALLED
#ifdef HAVE_OPT_SPLINE
    int spline_interpolate(int n_x, int n_y, int n_spline, double x[], double y[], double z[]);
#endif
#endif

#ifdef HAVE_MKL_INSTALLED
    int EigensolverFEAST_MKL_2D(double *v, int *nq, double ekin_param, double *stencil, int n_stencil, double e_min, double e_max, double *E, double *X);
#endif

#ifdef HAVE_ARMA_INSTALLED
    int EigensolverArmadillo_2D(double *v, int *nq, double ekin_param, double *stencil, int n_stencil, int n_out, double *E, double *X);
#endif


int main(int argc, char* argv[]){

// Conversion factors    1.0E20          Ang^2 / m^2
//                       1.0 / 4184.0    kcal / J
// Constants
    double lightspeed = 299792458;        // m/s
    double planck     = 6.626070040E-34;  // Js
    double avogadro   = 6.022140857E23;   // 1/mol

// Default values
    int dipole_flag = 0;
    int n_stencil   = 9;
    int n_spline    = 0;
    int analyse     = 0;
    int n_out       = 8;      // number of eigenstates (ARPACK)

    double ekin_factor = 1.0/4.184;     // (kcal/mol) / (kJ/mol)
    double epot_factor = 1.0;           // (output unit) / (input unit)
    double mass        = 1.0;           // g/mol
    double e_min       = 0.0;           // output energy unit
    double e_max       = 400.0;         // output energy unit
    double spacing_threshold = 1.0E-12; // abs(q[i] - q[i+1])
    double mu_factor = 1.0E20 * avogadro*avogadro * planck*planck/(4.0*M_PI*M_PI); // kJ/mol / (mol/g/angstrom^2)

// Which eigensolver to use
//  1 matches MKL FEAST
//  2 matches ARMADILLO ARPACK
    int Eigensolver = 1;

// file names
    char * input_file_name      = NULL;
    char * coriolis_input_file  = NULL;
    char * output_file_name     = "/dev/stdout";

    int control;
    int i, j, k, l;     // integers for loops
//------------------------------------------------------------------------------------------------------------------
//  FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS
//------------------------------------------------------------------------------------------------------------------
    if(argc == 1){ exit(Help(argv[0])); }
    // optstring contains a list of all short option indices,
    //  indices followed by a colon are options requiring an argument.
    const char         * optstring = "hm:k:v:n:l:u:N:s:ac:i:dPo:t:M:";
    const struct option longopts[] = {
    //  *name:      option name,
    //  has_arg:    if option requires argument,
    //  *flag:      if set to NULL getopt_long() returns val,
    //              else it returns 0 and flag points to a variable set to val
    //  val:        value to return
        {"help",                  no_argument, 0, 'h'},
        {"mass",            required_argument, 0, 'm'},
        {"fkin",            required_argument, 0, 'k'},
        {"fmu",             required_argument, 0, 'M'},
        {"fpot",            required_argument, 0, 'v'},
        {"n-stencil",       required_argument, 0, 'n'},
        {"lower-bound",     required_argument, 0, 'l'},
        {"upper-bound",     required_argument, 0, 'u'},
        {"nout",            required_argument, 0, 'N'},
        {"dq-threshold",    required_argument, 0, 't'},
        {"spline",          required_argument, 0, 's'},
        {"analyze",               no_argument, 0, 'a'},
        {"dipole",                no_argument, 0, 'd'},
        {"pipe",                  no_argument, 0, 'P'},
        {"input-file",      required_argument, 0, 'i'},
        {"coriolis-input",  required_argument, 0, 'c'},
        {"output-file",     required_argument, 0, 'o'},
        {"mkl",                   no_argument, &Eigensolver, 1},
        {"armadillo",             no_argument, &Eigensolver, 2},
        { 0, 0, 0, 0 }
    };

    optind = 1; // option index starting by 1, provided by <getopt.h>
    while(optind < argc){

    // control is the integer representation of the corresponding option, e.g. x = 120
    //  control = -1 corresponds to the end of the options
        control = getopt_long(argc, argv, optstring, longopts, &j);

    // iterate over options control
        switch(control){
            case 'h':
                control = Help(argv[0]);
                exit(control);

            case 'm':
                mass = atof(optarg);
                break;

            case 'k':
                ekin_factor = atof(optarg);
                break;

            case 'M':
                mu_factor = atof(optarg);
                break;

            case 'v':
                epot_factor = atof(optarg);
                break;

            case 'n':
                n_stencil = atoi(optarg);
                break;

            case 'l':
                e_min = atof(optarg);
                break;

            case 'u':
                e_max = atof(optarg);
                break;

            case 'N':
                n_out = atoi(optarg);
                break;

            case 'P':
                input_file_name = "/dev/stdin";
                break;

            case 's':
                n_spline = atoi(optarg);
                break;

            case 't':
                spacing_threshold = atof(optarg);
                break;

            case 'a':
                analyse = 1;
                break;

            case 'd':
                dipole_flag = 1;
                break;

            case 'i':
                input_file_name = optarg;
                break;

            case 'c':
                coriolis_input_file = optarg;
                break;

            case 'o':
                output_file_name = optarg;
                break;

        }
    }

//------------------------------------------------------------------------------------------------------------------
//       Check for usage of not compiled functionalities      Check for usage of not compiled functionalities
//------------------------------------------------------------------------------------------------------------------
#ifndef HAVE_GSL_INSTALLED

#ifndef HAVE_OPT_SPLINE
    if(n_spline != 0){
        fprintf(stderr,
            "\n (-) This version has been compiled without spline support"
            "\n     The setting will be ignored"
            "\n\n"
        );
    }
#endif

#endif

#ifndef HAVE_MKL_INSTALLED
#ifndef HAVE_ARMA_INSTALLED
    fprintf(stderr,
        "\n (-) The availability of a matrix eigensolver is a key requirement"
        "\n     of the Numerov procedure. Please make sure to compile the"
        "\n     program with at least one of the following defines:"
        "\n         -D HAVE_MKL_INSTALLED"
        "\n         -D HAVE_ARMA_INSTALLED"
        "\n\n"
    );
    exit(-1);
#endif
#endif

#ifndef HAVE_MKL_INSTALLED
    if(Eigensolver == 1){
        fprintf(stderr,
            "\n (-) MKL FEAST eigensolver not available."
            "\n     Please make sure to compile the -D HAVE_MKL_INSTALLED define"
            "\n     Trying next eigensolver of the list..."
            "\n\n"
        );
        Eigensolver++;
    }
#endif

#ifndef HAVE_ARMA_INSTALLED
    if(Eigensolver == 2){
        fprintf(stderr,
            "\n (-) Armadillo ARPACK eigensolver not available."
            "\n     Please make sure to compile the -D HAVE_ARMA_INSTALLED define"
            "\n     Trying next eigensolver of the list..."
            "\n\n"
        );
        Eigensolver++;
    }
#endif

//------------------------------------------------------------------------------------------------------------------
//   Declaration Declaration Declaration Declaration Declaration Declaration Declaration Declaration Declaration
//------------------------------------------------------------------------------------------------------------------
// Input
  // standard file
    int dimension = 2;
    int n_points  = 0;          // total number of entries per dimension
    int     * nq  = NULL;       // number of unique entries per dimension
    double ** q   = NULL;       // coordinate entries of all dimensions
    double  * v   = NULL;       // potential entries for each coordinate
    double ** dip = NULL;       // dipole moment for each coordinate
    double  * deltaq = NULL;    // delta q for each individual dimension (for check)
    double dq = 0;              // delta q after coordinate spacing check
    double v_min = 1.0E100;

  // Coriolis coefficients
    double **  q_coriolis = NULL;
    double *** zeta       = NULL;
    double *** mu         = NULL;

// kinetic energy factor:   - hbar^2/2 * 10^20          * 1000 * avogadro^2 / 1000 = -10^20 * hbar^2/2 * avogadro^2
//                            J kg m^2 * angstrom^2/m^2 * g/kg * (1/mol)^2  / kJ/J =  kJ/mol * g * angstrom^2 / mol
    double ekin_param = -1.0E20 * avogadro*avogadro * planck*planck/(8.0*M_PI*M_PI); // kJ/mol / (mol/g/angstrom^2)
    double kJmolToWavenumber = 10.0 / (avogadro*planck*lightspeed);              // cm^-1 / (kJ/mol)

// Eigenvalues and Eigenvectors
    double *E = NULL;   // eigenvalues
    double *X = NULL;   // eigenvectors
    int xsh, ysh;       // integers for applying stencil functions
    int element;

// Output
    int index;
    double freq;
    FILE * file_ptr = NULL;

    double * stencil   = NULL;
    double * integrand = NULL;
    double integral;

// dipole integration
    int n_ts_dip = 0;
    double ts_dip_square;
    double * ts_dip_x = NULL;
    double * ts_dip_y = NULL;
    double * ts_dip_z = NULL;


//------------------------------------------------------------------------------------------------------------------
// Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input
//------------------------------------------------------------------------------------------------------------------
// check input argument if the file is not present give a silly statement
    if(input_file_name == NULL){
        fprintf(stderr, "\n (-) Please specify input file...\n\n");
        exit (1);
    }

// create 2D array q
    q  = malloc(dimension * sizeof(double*));
    for(i = 0; i < dimension; ++i){
        q[i] = malloc(sizeof(double));
    }
    nq = calloc(dimension, sizeof(int));
    v  = malloc(sizeof(double));
    if(q == NULL || v  == NULL){
        fprintf(stderr,
            "\n (-) Error in memory allocation for q or v"
            "\n     Aborting...\n\n"
        );
        exit(1);
    }
    if(dipole_flag == 1){
        dip = malloc(3 * sizeof(double*));
        for(i = 0; i < 3; ++i){
            dip[i] = malloc(sizeof(double));

            if(dip[i] == NULL){
                fprintf(stderr,
                    "\n (-) Error in memory allocation for dipole array"
                    "\n     Aborting...\n\n"
                );
                exit(1);
            }
        }
        n_points = InputFunctionDipole(input_file_name, &q, nq, &v, &dip, dimension);
    }else{
        n_points = InputFunction(input_file_name, &q, nq, &v, dimension);
    }

// check if the "N nq[0] ... nq[dimension-1]" line in input file matches the number of data points
//  The following line is only needed in case of the 1D Numerov, since the N nq[0]...nq[dimension-1] flag
//  is not present in this particular case:
    if(dimension == 1){ nq[0] = n_points; }
    for(i = 0, control = 1; i < dimension; ++i){
        control *= nq[i];
    }
    if(n_points != control){
        fprintf(stderr, "\n (-) Error reading data from input-file: '%s'", input_file_name);
        fprintf(stderr, "\n     Number of Data points (\"%d\") doesn't match \"%d", n_points, nq[0]);
        for(i = 1; i < dimension; ++i){ fprintf(stderr, "*%d", nq[i]); } fprintf(stderr, "\"");
        fprintf(stderr, "\n     Aborting - please check your input...\n\n");
        exit(1);
    }

// there must be at least as many data points as stencil points (n_stencil ** dimension)
    control = 1;
    for(i = 0; i < dimension; ++i){ control *= n_stencil; }
    if(n_points < control){
        fprintf(stderr,
            "\n (-) Error reading data from input-file: '%s'"
            "\n     Insufficient number of data points (%d) "
                   "for stencil size %d (%d points)."
            "\n     Aborting - please check your input..."
            "\n\n"
            , input_file_name, n_points, n_stencil, control
        );
        exit(1);
    }

// input Coriolis coefficients file
    if(coriolis_input_file != NULL){

    // initialize q_coriolis 2D [D][data] double array
        q_coriolis = malloc(dimension * sizeof(double*));
        for(i = 0; i < dimension; ++i){
            q_coriolis[i] = malloc(sizeof(double));
        }

    // initialize zeta 3D [3][(D*D-D)/2][data] double array
        zeta = malloc(3 * sizeof(double**));
        for(i = 0; i < 3; ++i){
            zeta[i] = malloc((dimension*(dimension - 1))/2 * sizeof(double*));
            for(j = 0; j < ((dimension*(dimension - 1))/2); ++j){
                zeta[i][j] = malloc(sizeof(double));
            }
        }

    // initialize mu 3D [3][3][data] double array
        mu = malloc(3 * sizeof(double**));
        for(i = 0; i < 3; ++i){
            mu[i] = malloc(sizeof(double*));
            for(j = 0; j < 3; ++j){
                mu[i][j] = malloc(sizeof(double));
            }
        }

    // actual file input
        control = InputCoriolisCoefficients(coriolis_input_file, &q_coriolis, &zeta, &mu, dimension);

    // for every entry in Coriolis input file there must be exact one in the input file
        if(control != n_points){
            fprintf(stderr,
                "\n (-) Error reading data from input-file: '%s'"
                "\n     Number of Data points (%d) doesn't match number in '%s' (%d)"
                "\n     Aborting - please check your input..."
                "\n\n"
                , coriolis_input_file, control, input_file_name, n_points
            );
            exit(1);
        }

    // check if coordinates are the same
        for(i = 0; i < n_points; ++i){
            for(j = 0; j < dimension; ++j){
                if( (q[j][i] - q_coriolis[j][i])*(q[j][i] - q_coriolis[j][i]) > spacing_threshold*spacing_threshold ){
                    fprintf(stderr,
                        "\n (-) Error in Coriolis input file \"%s\"."
                        "\n     The coordinates do not match the ones in the"
                        "\n     standard input file \"%s\""
                        "\n     Aborting - please check your input..."
                        "\n\n"
                        , coriolis_input_file, input_file_name
                    );
                    exit(-1);

                }
            }
        }

    // free memory of q_coriolis
        for(i = 0; i < dimension; ++i){
            free(q_coriolis[i]);
            q_coriolis[i] = NULL;
        }
        free(q_coriolis); q_coriolis = NULL;

    }


//------------------------------------------------------------------------------------------------------------------
//  Check Coordinate Spacing    Check Coordinate Spacing    Check Coordinate Spacing    Check Coordinate Spacing
//------------------------------------------------------------------------------------------------------------------
    deltaq = calloc(dimension, sizeof(double));

// The Numerov method needs to be applied to an equispaced grid.
//  This means that the spacing within the particular coordinate axes (q[0] to q[dimension-1])
//  has to be constant and must be the same in all directions.

// Get the initial spacing:
//  for the last coordinate q[dimension-1] the next (different) value is the next one
//  for the second last coordinate q[dimension-2] the next (different) value is the nq[dimension-1]th, etc.
    for(i = (dimension - 1), j = 1; i >= 0; --i){
        deltaq[i] = q[i][j] - q[i][0];
        j *= nq[i];
    }

// check spacing within each dimension i
    for(i = (dimension - 1), k = 1; i >= 0; --i){

    // within dimension i check every subtraction of [j+k]th and [j]th values
    //  where j denotes the running index and k is the index jump (e.g. 1 for q[dimension-1],
    //  nq[dimension-1] for q[dimension-2] or nq[dimension-1]*nq[dimension-2] for q[dimension-3])
        for(j = 0; j < n_points-k; ++j){
        // don't check spacing on jump positions
            if( (j+1)%k*nq[i] != 0 ){
                if( (deltaq[i] - (q[i][j+k]-q[i][j]))*(deltaq[i] - (q[i][j+k]-q[i][j]) ) > spacing_threshold*spacing_threshold){
                    fprintf(stderr,
                        "\n (-) Error in input file."
                        "\n     Coordinate spacing not equivalent."
                        "\n     Aborting - please check your input..."
                        "\n\n"
                    );
                    exit(-1);
                }
            }
        }
        k *= nq[i];
    }

// check spacing between dimensions
    for(i = 1; i < dimension; ++i){
        if( (deltaq[i] - deltaq[i-1])*(deltaq[i] - deltaq[i-1]) > spacing_threshold*spacing_threshold){
            fprintf(stderr,
                "\n (-) Error in input file."
                "\n     Coordinate spacing not equivalent."
                "\n     Aborting - please check your input..."
                "\n\n"
            );
            exit(-1);
        }
    }
    dq = deltaq[0];
    free(deltaq); deltaq = NULL;


//------------------------------------------------------------------------------------------------------------------
//   Stencils  Stencils  Stencils  Stencils  Stencils  Stencils  Stencils  Stencils  Stencils  Stencils  Stencils
//------------------------------------------------------------------------------------------------------------------
// Allocate memory of n_stencil ** dimension for stencil
    for(i = 0, j = 1; i < dimension; ++i){ j *= n_stencil; }
    stencil = calloc(j, sizeof(double));
    if(stencil == NULL){
        fprintf(stderr,
            "\n (-) Error in memory allocation for stencil"
            "\n     Aborting..."
            "\n\n"
        );
        exit(1);
    }

// Get stencil through meta function
    control = MetaGetStencil(stencil, n_stencil, dimension);
    if(control != 0 ){
        fprintf(stderr,
            "\n (-) Error initialising stencil parameters."
            "\n     Aborting - please check your input..."
            "\n\n"
        );
        exit(1);
    }

// update kinetic energy pre-factor
    ekin_param *= (ekin_factor / dq / dq / mass);


//------------------------------------------------------------------------------------------------------------------
//  Shift potential    Shift potential    Shift potential    Shift potential    Shift potential    Shift potential
//------------------------------------------------------------------------------------------------------------------
// convert potential to output unit of energy
    for(i = 0; i < n_points; ++i){
        v[i] *= epot_factor;
    }

// if a Coriolis file has been provided as input:
//  add the second term of the Watson molecular Hamiltonian to the potential
    if(coriolis_input_file != NULL){
    // conversion of mu_factor from kJ/mol to desired output unit of energy
        mu_factor *= ekin_factor;
    // add -1/8 * hbar^2 * (mu_xx + mu_yy + mu_zz) to potential
        for(i = 0; i < n_points; ++i){
            v[i] -= (mu[0][0][i] + mu[1][1][i] + mu[2][2][i])*mu_factor/8.0;
        }
    }

// shift potential minimum to zero
    for(i = 0; i < n_points; ++i){
        if(v[i] < v_min){ v_min = v[i]; }
    }
    for(i = 0; i < n_points; ++i){
        v[i] -= v_min;
    }


#ifdef HAVE_OPT_SPLINE
//------------------------------------------------------------------------------------------------------------------
//  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline
//------------------------------------------------------------------------------------------------------------------
    if(n_spline > 0){
    // new number of points by splining
        n_points = ((nq[0] - 1) * (n_spline + 1) + 1) * ((nq[1] - 1) * (n_spline + 1) + 1);

    // reallocate memory and call spline function
        q[0] = realloc(q[0], n_points * sizeof(double));
        q[1] = realloc(q[1], n_points * sizeof(double));
        v  = realloc(v,  n_points * sizeof(double));
        if(q[0] == NULL || q[1] == NULL || v  == NULL){
            fprintf(stderr,
                "\n (-) Error in memory reallocation for q1, q2 or v"
                "\n     Aborting..."
                "\n\n"
            );
            exit(1);
        }
        control = spline_interpolate(nq[0], nq[1], n_spline, q[0], q[1], v);
        if(control != 0){
            fprintf(stderr,
                "\n (-) Error in execution of spline interpolation function."
                "\n     Aborting..."
                "\n\n"
            );
            exit(1);
        }
        if(dipole_flag == 1){
            dip[0] = realloc(dip[0], n_points * sizeof(double));
            dip[1] = realloc(dip[1], n_points * sizeof(double));
            dip[2] = realloc(dip[2], n_points * sizeof(double));
            if(dip[0] == NULL || dip[1] == NULL || dip[2] == NULL){
                fprintf(stderr,
                    "\n (-) Error in memory reallocation for dipole moment"
                    "\n     Aborting..."
                    "\n\n"
                );
                exit(1);
            }
            i = spline_interpolate(nq[0], nq[1], n_spline, q[0], q[1], dip[0]);
            j = spline_interpolate(nq[0], nq[1], n_spline, q[0], q[1], dip[1]);
            k = spline_interpolate(nq[0], nq[1], n_spline, q[0], q[1], dip[2]);
        }

    // set new values for number of points for q1 and q2 and new dq
        nq[0] = (nq[0] - 1) * (n_spline + 1) + 1;
        nq[1] = (nq[1] - 1) * (n_spline + 1) + 1;
        dq    = dq / (double) (n_spline + 1);

    // set new values for q1 and q2:
    //  add 1 step to q1 all "nq[1]"th iteration
    //  add 1 step to q2 every iteration, reset to 0 at every "nq[1]"th
        for(i = 0; i < n_points; i++){
            q[0][i] = q[0][0] + dq * (double) (i/nq[1]);
            q[1][i] = q[1][0] + dq * (double) (i%nq[1]);
        }
    }
#endif


//------------------------------------------------------------------------------------------------------------------
// eigenvalue solver  eigenvalue solver  eigenvalue solver  eigenvalue solver  eigenvalue solver  eigenvalue solver
//------------------------------------------------------------------------------------------------------------------
// allocate memory for eigenvector and eigenvalue arrays
    E  = calloc(n_points,          sizeof(double));   // Eigenvalues
    X  = calloc(n_points*n_points, sizeof(double));   // Eigenvectors
    if(E == NULL || X == NULL){
        fprintf(stderr,
            "\n (-) Error in memory allocation for eigenvalues E and/or eigenvectors X"
            "\n     Aborting...\n\n"
            "\n\n"
        );
        exit(1);
    }

// Eigensolver routines
    if(Eigensolver == 1){
        #ifdef HAVE_MKL_INSTALLED
            n_out = EigensolverFEAST_MKL_2D(v, nq, ekin_param, stencil, n_stencil, e_min, e_max, E, X);
        #endif
        #ifndef HAVE_MKL_INSTALLED
            n_out = -1;
        #endif
    }else if(Eigensolver == 2){
        #ifdef HAVE_ARMA_INSTALLED
            n_out = EigensolverArmadillo_2D(v, nq, ekin_param, stencil, n_stencil, n_out, E, X);
        #endif
        #ifndef HAVE_ARMA_INSTALLED
            n_out = -1;
        #endif
    }else{
        fprintf(stderr,
            "\n (-) Error in eigensolver execution: undefined eigensolver requested"
            "\n     Aborting..."
            "\n\n"
        );
        exit(1);
    }

    if(n_out <= 0){
        fprintf(stderr,
            "\n (-) Error in matrix diagonalisation, no eigenvalue found"
            "\n     Aborting..."
            "\n\n"
        );
        exit(1);
    }


// calculate norm
    integrand = calloc(n_points, sizeof(double));

    for(i = 0; i < n_out; i++){
        for (j = 0; j < n_points; j++){
            integrand[j] = X[j+i*n_points] * X[j+i*n_points];
        }

        integral = integrate_2d(nq[0], nq[1], dq, integrand);

        for (j = 0; j < n_points; j++){
            X[j+i*n_points] = X[j+i*n_points] / sqrt(integral);
        }
    }


//------------------------------------------------------------------------------------------------------------------
//  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output
//------------------------------------------------------------------------------------------------------------------
// open output file
    file_ptr = fopen(output_file_name, "w");
    if(file_ptr == NULL){
        printf("\n\n (-) Error opening output-file: '%s'", output_file_name);
        printf(  "\n     Exiting ... \n\n");
        exit(0);
    }

// output type of eigensolver
    if(Eigensolver == 1){
        fprintf(file_ptr, "# Eigensolver: MKL FEAST\n#\n");
    }else if(Eigensolver == 2){
        fprintf(file_ptr, "# Eigensolver: Armadillo ARPACK\n#\n");
    }

// output eigenvalues
    fprintf(file_ptr, "# Eigenvalues:");
    for(i = 0; i < n_out; i++){
        fprintf(file_ptr, " %24.16lf", E[i]);
    }

    fprintf(file_ptr, "\n# Mass:        %24.16lf", mass);

// and output frequencies
    fprintf(file_ptr, "\n#\n# Frequencies:\n#\n#");
    for(i = 0; i < (n_out - 1); i++){
        fprintf(file_ptr,"%11d   ", i);
    }
    for(i = 1; i < n_out; i++){
        fprintf(file_ptr, "\n#%3d",i);

        for(j = 0; j < i; j++){
            freq = (E[i] - E[j]) * kJmolToWavenumber / ekin_factor;
            fprintf(file_ptr, "  %12.5e", freq);
        }
    }
    fprintf(file_ptr, "\n#\n#");


//------------------------------------------------------------------------------------------------------------------
//  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output
//------------------------------------------------------------------------------------------------------------------
//    Analyze  Analyze  Analyze  Analyze  Analyze  Analyze  Analyze  Analyze  Analyze  Analyze  Analyze  Analyze
//------------------------------------------------------------------------------------------------------------------
// additional output information
    if(analyse == 1){
    // check for ortho-normality of evaluated wave functions:
    //  calculate int Psi_i*Psi_j dxdy (i.e. <X[i]|X[j]>)
        fprintf(file_ptr, "\n# Orthonormality:\n#\n#");

        for(i = 0; i < n_out; i++){
            fprintf(file_ptr, "%11d   ", i);
        }
        for(i = 0; i < n_out; i++){
            fprintf(file_ptr, "\n#%3d",i);

            for(j = 0; j < (i+1); j++){
                for(k = 0; k < n_points; k++){
                // generate integrand
                    integrand[k] = X[k + i*n_points]*X[k + j*n_points];
                }

                integral = integrate_2d(nq[0], nq[1], dq, integrand);
                fprintf(file_ptr, "  %12.5e", integral);
            }
        }
        fprintf(file_ptr, "\n#\n#");


    // Potential output
    //  calculate int Psi_i*V*Psi_j dxdy (i.e. <X[i]|V|X[j]>)
        fprintf(file_ptr, "\n# Potential:\n#\n#");

        for(i = 0; i < n_out; i++){
            fprintf(file_ptr, "%11d   ", i);
        }
        for(i = 0; i < n_out; i++){
            fprintf(file_ptr, "\n#%3d",i);

            for(j = 0; j < (i+1); j++){
                for(k = 0; k < n_points; k++){
                // generate integrand
                    integrand[k] = X[k + i*n_points]*X[k + j*n_points] * v[k];
                }

                integral = integrate_2d(nq[0], nq[1], dq, integrand);
                fprintf(file_ptr, "  %12.5e", integral);
            }
        }
        fprintf(file_ptr, "\n#\n#");


    // kinetic energy output
        fprintf(file_ptr, "\n# Kinetic Energy:\n#\n#");

        for(i = 0; i < n_out; i++){
            fprintf(file_ptr,"%11d   ", i);
        }

        for(i = 0; i < n_out; i++){
            fprintf(file_ptr, "\n#%3d",i);

            for(j = 0; j < (i+1); j++){
                for(k = 0; k < nq[0]; k++){
                    for(l = 0; l < nq[1]; l++){
                        index = k*nq[1] + l;
                        integrand[index]=0;
//------------------------------------------------------------------------------------------------------------------
                        for(xsh = -n_stencil/2; xsh < (n_stencil/2 + 1); xsh++){

                            if( (k+xsh > -1) && (k+xsh < nq[0]) ){
                                for(ysh = -n_stencil/2; ysh < (n_stencil/2 + 1); ysh++){

                                    if( (l+ysh > -1) && (l+ysh < nq[1]) ){
                                        element = (k + xsh)*nq[1] + l + ysh;

                                    // integrand has to be divided by d^2,
                                    //  but the division is already set in the "ekin_param" parameter
                                        integrand[index] = integrand[index] + X[element + i*n_points] * ekin_param * stencil[(xsh + n_stencil/2)*n_stencil + ysh + n_stencil/2]/2;
                                    }
                                }
                            }
                        }
//------------------------------------------------------------------------------------------------------------------
                        integrand[index] *= X[index + j*n_points];
                    }
                }

                integral = integrate_2d(nq[0], nq[1], dq, integrand);
                fprintf(file_ptr, "  %12.5e", integral);
            }
        }
        fprintf(file_ptr,"\n#\n#");


    // calculate coupling
        double *dichtematrix    = calloc(nq[0]*nq[0], sizeof(double));
        double *dichtematrix_sq = calloc(nq[0]*nq[0], sizeof(double));
        if(dichtematrix == NULL || dichtematrix_sq == NULL){
            fprintf(stderr, "\n (-) Error in memory allocation for dichtematrix or its square.");
            fprintf(stderr, "\n     Aborting...\n\n");
            exit(1);
        }
        double *dm_integrand    = calloc(nq[1],     sizeof(double));
        double *dm_integrand_sq = calloc(nq[0],     sizeof(double));
        if(dm_integrand == NULL || dm_integrand_sq == NULL){
            fprintf(stderr, "\n (-) Error in memory allocation for density matrix integrand or its square.");
            fprintf(stderr, "\n     Aborting...\n\n");
            exit(1);
        }
        int r1,r2;

        fprintf(file_ptr, "\n# Coupling:\n#");
        for(i = 0; i < n_out; i++){
        // calculate density-matrix for all wave functions
            for(r1 = 0; r1 < nq[0]; r1++){
                for(r2 = r1; r2 < nq[0]; r2++){
                    for(j = 0; j < nq[1]; j++){
                        dm_integrand[j] = X[i*n_points + r1*nq[1] + j]*X[i*n_points + r2*nq[1] + j];
                    }

                    integral = integrate_1d(nq[1], dq, dm_integrand);
                    dichtematrix[r1*nq[0] + r2] = integral;

                    if(r1 != r2){
                        dichtematrix[r2*nq[0] + r1] = integral;
                    }
                }
            }

        // calculate density-matrix square
        //  careful: density-matrix has dimension nq[0] times nq[0]!
            for(r1 = 0; r1 < nq[0]; r1++){
                for(r2 = r1; r2 < nq[0]; r2++){
                    for(j = 0; j < nq[0]; j++){
                        dm_integrand_sq[j] = dichtematrix[r1*nq[0] + j]*dichtematrix[j*nq[0] + r2];
                    }

                    integral = integrate_1d(nq[0], dq, dm_integrand_sq);
                    dichtematrix_sq[r1*nq[0] + r2] = integral;

                    if(r1 != r2){
                       dichtematrix_sq[r2*nq[0] + r1] = integral;
                    }
                }
            }

        // calculate the trace of density-matrix square
            for(j = 0; j < nq[0]; j++){
                dm_integrand_sq[j] = dichtematrix_sq[j*nq[0] + j];
            }
            integral = integrate_1d(nq[0], dq, dm_integrand_sq);

            fprintf(file_ptr, "\n# state %02d: %2.15lf", i, integral);
        }
        fprintf(file_ptr,"\n#\n#");

        free(dichtematrix);     dichtematrix    = NULL;
        free(dichtematrix_sq);  dichtematrix_sq = NULL;
        free(dm_integrand);     dm_integrand    = NULL;
        free(dm_integrand_sq);  dm_integrand_sq = NULL;
    }// end if(analyse == 1)
    free(stencil);  stencil = NULL;


//------------------------------------------------------------------------------------------------------------------
//  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output
//------------------------------------------------------------------------------------------------------------------
//    Dipole   Dipole   Dipole   Dipole   Dipole   Dipole   Dipole   Dipole   Dipole   Dipole   Dipole   Dipole
//------------------------------------------------------------------------------------------------------------------
    if(dipole_flag == 1){
        n_ts_dip = n_out * (n_out)/2;

        ts_dip_x = malloc(n_ts_dip * sizeof(double));
        ts_dip_y = malloc(n_ts_dip * sizeof(double));
        ts_dip_z = malloc(n_ts_dip * sizeof(double));
        if(ts_dip_x == NULL || ts_dip_y == NULL || ts_dip_z == NULL){
            fprintf(stderr, "\n (-) Error in memory allocation for ts_dip_x, ts_dip_y or ts_dip_z");
            fprintf(stderr, "\n     Aborting...\n\n");
            exit(1);
        }

    // Evaluate x-component of IR-intensity
    //  calculate int Psi_i * \mu_x * Psi_j dxdy (i.e. <X[i]|\mu_x|X[j]>
        fprintf(file_ptr, "\n# Dipole - x-component:\n#\n#");
        for(i = 0; i < n_out; i++){
            fprintf(file_ptr,"%11d   ", i);
        }

        element = 0;
        for(i = 0; i < n_out; i++){
            fprintf(file_ptr, "\n#%3d",i);

            for(j = 0; j < (i+1); j++){
                for(k = 0; k < n_points; k++){
                // generate integrand
                    integrand[k] = X[k + i*n_points]*X[k + j*n_points] * dip[0][k];
                }

                integral = integrate_2d(nq[0], nq[1], dq, integrand);
                fprintf(file_ptr, "  %12.5e", integral);

                if(i != j){
                    ts_dip_x[element] = integral;
                    element ++;
                }
            }
        }
        fprintf(file_ptr, "\n#\n#");

    // Evaluate y-component of IR-intensity
    //  calculate int Psi_i * \mu_y * Psi_j dxdy (i.e. <X[i]|\mu_y|X[j]>
        fprintf(file_ptr, "\n# Dipole - y-component:\n#\n#");
        for(i = 0; i < n_out; i++){
            fprintf(file_ptr,"%11d   ", i);
        }

        element = 0;
        for(i = 0; i < n_out; i++){
            fprintf(file_ptr, "\n#%3d",i);

            for(j = 0; j < (i+1); j++){
                for(k = 0; k < n_points; k++){
                // generate integrand
                    integrand[k] = X[k + i*n_points]*X[k + j*n_points] * dip[1][k];
                }

                integral = integrate_2d(nq[0], nq[1], dq, integrand);
                fprintf(file_ptr, "  %12.5e", integral);

                if(i != j){
                    ts_dip_y[element] = integral;
                    element ++;
                }
            }
        }
        fprintf(file_ptr, "\n#\n#");

    // Evaluate z-component of IR-intensity
    //  calculate int Psi_i * \mu_z * Psi_j dxdy (i.e. <X[i]|\mu_z|X[j]>
        fprintf(file_ptr, "\n# Dipole - z-component:\n#\n#");
        for(i = 0; i < n_out; i++){
            fprintf(file_ptr,"%11d   ", i);
        }

        element = 0;
        for(i = 0; i < n_out; i++){
            fprintf(file_ptr, "\n#%3d",i);

            for(j = 0; j < (i+1); j++){
                for(k = 0; k < n_points; k++){
                // generate integrand
                    integrand[k] = X[k + i*n_points]*X[k + j*n_points] * dip[2][k];
                }

                integral = integrate_2d(nq[0],nq[1], dq, integrand);
                fprintf(file_ptr, "  %12.5e", integral);

                if(i != j){
                    ts_dip_z[element] = integral;
                    element ++;
                }
            }
        }
        fprintf(file_ptr, "\n#\n#");

    // calculate oscillator strength:
    //  (4m\pi) / (3e^2\hbar) * (<Psi_i|\mu_x|Psi_j> + <Psi_i|\mu_y|Psi_j> + <Psi_i|\mu_z|Psi_j>) * (E_j - E_i)
        fprintf(file_ptr, "\n# Oscillator strength:\n#\n#");
        for(i = 0; i < (n_out - 1); i++){
            fprintf(file_ptr,"%11d   ", i);
        }

        element = 0;
        for(i = 1; i < n_out; i++){
            fprintf(file_ptr, "\n#%3d",i);

            for(j = 0; j < i; j++){
                ts_dip_square =   ts_dip_x[element]*ts_dip_x[element]
                                + ts_dip_y[element]*ts_dip_y[element]
                                + ts_dip_z[element]*ts_dip_z[element];
                freq = (E[i] - E[j]) * kJmolToWavenumber / ekin_factor;

                fprintf(file_ptr, "  %12.5e", 4.702E-7 * ts_dip_square * freq);
                element ++;
            }
        }
        fprintf(file_ptr, "\n#\n#");

    // free dipole 2D array:
        for(i = 0; i < 3; ++i){
            free(dip[i]); dip[i] = NULL;
        }
        free(dip);      dip      = NULL;
        free(ts_dip_x); ts_dip_x = NULL;
        free(ts_dip_y); ts_dip_y = NULL;
        free(ts_dip_z); ts_dip_z = NULL;
    }// end if(dipole_flag == 1)


//------------------------------------------------------------------------------------------------------------------
//  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output
//------------------------------------------------------------------------------------------------------------------
//  Eigenvectors  Eigenvectors  Eigenvectors  Eigenvectors  Eigenvectors  Eigenvectors  Eigenvectors  Eigenvectors
//------------------------------------------------------------------------------------------------------------------
// output eigenfunctions
    fprintf(file_ptr, "\n# Potential and Eigenfunctions: %d data-points", n_points);
    fprintf(file_ptr, "\n  N %2d %2d", nq[0], nq[1]);
    fprintf(file_ptr, "\n");

    for(i = 0; i < n_points; i++){
    // newline every "nq[1]"th line
        if(i%nq[1] == 0){
            fprintf(file_ptr, "\n");
        }

    // output coordinates q1 and q2 as well as potential
        fprintf(file_ptr,"%24.16lf    %24.16lf    %24.16lf", q[0][i], q[1][i], v[i]);

    // output wave functions
        for(j = 0; j < n_out; j++){
            fprintf(file_ptr, "    %24.16lf", X[i + j*n_points]);
        }

        fprintf(file_ptr,"\n");
    }

    fclose(file_ptr); file_ptr = NULL;
    for(i = 0; i < dimension; ++i){
        free(q[i]); q[i] = NULL;
    }
    free(q);   q = NULL;
    free(v);   v = NULL;
    free(X);   X = NULL;
    free(E);   E = NULL;

    return 0;
}


double integrate_2d(int nq1, int nq2, double dx, double integrand[]){

    int i, j,index;
    double integral = 0.0;

// inner part, without edges and "RAND"
    for(i = 1; i < (nq1-1); i++){
        for (j = 1; j < (nq2-1); j++){
            index = i*nq2 + j;
            integral = integral + integrand[index];
        }
    }

//randpunkte links rechts ohne ecken
    for(j = 1; j < (nq1-1); j++){
        index = j*nq2;
        integral = integral + 1/2*(integrand[index]+integrand[index + nq2 - 1]);
    }
// randpunkte oben unten
    for(i = 1; i < (nq2-1); i++){
        index = (nq1-1)*nq2 + i;
        integral = integral + 1/2*(integrand[i]+integrand[index]);
    }
    integral = integral + 1/4*(integrand[0]+integrand[nq1-1]+integrand[nq1*nq2-1]+ integrand[(nq1-1)*nq2]);

    return (integral * dx * dx);
}


double integrate_1d(int n, double dx, double integrand[]){

    int i;
    double integral = 0.0;

    integral = 0.5*(integrand[0] + integrand[n-1]);
    for(i = 1; i < (n-1); ++i){
        integral = integral + integrand[i];
    }
    return (integral*dx);
}
