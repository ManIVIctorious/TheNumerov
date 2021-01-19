#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "settings.h"

// settings
void usage(void);
settings SetDefaultSettings(void);
void GetSettingsControlFile(char* inputfile,  settings *defaults);
void GetSettingsGetopt(int argc, char** argv, settings *defaults);
void ValidateSettings(settings *set);
void PrintSettings(settings *prefs, FILE *fd);
void free_set_string_values(void);

//  data input
int InputDataFile(char* inputfile, double** *q, int* nq, double* *v, double** *mu, int dimension, char potential_true, char dipole_true);
int InputCoriolisCoefficients(char* inputfile, double** *q, double** zeta, double*** *mu, int dimension);
double CheckCoordinateSpacing(double** q, int* nq, double threshold, int dimension);

// meta functions
void MetaGetStencil(double* stencil, int n_stencil, int dimension);
int  MetaInterpolation(double* *v, int* nq, double dq, int dimension, int n_spline);
int  MetaEigensolver(settings *prefs, int* nq, double* v, double ekin_to_oue, double* stencil, double* *E, double* *X, double** q, double dq, double*** mu, double** zeta);
double Integrate(int dimension, int* nq, double dx, double* integrand);

// output functions
void    PrintFrequencies(FILE* fd, double kJpermol_to_oue, int n_out, double* E);
void PrintOrthonormality(FILE* fd, int dimensionality, int n_out, int n_points, int* nq, double dq, double* X);
void           PrintEPot(FILE* fd, int dimensionality, int n_out, int n_points, int* nq, double dq, double* X, double* v);
void           PrintEKin(FILE* fd, settings *prefs, int n_out, int n_points, int* nq, double dq, double* X, double* stencil, double ekin_to_oue);
void         PrintDipole(FILE* fd, settings *prefs, int n_out, int n_points, int* nq, double dq, double* E, double* X, double** dip);
void   PrintEigenvectors(FILE* fd, settings *prefs, int n_out, int n_points, int* nq, double** q, double* v, double*** mu, double* X, double** dip);


int main(int argc, char* argv[]){

//------------------------------------------------------------------------------------------------------------
//  Settings    Settings    Settings    Settings    Settings    Settings    Settings    Settings    Settings
//------------------------------------------------------------------------------------------------------------

    for(int i = 1; i < argc; ++i){
        if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0){
            usage();
            exit(EXIT_SUCCESS);
        }
    }
// get preferences from default values, control file and command line flags:
/*{{{
 *  1. Initialise settings struct with the default values
 *  2. Parse **argv for the literal string -C or --control-file
 *  3. If found pass the following argument to the GetSettingsControlFile() function updating prefs
 *  4. Call GetSettingsGetopt() function
 * This enables to set all settings in a control file (passed by the -C parameter)
 * and afterwards overwrite all settings with the corresponding command line flags
 *}}}*/
    settings prefs = SetDefaultSettings();
    for(int i = 1; i < (argc-1); ++i){
        if(strcmp(argv[i], "-C") == 0 || strcmp(argv[i], "--control-file") == 0){
            GetSettingsControlFile(argv[i+1], &prefs);
        }
    }
    GetSettingsGetopt(argc, argv, &prefs);

// check if at least one eigensolver is available
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
    exit(EXIT_FAILURE);
#endif
#endif

// validate and post process settings
    ValidateSettings(&prefs);

// print settings for verification
    PrintSettings(&prefs, NULL);


//------------------------------------------------------------------------------------------------------------
//  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input
//------------------------------------------------------------------------------------------------------------

// Input primary file
//------------------------------------------------------------------------------------------------------------
// initialise arrays:
//  2D array double q[dimension][entries] containing the coordinates of all dimensions
    double ** q = malloc(prefs.dimension * sizeof(double*));
    if(q == NULL){ perror("q"); exit(errno); }
    for(int i = 0; i < prefs.dimension; ++i){ q[i] = NULL; }

//  1D array int nq[dimension] containing the number of points for each dimension
    int * nq = calloc(prefs.dimension, sizeof(int));
    if(nq == NULL){ perror("nq"); exit(errno); }

//  1D array double v[entries] containing the potential energy values
    double * v = NULL;

//  2D array double dip[3][entries] containing the dipole moment's {x,y,z}-components for all coordinates
    double ** dip = NULL;
    if( prefs.dipole ){
        dip = malloc(3 * sizeof(double*));
        if(dip == NULL){ perror("dip"); exit(errno); }
        for(int i = 0; i < 3; ++i){ dip[i] = NULL; }
    }

// actual file input, the return value is the total number of valid entries in the input file
    int n_points = InputDataFile(prefs.input_file, &q, nq, &v, &dip, prefs.dimension, 1, prefs.dipole);

// check if the "N nq[0] ... nq[dimension-1]" line in input file matches the number of data points
    int control = 1;
    for(int i = 0; i < prefs.dimension; ++i){ control *= nq[i]; }
    if( control != n_points ){
        fprintf(stderr, "\n (-) Error reading data from input-file: '%s'", prefs.input_file);
        fprintf(stderr, "\n     Number of Data points (\"%d\") doesn't match \"%d", n_points, nq[0]);
        for(int i = 1; i < prefs.dimension; ++i){ fprintf(stderr, "*%d", nq[i]); }
        fprintf(stderr, "\"\n     Aborting - please check your input...\n\n");
        exit(EXIT_FAILURE);
    }


// Input of external dipole file
//------------------------------------------------------------------------------------------------------------
    if( prefs.ext_dip_file ){
    // initialise arrays:
    //  2D array double q_dip[D][entries] analogous to q
        double ** q_dip = malloc(prefs.dimension * sizeof(double*));
        if(q_dip == NULL){ perror("q_dip"); exit(errno); }
        for(int i = 0; i < prefs.dimension; ++i){ q_dip[i] = NULL; }

    //  2D array double dip[3][entries] containing the dipole moment's {x,y,z}-components for each coordinate
        dip = malloc(3 * sizeof(double*));
        if(dip == NULL){ perror("dip"); exit(errno); }
        for(int i = 0; i < 3; ++i){ dip[i] = NULL; }

    // read file
        control = InputDataFile(prefs.ext_dip_file, &q_dip, NULL, NULL, &dip, prefs.dimension, 0, 1);

    // for every entry in dipole input file there must be exactly one in the primary input file
        if( control != n_points ){
            fprintf(stderr,
                "\n (-) Error reading data from input-file: \"%s\""
                "\n     Number of Data points (%d) does not match number in \"%s\" (%d)"
                "\n     Aborting...\n\n"
                , prefs.ext_dip_file, control, prefs.input_file, n_points
            );
            exit(EXIT_FAILURE);
        }

    // check if coordinates are the same
        for(int i = 0; i < n_points; ++i){
            for(int j = 0; j < prefs.dimension; ++j){
                if( (q[j][i] - q_dip[j][i])*(q[j][i] - q_dip[j][i]) > prefs.threshold*prefs.threshold ){
                    fprintf(stderr,
                        "\n (-) Error in external dipole input-file: \"%s\", entry no %d."
                        "\n     The coordinates do not match the ones in \"%s\"."
                        "\n     Aborting...\n\n"
                        , prefs.ext_dip_file, i, prefs.input_file
                    );
                    exit(EXIT_FAILURE);
                }
            }
        }

    // when checks have passed free memory of q_dip
        for(int i = 0; i < prefs.dimension; ++i){ free(q_dip[i]); q_dip[i] = NULL; }
        free(q_dip); q_dip = NULL;

    // since dipole input was successful set prefs.dipole to true
        prefs.dipole = 1;
    }


// Input Coriolis coefficients file
//------------------------------------------------------------------------------------------------------------
    double **  zeta = NULL;   // Coriolis coefficients for all D*(D-1)/2 mode combinations
    double *** mu   = NULL;   // 3x3 "reciprocal moment of inertia tensor" for all coordinate entries

    if( prefs.coriolis_file ){
    // initialise arrays:
    //  2D array double q_coriolis[dimension][entries] analogue to q
        double ** q_coriolis = malloc( prefs.dimension * sizeof(double*) );
        if(q_coriolis == NULL){ perror("q_coriolis"); exit(errno); }
        for(int i = 0; i < prefs.dimension; ++i){ q_coriolis[i] = NULL; }

    //  2D array double zeta[3][(D*(D-1))/2] containing the Coriolis coefficients in {x,y,z}-direction
        zeta = malloc( 3 * sizeof(double*) );
        if(zeta == NULL){ perror("zeta"); exit(errno); }
        for(int i = 0; i < 3; ++i){
            zeta[i] = calloc( (prefs.dimension*(prefs.dimension - 1))/2 , sizeof(double));
            if(zeta[i] == NULL){ perror("zeta[i]"); exit(errno); }
        }

    //  3D array double mu[3][3][entries] containing the 3x3 "inverse moment of inertia tensor" for each configuration
        mu = malloc( 3 * sizeof(double**) );
        if(mu == NULL){ perror("mu"); exit(errno); }
        for(int i = 0; i < 3; ++i){
            mu[i] = malloc(3 * sizeof(double*));
            if(mu[i] == NULL){ perror("mu[i]"); exit(errno); }
            for(int j = 0; j < 3; ++j){ mu[i][j] = NULL; }
        }

    // actual file input
        control = InputCoriolisCoefficients(prefs.coriolis_file, &q_coriolis, zeta, &mu, prefs.dimension);

    // for every entry in Coriolis input file there must be exactly one in the input file
        if( control != n_points ){
            fprintf(stderr,
                "\n (-) Error reading data from input-file: \"%s\""
                "\n     Number of Data points (%d) does not match number in \"%s\" (%d)"
                "\n     Aborting...\n\n"
                , prefs.coriolis_file, control, prefs.input_file, n_points
            );
            exit(EXIT_FAILURE);
        }

    // check if coordinates are the same
        for(int i = 0; i < n_points; ++i){
            for(int j = 0; j < prefs.dimension; ++j){
                if( (q[j][i] - q_coriolis[j][i])*(q[j][i] - q_coriolis[j][i]) > prefs.threshold*prefs.threshold ){
                    fprintf(stderr,
                        "\n (-) Error in Coriolis input-file: \"%s\", entry no %d."
                        "\n     The coordinates do not match the ones in \"%s\"."
                        "\n     Aborting...\n\n"
                        , prefs.coriolis_file, i, prefs.input_file
                    );
                    exit(EXIT_FAILURE);
                }
            }
        }

    // free memory of q_coriolis
        for(int i = 0; i < prefs.dimension; ++i){ free(q_coriolis[i]); q_coriolis[i] = NULL; }
        free(q_coriolis); q_coriolis = NULL;


    // apply the watson threshold to the entries of the reciprocal moment of inertia tensor
        if( prefs.InvInertiaThreshold > 0.0 ){
            fprintf(stderr, "Applying reciprocal moment of inertia threshold\n");

            int count = 0;
            for(int a = 0; a < 3; ++a){
            for(int b = 0; b < 3; ++b){
                for(int i = 0; i < n_points; ++i){
                    if( mu[a][b][i] < -prefs.InvInertiaThreshold ){
                        mu[a][b][i] = -prefs.InvInertiaThreshold;
                        count++;
                    }
                    else if( mu[a][b][i] >  prefs.InvInertiaThreshold ){
                        mu[a][b][i] =  prefs.InvInertiaThreshold;
                        count++;
                    }
                }
            }
            }

            if( count ){
                fprintf(stderr, "%d values overwritten with threshold\n", count);
            }
        }
    }



//------------------------------------------------------------------------------------------------------------
//  Check Coordinate Spacing   Check Coordinate Spacing   Check Coordinate Spacing   Check Coordinate Spacing
//------------------------------------------------------------------------------------------------------------
// The Numerov method has to be applied on an equispaced, mass weighted grid.
//  This means that the spacing within each particular coordinate axis (q[0] to q[dimension-1])
//  has to be constant and must be the same in all directions.

// apply mass weighting to coordinates
    for(int i = 0; i < prefs.dimension; ++i){
        for(int j = 0; j < n_points; ++j){
            q[i][j] *= sqrt(prefs.masses[i]);
        }
    }

// perform spacing check
    double dq = 0;
    if( prefs.check_spacing ){
        dq = CheckCoordinateSpacing(q, nq, prefs.threshold, prefs.dimension);
    }else{
    // set default value if spacing check is omitted
        dq = q[prefs.dimension-1][1] - q[prefs.dimension-1][0];
    }


//------------------------------------------------------------------------------------------------------------
//  Interpolation  Interpolation  Interpolation  Interpolation  Interpolation  Interpolation  Interpolation
//------------------------------------------------------------------------------------------------------------
    if( prefs.n_spline ){

    // Interpolation of potential
        n_points = MetaInterpolation(&v, nq, dq, prefs.dimension, prefs.n_spline);

    // Calculate and verify expected return value
        int control = 1;
        for(int i = 0; i < prefs.dimension; ++i){
            control *= ((nq[i] - 1) * (prefs.n_spline + 1) + 1);
        }
        if( n_points != control ){
            fprintf(stderr,
                "\n (-) Error in execution of interpolation function."
                "\n     Aborting...\n\n"
            );
            exit(EXIT_FAILURE);
        }

    // update memory allocation to new n_points for all coordinates
        for(int i = 0; i < prefs.dimension; ++i){
            q[i] = realloc(q[i], n_points * sizeof(double));
            if(q[i] == NULL){ perror("reallocation of q[i]"); exit(errno); }
        }


    // if dipole moments are given interpolate them as well
        if( prefs.dipole ){
            int i = MetaInterpolation(&(dip[0]), nq, dq, prefs.dimension, prefs.n_spline);
            int j = MetaInterpolation(&(dip[1]), nq, dq, prefs.dimension, prefs.n_spline);
            int k = MetaInterpolation(&(dip[2]), nq, dq, prefs.dimension, prefs.n_spline);

        // check return values
            if(i != n_points || j != n_points || k != n_points){
                fprintf(stderr,
                    "\n (-) Error in execution of interpolation function."
                    "\n     Aborting..."
                    "\n\n"
                );
                exit(EXIT_FAILURE);
            }
        }


    // if Coriolis file is set also interpolate mu
        if( prefs.coriolis_file ){
        // mu is a symmetric 3 times 3 matrix. The lower triangle just points to the values of the
        // upper one, therefore an interpolation has only to be performed on the upper triangle
            int i = MetaInterpolation(&(mu[0][0]), nq, dq, prefs.dimension, prefs.n_spline);
            int j = MetaInterpolation(&(mu[0][1]), nq, dq, prefs.dimension, prefs.n_spline);
            int k = MetaInterpolation(&(mu[0][2]), nq, dq, prefs.dimension, prefs.n_spline);
            int l = MetaInterpolation(&(mu[1][1]), nq, dq, prefs.dimension, prefs.n_spline);
            int m = MetaInterpolation(&(mu[1][2]), nq, dq, prefs.dimension, prefs.n_spline);
            int n = MetaInterpolation(&(mu[2][2]), nq, dq, prefs.dimension, prefs.n_spline);

        // check return values
            if(    (i != n_points)
                || (j != n_points)
                || (k != n_points)
                || (l != n_points)
                || (m != n_points)
                || (n != n_points)
            ){
                fprintf(stderr,
                    "\n (-) Error in execution of interpolation function."
                    "\n     Aborting...\n\n"
                );
                exit(EXIT_FAILURE);
            }

        // ensure mu is symmetric again
        // (the mu returned by the interpolation function is at a different position in memory)
            mu[1][0] = mu[0][1];
            mu[2][0] = mu[0][2];
            mu[2][1] = mu[1][2];
        }


    // set new values for number of points for all dimensions and dq
        for(int i = 0; i < prefs.dimension; ++i){
            nq[i] = ((nq[i] - 1) * (prefs.n_spline + 1) + 1);
        }
        dq = dq / (double) (prefs.n_spline + 1);


    // set new values for all q entries
        int jump = 1;
        for(int i = (prefs.dimension - 1); i >= 0; --i){
            for(int j = 0; j < n_points/jump/nq[i]; ++j){
                for(int k = 0; k < nq[i]; ++k){
                    for(int l = 0; l < jump; ++l){

                        q[i][l + k*jump + j*nq[i]*jump] = q[i][0] + (double)k * dq;

                    }
                }
            }
            jump *= nq[i];
        }
    }


//------------------------------------------------------------------------------------------------------------
//  Stencils    Stencils    Stencils    Stencils    Stencils    Stencils    Stencils    Stencils    Stencils
//------------------------------------------------------------------------------------------------------------
// Allocate memory of n_stencil ** dimension for stencil
// Some of the filling routines only set the non zero values!
// Hence, the memory must be initialised to zero.
    control = 1;
    for(int i = 0; i < prefs.dimension; ++i){ control *= prefs.n_stencil; }
    double * stencil = calloc( control, sizeof(double) );
    if(stencil == NULL){ perror("Stencil"); exit(errno); }

// The number of data points in each dimension must be greater or equal
// to the number of stencil points in that direction
    for(int i = 0; i < prefs.dimension; ++i){

        if( nq[i] < prefs.n_stencil ){
            fprintf(stderr,
                "\n (-) Error: Insufficient number of data points."
                "\n     Dimension %d only contains %d points whereas %d are required."
                "\n     Aborting...\n\n"
                , i, nq[i], prefs.n_stencil
            );
            exit(EXIT_FAILURE);
        }

    }

// Get stencil through meta function
    MetaGetStencil(stencil, prefs.n_stencil, prefs.dimension);


//------------------------------------------------------------------------------------------------------------
//      Potential and kinetic energy      Potential and kinetic energy      Potential and kinetic energy
//------------------------------------------------------------------------------------------------------------
// convert potential to output unit of energy (oue)
    for(int i = 0; i < n_points; ++i){
        v[i] *= prefs.epot_to_oue;
    }

// shift potential minimum to zero
    double v_min = v[0];
    for(int i = 1; i < n_points; ++i){
        if(v[i] < v_min){ v_min = v[i]; }
    }
    for(int i = 0; i < n_points; ++i){
        v[i] -= v_min;
    }


// Conversion factor of kinetic energy to output unit of energy (oue)
/*
 *      E_kin = - hbar^2/2m d^2/dQ^2 Psi = - hbar^2/2 / (sqrt(m).dQ)^2 f.Psi
 *
 *  The mass and coordinates are expected to be in g/mol and Ångstrom, respectively
 *  which are then combined to mass weighted coordinates q = sqrt(m).Q resulting in
 *  [q] = sqrt(g/mol).Å
 *
 *        J.kg.m^2 * Å^2/m^2 * g/kg *  1/mol^2   * kJ/J / Å^2.g/mol
 *      - hbar^2/2 * 10^20   * 1000 * avogadro^2 / 1000 / dq / dq   =
 *
 *            = (-10^20 * avogadro^2 * hbar^2/2 / dq^2) kJ/mol
 */
// begin with ekin_to_kJ/mol and then convert it to output unit of energy
    double ekin_to_oue = -1.0E20 * 0.5*hbar*hbar * avogadro*avogadro / dq / dq;
    ekin_to_oue *= prefs.kJpermol_to_oue;


// if Coriolis file is set apply potential term of the Watson Hamiltonian
/*
 *                 V = V - 1/8 hbar^2 sum_{i=0}^2 mu_ii
 *
 *                 J.s.J.s <=> kJ.g.m^2    mol/(g*Å^2)
 *            - 0.125 * hbar*hbar  * (mu_xx + mu_yy + mu_zz)
 */
    if( prefs.coriolis_file ){
        double prefactor = 0.125 * 1.0E20*hbar*hbar*avogadro*avogadro //   kJ/mol . g.angstrom^2/mol
                         * prefs.InvInertia_to_molpergAasq            // * mol/(g.angstrom^2) / [mu]
                         * prefs.kJpermol_to_oue;                     // * oue / (kJ/mol)
        for(int i = 0; i < n_points; ++i){
            v[i] -= prefactor * (mu[0][0][i] + mu[1][1][i] + mu[2][2][i]);
        }
    }


//------------------------------------------------------------------------------------------------------------
//   eigen-value solver   eigen-value solver   eigen-value solver   eigen-value solver   eigen-value solver
//------------------------------------------------------------------------------------------------------------
    double * E          = NULL;  // eigenvalues
    double * X          = NULL;  // eigenvectors
    long long int n_out =   -1;  // Number of output eigenstates

// Solve the matrix eigenvalue problem
//  MetaEigensolver() forwards all arguments to the Solver<EigensolverName>() function.
//  This function calls an adequate matrix fill routine and then solves the matrix problem
    n_out = MetaEigensolver(&prefs, nq, v, ekin_to_oue, stencil, &E, &X, q, dq, mu, zeta);

// The x,y,z - Coriolis coefficients Zeta are not required anymore
    if( prefs.coriolis_file ){
        free(zeta[0]); zeta[0] = NULL;
        free(zeta[1]); zeta[1] = NULL;
        free(zeta[2]); zeta[2] = NULL;
        free(zeta);    zeta    = NULL;
    }


// Normalise eigenvectors
    double   integral;
    double * integrand = malloc(n_points * sizeof(double));
    if(integrand == NULL){ perror("Integrand"); exit(errno); }

    for(int i = 0; i < n_out; ++i){
        for(int j = 0; j < n_points; ++j){
            integrand[j] = X[i*n_points + j] * X[i*n_points + j];
        }

        integral = Integrate(prefs.dimension, nq, dq, integrand);

        for(int j = 0; j < n_points; ++j){
            X[i*n_points + j] = X[i*n_points + j] / sqrt(integral);
        }
    }
    free(integrand); integrand = NULL;


//------------------------------------------------------------------------------------------------------------
//   Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output
//------------------------------------------------------------------------------------------------------------
// open output file
    FILE * fd = fopen(prefs.output_file, "a");
    if(fd == NULL){ perror(prefs.output_file); exit(errno); }

// output eigenvalues in output unit of energy (oue) and frequencies in cm^-1
    if( prefs.frequencies ){
        PrintFrequencies(fd, prefs.kJpermol_to_oue, n_out, E);
    }

// analyse section
    if( prefs.analyse ){
    // output ortho-normality check (i.e. <X[i]|X[j]>)
        PrintOrthonormality(fd, prefs.dimension, n_out, n_points, nq, dq, X);
    // output potential energy (i.e. <X[i]|V|X[j]>)
        PrintEPot(fd, prefs.dimension, n_out, n_points, nq, dq, X, v);
    // output kinetic energy (i.e. <X[i]|ħ² * d²/dx²|X[j]>)
        if( (prefs.dimension == 2) && !prefs.coriolis_file ){
            PrintEKin(fd, &prefs, n_out, n_points, nq, dq, X, stencil, ekin_to_oue);
        }
    }
// last use of stencil in kinetic energy calculation
    free(stencil); stencil = NULL;

// Intensities and oscillator strengths
    if( prefs.dipole ){
        PrintDipole(fd, &prefs, n_out, n_points, nq, dq, E, X, dip);
    }


// Coordinate dependent
//------------------------------------------------------------------------------------------------------------
//  output coordinates, potential, dipole moments and wave functions
    PrintEigenvectors(fd, &prefs, n_out, n_points, nq, q, v, mu, X, dip);

// close output file
    fclose(fd); fd = NULL;


//------------------------------------------------------------------------------------------------------------
//  free memory    free memory    free memory    free memory    free memory    free memory    free memory
//------------------------------------------------------------------------------------------------------------
    free_set_string_values();
// coordinates q and dimensions nq
    for(int i = 0; i < prefs.dimension; ++i){
        free(q[i]); q[i] = NULL;
    }
    free( q);  q = NULL;
    free(nq); nq = NULL;

// potential v, eigenvectors X and eigenvalues E
    free(v);  v  = NULL;
    free(X);  X  = NULL;
    free(E);  E  = NULL;

// dipole moments in x,y and z-direction
    if(prefs.dipole){
        free(dip[0]); dip[0] = NULL;
        free(dip[1]); dip[1] = NULL;
        free(dip[2]); dip[2] = NULL;
        free(dip);    dip    = NULL;
    }

// mu (inverse moment of inertia 3x3 tensor):
    if( prefs.coriolis_file ){
    // If interpolation is active the lower triangle is already freed.
    //  lower and upper triangle point to the same memory address causing a double call to free()
    //  according to free() specification: "If ptr is NULL, no operation is performed"
    //  => pointing the redundant pointers to NULL is obligatory to prevent undefined behaviour
        for(int i = 0; i < 3; ++i){
            for(int j = i; j < 3; ++j){
                free(mu[i][j]); mu[i][j] = NULL;
            }
            free(mu[i]); mu[i] = NULL;
        }
        free(mu); mu = NULL;
    }

    return 0;
}
