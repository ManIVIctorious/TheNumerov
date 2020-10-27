#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "typedefinitions.h"

// input functions
//  settings input
settings SetDefaultSettings();
void GetSettingsControlFile(char* inputfile,  settings *defaults);
void GetSettingsGetopt(int argc, char** argv, settings *defaults);
//  data input
int InputDataFile(char* inputfile, double** *q, int* nq, double* *v, double** *mu, int dimension, char potential_true, char dipole_true);
int InputCoriolisCoefficients(char* inputfile, double** *q, double** zeta, double*** *mu, int dimension);
double CheckCoordinateSpacing(double** q, int* nq, double threshold, int dimension);

// meta functions
int MetaGetStencil(double* stencil, int n_stencil, int dimension);
int MetaInterpolation(double* *v, int* nq, double dq, int dimension, int n_spline);
int MetaEigensolver(settings prefs, int* nq, double* v, double ekin_param, double* stencil, double* *E, double* *X, double** q, double dq, double*** mu, double** zeta);
double Integrate(int dimension, int* nq, double dx, double* integrand);

// output functions
int Help();
int OutputSettings(FILE* fd, settings *preferences);
int    TextOut_Frequencies(FILE* fd, double ekin_factor, int n_out, double* E);
int TextOut_Orthonormality(FILE* fd, int dimensionality, int n_out, int n_points, int* nq, double* integrand, double dq, double* X);
int      TextOut_Potential(FILE* fd, int dimensionality, int n_out, int n_points, int* nq, double* integrand, double dq, double* X, double* v);
int       TextOut_EKinetic(FILE* fd, settings prefs, int n_out, int n_points, int* nq, double* integrand, double dq, double* X, double* stencil, double ekin_param);
int         TextOut_Dipole(FILE* fd, settings prefs, int n_out, int n_points, int* nq, double* integrand, double dq, double* E, double* X, double** dip);
int   TextOut_Eigenvectors(FILE* fd, settings *prefs, int n_out, int n_points, int* nq, double** q, double* v, double*** mu, double* X, double** dip);


int main(int argc, char* argv[]){

    int i, j, k, l, m, control;

// get preferences out of command line flags, control file and defaults struct:
//  1. Initialize struct "prefs" with the default values
//  2. Parse **argv for the literal string -C or --control-file
//  3. If found pass the following argument to the GetSettingsControlFile() function updating prefs
//  4. Call GetSettingsGetopt() function
// This allows to set all settings in a control file (passed by the -C parameter)
// and afterwards overwrite all settings with the corresponding command line flags
    settings prefs = SetDefaultSettings();
    for(i = 1; i < (argc-1); ++i){
        if(strncmp(argv[i], "-C", 2) == 0 || strncmp(argv[i], "--control-file", 14) == 0){
            GetSettingsControlFile(argv[i+1], &prefs);
        }
    }
    GetSettingsGetopt(argc, argv, &prefs);


//------------------------------------------------------------------------------------------------------------
//    Check for usage of not compiled functionalities      Check for usage of not compiled functionalities
//------------------------------------------------------------------------------------------------------------
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


//------------------------------------------------------------------------------------------------------------
//       Post process settings input       Post process settings input       Post process settings input
//------------------------------------------------------------------------------------------------------------
//Reduced masses:
// The reduced masses are given as a colon separated string array with the number of
//  entries being the dimensionality. To ensure the dimensionality is set beforehand
//  the string is saved temporarily and the conversion to a double array is
//  performed after the initial input

// allocate memory for prefs.masses array and initialize all entries to 1.0 g/mol
    prefs.masses = malloc(prefs.dimension * sizeof(double));
    if(prefs.masses == NULL){ perror("prefs.masses"); exit(errno); }
    for(i = 0; i < prefs.dimension; ++i){
        prefs.masses[i] = 1.0;
    }

// get reduced masses from masses_string
    if(prefs.masses_string_set){

        char * stringp = prefs.masses_string;
        char * token   = NULL;
        char * endptr  = NULL;

        for(i = 0; i < prefs.dimension; ++i){

        // masses string is separated by colon => split at colon position
            token = strsep(&stringp, ":");

            if(token){

                prefs.masses[i] = strtod(token, &endptr);
                if(errno){ perror("strtod()"); exit(errno); }
                if(endptr == token){
                    fprintf(stderr,
                        "\n (-) Error: No digits were found in masses string"
                        "\n     The reduced masses have to be passed as a colon separated array,"
                        "\n     which has the same number of entries as the given dimensionality."
                        "\n     e.g. 1D: -m 1.070908503521"
                        "\n          2D: -m 1.070908503521:1.262900145313, etc."
                        "\n     Aborting...\n\n"
                    );
                    exit(EXIT_FAILURE);
                }

            }
        }
    }


// output settings for verification purpose
    FILE * fd = fopen(prefs.output_file, "w");
    if(fd == NULL){ perror(prefs.output_file); exit(errno); }
    OutputSettings(fd, &prefs);
    fclose(fd); fd = NULL;


//------------------------------------------------------------------------------------------------------------
//  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input
//------------------------------------------------------------------------------------------------------------
    if(!prefs.input_file_set){
        fprintf(stderr, "\n (-) Error: No input file specified.\n     Aborting...\n\n");
        exit(EXIT_FAILURE);
    }

// Input primary file
//------------------------------------------------------------------------------------------------------------
// disable reading dipole file from primary input file when external dipole file is set to be used
    if(prefs.ext_dip_file_set){ prefs.dipole = 0; }

// initialize arrays:
//  2D array double q[dimension][entries] containing the coordinates of all dimensions
    double ** q = malloc(prefs.dimension * sizeof(double*));
    if(q == NULL){ perror("q"); exit(errno); }
    for(i = 0; i < prefs.dimension; ++i){ q[i] = NULL; }

//  1D array int nq[dimension] containing the number of points for each dimension
    int * nq = calloc(prefs.dimension, sizeof(int));
    if(nq == NULL){ perror("nq"); exit(errno); }

//  1D array double v[entries] containing the potential energy values
    double * v = NULL;

//  2D array double dip[3][entries] containing the dipole moment's {x,y,z}-components for all coordinates
    double ** dip = NULL;
    if(prefs.dipole){
        dip = malloc(3 * sizeof(double*));
        if(dip == NULL){ perror("dip"); exit(errno); }
        for(i = 0; i < 3; ++i){ dip[i] = NULL; }
    }

// actual file input, the return value is the total number of valid entries in the input file
    int n_points = InputDataFile(prefs.input_file, &q, nq, &v, &dip, prefs.dimension, 1, prefs.dipole);

// check if the "N nq[0] ... nq[dimension-1]" line in input file matches the number of data points
    for(i = 0, control = 1; i < prefs.dimension; ++i){ control *= nq[i]; }
    if(n_points != control){
        fprintf(stderr, "\n (-) Error reading data from input-file: '%s'", prefs.input_file);
        fprintf(stderr, "\n     Number of Data points (\"%d\") doesn't match \"%d", n_points, nq[0]);
        for(i = 1; i < prefs.dimension; ++i){ fprintf(stderr, "*%d", nq[i]); } fprintf(stderr, "\"");
        fprintf(stderr, "\n     Aborting - please check your input...\n\n");
        exit(EXIT_FAILURE);
    }


// Input of external dipole file
//------------------------------------------------------------------------------------------------------------
    if(prefs.ext_dip_file_set){
    // initialize arrays:
    //  2D array double q_dip[D][entries] analogous to q
        double ** q_dip = malloc(prefs.dimension * sizeof(double*));
        if(q_dip == NULL){ perror("q_dip"); exit(errno); }
        for(i = 0; i < prefs.dimension; ++i){ q_dip[i] = NULL; }

    //  2D array double dip[3][entries] containing the dipole moment's {x,y,z}-components for all coordinates
        dip = malloc(3 * sizeof(double*));
        if(dip == NULL){ perror("dip"); exit(errno); }
        for(i = 0; i < 3; ++i){ dip[i] = NULL; }

    // read file
        control = InputDataFile(prefs.ext_dip_file, &q_dip, NULL, NULL, &dip, prefs.dimension, 0, 1);

    // for every entry in dipole input file there must be exactly one in the primary input file
        if(control != n_points){
            fprintf(stderr,
                "\n (-) Error reading data from input-file: \"%s\""
                "\n     Number of Data points (%d) does not match number in \"%s\" (%d)"
                "\n     Aborting...\n\n"
                , prefs.ext_dip_file, control, prefs.input_file, n_points
            );
            exit(EXIT_FAILURE);
        }

    // check if coordinates are the same
        for(i = 0; i < n_points; ++i){
            for(j = 0; j < prefs.dimension; ++j){
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
        for(i = 0; i < prefs.dimension; ++i){
            free(q_dip[i]); q_dip[i] = NULL;
        }
        free(q_dip); q_dip = NULL;

    // since dipole input was successful set prefs.dipole flag to true
        prefs.dipole = 1;
    }


// Input Coriolis coefficients file
//------------------------------------------------------------------------------------------------------------
    double **  zeta = NULL;   // Coriolis coefficients for all D*(D-1)/2 mode combinations
    double *** mu   = NULL;   // 3*3 "reciprocal moment of inertia tensor" for all coordinate entries

    if(prefs.coriolis_file_set){
    // initialize arrays:
    //  2D array double q_coriolis[dimension][entries] analogue to q
        double **  q_coriolis = malloc(prefs.dimension * sizeof(double*));
        if(q_coriolis == NULL){ perror("q_coriolis"); exit(errno); }
        for(i = 0; i < prefs.dimension; ++i){ q_coriolis[i] = NULL; }

    //  2D array double zeta[3][(D*(D-1))/2] containing the Coriolis coefficients in {x,y,z}-direction
        zeta = malloc(3 * sizeof(double*));
        if(zeta == NULL){ perror("zeta"); exit(errno); }
        for(i = 0; i < 3; ++i){
            zeta[i] = calloc((prefs.dimension*(prefs.dimension - 1))/2, sizeof(double));
            if(zeta[i] == NULL){ perror("zeta[i]"); exit(errno); }
        }

    //  3D array double mu[3][3][entries] containing the 3*3 "inverse moment of inertia tensor" for each configuration
        mu = malloc(3 * sizeof(double**));
        if(mu == NULL){ perror("mu"); exit(errno); }
        for(i = 0; i < 3; ++i){
            mu[i] = malloc(3 * sizeof(double*));
            if(mu[i] == NULL){ perror("mu[i]"); exit(errno); }
            for(j = 0; j < 3; ++j){ mu[i][j] = NULL; }
        }

    // actual file input
        control = InputCoriolisCoefficients(prefs.coriolis_file, &q_coriolis, zeta, &mu, prefs.dimension);

    // for every entry in Coriolis input file there must be exactly one in the input file
        if(control != n_points){
            fprintf(stderr,
                "\n (-) Error reading data from input-file: \"%s\""
                "\n     Number of Data points (%d) does not match number in \"%s\" (%d)"
                "\n     Aborting...\n\n"
                , prefs.coriolis_file, control, prefs.input_file, n_points
            );
            exit(EXIT_FAILURE);
        }

    // check if coordinates are the same
        for(i = 0; i < n_points; ++i){
            for(j = 0; j < prefs.dimension; ++j){
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
        for(i = 0; i < prefs.dimension; ++i){
            free(q_coriolis[i]); q_coriolis[i] = NULL;
        }
        free(q_coriolis); q_coriolis = NULL;
    }


//------------------------------------------------------------------------------------------------------------
//  Check Coordinate Spacing   Check Coordinate Spacing   Check Coordinate Spacing   Check Coordinate Spacing
//------------------------------------------------------------------------------------------------------------
// The Numerov method has to be applied on an equispaced, mass weighted grid.
//  This means that the spacing within each particular coordinate axis (q[0] to q[dimension-1])
//  has to be constant and must be the same in all directions.

// apply mass weighting to coordinates
    for(i = 0; i < prefs.dimension; ++i){
        for(j = 0; j < n_points; ++j){
            q[i][j] *= sqrt(prefs.masses[i]);
        }
    }

// perform spacing check
    double dq = 0;
    if(prefs.check_spacing){
        dq = CheckCoordinateSpacing(q, nq, prefs.threshold, prefs.dimension);
    }else{
    // set default value if no spacing check is performed
        dq = q[prefs.dimension-1][1] - q[prefs.dimension-1][0];
    }


//------------------------------------------------------------------------------------------------------------
//  Interpolation  Interpolation  Interpolation  Interpolation  Interpolation  Interpolation  Interpolation
//------------------------------------------------------------------------------------------------------------
    if(prefs.n_spline){

    // Interpolation of potential
        n_points = MetaInterpolation(&v, nq, dq, prefs.dimension, prefs.n_spline);

    // Calculate and verify expected return value
        for(i = 0, control = 1; i < prefs.dimension; ++i){
            control *= ((nq[i] - 1) * (prefs.n_spline + 1) + 1);
        }
        if(n_points != control){
            fprintf(stderr,
                "\n (-) Error in execution of interpolation function."
                "\n     Aborting...\n\n"
            );
            exit(EXIT_FAILURE);
        }

    // update memory allocation to new n_points for all coordinates
        for(i = 0; i < prefs.dimension; ++i){
            q[i] = realloc(q[i], n_points * sizeof(double));
            if(q[i] == NULL){ perror("reallocation of q[i]"); exit(errno); }
        }


    // if dipole moments are given interpolate them as well
        if(prefs.dipole){
            i = MetaInterpolation(&(dip[0]), nq, dq, prefs.dimension, prefs.n_spline);
            j = MetaInterpolation(&(dip[1]), nq, dq, prefs.dimension, prefs.n_spline);
            k = MetaInterpolation(&(dip[2]), nq, dq, prefs.dimension, prefs.n_spline);

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
        if(prefs.coriolis_file_set){
        // mu is a symmetric 3 times 3 matrix, so an interpolation of the upper or
        //  lower triangle should be identical to a full interpolation
            control = MetaInterpolation(&(mu[0][0]), nq, dq, prefs.dimension, prefs.n_spline);
               i    = MetaInterpolation(&(mu[0][1]), nq, dq, prefs.dimension, prefs.n_spline);
               j    = MetaInterpolation(&(mu[0][2]), nq, dq, prefs.dimension, prefs.n_spline);
               k    = MetaInterpolation(&(mu[1][1]), nq, dq, prefs.dimension, prefs.n_spline);
               l    = MetaInterpolation(&(mu[1][2]), nq, dq, prefs.dimension, prefs.n_spline);
               m    = MetaInterpolation(&(mu[2][2]), nq, dq, prefs.dimension, prefs.n_spline);

        // check return values
            if(control != n_points || i != n_points || j != n_points || k != n_points || l != n_points || m != n_points){
                fprintf(stderr,
                    "\n (-) Error in execution of interpolation function."
                    "\n     Aborting...\n\n"
                );
                exit(EXIT_FAILURE);
            }

        // make mu symmetric again:
        //  free memory of lower triangle and point it to the upper one
            free(mu[1][0]); mu[1][0] = mu[0][1];
            free(mu[2][0]); mu[2][0] = mu[0][2];
            free(mu[2][1]); mu[2][1] = mu[1][2];
        }


    // set new values for number of points for all dimensions and dq
        for(i = 0; i < prefs.dimension; ++i){
            nq[i] = ((nq[i] - 1) * (prefs.n_spline + 1) + 1);
        }
        dq = dq / (double) (prefs.n_spline + 1);


    // set new values for all q entries
        int jump = 1;
        for(i = (prefs.dimension - 1); i >= 0; --i){
            for(j = 0; j < n_points/jump/nq[i]; ++j){
                for(k = 0; k < nq[i]; ++k){
                    for(l = 0; l < jump; ++l){

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
    double * stencil = NULL;

// Allocate memory of n_stencil ** dimension for stencil
    for(i = 0, control = 1; i < prefs.dimension; ++i){ control *= prefs.n_stencil; }
    stencil = calloc(control, sizeof(double));
    if(stencil == NULL){ perror("Stencil"); exit(errno); }

// The number of data points (n_points) must be at least as big as the number of
//  stencil points (n_stencil ** dimension)
    if(n_points < control){
        fprintf(stderr,
            "\n (-) Error reading data from input-file: '%s'"
            "\n     Insufficient number of data points (%d) for stencil size %d (%d points)."
            "\n     Aborting...\n\n"
            , prefs.input_file, n_points, prefs.n_stencil, control
        );
        exit(EXIT_FAILURE);
    }

// Get stencil through meta function
    control = MetaGetStencil(stencil, prefs.n_stencil, prefs.dimension);


//------------------------------------------------------------------------------------------------------------
//      Potential and kinetic energy      Potential and kinetic energy      Potential and kinetic energy
//------------------------------------------------------------------------------------------------------------
// convert potential to output unit of energy
    for(i = 0; i < n_points; ++i){
        v[i] *= prefs.epot_factor;
    }

// shift potential minimum to zero
    double v_min = v[0];
    for(i = 1; i < n_points; ++i){
        if(v[i] < v_min){ v_min = v[i]; }
    }
    for(i = 0; i < n_points; ++i){
        v[i] -= v_min;
    }


// kinetic energy parameter:    [ekin_param] = kJ/mol / (mol/g/angstrom^2)
//  - hbar^2/2 * 10^20          * 1000 * avogadro^2 / 1000 = -10^20 * hbar^2/2 * avogadro^2
//    J kg m^2 * angstrom^2/m^2 * g/kg * (1/mol)^2  / kJ/J =  kJ/mol * g / mol * angstrom^2
    double ekin_param = -1.0E20 * avogadro*avogadro * planck*planck/(8.0*M_PI*M_PI);

// convert kinetic energy parameter to output unit of energy
//  ekin_param          is given in     kJ/mol / (mol/g/angstrom^2)
//  dq                  is given in     angstrom * sqrt(g/mol)
//  prefs.ekin_factor   is given in     (output unit of energy) / (kJ/mol)
    ekin_param *= (prefs.ekin_factor / dq / dq);


// if Coriolis file is set apply potential term of the Watson Hamiltonian ( -1/8 sum_{i=0}^2 mu_ii )
//  mu                  is given in     g/mol/angstrom^2
//  prefs.mu_factor     is given in     kJ/mol / [mu]
//  prefs.ekin_factor   is given in     (output unit of energy) / (kJ/mol)
    if(prefs.coriolis_file_set){
        for(i = 0; i < n_points; ++i){
            v[i] -= ((mu[0][0][i] + mu[1][1][i] + mu[2][2][i]) / 8.0 * (prefs.mu_factor * prefs.ekin_factor));
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
    n_out = MetaEigensolver(prefs, nq, v, ekin_param, stencil, &E, &X, q, dq, mu, zeta);

// The x,y,z - Coriolis coefficients Zeta are not required anymore
    if(prefs.coriolis_file_set){
        free(zeta[0]); zeta[0] = NULL;
        free(zeta[1]); zeta[1] = NULL;
        free(zeta[2]); zeta[2] = NULL;
        free(zeta);    zeta    = NULL;
    }


// Normalize eigen vectors
    double   integral;
    double * integrand = malloc(n_points * sizeof(double));
    if(integrand == NULL){ perror("Integrand"); exit(errno); }

    for(i = 0; i < n_out; ++i){
        for (j = 0; j < n_points; ++j){
            integrand[j] = X[i*n_points + j] * X[i*n_points + j];
        }

        integral = Integrate(prefs.dimension, nq, dq, integrand);

        for (j = 0; j < n_points; ++j){
            X[i*n_points + j] = X[i*n_points + j] / sqrt(integral);
        }
    }


//------------------------------------------------------------------------------------------------------------
//   Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output
//------------------------------------------------------------------------------------------------------------
// open output file
    fd = fopen(prefs.output_file, "a");
    if(fd == NULL){ perror(prefs.output_file); exit(errno); }

// output eigenvalues and frequencies (in cm^-1)
    TextOut_Frequencies(fd, prefs.ekin_factor, n_out, E);

// Analyze section
//------------------------------------------------------------------------------------------------------------
    if(prefs.analyze){

    // output ortho-normality check (i.e. <X[i]|X[j]>)
        TextOut_Orthonormality(fd, prefs.dimension, n_out, n_points, nq, integrand, dq, X);

    // output potential energy (i.e. <X[i]|V|X[j]>)
        TextOut_Potential(fd, prefs.dimension, n_out, n_points, nq, integrand, dq, X, v);

    // output kinetic energy (i.e. <X[i]|ħ² * d²/dx²|X[j]>)
        if(prefs.dimension == 2 && !prefs.coriolis_file_set){
            TextOut_EKinetic(fd, prefs, n_out, n_points, nq, integrand, dq, X, stencil, ekin_param);
        }

    }
// last use of stencil in kinetic energy part
    free(stencil); stencil = NULL;


// dipole section
//------------------------------------------------------------------------------------------------------------
    if(prefs.dipole){
    // Intensities and oscillator strength
        TextOut_Dipole(fd, prefs, n_out, n_points, nq, integrand, dq, E, X, dip);
    }
// last use of integrand in dipole moment part
    free(integrand); integrand = NULL;


// Coordinate dependent
//------------------------------------------------------------------------------------------------------------
//  output coordinates, potential, dipole moments and wave functions
    TextOut_Eigenvectors(fd, &prefs, n_out, n_points, nq, q, v, mu, X, dip);

// close output file
    fclose(fd); fd = NULL;


//------------------------------------------------------------------------------------------------------------
//  free memory    free memory    free memory    free memory    free memory    free memory    free memory
//------------------------------------------------------------------------------------------------------------
// coordinates q and dimensions nq
    for(i = 0; i < prefs.dimension; ++i){
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
    if(prefs.coriolis_file_set){
    // If interpolation is active the lower triangle is already freed.
    //  lower and upper triangle point to the same memory address causing a double call to free()
    //  according to free() specification: "If ptr is NULL, no operation is performed"
    //  => pointing the redundant pointers to NULL is obligatory to prevent undefined behaviour
        if(prefs.n_spline){
            mu[1][0] = NULL;
            mu[2][0] = NULL;
            mu[2][1] = NULL;
        }
        for(i = 0; i < 3; ++i){
            for(j = i; j < 3; ++j){
                free(mu[i][j]); mu[i][j] = NULL;
            }
            free(mu[i]); mu[i] = NULL;
        }
        free(mu); mu = NULL;
    }

    return 0;
}
