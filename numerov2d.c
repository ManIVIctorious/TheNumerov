#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include "mkl_solvers_ee.h"

#include "typedefinitions.h"

// input/output functions
settings GetSettingsGetopt(settings defaults, int argc, char **argv);
int InputFunction(char *inputfile, double ***q, int *nq, double **V, int dimension);
int InputFunctionDipole(char *inputfile, double ***q, int *nq, double **V, double ***mu, int dimension);
int InputCoriolisCoefficients(char *inputfile, double ***q, double ****zeta, double ****mu, int dimension);
int OutputSettings(FILE *fd, settings preferences);
double CheckCoordinateSpacing(double **q, int *nq, double threshold, int dimension);

// meta functions
int MetaGetStencil(double *stencil, int n_stencil, int dimension);
int MetaInterpolation(double ** v, int * nq, double dq, int dimension, int n_spline);
int MetaEigensolver(settings preferences, double *v, int *nq, double ekin_param, double *stencil, double *E, double *X);

// other
double integrate_1d(int n, double dx, double integrand[]);
double integrate_2d(int nx, int ny, double dx, double integrand[]);

// functions requiring compile time flags
#ifdef HAVE_MKL_INSTALLED
    int EigensolverFEAST_MKL_2D(double *v, int *nq, double ekin_param, double *stencil, int n_stencil, double e_min, double e_max, double *E, double *X);
#endif

#ifdef HAVE_ARMA_INSTALLED
    int EigensolverArmadillo_2D(double *v, int *nq, double ekin_param, double *stencil, int n_stencil, int n_out, double *E, double *X);
#endif


int main(int argc, char* argv[]){

// Conversion factors
//      1.0E20       Angstrom^2 / m^2
//      1.0 / 4.184         cal / J
//      2625.49962     (kJ/mol) / Hartree
//      627.509469   (kcal/mol) / Hartree
// Constants
    double lightspeed = 299792458;        // m/s
    double planck     = 6.62607004E-34;   // Js
    double avogadro   = 6.022140857E23;   // 1/mol

    settings defaults, prefs;

// fill defaults struct with default values
    defaults = (struct settings) {
    // integer values
        .dimension = 2,     // dimension of the problem
        .n_stencil = 9,     // 1D stencil size
        .n_spline  = 0,     // number of interpolation points

    // double precision values
        .mass   = 1.0,      // g/mol

        .ekin_factor = 1.0/4.184,     // (kcal/mol) / (kJ/mol)
        .epot_factor = 1.0,           // (output unit) / (input unit)
        .mu_factor   = 1.0E20 * avogadro*avogadro * planck*planck/(4.0*M_PI*M_PI), // kJ/mol / (mol/g/angstrom^2)
        .threshold = 1.0E-10, // abs(q[i] - q[i+1])

    // Flags
        .analyze = 0,
        .dipole  = 0,
        .check_spacing = 1,

    // Eigensolver specific values
        .Eigensolver = 1,   // 1 = MKL FEAST; 2 = ARMADILLO ARPACK
        .n_out  = 8,        // Number of output eigenstates (ARPACK)
        .e_min  = 0.0,      // minimal energy in output energy unit (FEAST)
        .e_max  = 400.0,    // maximal energy in output energy unit (FEAST)

    // file names
        .input_file    = NULL,
        .coriolis_file = NULL,
        .output_file   = "/dev/stdout",
    };

// get preferences out of command line flags and defaults struct
    prefs = GetSettingsGetopt(defaults, argc, argv);


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
    exit(-1);
#endif
#endif


//------------------------------------------------------------------------------------------------------------
//   Declaration  Declaration  Declaration  Declaration  Declaration  Declaration  Declaration  Declaration
//------------------------------------------------------------------------------------------------------------
// Input
  // standard file
    int control;
    int i, j, k, l;
    long long int n_out = -1;
    int n_points  = 0;          // total number of entries per dimension
    int     * nq  = NULL;       // number of unique entries per dimension
    double ** q   = NULL;       // coordinate entries of all dimensions
    double  * v   = NULL;       // potential entries for each coordinate
    double ** dip = NULL;       // dipole moment for each coordinate
    double dq = 0;              // delta q after coordinate spacing check
    double v_min;


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

// jump parameter for generation of dimension-generalized
//  checks and sets
    int jump;


//------------------------------------------------------------------------------------------------------------
//  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input  Input
//------------------------------------------------------------------------------------------------------------
// check input argument if the file is not present give a silly statement
    if(prefs.input_file == NULL){
        fprintf(stderr, "\n (-) Please specify input file...\n\n");
        exit (1);
    }

// create 2D array q[dimension][entries] containing coordinates
//  array nq, containing number of entries for each dimension
//  and array v, containing all energy values
    q  = malloc(prefs.dimension * sizeof(double*));
    for(i = 0; i < prefs.dimension; ++i){
        q[i] = malloc(sizeof(double));
    }
    nq = calloc(prefs.dimension, sizeof(int));
    v  = malloc(sizeof(double));
    if(q == NULL || v  == NULL){
        fprintf(stderr,
            "\n (-) Error in memory allocation for q or v"
            "\n     Aborting...\n\n"
        );
        exit(1);
    }

// if dipole flag is set use different input routine,
//  expecting 3 additional columns (dip_x, dip_y and dip_z)
    if(prefs.dipole == 1){
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
        n_points = InputFunctionDipole(prefs.input_file, &q, nq, &v, &dip, prefs.dimension);
    }else{
        n_points = InputFunction(prefs.input_file, &q, nq, &v, prefs.dimension);
    }

// check if the "N nq[0] ... nq[dimension-1]" line in input file matches the number of data points
    for(i = 0, control = 1; i < prefs.dimension; ++i){
        control *= nq[i];
    }
    if(n_points != control){
        fprintf(stderr, "\n (-) Error reading data from input-file: '%s'", prefs.input_file);
        fprintf(stderr, "\n     Number of Data points (\"%d\") doesn't match \"%d", n_points, nq[0]);
        for(i = 1; i < prefs.dimension; ++i){ fprintf(stderr, "*%d", nq[i]); } fprintf(stderr, "\"");
        fprintf(stderr, "\n     Aborting - please check your input...\n\n");
        exit(1);
    }

// there must be at least as many data points as stencil points (n_stencil ** dimension)
    for(i = 0, control = 1; i < prefs.dimension; ++i){
        control *= prefs.n_stencil;
    }
    if(n_points < control){
        fprintf(stderr,
            "\n (-) Error reading data from input-file: '%s'"
            "\n     Insufficient number of data points (%d) "
                   "for stencil size %d (%d points)."
            "\n     Aborting - please check your input..."
            "\n\n"
            , prefs.input_file, n_points, prefs.n_stencil, control
        );
        exit(1);
    }

// input Coriolis coefficients file
    if(prefs.coriolis_file != NULL){

    // initialize q_coriolis 2D [D][data] double array
        q_coriolis = malloc(prefs.dimension * sizeof(double*));
        for(i = 0; i < prefs.dimension; ++i){
            q_coriolis[i] = malloc(sizeof(double));
        }

    // initialize zeta 3D [3][(D*D-D)/2][data] double array
        zeta = malloc(3 * sizeof(double**));
        for(i = 0; i < 3; ++i){
            zeta[i] = malloc((prefs.dimension*(prefs.dimension - 1))/2 * sizeof(double*));
            for(j = 0; j < ((prefs.dimension*(prefs.dimension - 1))/2); ++j){
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
        control = InputCoriolisCoefficients(prefs.coriolis_file, &q_coriolis, &zeta, &mu, prefs.dimension);

    // for every entry in Coriolis input file there must be exact one in the input file
        if(control != n_points){
            fprintf(stderr,
                "\n (-) Error reading data from input-file: '%s'"
                "\n     Number of Data points (%d) doesn't match number in '%s' (%d)"
                "\n     Aborting - please check your input..."
                "\n\n"
                , prefs.coriolis_file, control, prefs.input_file, n_points
            );
            exit(1);
        }

    // check if coordinates are the same
        for(i = 0; i < n_points; ++i){
            for(j = 0; j < prefs.dimension; ++j){
                if( (q[j][i] - q_coriolis[j][i])*(q[j][i] - q_coriolis[j][i]) > prefs.threshold*prefs.threshold ){
                    fprintf(stderr,
                        "\n (-) Error in Coriolis input file \"%s\"."
                        "\n     The coordinates do not match the ones in the"
                        "\n     standard input file \"%s\""
                        "\n     Aborting - please check your input..."
                        "\n\n"
                        , prefs.coriolis_file, prefs.input_file
                    );
                    exit(-1);

                }
            }
        }

    // free memory of q_coriolis
        for(i = 0; i < prefs.dimension; ++i){
            free(q_coriolis[i]);
            q_coriolis[i] = NULL;
        }
        free(q_coriolis); q_coriolis = NULL;
    }


//------------------------------------------------------------------------------------------------------------
//  Check Coordinate Spacing   Check Coordinate Spacing   Check Coordinate Spacing   Check Coordinate Spacing
//------------------------------------------------------------------------------------------------------------
// The Numerov method has to be applied to an equispaced grid.
//  This means that the spacing within each particular coordinate axis (q[0] to q[dimension-1])
//  has to be constant and must be the same in all directions.

    if(prefs.check_spacing == 1){
        dq = CheckCoordinateSpacing(q, nq, prefs.threshold, prefs.dimension);
    }

//------------------------------------------------------------------------------------------------------------
//  Stencils    Stencils    Stencils    Stencils    Stencils    Stencils    Stencils    Stencils    Stencils
//------------------------------------------------------------------------------------------------------------
// Allocate memory of n_stencil ** dimension for stencil
    for(i = 0, j = 1; i < prefs.dimension; ++i){ j *= prefs.n_stencil; }
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
    control = MetaGetStencil(stencil, prefs.n_stencil, prefs.dimension);
    if(control != 0 ){
        fprintf(stderr,
            "\n (-) Error initialising stencil parameters."
            "\n     Aborting - please check your input..."
            "\n\n"
        );
        exit(1);
    }

// update kinetic energy pre-factor
    ekin_param *= (prefs.ekin_factor / dq / dq / prefs.mass);


//------------------------------------------------------------------------------------------------------------
//  Shift potential   Shift potential   Shift potential   Shift potential   Shift potential   Shift potential
//------------------------------------------------------------------------------------------------------------
// convert potential to output unit of energy
    for(i = 0; i < n_points; ++i){
        v[i] *= prefs.epot_factor;
    }

// shift potential minimum to zero
    for(i = 1, v_min = v[0]; i < n_points; ++i){
        if(v[i] < v_min){ v_min = v[i]; }
    }
    for(i = 0; i < n_points; ++i){
        v[i] -= v_min;
    }


//------------------------------------------------------------------------------------------------------------
//  Interpolation  Interpolation  Interpolation  Interpolation  Interpolation  Interpolation  Interpolation
//------------------------------------------------------------------------------------------------------------
    if(prefs.n_spline > 0){

    // Interpolation
        n_points = MetaInterpolation(&v, nq, dq, prefs.dimension, prefs.n_spline);

    // Calculate expected return value and check it against n_points
        for(i = 0, control = 1; i < prefs.dimension; ++i){
            control *= ((nq[i] - 1) * (prefs.n_spline + 1) + 1);
        }
        if(n_points != control){
            fprintf(stderr,
                "\n (-) Error in execution of interpolation function."
                "\n     Aborting..."
                "\n\n"
            );
            exit(1);
        }

    // update memory allocation to new n_points for all coordinates
        for(i = 0; i < prefs.dimension; ++i){
            q[i] = realloc(q[i], n_points * sizeof(double));

            if(q[i] == NULL){
                fprintf(stderr,
                    "\n (-) Error in memory reallocation of q%d"
                    "\n     Aborting..."
                    "\n\n"
                    , i
                );
                exit(1);
            }
        }


    // if prefs.dipole is set also interpolate dipole moments
        if(prefs.dipole == 1){
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
                exit(1);
            }
        }


    // set new values for number of points for all dimensions and dq
        for(i = 0; i < prefs.dimension; ++i){
            nq[i] = ((nq[i] - 1) * (prefs.n_spline + 1) + 1);
        }
        dq = dq / (double) (prefs.n_spline + 1);


    // set new values for all q entries
        for(i = prefs.dimension-1, jump = 1; i >= 0; --i){
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
//   eigen-value solver   eigen-value solver   eigen-value solver   eigen-value solver   eigen-value solver
//------------------------------------------------------------------------------------------------------------
//########################################################################################################################

double onestencil[52]={0.0, 0.0, 0.0, 0.0, 0.0, -1.0/2.0, 0.0, 1.0/2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0/60.0, 3.0/20.0, -3.0/4.0, 0.0, 3.0/4.0, -3.0/20.0, 1.0/60.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/280.0, -4.0/105.0, 1.0/5.0, -4.0/5.0, 0.0, 4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0, 0.0, 0.0};
//double onestencil[52]={0.0,0.0,0.0,0.0,0.0,-1/2,0.0,1/2,0.0,0.0,0.0,0.0,0.0};
double second_der[78]={0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, -1.0/12.0, 4.0/3.0, -5.0/2.0, 4.0/3.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 1.0/90.0, -3.0/20.0, 3.0/2.0, -49.0/18.0, 3.0/2.0, -3.0/20.0, 1.0/90.0, 0.0, 0.0, 0.0,0.0, 0.0, -1.0/560.0, 8.0/315.0, -1.0/5.0,  8.0/5.0,-205.0/72.0, 8.0/5.0, -1.0/5.0, 8.0/315.0, -1.0/560.0, 0.0, 0.0,0.0, 1.0/3150.0, -5.0/1008.0, 5.0/126.0, -5.0/21.0, 5.0/3.0, -5296.0/1800.0, 5.0/3.0, -5.0/21.0, 5.0/126.0, -5.0/1008.0, 1.0/3150.0, 0.0,-1.0/16632.0, 2.0/1925.0, -1.0/112.0, 10.0/189.0, -15.0/56.0, 12.0/7.0, -5369.0/1800.0, 12.0/7.0, -15.0/56.0,10.0/189.0, -1.0/112.0, 2.0/1925.0, -1.0/16632.0};

int sec_st;
if(prefs.n_stencil<10)
{sec_st=(prefs.n_stencil-1)/2-1;}
else
{sec_st=3;}

    char  UPLO = 'F';
    const MKL_INT N = n_points;
    int   n_entries = 0;

    // calculate max_entries
    int max_entries   = 0;
    int sum_q1 = nq[0];
    int sum_q2 = nq[1];

    for(i = 1; i < (prefs.n_stencil/2 + 1); i++){
        sum_q1=sum_q1 + 2*(nq[0]-i);
        sum_q2=sum_q2 + 2*(nq[1]-i);
    }
    max_entries = sum_q1*sum_q2; // upper estimation for nnz entries in the matrix, but the easy way to code.

    MKL_INT   rows_A[n_points+1];
    MKL_INT * cols_A = malloc(max_entries * sizeof(MKL_INT));
    if(cols_A == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for cols_A");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }
    double  * vals_A = malloc(max_entries * sizeof(double));
    if(vals_A == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for vals_A");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }


int m,n;
double vorfaktor=0.0;
double nothing[13]={0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0};// für die einzelnen ableitungen.

    for(i = 0; i < nq[0]; i++)
    {
        for(j = 0; j < nq[1]; j++)
       {

            for(xsh = -prefs.n_stencil/2; xsh < prefs.n_stencil/2 + 1; xsh++)
        {
                if( (i+xsh > -1) && (i+xsh < nq[0]) )
            {

                    for(ysh = -prefs.n_stencil/2; ysh < prefs.n_stencil/2 + 1; ysh++)
            {
                        if( (j+ysh > -1) && ( j+ysh < nq[1]) )
            {

                            element = (i + xsh)*nq[1] + j+ysh;
                            cols_A[n_entries] = element+1; // wieso +1? weil intel!!
                        

                        // stencil entries have to be divided by 2 to get the right result.
                        //  in three dimensions it should be a division by 4
                          if(prefs.coriolis_file != NULL)
                          {   vorfaktor=0.0;

                              for (n=0;n<3;n++)
                                 {
                                 for(m=0;m<3;m++) 
                                 {

                                  vorfaktor=vorfaktor -   zeta[n][0][i*nq[1]+j]*zeta[m][0][i*nq[1]+j]*mu[n][m][i*nq[1]+j];
                                 
                                 }// for m
                                }// for n


                              vals_A[n_entries] = ekin_param * stencil[(xsh+prefs.n_stencil/2)*prefs.n_stencil+ysh+prefs.n_stencil/2]/2.0;

     vals_A[n_entries] -= ((prefs.mu_factor * prefs.ekin_factor)/2.0) * vorfaktor * (q[0][i*nq[1]+j] *                   onestencil[sec_st*13+6+xsh] * nothing[6+ysh]              / dq
                                                                                   + q[1][i*nq[1]+j] *                   onestencil[sec_st*13+6+ysh] * nothing[6+xsh]              / dq
                                                                                   - q[0][i*nq[1]+j] * q[0][i*nq[1]+j] * second_der[sec_st*13+6+ysh] * nothing[6+xsh]              / dq / dq
                                                                                   - q[1][i*nq[1]+j] * q[1][i*nq[1]+j] * second_der[sec_st*13+6+xsh] * nothing[6+ysh]              / dq / dq
                                                                                +2.0*q[0][i*nq[1]+j] * q[1][i*nq[1]+j] * onestencil[sec_st*13+6+xsh] * onestencil[sec_st*13+6+ysh] / dq / dq);
         

                                if(xsh == 0 && ysh ==0)                        // add potential to diagonal element
                                {
                                 vals_A[n_entries] += v[i*nq[1]+j] - ((mu[0][0][i*nq[1]+j] + mu[1][1][i*nq[1]+j] + mu[2][2][i*nq[1]+j]) * (prefs.mu_factor * prefs.ekin_factor) / 8.0);// *dq*dq*mass*1.0/4.0; // watson pot 

                                 }// end if xsh=ysh=0

                            }// end if coriolis
                            else
                            {
                               vals_A[n_entries] = ekin_param * stencil[(xsh+prefs.n_stencil/2)*prefs.n_stencil+ysh+prefs.n_stencil/2]/2.0;
                               if(xsh == 0 && ysh ==0)                        // add potential to diagonal element
                                {
                                vals_A[n_entries] = vals_A[n_entries] + v[i*nq[1]+j];
                                }// end if xsh=ysh=0
                            }//end else coriolis

                            n_entries ++;
                        }// end if 0<=ysh<nq[1]
                    }// end for ysh
                } // end if 0<=xsh<nq[0]  
            }// end for xsh
      // after inserting all entries in a row the total number of entries is inserted in the CSR format.
        rows_A[i*nq[1]+j+1]=n_entries+1;
        } // end for j
    }//end for i
    rows_A[0] = 1;

//########################################################################################################################


///// START EIGENVALUE CALCULATION

    MKL_INT      fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */

    double       epsout;        /* Relative error on the trace */
    MKL_INT      loop;          /* Number of refinement loop */
    MKL_INT      L = n_points/2;
    MKL_INT      M0;            /* Initial guess for subspace dimension to be used */


    E = calloc(n_points, sizeof(double));
    X = calloc(n_points*n_points, sizeof(double));
    if(X == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for Eigenvectors X");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }

    double       res[n_points];       /* Residual */

    MKL_INT      info;          /* Errors */
    M0    = L;
    n_out = n_points/2;;
    loop  = 0;
    info  = 0;
    epsout = 0.0;

    feastinit(fpm); /* OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */

    dfeast_scsrev(
        &UPLO,      // IN: UPLO = 'F', stores the full matrix
        &N,         // IN: Size of the problem
        vals_A,     // IN: CSR matrix A, values of non-zero elements
        rows_A,     // IN: CSR matrix A, index of the first non-zero element in row
        cols_A,     // IN: CSR matrix A, columns indices for each non-zero element
        fpm,        // IN/OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers
        &epsout,    // OUT: Relative error of on the trace
        &loop,      // OUT: Contains the number of refinement loop executed
        &prefs.e_min,     // IN: Lower bound of search interval
        &prefs.e_max,     // IN: Upper bound of search interval
        &M0,        // IN: The initial guess for subspace dimension to be used.
        E,          // OUT: The first M entries of Eigenvalues
        X,          // IN/OUT: The first M entries of Eigenvectors
        &n_out,     // OUT: The total number of eigenvalues found in the interval
        res,        // OUT: The first n_out components contain the relative residual vector
        &info       // OUT: Error code
        );

    // Error output
    if ( info != 0 ){
        printf("\n(-) Routine sfeast_scsrev returns code of ERROR: %i\n\n", (int)info);
        exit((int)info);
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
    file_ptr = fopen(prefs.output_file, "w");
    if(file_ptr == NULL){
        printf("\n\n (-) Error opening output-file: '%s'", prefs.output_file);
        printf(  "\n     Exiting ... \n\n");
        exit(0);
    }

// output all available settings
    OutputSettings(file_ptr, prefs);

// output eigenvalues
    fprintf(file_ptr, "# Eigenvalues:");
    for(i = 0; i < n_out; i++){
        fprintf(file_ptr, " %24.16lf", E[i]);
    }

// and output frequencies
    fprintf(file_ptr, "\n#\n# Frequencies:\n#\n#");
    for(i = 0; i < (n_out - 1); i++){
        fprintf(file_ptr,"%11d   ", i);
    }
    for(i = 1; i < n_out; i++){
        fprintf(file_ptr, "\n#%3d",i);

        for(j = 0; j < i; j++){
            freq = (E[i] - E[j]) * kJmolToWavenumber / prefs.ekin_factor;
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
    if(prefs.analyze == 1){
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
                        for(xsh = -prefs.n_stencil/2; xsh < (prefs.n_stencil/2 + 1); xsh++){

                            if( (k+xsh > -1) && (k+xsh < nq[0]) ){
                                for(ysh = -prefs.n_stencil/2; ysh < (prefs.n_stencil/2 + 1); ysh++){

                                    if( (l+ysh > -1) && (l+ysh < nq[1]) ){
                                        element = (k + xsh)*nq[1] + l + ysh;

                                    // integrand has to be divided by d^2,
                                    //  but the division is already set in the "ekin_param" parameter
                                        integrand[index] = integrand[index] + X[element + i*n_points] * ekin_param * stencil[(xsh + prefs.n_stencil/2)*prefs.n_stencil + ysh + prefs.n_stencil/2]/2;
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
    }// end if(prefs.analyze == 1)
    free(stencil);  stencil = NULL;


//------------------------------------------------------------------------------------------------------------------
//  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output
//------------------------------------------------------------------------------------------------------------------
//    Dipole   Dipole   Dipole   Dipole   Dipole   Dipole   Dipole   Dipole   Dipole   Dipole   Dipole   Dipole
//------------------------------------------------------------------------------------------------------------------
    if(prefs.dipole == 1){
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
                freq = (E[i] - E[j]) * kJmolToWavenumber / prefs.ekin_factor;

                fprintf(file_ptr, "  %12.5e", 4.702E-7 * ts_dip_square * freq);
                element ++;
            }
        }
        fprintf(file_ptr, "\n#\n#");

        free(ts_dip_x); ts_dip_x = NULL;
        free(ts_dip_y); ts_dip_y = NULL;
        free(ts_dip_z); ts_dip_z = NULL;
    }// end if(prefs.dipole == 1)


//------------------------------------------------------------------------------------------------------------------
//  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output  Output
//------------------------------------------------------------------------------------------------------------------
//  Eigenvectors  Eigenvectors  Eigenvectors  Eigenvectors  Eigenvectors  Eigenvectors  Eigenvectors  Eigenvectors
//------------------------------------------------------------------------------------------------------------------
    fprintf(file_ptr, "\n# Potential and Eigenfunctions: %d data-points", n_points);

// output size per dimension
    fprintf(file_ptr, "\n  N");
    for(i = 0; i < prefs.dimension; ++i){
        fprintf(file_ptr, " %4d", nq[i]);
    }
    fprintf(file_ptr, "\n");

// output key
    fprintf(file_ptr, "\n#");
    for(i = 0; i < prefs.dimension; ++i){
        fprintf(file_ptr, "\t          q[%d]          ", i);
    }
        fprintf(file_ptr, "\t           v(q)          ");

    if(prefs.dipole != 0){
        fprintf(file_ptr, "\t          dip_x          ");
        fprintf(file_ptr, "\t          dip_y          ");
        fprintf(file_ptr, "\t          dip_z          ");
    }

    if(prefs.coriolis_file != NULL){
        fprintf(file_ptr, "\tv(q) - sum_i(mu[i][i])/8 ");
    }

    for(i = 0; i < n_out; ++i){
        fprintf(file_ptr, "\t        Psi[%d]          ", i);
    }


// output data
    for(i = 0; i < n_points; ++i){
    // add a newline every time a index jumps (from max to min)
        for(j = (prefs.dimension - 1), k = 1; j >= 1; --j){
            k *= nq[j];
            if(i%k == 0){
                fprintf(file_ptr, "\n");
            }
        }

    // output coordinates q[j][i] and potential v[i]
        for(j = 0; j < prefs.dimension; ++j){
            fprintf(file_ptr, "\t% 24.16lf", q[j][i]);
        }
        fprintf(file_ptr, "\t% 24.16lf", v[i]);

    // output dipole moment components
        if(prefs.dipole != 0){
            fprintf(file_ptr, "\t% 24.16lf", dip[0][i]);
            fprintf(file_ptr, "\t% 24.16lf", dip[1][i]);
            fprintf(file_ptr, "\t% 24.16lf", dip[2][i]);
        }


    // output potential after addition of Watson potential term
        if(prefs.coriolis_file != NULL){
            fprintf(file_ptr, "\t% 24.16lf", v[i] - ((mu[0][0][i] + mu[1][1][i] + mu[2][2][i])/8.0 * (prefs.mu_factor * prefs.ekin_factor)));
        }

    // output wave functions
        for(j = 0; j < n_out; ++j){
            fprintf(file_ptr, "\t% 24.16lf", X[i + j*n_points]);
        }

        fprintf(file_ptr, "\n");
    }


// free memory
    fclose(file_ptr); file_ptr = NULL;
    for(i = 0; i < prefs.dimension; ++i){
        free(q[i]); q[i] = NULL;
    }
    free(q); q = NULL;
    for(i = 0; i < 3; ++i){
        free(dip[i]); dip[i] = NULL;
    }
    free(dip); dip = NULL;
    free(q);    q  = NULL;
    free(v);    v  = NULL;
    free(X);    X  = NULL;
    free(E);    E  = NULL;

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
