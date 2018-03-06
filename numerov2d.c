#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include "mkl_solvers_ee.h"


// input functions
int InputFunction(char *inputfile, double ***q, int *nq, double **V, int dimension);
int InputFunctionDipole(char *inputfile, double ***q, int *nq, double **V, double ***mu, int dimension);
int InputCoriolisCoefficients(char *inputfile, double ***q, double ****zeta, double ****mu, int dimension);

// output functions
int Help(char *filename);

// other
int get_stencil(double stencil[], int n_stencil);
int spline_interpolate(int n_x, int n_y, int n_spline, double x[], double y[], double z[]);
double integrate_1d(int n, double dx, double integrand[]);
double integrate_2d(int nx, int ny, double dx, double integrand[]);


int main(int argc, char* argv[]){

// Conversion factors    1.0E20          Ang^2 / m^2
//                       1.0 / 4184.0    kcal / J
//                       6.275101074905574e+02   kcal/mol/hartree
// Constants
    double lightspeed = 299792458;        // m/s
    double planck     = 6.626070040E-34;  // Js
    double avogadro   = 6.022140857E23;   // 1/mol

// Default values
    int dipole_flag = 0;
    int n_stencil   = 9;
    int n_spline    = 0;
    int coriolis    = 0;
    int analyse     = 0;


// file names
    char * input_file_name      = NULL;
    char * coriolis_input_file  = NULL;
    char * output_file_name     = "/dev/stdout";

    double ekin_factor = 1.0/4.184;     // (kcal/mol) / (kJ/mol)
    double epot_factor = 1.0;           // (output unit) / (input unit)
    double mass        = 1.0;           // g/mol
    double e_min       = 0.0;           // output energy unit
    double e_max       = 400.0;         // output energy unit
    double spacing_threshold = 1.0E-12; // abs(q[i] - q[i+1])


    int control;
    int i, j, k, l;     // integers for loops
//------------------------------------------------------------------------------------------------------------------
//  FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS
//------------------------------------------------------------------------------------------------------------------
    if(argc == 1){ exit(Help(argv[0])); }
    // optstring contains a list of all short option indices,
    //  indices followed by a colon are options requiring an argument.
    const char         * optstring = "hm:k:v:n:l:u:s:ac:i:dPo:t:";
    const struct option longopts[] = {
    //  *name:      option name,
    //  has_arg:    if option requires argument,
    //  *flag:      if set to NULL getopt_long() returns val,
    //              else it returns 0 and flag points to a variable set to val
    //  val:        value to return
        {"help",                  no_argument, 0, 'h'},
        {"mass",            required_argument, 0, 'm'},
        {"fkin",            required_argument, 0, 'k'},
        {"fpot",            required_argument, 0, 'v'},
        {"n-stencil",       required_argument, 0, 'n'},
        {"lower-bound",     required_argument, 0, 'l'},
        {"upper-bound",     required_argument, 0, 'u'},
        {"dq-threshold",    required_argument, 0, 't'},
        {"spline",          required_argument, 0, 's'},
        {"analyze",               no_argument, 0, 'a'},
        {"dipole",                no_argument, 0, 'd'},
        {"pipe",                  no_argument, 0, 'P'},
        {"input-file",      required_argument, 0, 'i'},
        {"coriolis-input",  required_argument, 0, 'c'},
        {"output-file",     required_argument, 0, 'o'},
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
                coriolis = 1;
                break;

            case 'o':
                output_file_name = optarg;
                break;

            default:
                control = Help(argv[0]);
                exit(control);
        }
    }

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
    double watson_deriv_param = -1.05457180013E-34 * 1.05457180013E-34 / 2.0 / 1.66053904020E-27 * 1.0E20 * 6.02214085774E23 / 4184.0;
    double watson_pot_param = -1.05457180013E-34 * 1.05457180013E-34 / 8.0 / 1.66053904020E-27 * 1.0E20 * 6.02214085774E23 / 4184.0;
//  double ekin_param = -1.05457180013E-34 * 1.05457180013E-34 / 2.0 / 1.66053904020E-27 * 1.0E20 * 6.02214085774E23 / 4184.0;
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

double onestencil[52]={0.0, 0.0, 0.0, 0.0, 0.0, -1.0/2.0, 0.0, 1.0/2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0/60.0, 3.0/20.0, -3.0/4.0, 0.0, 3.0/4.0, -3.0/20.0, 1.0/60.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/280.0, -4.0/105.0, 1.0/5.0, -4.0/5.0, 0.0, 4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0, 0.0, 0.0};
//double onestencil[52]={0.0,0.0,0.0,0.0,0.0,-1/2,0.0,1/2,0.0,0.0,0.0,0.0,0.0};
double second_der[78]={0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, -1.0/12.0, 4.0/3.0, -5.0/2.0, 4.0/3.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 1.0/90.0, -3.0/20.0, 3.0/2.0, -49.0/18.0, 3.0/2.0, -3.0/20.0, 1.0/90.0, 0.0, 0.0, 0.0,0.0, 0.0, -1.0/560.0, 8.0/315.0, -1.0/5.0,  8.0/5.0,-205.0/72.0, 8.0/5.0, -1.0/5.0, 8.0/315.0, -1.0/560.0, 0.0, 0.0,0.0, 1.0/3150.0, -5.0/1008.0, 5.0/126.0, -5.0/21.0, 5.0/3.0, -5296.0/1800.0, 5.0/3.0, -5.0/21.0, 5.0/126.0, -5.0/1008.0, 1.0/3150.0, 0.0,-1.0/16632.0, 2.0/1925.0, -1.0/112.0, 10.0/189.0, -15.0/56.0, 12.0/7.0, -5369.0/1800.0, 12.0/7.0, -15.0/56.0,10.0/189.0, -1.0/112.0, 2.0/1925.0, -1.0/16632.0};
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
        fprintf(stderr, "\n (-) Error in memory allocation for q or v");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }
    if(dipole_flag == 1){
        dip = malloc(3 * sizeof(double*));
        for(i = 0; i < 3; ++i){
            dip[i] = malloc(sizeof(double));

            if(dip[i] == NULL){
                fprintf(stderr, "\n (-) Error in memory allocation for dipole array");
                fprintf(stderr, "\n     Aborting...\n\n");
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

// there must be at least as many data points as the stencil size
    if(n_points < n_stencil){
        fprintf(stderr, "\n (-) Error reading data from input-file: '%s'", input_file_name);
        fprintf(stderr, "\n     Insufficient number of data points %d for stencil size %d.", n_points, n_stencil);
        fprintf(stderr, "\n     Aborting - please check your input...\n\n");
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
            fprintf(stderr, "\n (-) Error reading data from input-file: '%s'", coriolis_input_file);
            fprintf(stderr, "\n     Number of Data points (%d) doesn't match number in '%s' (%d)", control, input_file_name, n_points);
            fprintf(stderr, "\n     Aborting - please check your input...\n\n");
            exit(1);
        }

    // check if coordinates are the same
        for(i = 0; i < n_points; ++i){
            for(j = 0; j < dimension; ++j){
                if( (q[j][i] - q_coriolis[j][i])*(q[j][i] - q_coriolis[j][i]) > spacing_threshold*spacing_threshold ){
                    fprintf(stderr, "\n (-) Error in Coriolis input file \"%s\".", coriolis_input_file);
                    fprintf(stderr, "\n     The coordinates do not match the ones in the");
                    fprintf(stderr, "\n     standard input file \"%s\"", input_file_name);
                    fprintf(stderr, "\n     Aborting - please check your input...\n\n");
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
                    fprintf(stderr, "\n (-) Error in input file.");
                    fprintf(stderr, "\n     Coordinate spacing not equivalent. (1st test)");
                    fprintf(stderr, "\n     Aborting - please check your input...\n\n");
                    exit(-1);
                }
            }
        }
        k *= nq[i];
    }

// check spacing between dimensions
    for(i = 1; i < dimension; ++i){
        if( (deltaq[i] - deltaq[i-1])*(deltaq[i] - deltaq[i-1]) > spacing_threshold*spacing_threshold){
            fprintf(stderr, "\n (-) Error in input file.");
            fprintf(stderr, "\n     Coordinate spacing not equivalent. (3rd test)");
            fprintf(stderr, "\n     Aborting - please check your input...\n\n");
            exit(-1);
        }
    }
    dq = deltaq[0];
    free(deltaq); deltaq = NULL;


//------------------------------------------------------------------------------------------------------------------
//  Shift potential    Shift potential    Shift potential    Shift potential    Shift potential    Shift potential
//------------------------------------------------------------------------------------------------------------------
// get potential minimum, subtract it from the potential
//  and apply potential energy factor to convert to desired energy output
    for(i = 0; i < n_points; ++i){
        if(v[i] < v_min){ v_min = v[i]; }
    }
    for(i = 0; i < n_points; ++i){
        v[i] = (v[i] - v_min) * epot_factor;
    }


//------------------------------------------------------------------------------------------------------------------
//   Stencils  Stencils  Stencils  Stencils  Stencils  Stencils  Stencils  Stencils  Stencils  Stencils  Stencils
//------------------------------------------------------------------------------------------------------------------
// stencil has to have an odd number of entries
    if(n_stencil%2 == 0){
        fprintf(stderr, "\n (-) Stencil size is given as even (%d), but must be an odd number.", n_stencil);
        fprintf(stderr, "\n     Aborting - please check your input...\n\n");
        exit(1);
    }

// get stencil, in two dimensions the size is n_stencil * n_stencil.
    stencil = malloc(n_stencil * n_stencil * sizeof(double));
    if(stencil == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for stencil");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }
    control = get_stencil(stencil, n_stencil);

    if(control != 0 ){
        fprintf(stderr, "\n (-) Error initialising stencil parameters.");
        fprintf(stderr, "\n     Aborting - please check your input...\n\n");
        exit(1);
    }


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
            fprintf(stderr, "\n (-) Error in memory reallocation for q1, q2 or v");
            fprintf(stderr, "\n     Aborting...\n\n");
            exit(1);
        }
        control = spline_interpolate(nq[0], nq[1], n_spline, q[0], q[1], v);
        if(control != 0){
            fprintf(stderr, "\n (-) Error in execution of spline interpolation function.");
            fprintf(stderr, "\n     Aborting...\n\n");
            exit(1);
        }
        if(dipole_flag == 1){
            dip[0] = realloc(dip[0], n_points * sizeof(double));
            dip[1] = realloc(dip[1], n_points * sizeof(double));
            dip[2] = realloc(dip[2], n_points * sizeof(double));
            if(dip[0] == NULL || dip[1] == NULL || dip[2] == NULL){
                fprintf(stderr, "\n (-) Error in memory reallocation for dipole moment");
                fprintf(stderr, "\n     Aborting...\n\n");
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

// apply kinetic energy factor and spacing to ekin_param
    ekin_param = ekin_param / dq / dq / mass;
    ekin_param *= ekin_factor;

int sec_st;
if(n_stencil<10)
{sec_st=(n_stencil-1)/2-1;}
else
{sec_st=3;}

//------------------------------------------------------------------------------------------------------------------
// MKL FEAST eigenvalue solver  MKL FEAST eigenvalue solver  MKL FEAST eigenvalue solver MKL FEAST eigenvalue solver
//------------------------------------------------------------------------------------------------------------------
    char  UPLO = 'F';
    const MKL_INT N = n_points;
    int   n_entries = 0;
    int   xsh, ysh;
    int   element;

    // calculate max_entries
    int max_entries   = 0;
    int sum_q1 = nq[0];
    int sum_q2 = nq[1];

    for(i = 1; i < (n_stencil/2 + 1); i++){
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
//#############################################################################
// TEST   TEST   TEST TEST TES TEST TEST TEST
//##############################################################################
// nothing to test atm.
//#############################################################################
// TEST   TEST   TEST TEST TES TEST TEST TEST
//##############################################################################

    for(i = 0; i < nq[0]; i++)
	{
        for(j = 0; j < nq[1]; j++)
	   {

            for(xsh = -n_stencil/2; xsh < n_stencil/2 + 1; xsh++)
		{
                if( (i+xsh > -1) && (i+xsh < nq[0]) )
		    {

                    for(ysh = -n_stencil/2; ysh < n_stencil/2 + 1; ysh++)
			{
                        if( (j+ysh > -1) && ( j+ysh < nq[1]) )
			{

                            element = (i + xsh)*nq[1] + j+ysh;
                            cols_A[n_entries] = element+1; // wieso +1? weil intel!!
                        

                        // stencil entries have to be divided by 2 to get the right result.
                        //  in three dimensions it should be a division by 4
                          if(coriolis==1)
                          {   vorfaktor=0.0;

                              for (n=0;n<3;n++)
                                 {
                                 for(m=0;m<3;m++) 
                                 {

                                  vorfaktor=vorfaktor -   zeta[n][0][i*nq[1]+j]*zeta[m][0][i*nq[1]+j]*mu[n][m][i*nq[1]+j];
                                 
                                 }// for m
                                }// for n


                              vals_A[n_entries] = ekin_param * stencil[(xsh+n_stencil/2)*n_stencil+ysh+n_stencil/2]/2.0;

     vals_A[n_entries] =vals_A[n_entries]+watson_deriv_param * vorfaktor * (q[0][i*nq[1]+j] *                   onestencil[sec_st*13+6+xsh] * nothing[6+ysh]              / dq
                                                                          + q[1][i*nq[1]+j] *                   onestencil[sec_st*13+6+ysh] * nothing[6+xsh]              / dq
                                                                          - q[0][i*nq[1]+j] * q[0][i*nq[1]+j] * second_der[sec_st*13+6+ysh] * nothing[6+xsh]              / dq / dq
                                                                          - q[1][i*nq[1]+j] * q[1][i*nq[1]+j] * second_der[sec_st*13+6+xsh] * nothing[6+ysh]              / dq / dq
                                                                       +2.0*q[0][i*nq[1]+j] * q[1][i*nq[1]+j] * onestencil[sec_st*13+6+xsh] * onestencil[sec_st*13+6+ysh] / dq / dq);


                                if(xsh == 0 && ysh ==0)                        // add potential to diagonal element
                                {
                                 vals_A[n_entries] = vals_A[n_entries] + v[i*nq[1]+j] +  watson_pot_param * (mu[0][0][i*nq[1]+j]+mu[1][1][i*nq[1]+j]+mu[2][2][i*nq[1]+j]);// *dq*dq*mass*1.0/4.0; // watson pot 

                                 }// end if xsh=ysh=0

                            }// end if coriolis
                            else
                            {
                               vals_A[n_entries] = ekin_param * stencil[(xsh+n_stencil/2)*n_stencil+ysh+n_stencil/2]/2.0;
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



///// START EIGENVALUE CALCULATION

    MKL_INT      fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */

    double       epsout;        /* Relative error on the trace */
    MKL_INT      loop;          /* Number of refinement loop */
    MKL_INT      L = n_points/2;
    MKL_INT      M0;            /* Initial guess for subspace dimension to be used */
    MKL_INT      n_out = 5;


    double       E[n_points];         /* Eigenvalues */

    double * X = calloc (n_points*n_points, sizeof (double));
    if(X == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for Eigenvectors X");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }

    double       res[n_points];       /* Residual */

    MKL_INT      info;          /* Errors */
    M0    = L;
    n_out = L;
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
        &e_min,     // IN: Lower bound of search interval
        &e_max,     // IN: Upper bound of search interval
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
//'#################################################################################################################
//'#################################################################################################################
//'#################################################################################################################
    // get norm  <-------------- sollte nun auch 2d lafn
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

    free(cols_A);   cols_A  = NULL;
    free(vals_A);   vals_A  = NULL;

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
                                        integrand[index] = integrand[index] + X[element + i*n_points] * stencil[(xsh + n_stencil/2)*n_stencil + ysh + n_stencil/2]/2;
                                    }
                                }
                            }
                        }
//------------------------------------------------------------------------------------------------------------------
                        integrand[index]=integrand[index]*ekin_param*X[index + j*n_points];
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
//####################################################################################################################################################<----------here is something to delete
         // summe_mu=0;
         // for (n=0;n<3;n++)
         //  {summe_mu+=mu[n][n][i];}
//        fprintf(file_ptr,"%24.16lf    %24.16lf    %24.16lf   %24.16lf", q[0][i], q[1][i], v[i],v[i]-1/4*ekin_param*summe_mu*dq*dq*mass);
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
