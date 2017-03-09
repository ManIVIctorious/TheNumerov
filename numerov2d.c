#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include "mkl_solvers_ee.h"

int Help(char *filename);
int get_stencil(double stencil[], int n_stencil);
double integrate_1d(int n, double dx, double integrand[]);
double integrate_2d(int nx, int ny, double dx, double integrand[]);
int InputFunction(char *inputfile, double **q1, double **q2, double **V, int *nq1, int *nq2);
int InputFunctionDipole(char *inputfile, double **q1, double **q2, double **V, double **mux, double **muy, double **muz, int *nq1, int *nq2);
int spline_interpolate(int n_x, int n_y, int n_spline, double x[], double y[], double z[]);


int main(int argc, char* argv[]){

// Conversion factors    1.0E20          Ang^2 / m^2
//                       1.0 / 4184.0    kcal / J
// Constants
    double pi         = 3.1415926535897932384626433832795;
    double lightspeed = 299792458;        // m/s
    double planck     = 6.626070040E-34;  // Js
    double avogadro   = 6.022140857E23;   // 1/mol

// Default values
    int dipole_flag = 0;
    int n_stencil   = 9;
    int n_spline    = 0;
    int analyse     = 0;

    char *input_file_name  = NULL;
    char *output_file_name = "/dev/stdout";

    double ekin_factor = 1.0/4.184;     // (kcal/mol) / (kJ/mol)
    double epot_factor = 1.0;           // (output unit) / (input unit)
    double mass        = 1.0;           // g/mol
    double e_min       = 0.0;           // output energy unit
    double e_max       = 100.0;         // output energy unit
    double spacing_threshold = 1.0E-12; // abs(q[i] - q[i+1])


    int control;
    if(argc == 1){
       control = Help(argv[0]);
       exit(control);
    }
//------------------------------------------------------------------------------------------------------------------
//  FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS   FLAGS
//------------------------------------------------------------------------------------------------------------------
    while(1){
        static struct option long_options[] = {
            {"help",               no_argument, 0, 'h'},
            {"mass",         required_argument, 0, 'm'},
            {"fkin",         required_argument, 0, 'k'},
            {"fpot",         required_argument, 0, 'v'},
            {"n-stencil",    required_argument, 0, 'n'},
            {"lower-bound",  required_argument, 0, 'l'},
            {"upper-bound",  required_argument, 0, 'u'},
            {"dq-threshold", required_argument, 0, 't'},
            {"spline",       required_argument, 0, 's'},
            {"analyze",            no_argument, 0, 'a'},
            {"dipole",             no_argument, 0, 'd'},
            {"pipe",               no_argument, 0, 'P'},
            {"input-file",   required_argument, 0, 'i'},
            {"output-file",  required_argument, 0, 'o'},
            { 0, 0, 0, 0 }
        };
    // getopt_long stores the option index here.
        int option_index  = 0;

        control = getopt_long(argc, argv, "hm:k:v:n:l:u:s:adi:Po:", long_options, &option_index);

    // Detect the end of the options.
        if(control == -1)
            break;

        switch(control){
            case 0:
                /* If this option set a flag, do nothing else now. */
                if(long_options[option_index].flag != 0)    break;

                printf("option %s", long_options[option_index].name);
                if(optarg)  printf(" with arg %s", optarg);

                printf("\n");
                break;

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

    int i, j, k, l, index;
    int nq1, nq2, n_points;

// Input
    double * q1 = NULL;
    double * q2 = NULL;
    double * v  = NULL;
    double * dip_x = NULL;
    double * dip_y = NULL;
    double * dip_z = NULL;
    double v_min = 1.0E100;
    double dq, dq1, dq2;

    double * stencil   = NULL;
    double * integrand = NULL;
    double integral;

// kinetic energy factor:   - hbar^2/2 * 10^20          * 1000 * avogadro^2 / 1000 = -10^20 * hbar^2/2 * avogadro^2
//                            J kg m^2 * angstrom^2/m^2 * g/kg * (1/mol)^2  / kJ/J =  kJ/mol * g * angstrom^2 / mol
    double ekin_param = -1.0E20 * avogadro*avogadro * planck*planck/(8.0*pi*pi); // kJ/mol / (mol/g/angstrom^2)
    double kJmolToWavenumber = 10.0 / (avogadro*planck*lightspeed);              // cm^-1 / (kJ/mol)

// dipole integration
    int n_ts_dip = 0;
    double ts_dip_square;
    double * ts_dip_x = NULL;
    double * ts_dip_y = NULL;
    double * ts_dip_z = NULL;

// Output
    double freq;
    FILE * file_ptr = NULL;


//------------------------------------------------------------------------------------------------------------------
// Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input
//------------------------------------------------------------------------------------------------------------------
// check input argument if the file is not present give a silly statement
    if(input_file_name == NULL){
        fprintf(stderr, "\n (-) Please specify input file...\n\n");
        exit (1);
    }

    q1 = malloc(sizeof(double));
    q2 = malloc(sizeof(double));
    v  = malloc(sizeof(double));
    if(q1 == NULL || q2 == NULL || v  == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for q1, q2 or v");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }
    if(dipole_flag == 1){
        dip_x = malloc(sizeof(double));
        dip_y = malloc(sizeof(double));
        dip_z = malloc(sizeof(double));
        if(dip_x == NULL || dip_y == NULL || dip_z == NULL){
            fprintf(stderr, "\n (-) Error in memory allocation for dip_x, dip_y or dip_z");
            fprintf(stderr, "\n     Aborting...\n\n");
            exit(1);
        }
        n_points = InputFunctionDipole(input_file_name, &q1, &q2, &v, &dip_x, &dip_y, &dip_z, &nq1, &nq2);
    }
    else{
        n_points = InputFunction(input_file_name, &q1, &q2, &v, &nq1, &nq2);
    }

// check if the "N nq1 nq2" line in input file matches the number of data points
    if(n_points != nq1*nq2){
        fprintf(stderr, "\n (-) Error reading data from input-file: '%s'", input_file_name);
        fprintf(stderr, "\n     Number of Data points (\"%d\") doesn't match \"%d*%d\"", n_points, nq1, nq2);
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
        fprintf(stderr, "\n (-) Error in memory allocation for dip_x");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }
    control = get_stencil(stencil, n_stencil);

    if(control != 0 ){
        fprintf(stderr, "\n (-) Error initialising stencil parameters.");
        fprintf(stderr, "\n     Aborting - please check your input...\n\n");
        exit(1);
    }


// get potential minimum, subtract it from the potential
//  and apply potential energy factor to convert to desired energy output
    for(i = 0; i < n_points; ++i){
        if(v[i] < v_min)   v_min = v[i];
    }
    for(i = 0; i < n_points; ++i){
        v[i] = (v[i] - v_min) * epot_factor;
    }


//------------------------------------------------------------------------------------------------------------------
//  Check Coordinate Spacing    Check Coordinate Spacing    Check Coordinate Spacing    Check Coordinate Spacing
//------------------------------------------------------------------------------------------------------------------
    dq1 = q1[nq2] - q1[0];
    dq2 = q2[ 1 ] - q2[0];
    for(i = 1; i < (n_points-nq2); ++i){

    // spacing between q1[i] and q1[i-1]
        if( (dq1 - (q1[i+nq2] - q1[i]))*(dq1 - (q1[i+nq2] - q1[i])) > spacing_threshold*spacing_threshold ){
            fprintf(stderr, "\n (-) Error in input file.");
            fprintf(stderr, "\n     Coordinate spacing not equivalent.");
            fprintf(stderr, "\n     Aborting - please check your input...\n\n");
            exit(-1);
        }
        dq1 = q1[i+nq2] - q1[i];

    // spacing between q2[i] and q2[i-1]
        if( (i+1)%nq2 != 0 ){
            if( (dq2 - (q2[i+1] - q2[i]))*(dq2 - (q2[i+1] - q2[i])) > spacing_threshold*spacing_threshold ){
                fprintf(stderr, "\n (-) Error in input file.");
                fprintf(stderr, "\n     Coordinate spacing not equivalent.");
                fprintf(stderr, "\n     Aborting - please check your input...\n\n");
                exit(-1);
            }
            dq2 = q2[i+1] - q2[i];
        }

    // spacing between q1[i] and q2[i]
        if( (dq1 - dq2)*(dq1 - dq2) > spacing_threshold*spacing_threshold ){
            fprintf(stderr, "\n (-) Error in input file.");
            fprintf(stderr, "\n     Coordinate spacing not equivalent.");
            fprintf(stderr, "\n     Aborting - please check your input...\n\n");
            exit(-1);
        }
    }
    dq = dq1;


//------------------------------------------------------------------------------------------------------------------
//  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline
//------------------------------------------------------------------------------------------------------------------
    if(n_spline > 0){
    // new number of points by splining
        n_points = ((nq1 - 1) * (n_spline + 1) + 1) * ((nq2 - 1) * (n_spline + 1) + 1);

    // reallocate memory and call spline function
        q1 = realloc(q1, n_points * sizeof(double));
        q2 = realloc(q2, n_points * sizeof(double));
        v  = realloc(v,  n_points * sizeof(double));
        if(q1 == NULL || q2 == NULL || v  == NULL){
            fprintf(stderr, "\n (-) Error in memory reallocation for q1, q2 or v");
            fprintf(stderr, "\n     Aborting...\n\n");
            exit(1);
        }
        control = spline_interpolate(nq1, nq2, n_spline, q1, q2, v);
        if(control != 0){
            fprintf(stderr, "\n (-) Error in execution of spline interpolation function.");
            fprintf(stderr, "\n     Aborting...\n\n");
            exit(1);
        }
        if(dipole_flag == 1){
            dip_x = realloc(dip_x, n_points * sizeof(double));
            dip_y = realloc(dip_y, n_points * sizeof(double));
            dip_z = realloc(dip_z, n_points * sizeof(double));
            if(dip_x == NULL || dip_y == NULL || dip_z == NULL){
                fprintf(stderr, "\n (-) Error in memory reallocation for dip_x, dip_y or dip_z");
                fprintf(stderr, "\n     Aborting...\n\n");
                exit(1);
            }
            i = spline_interpolate(nq1, nq2, n_spline, q1, q2, dip_x);
            j = spline_interpolate(nq1, nq2, n_spline, q1, q2, dip_y);
            k = spline_interpolate(nq1, nq2, n_spline, q1, q2, dip_z);
        }

    // set new values for number of points for q1 and q2 and new dq
        nq1 = (nq1 - 1) * (n_spline + 1) + 1;
        nq2 = (nq2 - 1) * (n_spline + 1) + 1;
        dq  = dq / (double) (n_spline + 1);

    // set new values for q1 and q2:
    //  add 1 step to q1 all "nq2"th iteration
    //  add 1 step to q2 every iteration, reset to 0 at every "nq2"th
        for(i = 0; i < n_points; i++){
            q1[i] = q1[0] + dq * (double) (i/nq2);
            q2[i] = q2[0] + dq * (double) (i%nq2);
        }
    }

// apply kinetic energy factor and spacing to ekin_param
    ekin_param = ekin_param / dq / dq / mass;
    ekin_param *= ekin_factor;


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
    int sum_q1 = nq1;
    int sum_q2 = nq2;

    for(i = 1; i < (n_stencil/2 + 1); i++){
        sum_q1=sum_q1 + 2*(nq1-i);
        sum_q2=sum_q2 + 2*(nq2-i);
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

    for(i = 0; i < nq1; i++){
        for(j = 0; j < nq2; j++){
            for(xsh = -n_stencil/2; xsh < n_stencil/2 + 1; xsh++){

                if( (i+xsh > -1) && (i+xsh < nq1) ){
                    for(ysh = -n_stencil/2; ysh < n_stencil/2 + 1; ysh++){

                        if( (j+ysh > -1) && ( j+ysh < nq2) ){
                            element = (i + xsh)*nq2 + j+ysh;
                            cols_A[n_entries] = element+1; // wieso +1? weil intel!!

                        // stencil entries have to be divided by 2 to get the right result.
                        //  in three dimensions it should be a division by 4
                            vals_A[n_entries] = ekin_param * stencil[(xsh+n_stencil/2)*n_stencil+ysh+n_stencil/2]/2;

                        // add potential to diagonal element
                            if(xsh == 0 && ysh ==0){
                                vals_A[n_entries] = vals_A[n_entries] + v[i*nq2+j];
                            }

                            n_entries ++;
                        }
                    }
                }
            }
      // after inserting all entries in a row the total number of entries is inserted in the CSR format.
        rows_A[i*nq2+j+1]=n_entries+1;
        }
    }
    rows_A[0] = 1;


///// START EIGENVALUE CALCULATION

    MKL_INT      fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */

    double       epsout;        /* Relative error on the trace */
    MKL_INT      loop;          /* Number of refinement loop */
    MKL_INT      L = n_points/2;
    MKL_INT      M0;            /* Initial guess for subspace dimension to be used */
    MKL_INT      n_out = 5;


    double       E[n_points];         /* Eigenvalues */

    double * X = calloc (n_points*(n_points-1)/2, sizeof (double));
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

        integral = integrate_2d(nq1,nq2, dq, integrand);

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

                integral = integrate_2d(nq1, nq2, dq, integrand);
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

                integral = integrate_2d(nq1, nq2, dq, integrand);
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
                for(k = 0; k < nq1; k++){
                    for(l = 0; l < nq2; l++){
                        index = k*nq2 + l;
                        integrand[index]=0;
//------------------------------------------------------------------------------------------------------------------
                        for(xsh = -n_stencil/2; xsh < (n_stencil/2 + 1); xsh++){

                            if( (k+xsh > -1) && (k+xsh < nq1) ){
                                for(ysh = -n_stencil/2; ysh < (n_stencil/2 + 1); ysh++){

                                    if( (l+ysh > -1) && (l+ysh < nq2) ){
                                        element = (k + xsh)*nq2 + l + ysh;

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

                integral = integrate_2d(nq1, nq2, dq, integrand);
                fprintf(file_ptr, "  %12.5e", integral);
            }
        }
        fprintf(file_ptr,"\n#\n#");


    // calculate coupling
        double *dichtematrix    = calloc(nq1*nq1, sizeof(double));
        double *dichtematrix_sq = calloc(nq1*nq1, sizeof(double));
        if(dichtematrix == NULL || dichtematrix_sq == NULL){
            fprintf(stderr, "\n (-) Error in memory allocation for dichtematrix or its square.");
            fprintf(stderr, "\n     Aborting...\n\n");
            exit(1);
        }
        double *dm_integrand    = calloc(nq2,     sizeof(double));
        double *dm_integrand_sq = calloc(nq1,     sizeof(double));
        if(dm_integrand == NULL || dm_integrand_sq == NULL){
            fprintf(stderr, "\n (-) Error in memory allocation for density matrix integrand or its square.");
            fprintf(stderr, "\n     Aborting...\n\n");
            exit(1);
        }
        int r1,r2;

        fprintf(file_ptr, "\n# Coupling:\n#");
        for(i = 0; i < n_out; i++){
        // calculate density-matrix for all wave functions
            for(r1 = 0; r1 < nq1; r1++){
                for(r2 = r1; r2 < nq1; r2++){
                    for(j = 0; j < nq2; j++){
                        dm_integrand[j] = X[i*n_points + r1*nq2 + j]*X[i*n_points + r2*nq2 + j];
                    }

                    integral = integrate_1d(nq2, dq, dm_integrand);
                    dichtematrix[r1*nq1 + r2] = integral;

                    if(r1 != r2){
                        dichtematrix[r2*nq1 + r1] = integral;
                    }
                }
            }

        // calculate density-matrix square
        //  careful: density-matrix has dimension nq1 times nq1!
            for(r1 = 0; r1 < nq1; r1++){
                for(r2 = r1; r2 < nq1; r2++){
                    for(j = 0; j < nq1; j++){
                        dm_integrand_sq[j] = dichtematrix[r1*nq1 + j]*dichtematrix[j*nq1 + r2];
                    }

                    integral = integrate_1d(nq1, dq, dm_integrand_sq);
                    dichtematrix_sq[r1*nq1 + r2] = integral;

                    if(r1 != r2){
                       dichtematrix_sq[r2*nq1 + r1] = integral;
                    }
                }
            }

        // calculate the trace of density-matrix square
            for(j = 0; j < nq1; j++){
                dm_integrand_sq[j] = dichtematrix_sq[j*nq1 + j];
            }
            integral = integrate_1d(nq1, dq, dm_integrand_sq);

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
                    integrand[k] = X[k + i*n_points]*X[k + j*n_points] * dip_x[k];
                }

                integral = integrate_2d(nq1, nq2, dq, integrand);
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
                    integrand[k] = X[k + i*n_points]*X[k + j*n_points] * dip_y[k];
                }

                integral = integrate_2d(nq1, nq2, dq, integrand);
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
                    integrand[k] = X[k + i*n_points]*X[k + j*n_points] * dip_z[k];
                }

                integral = integrate_2d(nq1,nq2, dq, integrand);
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

        free(dip_x);    dip_x    = NULL;
        free(dip_y);    dip_y    = NULL;
        free(dip_z);    dip_z    = NULL;
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
    fprintf(file_ptr, "\n# N %2d %2d", nq1, nq2);
    fprintf(file_ptr, "\n");

    for(i = 0; i < n_points; i++){
    // newline every "nq2"th line
        if(i%nq2 == 0){
            fprintf(file_ptr, "\n");
        }

    // output coordinates q1 and q2 as well as potential
        fprintf(file_ptr,"%24.16lf    %24.16lf    %24.16lf", q1[i], q2[i], v[i]);

    // output wave functions
        for(j = 0; j < n_out; j++){
            fprintf(file_ptr, "    %24.16lf", X[i + j*n_points]);
        }

        fprintf(file_ptr,"\n");
    }

    fclose(file_ptr);
    free(q1);   q1 = NULL;
    free(q2);   q2 = NULL;
    free(v);    v  = NULL;
    free(X);    X  = NULL;

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
