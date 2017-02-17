#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include "spline_interpolate.h"
#include "mkl_solvers_ee.h"

int get_stencil(double stencil[], int n_stencil);
double integrate_1d(int n, double dx, double integrand[]);
double integrate_2d(int nx, int ny, double dx, double integrand[]);
int InputFunction(char *inputfile, double **q1, double **q2, double **V, int *nq1, int *nq2);
int InputFunctionDipole(char *inputfile, double **q1, double **q2, double **V, double **mux, double **muy, double **muz, int *nq1, int *nq2);

// used constants
//
// factor for -h_bar^2/2m for reduced mass in g
// h_bar                 1.05457180013E-34 Js
// Avogadro's constant   6.02214085774E23 1/mol
// atomic mass unit      1.66053904020E-27 kg
//
// Conversion factors    1.0E20          Ang^2 / m^2
//                       1.0 / 4184.0    kcal / J

int main(int argc, char* argv[]){

  int i, j, k, l, xsh, ysh;
  int nq1, nq2;

  int control;
  int element;

  int index;
  int option_index  = 0;
  int n_points      = 0;
  int n_ts_dip      = 0;
  int n_stencil     = 9;
  int n_entries     = 0;
  int max_entries   = 0;

  int n_spline      = 0;
  int analyse       = 0;

  double dq,dx,dy;
  double mass = 1.0;
  double ekin_param = -1.05457180013E-34 * 1.05457180013E-34 / 2.0 / 1.66053904020E-27 * 1.0E20 * 6.02214085774E23 / 4184.0;
  double kcal_per_mol_to_inv_cm = 219474.6313705/627.509469;

  double ekin_factor = 1.0;
  double epot_factor = 1.0;
  double e_min = 0.0;
  double e_max = 100.0;
  double spacing_threshold = 1.0E-12;;

  double v_min = 1.0E100;

  double freq;
  double integral;
  double ts_dip_square;

  double * q1 = NULL;
  double * q2 = NULL;
  double * v = NULL;

  double * dip_x = NULL;
  double * dip_y = NULL;
  double * dip_z = NULL;

  double * ts_dip_x = NULL;
  double * ts_dip_y = NULL;
  double * ts_dip_z = NULL;

  double * stencil;
  double * integrand;

  int dipole_flag = 0;
  char *input_file_name  = NULL;
  char *output_file_name = "/dev/stdout";

  FILE * file_ptr;

/// FLAGS // FLAGS /// FLAGS // FLAGS /// FLAGS // FLAGS /// FLAGS // FLAGS /// FLAGS // FLAGS
    while(1){
        static struct option long_options[] = {
            {"help",             no_argument, 0, 'h'},
            {"mass",       required_argument, 0, 'm'},
            {"fkin",       required_argument, 0, 'k'},
            {"fpot",       required_argument, 0, 'v'},
            {"n_stencil",  required_argument, 0, 'n'},
            {"e_min",      required_argument, 0, 'l'},
            {"e_max",      required_argument, 0, 'u'},
            {"spline",     required_argument, 0, 's'},
            {"pipe",             no_argument, 0, 'P'},
            {"analye",           no_argument, 0, 'a'},
            {"dipole",           no_argument, 0, 'd'},
            {"in-file",    required_argument, 0, 'i'},
            {"out-file",   required_argument, 0, 'o'},
            { 0, 0, 0, 0 }
        };
    // getopt_long stores the option index here.
        i = getopt_long(argc, argv, "hm:k:v:n:l:u:s:adi:Po:", long_options, &option_index);

    // Detect the end of the options.
        if(i == -1)
            break;

        switch(i){
            case 0:
                /* If this option set a flag, do nothing else now. */
                if(long_options[option_index].flag != 0)    break;

                printf("option %s", long_options[option_index].name);
                if(optarg)  printf(" with arg %s", optarg);

                printf("\n");
                break;

            case 'h':
//                Help(argv[0]);
                exit (0);

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
//               Help(argv[0]);
                exit(0);
        }
    }

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
    if(dipole_flag == 1){
        dip_x = malloc(sizeof(double));
        dip_y = malloc(sizeof(double));
        dip_z = malloc(sizeof(double));
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
    control = get_stencil(stencil, n_stencil);

    if(control != 0 ){
        fprintf(stderr, "\n (-) Error initialising stencil parameters.");
        fprintf(stderr, "\n     Aborting - please check your input...\n\n");
        exit(1);
    }


//------------------------------------------------------------------------------------------------------------------
// Body  Body  Body  Body  Body  Body  Body  Body  Body  Body  Body  Body  Body  Body  Body  Body  Body  Body  Body
//------------------------------------------------------------------------------------------------------------------
// get potential minimum and subtract it from the potential
    for(i = 0; i < n_points; ++i){
        if(v[i] < v_min)   v_min = v[i];
    }
    for(i = 0; i < n_points; ++i){
        v[i] -= v_min;
    }

// get spacing interval dx and check for uniform spacing
    dq = q2[1] - q2[0];

// testing equal spacing is more complicated in two dimensions.
//  for (i = 1; i < n_points-1; i++){
//    if ( fabs(q1[i+1] - q1[i] - dq) > spacing_threshold){
//      printf("\n\n (-) Error reading data from input-file: '%s'", input_file_name);
//      printf(  "\n     Data not uniformly spaced. Exiting ... \n\n");
//      exit(0);
//    }
//  }


//------------------------------------------------------------------------------------------------------------------
//  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline  Spline
//------------------------------------------------------------------------------------------------------------------
    if(n_spline > 0){
    // new number of points by splining
        n_points = ((nq1 - 1) * (n_spline + 1) + 1) * ((nq2 - 1) * (n_spline + 1) + 1);

    // reallocate memory and call spline function
        q1  = realloc(q1, n_points * sizeof(double));
        q2  = realloc(q2, n_points * sizeof(double));
        v   = realloc(v,  n_points * sizeof(double));
        spline_interpolate(nq1, nq2, n_spline, q1, q2, v);

    // set new values for number of points for q1 and q2 and new dq
        nq1 = (nq1 - 1) * (n_spline + 1) + 1;
        nq2 = (nq2 - 1) * (n_spline + 1) + 1;
        dq = dq / (double) (n_spline + 1);

    // set new values for q1 and q2:
    //  add 1 step to q1 all "nq2"th iteration
    //  add 1 step to q2 every iteration, reset to 0 at every "nq2"th
        for(i = 0; i < n_points; i++){
            q1[i] = q1[0] + dq * (double) (i/nq2);
            q2[i] = q2[0] + dq * (double) (i%nq2);
        }
    }

// apply kinetic energy factor and spacing to ekin_param
    ekin_param = ekin_param * ekin_factor / dq / dq / mass;
//printf("%lf, x0 = %lf, x1 = %lf \n",ekin_param,q2[0],q2[1]);


//------------------------------------------------------------------------------------------------------------------
// MKL FEAST eigenvalue solver  MKL FEAST eigenvalue solver  MKL FEAST eigenvalue solver MKL FEAST eigenvalue solver
//------------------------------------------------------------------------------------------------------------------
    char          UPLO = 'F';
    const MKL_INT N = n_points;
    MKL_INT       rows_A[n_points+1];

    rows_A[0] = 1;

    // calculate max_entries
    int sum_q1 = nq1;
    int sum_q2 = nq2;

    for(i = 1; i < (n_stencil/2+1); i++){
        sum_q1=sum_q1 + 2*(nq1-i);
        sum_q2=sum_q2 + 2*(nq2-i);
    }
    max_entries = sum_q1*sum_q2; // upper estimation for nnz entries in the matrix, but the easy way to code.


    MKL_INT     *cols_A = (MKL_INT *) malloc( max_entries * sizeof(MKL_INT));
    double      *vals_A = (double * ) malloc( max_entries * sizeof(double) );

    for(i = 0; i < nq1; i++){
        for(j = 0; j < nq2; j++){
            for(xsh = -n_stencil/2; xsh < n_stencil/2 + 1; xsh++){
                
                if( (i+xsh > -1) && (i+xsh < nq1) ){
                    for(ysh = -n_stencil/2; ysh < n_stencil/2 + 1; ysh++){

                        if( (j+ysh > -1) && ( j+ysh < nq2) ){
                            element = (i + xsh)*nq2 + j+ysh;
                            cols_A[n_entries] = element+1; // wieso +1? weil intel!!

                        // stencil entries have to be divised by 2 to get the right result.
                        //  in three dimensions it should be a division by 4
                            vals_A[n_entries] = ekin_param * stencil[(xsh+n_stencil/2)*n_stencil+ysh+n_stencil/2]/2;

                        // add potential to diagonal element
                            if(xsh == 0 && ysh ==0){
                                vals_A[n_entries] = vals_A[n_entries] + v[i*nq2+j] * epot_factor;
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



///// START EIGENVALUE CALCULATION

    MKL_INT      fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */

    double       epsout;        /* Relative error on the trace */
    MKL_INT      loop;          /* Number of refinement loop */
    MKL_INT      L = n_points/2;
    MKL_INT      M0;            /* Initial guess for subspace dimension to be used */
    MKL_INT      n_out = 5;


    double       E[n_points];         /* Eigenvalues */

    double *X;

    X = calloc (n_points*(n_points-1)/2, sizeof (double));  /* Eigenvectors */

    double       res[n_points];       /* Residual */

    MKL_INT      info;          /* Errors */

    char         SGEMMC = 'T';   /* Character for GEMM routine, transposed case */
    char         SGEMMN = 'N';   /* Character for GEMM routine, non-transposed case */
    double       one = 1.0;      /* alpha parameter for GEMM */
    double       zero = 0.0;     /* beta  parameter for GEMM */
    MKL_INT      ldx = n_points; /* Leading dimension for source arrays in GEMM */
    MKL_INT      ldy = n_points; /* Leading dimension for destination array in GEMM */

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
    if ( info != 0 )
    {
        printf(" (-) Routine sfeast_scsrev returns code of ERROR: %i\n\n\n", (int)info);
        return 1;
    }
//'#################################################################################################################
//'#################################################################################################################
//'#################################################################################################################
    // get norm  <-------------- sollte nun auch 2d lafn
    integrand = (double *) calloc(n_points, sizeof(double) );

    for (i = 0; i < n_out; i++)
    {

      for (j = 0; j < n_points; j++)
      {
        integrand[j] = X[j+i*n_points] * X[j+i*n_points];
      }

      integral = integrate_2d(nq1,nq2, dq, integrand);

      for (j = 0; j < n_points; j++)
      {
    X[j+i*n_points] = X[j+i*n_points] / sqrt(integral);
      }
  }

  if ( (file_ptr = fopen(output_file_name,"w")) == NULL)
  {
    printf("\n\n (-) Error opening output-file: '%s'", output_file_name);
    printf(  "\n     Exiting ... \n\n");

    exit(0);
  }
//'#################################################################################################################
//'#################################################################################################################
//'#################################################################################################################
  // Output eigenvalues
  fprintf(file_ptr, "# Eigenvalues");

  for (i=0; i < n_out; i++)
  {
    fprintf(file_ptr, "  %24.16lf", E[i]);
  }

  fprintf(file_ptr, "\n# Mass         %24.16lf", mass);

  fprintf(file_ptr,"\n#");

    fprintf(file_ptr, "\n# Frequencies:\n#");

    fprintf(file_ptr, "\n# ");

    for (i = 0; i < (n_out - 1); i++)
      fprintf(file_ptr,"%11d   ", i);

    element = 0;
    for (i = 1; i < n_out; i++)
    {
      fprintf(file_ptr, "\n#%3d",i);

      for (j = 0; j < i; j++)
      {
    freq = 219474.6313705/627.509469 * epot_factor * (E[i] - E[j]);
    fprintf(file_ptr, "  %12.5e", freq);
      }
  }
//'#################################################################################################################
//'#################################################################################################################
//'#################################################################################################################

  if (analyse == 1)
  {
    // Orthogonality output

    fprintf(file_ptr, "\n#");

    fprintf(file_ptr, "\n# Orthonormality:\n#");

    fprintf(file_ptr, "\n# ");

    for (i=0; i < n_out; i++)
      fprintf(file_ptr,"%11d   ", i);

    for (i = 0; i < n_out; i++)
    {
      fprintf(file_ptr, "\n#%3d",i);

      for (j = 0; j < (i+1); j++)
      {
    for (k = 0; k < n_points; k++)
        {
      // generate integrand
      integrand[k] = X[k+i*n_points]*X[k+j*n_points];
    }

    integral = integrate_2d(nq1, nq2, dq, integrand);

    fprintf(file_ptr, "  %12.5e", integral);
      }
    }
    fprintf(file_ptr, "\n#");

    // Potential output

    fprintf(file_ptr, "\n#");

    fprintf(file_ptr, "\n# Potential:\n#");

    fprintf(file_ptr, "\n# ");

    for (i=0; i < n_out; i++)
      fprintf(file_ptr,"%11d   ", i);

    for (i = 0; i < n_out; i++)
    {
      fprintf(file_ptr, "\n#%3d",i);

      for (j = 0; j < (i+1); j++)
      {
    for (k = 0; k < n_points; k++)
        {
      // generate integrand
      integrand[k] = X[k+i*n_points]*X[k+j*n_points] * v[k];
    }

    integral = integrate_2d(nq1, nq2, dq, integrand);

    fprintf(file_ptr, "  %12.5e", integral);
      }
    }
    fprintf(file_ptr, "\n#\n#");
//'#################################################################################################################
//'#################################################################################################################
// STENCIL TO BE APPLIED
//'#################################################################################################################
    // kinetic energy output

    fprintf(file_ptr, "\n# Kinetic Energy:\n#");

    fprintf(file_ptr, "\n# ");

    for (i=0; i < n_out; i++)
      fprintf(file_ptr,"%11d   ", i);

    for (i = 0; i < n_out; i++)
    {
      fprintf(file_ptr, "\n#%3d",i);

      for (j = 0; j < (i+1); j++)
      {
//-----------------------------------------------------------------------------------------------------
        for (k = 0; k < nq1; k++)
        {
           for (l = 0; l < nq2; l++)
           {
              index = k*nq2+l;
              integrand[index]=0;
              for( xsh = -n_stencil/2; xsh < n_stencil/2+1; xsh++)
              {
             if ( (k+xsh > -1) && ( k+xsh < nq1))
             {
                   for (ysh = -n_stencil/2; ysh < n_stencil/2 + 1; ysh++)
                   {
                      if ( (l+ysh > -1) && ( l+ysh < nq2))
                  {
                         element = (k + xsh)*nq2 + l+ysh;
                        integrand[index] = integrand[index] + X[element+i*n_points] * stencil[(xsh+n_stencil/2)*n_stencil+ysh+n_stencil/2]/2; // durch d^2 ist dringend nötig aber schon in ekin_param enthalten.
                      }
                   }
                 }
              }
              integrand[index]=integrand[index]*ekin_param*X[index+j*n_points];
           }
        }
//--------------------------------------------------------------------------------------------------
/*  for (k = 0; k < n_points; k++)
        {
     // reset integrand because
     // ekin-operator has to be applyed via a loop
     integrand[k] = 0.0;

     // apply stencil
         for (l = 0; l < n_stencil + 1; l++)
     {
       element = k+l-n_stencil/2;

       if ( (element > -1) && (element < (n_points +1) ) )
       {
         integrand[k] = integrand[k] + X[element+i*n_points] * stencil[l];
       }
     }

     integrand[k] = ekin_param * integrand[k] * X[k+j*n_points];

    }*/
//-------------------------------------------------------------------------------------------------------

    integral = integrate_2d(nq1, nq2, dq, integrand);

    fprintf(file_ptr, "  %12.5e", integral);
      }
    }

    fprintf(file_ptr,"\n#");

  // calculate coupling
//#####################################################################################
///*
double *dichtematrix    = (double *) calloc (nq1*nq1, sizeof (double));
double *dichtematrix_sq = (double *) calloc (nq1*nq1, sizeof (double));
//printf("im start\n");
double *dm_integrand    = (double *) calloc (nq2, sizeof (double));
double *dm_integrand_sq = (double *) calloc (nq1, sizeof (double));
int r1,r2;
     fprintf(file_ptr, "\n# COUPLING\n#");
for (i = 0; i < n_out; i++)// for all psi
    {
    // calculate dichtematrix
    for (r1=0;r1<nq1;r1++)
        {
        for (r2=r1; r2<nq1;r2++)
            {
             for(j=0;j<nq2;j++)
                {
                  dm_integrand[j]=X[i*n_points + r1*nq2 + j]*X[i*n_points + r2*nq2 + j];

                }//endfor j
             integral = integrate_1d(nq2, dq, dm_integrand);
             dichtematrix[r1*nq1+r2]= integral;
             if(r1!=r2)
             dichtematrix[r2*nq1+r1]= integral;
            }//endfor r2
        }//endfor r1

      //calculate dichtematrix quadrat
    for (r1=0;r1<nq1;r1++)
        {
        for (r2=r1; r2<nq1;r2++)
            {
             for(j=0;j<nq1;j++)// because dichtematrix is nq1 times nq1
                {
                  dm_integrand_sq[j]=dichtematrix[r1*nq1+j]*dichtematrix[j*nq1+r2];

                }//endfor j
             integral = integrate_1d(nq1, dq, dm_integrand_sq);
             dichtematrix_sq[r1*nq1+r2]= integral;
             if(r1!=r2)
             dichtematrix_sq[r2*nq1+r1]= integral;
            }//endfor r2
        }//endfor r1

        // jetzt sollte noch die spur berechnet werden.
     for(j=0;j<nq1;j++)
     {
      dm_integrand_sq[j] = dichtematrix_sq[j*nq1+j];
     }
     integral =  integrate_1d(nq1, dq, dm_integrand_sq);

     fprintf(file_ptr, "# state %d: %2.15lf \n",i,integral);
     //free(dm_integrand_sq);
     //free(dm_integrand);
     //free(dichtematrix);
     //free(dichtematrix_sq);
    }// endfor i
//#####################################################################################
//printf("ENDE\n");//*/
  }// if (analyse == 1)
//'#################################################################################################################
//'#################################################################################################################
//'#################################################################################################################
// NOCH ZU ERLEDIGEN
  if (dipole_flag == 1)
  {

    n_ts_dip = n_out * (n_out)/2;

    ts_dip_x = (double *) malloc( n_ts_dip * sizeof(double) );
    ts_dip_y = (double *) malloc( n_ts_dip * sizeof(double) );
    ts_dip_z = (double *) malloc( n_ts_dip * sizeof(double) );

    fprintf(file_ptr, "\n#");

    fprintf(file_ptr, "\n# Dipole - x-component\n#");

    fprintf(file_ptr, "\n# ");

    for (i=0; i < n_out; i++)
      fprintf(file_ptr,"%11d   ", i);

    element = 0;

    for (i = 0; i < n_out; i++)
    {
      fprintf(file_ptr, "\n#%3d",i);

      for (j = 0; j < (i+1); j++)
      {
    for (k = 0; k < n_points; k++)
        {
      // generate integrand
      integrand[k] = X[k+i*n_points]*X[k+j*n_points] * dip_x[k];
    }

    integral = integrate_2d(nq1, nq2, dq, integrand);

    fprintf(file_ptr, "  %12.5e", integral);

    if (i != j)
    {
      ts_dip_x[element] = integral;
      element ++;
    }
      }
    }
    fprintf(file_ptr, "\n#\n#");

    fprintf(file_ptr, "\n# Dipole - y-component\n#");

    fprintf(file_ptr, "\n# ");

    for (i = 0; i < n_out; i++)
      fprintf(file_ptr,"%11d   ", i);

    element = 0;

    for (i = 0; i < n_out; i++)
    {
      fprintf(file_ptr, "\n#%3d",i);

      for (j = 0; j < (i+1); j++)
      {
    for (k = 0; k < n_points; k++)
        {
      // generate integrand
      integrand[k] = X[k+i*n_points]*X[k+j*n_points] * dip_y[k];
    }

    integral = integrate_2d(nq1, nq2, dq, integrand);

    fprintf(file_ptr, "  %12.5e", integral);

    if (i != j)
    {
      ts_dip_y[element] = integral;
      element ++;
    }
      }
    }
    fprintf(file_ptr, "\n#\n#");

    fprintf(file_ptr, "\n# Dipole - z-component\n#");

    fprintf(file_ptr, "\n# ");

    for (i = 0; i < n_out; i++)
      fprintf(file_ptr,"%11d   ", i);

    element = 0;

    for (i = 0; i < n_out; i++)
    {
      fprintf(file_ptr, "\n#%3d",i);

      for (j = 0; j < (i+1); j++)
      {
    for (k = 0; k < n_points; k++)
        {
      // generate integrand
      integrand[k] = X[k+i*n_points]*X[k+j*n_points] * dip_z[k];
    }

    integral = integrate_2d(nq1,nq2, dq, integrand);

    fprintf(file_ptr, "  %12.5e", integral);

    if (i != j)
    {
      ts_dip_z[element] = integral;
      element ++;
    }
      }
    }
//<------------------------------------------------

    fprintf(file_ptr, "\n#\n#");

    fprintf(file_ptr, "\n# Oscillator strength\n#");

    fprintf(file_ptr, "\n# ");

    for (i = 0; i < (n_out - 1); i++)
      fprintf(file_ptr,"%11d   ", i);

    element = 0;
    for (i = 1; i < n_out; i++)
    {
      fprintf(file_ptr, "\n#%3d",i);

      for (j = 0; j < i; j++)
      {
    ts_dip_square =  ts_dip_x[element] * ts_dip_x[element] + ts_dip_y[element] * ts_dip_y[element]+ ts_dip_z[element] * ts_dip_z[element];

    freq = 219474.6313705/627.509469 * epot_factor * (E[i] - E[j]);

    fprintf(file_ptr, "  %12.5e", 4.702E-7 * ts_dip_square * freq);
    //fprintf(file_ptr, "  %12.5e", ts_dip_z[element]);

    element ++;
      }
    ///////////////////

    }
    fprintf(file_ptr, "\n#\n#");
  }
//}
//'#################################################################################################################
//'#################################################################################################################
  fprintf(file_ptr,"\n# Potential and Eigenfunctions: %d datapoints\n", n_points);

  // Output eigenfunctions
  for (i = 0; i < n_points; i++)
  {
    if(i%nq2==0)
      fprintf(file_ptr, "\n");

    fprintf(file_ptr,"%24.16lf %24.16lf    %24.16lf", q1[i], q2[i], v[i]);

    for (j=0; j<n_out; j++)
    {
      fprintf(file_ptr, "  %24.16lf", X[i+j*n_points]);
    }

    fprintf(file_ptr,"\n");
  }

  fprintf(file_ptr,"\n\n");

  fclose(file_ptr);


  printf("\n\n");

  exit (0);
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
