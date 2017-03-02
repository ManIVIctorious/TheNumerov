#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

int cubic_spline(double x[], double y[], double b[], double c [], double d[], int n_points);
int spline_interpolate(int n_x, int n_y, int n_spline, double x[], double y[], double z[]);


int spline_interpolate(int n_x, int n_y, int n_spline, double x[], double y[], double z[]){

  int i, j, k, l, m, n;

  int max_dim;
  int count = 0;

  int n_x_new = (n_x - 1) * (n_spline +1) + 1;
  int n_y_new = (n_y - 1) * (n_spline +1) + 1;

  int n_data = n_x * n_y;
  int n_points = n_x_new * n_y_new;

  int signum;

  int idx, idx_new, idx_mm, idx_mp, idx_pm, idx_pp;

  double * mixed_deriv;

  if(n_x < n_y){
    max_dim = n_y;
  }
  else{
    max_dim = n_x;
  }

  double spline_x[max_dim];
  double spline_y[max_dim];

  double b[max_dim];
  double c[max_dim];
  double d[max_dim];

  double alpha[16];
  double beta[16];

  double x_new, x_int;
  double y_int;
  double z_int;

  double first_deriv_x[n_data];
  double first_deriv_y[n_data];

  double *z_new;
  double dx, dy;

  z_new = (double *) calloc(n_points, sizeof(double) );

  gsl_matrix * spline_matrix = gsl_matrix_alloc (16, 16);
  gsl_matrix * inverse_matrix = gsl_matrix_alloc (16, 16);

  gsl_permutation * permute = gsl_permutation_alloc (16);

  double M_inv[16][16];

  double **M;

  M = (double **) calloc(16, sizeof(double*) );

  for (i=0; i< 16; i++)
  {
    M[i] = (double *) calloc(16, sizeof(double) );
  }



  // work on y-part first
  for (i = 0; i < n_x; i++)
  {
    for (j = 0; j < n_y; j++)
    {
      spline_x[j] = y[j+i*n_y];
      spline_y[j] = z[j+i*n_y];
    }

    // execute spline interpolation
    cubic_spline(spline_x, spline_y, b, c, d, n_y);

    // get y-derivatives at grid points
    for (j = 0; j < n_y; j++)
    {
      first_deriv_y[j + i * n_y]   = b[j];
    }

    // Interpolate along spline
    count = 0;

    for (j = 0; j < n_y - 1; j++)
    {
      dx = (spline_x[j + 1] - spline_x[j]) / (double) (1 + n_spline);

      for (x_new = spline_x[j]; x_new < spline_x[j + 1]-dx/2; x_new = x_new + dx)
      {
	x_int = x_new - spline_x[j];

	z_new[count + i * n_y_new*(1+n_spline)] = spline_y[j] + x_int * (b[j] + x_int * (c[j] + x_int * d[j]));

	count ++;
      }

    }// Interpolate for (j=0; j < n_points - 1; j++)

    // last point manually
    z_new[count + i * n_y_new*(1+n_spline)] = spline_y[j];

  } // work on y-part

  // printf("\nX-Part\n");

  // work on x-part
  for (i = 0; i < n_y; i++)
  {
    for (j = 0; j < n_x; j++)
    {
      spline_x[j] = x[i+j*n_y];
      spline_y[j] = z[i+j*n_y];
    }

    // execute spline interpolation
    cubic_spline(spline_x, spline_y, b, c, d, n_x);

    // get x-derivatives at grid points
    for (j = 0; j < n_x; j++)
    {
      first_deriv_x[i + j * n_y]   = b[j];
    }

    // Interpolate along spline
    count = 0;

    for (j = 0; j < n_x - 1; j++)
    {
      dx = (spline_x[j + 1] - spline_x[j]) / (double) (1 + n_spline);

      for (x_new = spline_x[j]; x_new < spline_x[j + 1]-dx/2; x_new = x_new + dx)
      {
	x_int = x_new - spline_x[j];

	z_new[count * n_y_new + i*(n_spline+1)] = spline_y[j] + x_int * (b[j] + x_int * (c[j] + x_int * d[j]));

	// printf("%d  %d x_new  %20.14lf   %20.14lf   %20.14lf  %d\n", j, count, x_new, x_int, z_new[count * n_y_new + i*(n_spline+1)], count * n_y_new + i*(n_spline+1) );

	count ++;
      }

    }// Interpolate for (j=0; j < n_points - 1; j++)

    // last point manually
    z_new[count * n_y_new + i*(n_spline+1)] = spline_y[j];

    // printf("%d  %d x_new  %20.14lf   %20.14lf   %20.14lf  %d\n", j, count, x_new, x_int, z_new[count * n_y_new + i*(n_spline+1)], count * n_y_new + i*(n_spline+1) );


  } // work on x-part


  ///////////////




  mixed_deriv = (double *) calloc(n_points, sizeof(double) );

  // get mixed derivatives

  for (i = 1; i < n_x - 1; i++)
  {
    for (j = 1; j < n_y - 1; j++)
    {
      idx = i * n_y + j;

      idx_mm = (i-1) * n_y + (j-1);
      idx_mp = (i-1) * n_y + (j+1);
      idx_pm = (i+1) * n_y + (j-1);
      idx_pp = (i+1) * n_y + (j+1);

      mixed_deriv[idx] = 0.25 * ( z[idx_mm] - z[idx_mp] - z[idx_pm] + z[idx_pp] ) / (x[idx_pm] - x[idx_mm]) /(y[idx_mp] - y[idx_mm]);
    }
  }

  // get cross-derivatives at the edge via cubic spline extrapolation


  ////////////////////
  // work on y-part first
  // careful range exclude edges
  for (i = 1; i < n_x - 1 ; i++)
  {
    for (j = 1; j < n_y - 1; j++)
    {
      spline_x[j-1] = y[j+i*n_y];
      spline_y[j-1] = mixed_deriv[j+i*n_y];
    }

    // execute spline interpolation
    cubic_spline(spline_x, spline_y, b, c, d, n_y-2);

    // extrapolate left side
    x_int = y[i*n_y] - y[i*n_y+1];

    mixed_deriv[i*n_y] = spline_y[0] + x_int * (b[0] + x_int * (c[0] + x_int * d[0]));

    // extrapolate right side -- carefule: ny-2 points, ny-3 intervals!!!
    x_int = y[(i+1)*n_y - 1] - y[(i+1)*n_y-2];

    mixed_deriv[(i+1)*n_y - 1] = spline_y[n_y - 3] + x_int * (b[n_y - 3] + x_int * (c[n_y - 3] + x_int * d[n_y - 3]));

  } // work on y-part

  // work on x-part
  // careful range exclude edges

  for (i = 1; i < n_y-1; i++)
  {
    for (j = 1; j < n_x-1; j++)
    {
      spline_x[j-1] = x[i+j*n_y];
      spline_y[j-1] = mixed_deriv[i+j*n_y];
    }

    // execute spline interpolation
    cubic_spline(spline_x, spline_y, b, c, d, n_x-2);

    // extrapolate left side
    x_int = x[i] - x[i+n_y];

    mixed_deriv[i] = spline_y[0] + x_int * (b[0] + x_int * (c[0] + x_int * d[0]));

     // extrapolate right side -- careful: ny-2 points, ny-3 intervals!!!
    x_int = x[i+n_y*(n_x-1)] - x[i+n_y*(n_x-2)];

    mixed_deriv[i+n_y*(n_x-1)] = spline_y[n_y - 3] + x_int * (b[n_y - 3] + x_int * (c[n_y - 3] + x_int * d[n_y - 3]));
  }

  // finally extrapolate cross-derivatives at corners via diagonal
/*
  // diagonal #1
    for (i = 1; i < n_x-1; i++)
    {
      spline_x[i-1] = x[i+i*n_y]/fabs((x[i+i*n_y]+0.0000000000000000001))*sqrt(x[i+i*n_y] * x[i+i*n_y] + y[i+i*n_y]*y[i+i*n_y]);
      spline_y[i-1] = mixed_deriv[i+i*n_y];
    }

    cubic_spline(spline_x, spline_y, b, c, d, n_x-2);

    // extrapolate left side
    x_int = sqrt( (y[0] - y[n_y+1]) * (y[0] - y[n_y+1]) + (x[0]-x[n_y+1])*(x[0]-x[n_y+1])) * (x[0]-x[n_y])/fabs(x[0]-x[n_y+1]) ;

    mixed_deriv[0] = spline_y[0] + x_int * (b[0] + x_int * (c[0] + x_int * d[0]));

    // extrapolate right side -- carefule: ny-2 points, ny-3 intervals!!!
    x_int = sqrt( (y[n_x*n_y-1] - y[n_x*(n_y-1) - 2]) * (y[n_x*n_y-1] - y[n_x*(n_y-1) - 2]) +
                  (x[n_x*n_y-1] - x[n_x*(n_y-1) - 2]) * (x[n_x*n_y-1] - x[n_x*(n_y-1) - 2]) ) *
      (y[n_x*n_y-1] - y[(n_x-1)*(n_y-1) - 1])/fabs(y[n_x*n_y-1] - y[(n_x-1)*(n_y-1) - 1]);

    mixed_deriv[n_x*n_y-1] = spline_y[n_y - 3] + x_int * (b[n_y - 3] + x_int * (c[n_y - 3] + x_int * d[n_y - 3]));

    // diagonal #2
    for (i = 1; i < n_x-1; i++)
    {
      spline_x[i-1] = x[n_y * (n_x-i-1)+i]/fabs((x[n_y * (n_x-i-1)+i]+0.0000000000000000001))*sqrt(x[n_y * (n_x-i-1)+i] * x[n_y * (n_x-i-1)+i] + y[n_y * (n_x-i-1)+i]*y[n_y * (n_x-i-1)+i]);
      spline_y[i-1] = mixed_deriv[n_y * (n_x-i-1)+i];

    }

    cubic_spline(spline_x, spline_y, b, c, d, n_x-2);

    // extrapolate left side
    x_int = sqrt( (y[n_y*(n_x-2) + 1] - y[n_y*(n_x-1)]) * (y[n_y*(n_x-2) + 1] - y[n_y*(n_x-1)]) +
                  (x[n_y*(n_x-2) + 1] - x[n_y*(n_x-1)]) * (x[n_y*(n_x-2) + 1] - x[n_y*(n_x-1)]) ) *
                  (x[n_y*(n_x-2) + 1] - x[n_y*(n_x-1)]) / fabs(x[n_y*(n_x-2) + 1] - x[n_y*(n_x-1)]);

    mixed_deriv[n_y*(n_x-1)] = spline_y[n_y - 3] + x_int * (b[n_y - 3] + x_int * (c[n_y - 3] + x_int * d[n_y - 3]));

    // extrapolate right side -- carefule: ny-2 points, ny-3 intervals!!!
    x_int = sqrt( (y[n_y-1] - y[2*(n_y-1)]) * (y[n_y-1] - y[2*(n_y-1)]) +
                  (x[n_y-1] - x[2*(n_y-1)]) * (x[n_y-1] - x[2*(n_y-1)]) ) *
                  (y[n_y-1] - y[2*(n_y-1)])/fabs(y[n_y-1] - y[2*(n_y-1)]);

    mixed_deriv[n_y-1] = spline_y[n_y - 3] + x_int * (b[n_y - 3] + x_int * (c[n_y - 3] + x_int * d[n_y - 3]));
*/

  mixed_deriv[n_x*n_y-1]   = 0.0;
  mixed_deriv[n_y*(n_x-1)] = 0.0;
  mixed_deriv[n_y-1]       = 0.0;
  mixed_deriv[n_x*n_y-1]   = 0.0;

    // Finally interpolate 2d
    dx = x[n_y] - x[0];
    dy = y[1]   - y[0];

    M[0][0]  = 1.0;

    M[1][0]  = 1.0;
    M[1][1]  = dx;
    M[1][2]  = dx * dx;
    M[1][3]  = dx * dx * dx;

    M[2][0]  = 1.0;
    M[2][4]  = dy;
    M[2][8]  = dy * dy;
    M[2][12] = dy * dy * dy;

    M[3][0]  = 1.0;
    M[3][1]  =  dx;
    M[3][2]  =  dx * dx;
    M[3][3]  =  dx * dx * dx;
    M[3][4]  =  dy;
    M[3][5]  =  dx * dy;
    M[3][6]  =  dx * dx * dy;
    M[3][7]  =  dx * dx * dx * dy;
    M[3][8]  =  dy * dy;
    M[3][9]  =  dx * dy * dy;
    M[3][10] =  dx * dx * dy * dy;
    M[3][11] =  dx * dx * dx * dy * dy;
    M[3][12] =  dy * dy * dy;
    M[3][13] =  dx * dy * dy * dy;
    M[3][14] =  dx * dx * dy * dy * dy;
    M[3][15] =  dx * dx * dx * dy * dy * dy;

    M[4][1]  = 1.0;

    M[5][1]  = 1.0;
    M[5][2]  = 2.0 * dx;
    M[5][3]  = 3.0 * dx * dx;

    M[6][1]  = 1.0;
    M[6][5]  =  dy;
    M[6][9]  =  dy * dy;
    M[6][13] =  dy * dy * dy;

    M[7][1]  = 1.0;
    M[7][2]  = 2.0 * dx;
    M[7][3]  = 3.0 * dx * dx;
    M[7][5]  =  dy;
    M[7][6]  = 2.0 * dx * dy;
    M[7][7]  = 3.0 * dx * dx * dy;
    M[7][9]  =  dy * dy;
    M[7][10] = 2.0 * dx * dy * dy;
    M[7][11] = 3.0 * dx * dx * dy * dy;
    M[7][13] =  dy * dy * dy;
    M[7][14] = 2.0 * dx * dy * dy * dy;
    M[7][15] = 3.0 * dx * dx * dy * dy * dy;

    M[8][4]  = 1.0;

    M[9][4]  = 1.0;
    M[9][5]  =  dx;
    M[9][6]  =  dx * dx;
    M[9][7]  =  dx * dx * dx;

    M[10][4]  = 1.0;
    M[10][8]  = 2.0 * dy;
    M[10][12] = 3.0 * dy * dy;

    M[11][4]  = 1.0;
    M[11][5]  =  dx;
    M[11][6]  =  dx * dx;
    M[11][7]  =  dx * dx * dx;
    M[11][8]  = 2.0 * dy;
    M[11][9]  = 2.0 * dx * dy;
    M[11][10] = 2.0 * dx * dx * dy;
    M[11][11] = 2.0 * dx * dx * dx * dy;
    M[11][12] = 3.0 * dy * dy;
    M[11][13] = 3.0 * dx * dy * dy;
    M[11][14] = 3.0 * dx * dx * dy * dy;
    M[11][15] = 3.0 * dx * dx * dx * dy * dy;

    M[12][5]  = 1.0;

    M[13][5]  = 1.0;
    M[13][6]  = 2.0 * dx;
    M[13][7]  = 3.0 * dx * dx;

    M[14][5]  = 1.0;
    M[14][9]  = 2.0 * dy;
    M[14][13] = 3.0 * dy * dy;

    M[15][5]  = 1.0;
    M[15][6]  = 2.0 * dx;
    M[15][7]  = 3.0 * dx * dx;
    M[15][9]  = 2.0 * dy;
    M[15][10] = 4.0 * dx * dy;
    M[15][11] = 6.0 * dx * dx * dy;
    M[15][13] = 3.0 * dy * dy;
    M[15][14] = 6.0 * dx * dy * dy;
    M[15][15] = 9.0 * dx * dx * dy * dy;

    for (i = 0; i < 16; i++)
      for (j = 0; j < 16; j++)
	gsl_matrix_set (spline_matrix, i, j, M[i][j]);

    //  inverse matrix with the help of GSL
    gsl_linalg_LU_decomp (spline_matrix, permute, &signum);
    gsl_linalg_LU_invert (spline_matrix, permute, inverse_matrix);

    // extract inverse matrix from gsl
    for (i = 0; i < 16; i++)
    {
      for (j = 0; j < 16; j++)
      {
	M_inv[i][j] = gsl_matrix_get (inverse_matrix, i, j);
      }
    }

    // loop over grid
    for (i = 0; i < n_x - 1; i++)
    {
      for (j = 0; j < n_y - 1; j++)
      {
	// generate vector beta

	beta[0]  = z[i     * n_y + j    ];
	beta[1]  = z[i     * n_y + j + 1];
	beta[2]  = z[(i+1) * n_y + j    ];
	beta[3]  = z[(i+1) * n_y + j + 1];

	beta[4]  = first_deriv_y[i     * n_y + j    ] ;
	beta[5]  = first_deriv_y[i     * n_y + j + 1] ;
	beta[6]  = first_deriv_y[(i+1) * n_y + j    ] ;
	beta[7]  = first_deriv_y[(i+1) * n_y + j + 1] ;

	beta[8]  = first_deriv_x[i     * n_y + j    ] ;
	beta[9]  = first_deriv_x[i     * n_y + j + 1] ;
	beta[10] = first_deriv_x[(i+1) * n_y + j    ] ;
	beta[11] = first_deriv_x[(i+1) * n_y + j + 1] ;

	beta[12] = mixed_deriv[i     * n_y + j    ] ;
	beta[13] = mixed_deriv[i     * n_y + j + 1] ;
	beta[14] = mixed_deriv[(i+1) * n_y + j    ] ;
	beta[15] = mixed_deriv[(i+1) * n_y + j + 1] ;

	// multply M_inv with beta
	for (k = 0; k < 16; k++)
	{
          alpha[k] = 0.0;

          for (l  = 0; l < 16; l++)
	  {
	    alpha[k] = alpha[k] + M_inv[k][l] * beta[l];
	  }
	}

	for ( m = 0; m < n_spline; m = m + 1)
	{
	  for ( n = 0; n < n_spline; n = n +1)
	    {

	      x_int = dx * (double) (m + 1) / (double) (n_spline + 1);
	      y_int = dy * (double) (n + 1) / (double) (n_spline + 1);

           z_int =  alpha[0]  +
                    alpha[1]  * pow(y_int, 1.0) +
	                alpha[2]  * pow(y_int, 2.0) +
	                alpha[3]  * pow(y_int, 3.0) +
	                alpha[4]                    * pow(x_int, 1.0) +
	                alpha[5]  * pow(y_int, 1.0) * pow(x_int, 1.0) +
	                alpha[6]  * pow(y_int, 2.0) * pow(x_int, 1.0) +
	                alpha[7]  * pow(y_int, 3.0) * pow(x_int, 1.0) +
	                alpha[8]                    * pow(x_int, 2.0) +
	                alpha[9]  * pow(y_int, 1.0) * pow(x_int, 2.0) +
	                alpha[10] * pow(y_int, 2.0) * pow(x_int, 2.0) +
	                alpha[11] * pow(y_int, 3.0) * pow(x_int, 2.0) +
	                alpha[12]                   * pow(x_int, 3.0) +
	                alpha[13] * pow(y_int, 1.0) * pow(x_int, 3.0) +
	                alpha[14] * pow(y_int, 2.0) * pow(x_int, 3.0) +
	                alpha[15] * pow(y_int, 3.0) * pow(x_int, 3.0);

	   idx = i * n_y + j;

	   idx_new = (i * (n_spline + 1) + m + 1) * n_y_new + j * (n_spline + 1) + n  + 1;

          x[idx_new] = x[i * n_y + j] + x_int;
          y[idx_new] = y[i * n_y + j] + y_int;
	   z_new[idx_new] = z_int;
	  }  // loop over interpolation
	}  // loop over interpolation
      } // loop over grid
    } // loop over grid



    /*
    printf("N %d %d\n\n", n_x_new , n_y_new);



    double y_new;

    for (i=0; i < n_points; i++)
    {
      j=i/n_y_new;
      k=i%n_y_new;

      x_new = x[0] + dx * (double) (i/n_y_new) / (double) (n_spline + 1);

      y_new = y[0] + dy * (double) (i%n_y_new) / (double) (n_spline + 1); // careful - change to dy

      // printf("%4d  %3d  %4d   -  %12.8lf   %12.8lf   %12.8lf\n", i, j, k, x_new, y_new, z_new[i]);
      fprintf(stderr,"Debug %12.8lf   %12.8lf   %12.8lf\n", x_new, y_new, z_new[i]);

      if ( k == n_y_new-1)
	    fprintf(stderr, "\n");
    }
    */


    // transfer z back
    // determine x and y in previuos subroutine
    for (i=0; i < n_points; i++)
    {
      z[i] = z_new[i];
    }

    free(M);                         M              = NULL;
    free(z_new);                     z_new          = NULL;
    free(mixed_deriv);               mixed_deriv    = NULL;
    gsl_matrix_free(spline_matrix);  spline_matrix  = NULL;
    gsl_matrix_free(inverse_matrix); inverse_matrix = NULL;
    gsl_permutation_free(permute);   permute        = NULL;

    return(0);
}

int spline_equalise(int n_points, int n_equal, double x[], double y[])
{
  int cubic_spline(double x[], double y[], double b[], double c [], double d[], int n_points);

  int i;
  int j = 0;
  int k;
  int l;

  int count = 0;

  double b[n_points];
  double c[n_points];
  double d[n_points];

  double x_new, x_int;
  double y_new[n_equal];

  double dx;

  // execute spline interpolation
  cubic_spline(x, y, b, c, d, n_points);

  dx = (x[n_points-1] - x[0]) / (double) (n_equal - 1);

  // Equalise grid
  for (i=0; i < n_equal - 1; i++)
  { 
    x_new = x[0] + dx * (double) i;

      // identfy interval
      if ((x[j] > x_new) || (x[j+1] < x_new))
      {
	j = 0;
	k = n_points;

	do
	{
	  l = (j + k) / 2;
	      
	  if (x_new < x[l])
	    k = l;
	      
	  if (x_new >= x[l])
	    j = l;

	} while (k > j+1);
      }

      x_int = x_new - x[j];

      y_new[count] = y[j] + x_int * (b[j] + x_int * (c[j] + x_int * d[j]));

      //printf("%d  %20.14lf   %20.14lf   %20.14lf  %20.14lf  %d   %d\n", i, x_new, x[k], y_new[count], y[k], j, k );

      count ++;
  }

  y_new[count] = y[n_points-1];

  for (i=0; i < n_equal; i++)
  { 
    y[i] = y_new[i];
  }

  return(1);
}


int cubic_spline(double x[], double y[], double b[], double c [], double d[], int n_points){

  int i, j;

  double ratio;

  d[0] = x[1] - x[0];
  c[1] = (y[1] - y[0]) / d[0];

  // setup - loop over n-1 intervals
  for (i = 1; i < (n_points - 1); i++)
  {
    d[i]   = x[i+1] - x[i];
    b[i]   = 2.0 * (d[i-1] + d[i]);
    c[i+1] = (y[i+1] - y[i]) / d[i];
    c[i]   = c[i+1] - c[i];
  }

  // standard end conditions
  b[0]            = -d[0];
  b[n_points - 1] = -d[n_points-2];
  c[0]            = 0.0;
  c[n_points - 1] = 0.0;

  if (n_points > 3)
  {
    c[0]          = c[2] / (x[3] - x[1]) - c[1] / (x[2] - x[0]);
    c[n_points-1] = c[n_points-2] / (x[n_points-1] - x[n_points-3]) - c[n_points-3] / (x[n_points-2] - x[n_points-4]);
    c[0]          = c[0] * d[0] * d[0] / (x[3] - x[0]);
    c[n_points-1] = -c[n_points-1] * d[n_points-2] * d[n_points-2] / (x[n_points-1] - x[n_points-4]);
  }

  // Forward elimination
  for (i = 1; i < n_points; i++)
  {
    ratio = d[i-1] / b[i-1];
    b[i]  = b[i] - ratio * d[i-1];
    c[i]  = c[i] - ratio * c[i-1];
  }

  // Back substitution
  c[n_points-1] = c[n_points-1] / b[n_points-1];

  for (i = 0; i < (n_points - 1); i++)
  {
    j    = n_points - i - 2;
    c[j] = (c[j] - d[j] * c[j+1]) / b[j];
  }


  // Evaluate spline coefficient b[], c[] and d[]
  b[n_points-1] = ( y[n_points-1] - y[n_points-2] ) / d[n_points-2] + d[n_points-2] * ( c[n_points-2] + 2.0 * c[n_points-1] );

  for (i = 0; i < (n_points - 1); i++)
  {
    b[i] = (y[i+1] - y[i]) / d[i] - d[i] * (c[i+1] + 2.0 * c[i]);
    d[i] = (c[i+1] - c[i]) / d[i];
    c[i] = 3.0 * c[i];
  }
  c[n_points-1] = 3.0 * c[n_points-1];
  d[n_points-1] = d[n_points-2];

    return 0;
}

