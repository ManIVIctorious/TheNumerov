#include <stdio.h>
#include <stdlib.h>



int spline1d(int n_points, int n_spline, double x[], double y[])
{
  int i;
  int j = 0;
  int k;
  int l;

  int count = 0;
  int n_new = (n_points - 1) * (n_spline + 1) + 1;

  double b[n_points];
  double c[n_points];
  double d[n_points];

  double x_new, x_int;
  double y_new[n_new];

  double dx;
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


  // Interpolate
  for (i=0; i < n_points - 1; i++)
  { 
    dx = (x[i+1] - x[i]) / (double) (1 + n_spline);
 
    for (x_new = x[i]; x_new < x[i+1]-dx/2; x_new = x_new + dx)
    {
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

      x_int = x_new - x[i];

      y_new[count] = y[i] + x_int * (b[i] + x_int * (c[i] + x_int * d[i]));

      //printf("%d   x_new  %20.14lf      %20.14lf  %d\n", count, x_new, y_new[count], i);

      count ++;
    }
  }

  y_new[count] = y[i];

  x_new = x[0];

  for (i=0; i < n_new; i++)
  { 
    //x[i] = x_new + dx * (double) i;
    y[i] = y_new[i];
    
    //printf("%12.8lf   %12.8lf\n", x[i], y[i]);
  }

    //printf("\n\n");
    //exit(0);
    return(1);
}

