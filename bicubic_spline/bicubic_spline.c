#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>


int main(int argc, char* argv[])
{
  int spline_interpolate(int n_x, int n_y, int n_spline, double x[], double y[], double z[]);

  int i, j, k;
  int control;
  int n_error = 0;
  int line_number = 0;

  int n_pot    = 0;
  int n_points = 0;
  int n_spline;


  int nx = 0;
  int ny = 0;
  int n_x_new, n_y_new;

  double dx, dy;
  double x_new, y_new;

  int max;

  char c;
  char dummy[128];
  char line[2048];

  char * pEnd;
  char * inputfile_name = NULL;

  double * x = NULL;
  double * y = NULL;
  double * v = NULL;

  FILE * infile;

  // check input argument if the file is not present give a silly statement
  if (argv[1] == NULL) 
  { 
    printf("\n\n (-) Please specify an input file ... \n\n");
    exit (1);
    // usage
  }

  // define string for input filename
  inputfile_name = (char*) calloc(1+(strlen(argv[1])),sizeof(char));
  strcpy(inputfile_name, argv[1]);   
  strcat(inputfile_name, "\0");

  // check input argument if the file is not present give a silly statement
  if (argv[2] == NULL) 
  { 
    printf("\n\n (-) Please specify the number of spline points ... \n\n");
    exit (1);
    // usage
  }

  n_spline = atoi(argv[2]);

  // get memory and read input data 
  x  = (double *) malloc( sizeof(double) );
  v  = (double *) malloc( sizeof(double) );

  if ( (infile = fopen(inputfile_name,"r")) == NULL)
  {  
    printf("\n\n\t (-) Error opening input-file: '%s'", inputfile_name);
    printf(  "\n\t     Exiting ... \n\n");

    exit(0);
  } 

  while ( (fgets(line, 2048, infile) ) != NULL )
  {
    // remove leading white space
    // remove leading blanks and tabs  
    i=0;

    line_number ++;

    do 
    {
      c=line[i];
      if ((c == ' ') || (c == '\t'))  
	for (j = i; j < strlen(line) ; j++)        
	  line[j] = line[j+1];
      else 
	i++;        
    }
    while (c == ' '); 

    // if line is  blank read next line 
    if (line[0] == '\n') continue;

    // if line is a comment read next line 
    if (line[0] == '#') continue;

    // check for n-entry 
    control=sscanf(line,"%s", dummy);

    i=0;
    while (dummy[i])
    {
      dummy[i] = tolower(dummy[i]);
      ++i;
    }
    dummy[i]='\0';

    // assign stepcounter
    if (strcmp(dummy,"n") == 0) 
    {
      control=sscanf(line,"%s   %d   %d", dummy, &nx, &ny);
      continue;
    }

    x  = (double *) realloc(x, (n_pot + 1) * sizeof(double) );
    y  = (double *) realloc(y, (n_pot + 1) * sizeof(double) );
    v  = (double *) realloc(v, (n_pot + 1) * sizeof(double) );
     
    control=sscanf(line,"%lf  %lf %lf", &x[n_pot], &y[n_pot], &v[n_pot]);

    if (control != 3)
    {
       printf("\n\n\t (-) Error reading data in line %d from input-file: '%s'",  line_number, inputfile_name);
       printf(  "\n\t     Exiting ... \n\n");

      exit(0);
    }

    n_pot ++;
  }

  fclose(infile);

  if (nx < 1)
  {
    printf("\n\n\t (-) Error reading data from input-file: '%s'", inputfile_name);
    printf(  "\n\t     Number of points in x-direction lower or equal zero.");
    printf(  "\n\t     Exiting ... \n\n");

    ++ n_error;
  }

  if (ny < 1)
  {
    printf("\n\n\t (-) Error reading data from input-file: '%s'", inputfile_name);
    printf(  "\n\t     Number of points in y-direction lower or equal zero.");
    printf(  "\n\t     Exiting ... \n\n");

    ++ n_error;
  }

  if ( (nx*ny) != n_pot)
  {
      printf("\n\n\t (-) Error reading data from input-file: '%s'", inputfile_name);
      printf(  "\n\t     Number of points %d does not form a %dx%d-grid.", n_pot, nx, ny);
      printf(  "\n\t     Exiting ... \n\n");

      ++ n_error;
  }


  if (n_error !=  0)
  {

    printf( "\n\n (-) %d Error(s) encoutered while reading the input.",n_error);
    printf(   "\n     Exiting - please correct the data file ...\n\n\n\n");

    exit (1);
  }

    // Cubic spline interpolation of the potential 
    if (n_spline > 0)
    {
      n_x_new = (nx - 1) * (n_spline +1) + 1;
      n_y_new = (ny - 1) * (n_spline +1) + 1;

      n_points = n_x_new * n_y_new;

      x  = (double *) realloc(x, (n_points + 1) * sizeof(double) );
      y  = (double *) realloc(y, (n_points + 1) * sizeof(double) );
      v  = (double *) realloc(v, (n_points + 1) * sizeof(double) );

      spline_interpolate(nx, ny, n_spline, x, y, v);
    


     // Finally interpolate 2d
      dx = x[ny]  - x[0];
      dy = y[1]   - y[0];

      for (i=0; i < n_points; i++)
      {
	j=i/n_y_new;
	k=i%n_y_new;

	x_new = x[0] + dx * (double) (i/n_y_new) / (double) (n_spline + 1);
	
	y_new = y[0] + dy * (double) (i%n_y_new) / (double) (n_spline + 1); // careful - change to dy

	// printf("%4d  %3d  %4d   -  %12.8lf   %12.8lf   %12.8lf\n", i, j, k, x_new, y_new, z_new[i]);
	printf(" %12.8lf   %12.8lf   %12.8lf\n", x_new, y_new, v[i]);

	if ( k == n_y_new-1)
	  printf("\n");
      }

    } // if (n_spline > 0)
}



