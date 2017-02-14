#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

int main(int argc, char* argv[])
{
  int i, j, k;

  int nx = 0;
  int ny = 0;

  int idx;
  int idx_x;
  int idx_y;

  int zero_x = -1;
  int zero_y = -1;

  int n_points = 0;
  int n_entries = -1;
  int number;

  int count;
  int control;
  int line_number = 0;

  int marker;
  int position;

  int n_min;
  int nx_min;
  int ny_min;

  char c;

  char * inputfile_name = NULL;
  char keyword[128];
  char line[2048];
  char string[2048] = {0};

  double dr_x = -0.1;
  double dr_y = -0.1;

  double red_m_x = -0.1;
  double red_m_y = -0.1;

  double e_min = 1.0E100;
  double integral;

  double * pot_x;
  double * pot_y;

  double ** data;
  double * integrand;

  FILE * infile;

  double integral_2d(double * integrand, int nx, int ny, double dx, double dy);

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

  if ( (infile = fopen(inputfile_name,"r")) == NULL)
  {  
    printf("\n\n\t (-) Error opening input-file: '%s'", inputfile_name);
    printf(  "\n\t     Exiting ... \n\n");

    exit(0);
  } 

  // first read of infile for initialisation
  while ( (fgets(line, 2048, infile) ) != NULL )
  {
    line_number ++;

    // remove leading white space
    // remove leading blanks and tabs  
    i=0;
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

    marker = 0;
   // check if an ENTER is somewhere in the string
    for (i = 0; i < strlen(line)  ; i++)
    {
       c=line[i];
      if (c == '\n') marker = 1;
      if (marker == 1) line[i] = 0;  
    }       

   // terminate if no ENTER was found 
    if (marker == 0)
    {
       printf("\n\n (-) Invalid instruction line %d - the line is too long.",line_number);
       printf(  "\n     Exiting - please correct the input file ...\n\n\n\n");
       exit (1); 
    }

    // check for keyword-entries
    control=sscanf(line,"%s", keyword);

    i=0;
    while (keyword[i])
    {
      keyword[i] = tolower(keyword[i]);
      ++i;
    }
    keyword[i]='\0';

    // assign stepcounter
    if (strcmp(keyword,"n") == 0) 
    {
      control=sscanf(line, "%s   %d   %d", keyword, &nx, &ny);
      continue;
    }
    else if (strcmp(keyword, "zero") == 0)
    {
      control=sscanf(line,"%s   %d   %d", keyword, &zero_x, &zero_y);
      continue;
    }
    else if (strcmp(keyword, "dr") == 0)
    {
      control=sscanf(line,"%s   %lf   %lf", keyword, &dr_x, &dr_y);
      continue;
    } 
    else if (strcmp(keyword, "red_mass") == 0)
    {
      control=sscanf(line,"%s   %lf   %lf", keyword, &dr_x, &dr_y);
      continue;
    } 

    // read entries
    count = 0;
    while (sscanf(line,"%s", &string) != -1)
    { 
      // eliminate recent number from line
      position = 0;
      while( (line[position] != ' ') && (line[position] != '\t') )
	position ++;

      j = 0;
      for (i = position; i < strlen(line); i++)
	{
	  line[j] = line[i];
	  j++;
	}
      line[j] = '\0';
    
      // remove additional leading whitespace characters
      i=0;
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

      count ++;

    }
    
    // check consistency in number of entries
    if (n_points == 0)
    {

      if( count < 5)
      {
       printf("\n\n\t (-) Error reading data in line %d from input-file: '%s'",  line_number, inputfile_name);
       printf(  "\n\t     Less than five entries per line detected.");
       printf(  "\n\t     Exiting ... \n\n");
      }

      n_entries = count;
    }
    else if( count != n_entries)
    {
       printf("\n\n\t (-) Error reading data in line %d from input-file: '%s'",  line_number, inputfile_name);
       printf(  "\n\t     Inconsisten number of entries per line.");
       printf(  "\n\t     Exiting ... \n\n");
    }

    n_points ++;

  } // end read of infile for initialisation

  rewind(infile);

  // initialize data arras
  data = (double **) calloc(n_entries + 1, sizeof(double *) );

  for (i = 0; i < (n_entries+1); i++)
    data[i] = (double *) calloc(n_points, sizeof(double) );

  integrand = (double *) calloc(n_points, sizeof(double) );

  // second read to store data
  line_number = 0;
  n_points = 0;

  while ( (fgets(line, 2048, infile) ) != NULL )
  {
    line_number ++;

    // remove leading white space
    // remove leading blanks and tabs  
    i=0;
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

    marker = 0;
   // check if an ENTER is somewhere in the string
    for (i = 0; i < strlen(line)  ; i++)
    {
       c=line[i];
      if (c == '\n') marker = 1;
      if (marker == 1) line[i] = 0;  
    }       

   // terminate if no ENTER was found 
    if (marker == 0)
    {
       printf("\n\n (-) Invalid instruction line %d - the line is too long.",line_number);
       printf(  "\n     Exiting - please correct the input file ...\n\n\n\n");
       exit (1); 
    }

    // check for keyword-entries
    control=sscanf(line,"%s", keyword);

    i=0;
    while (keyword[i])
    {
      keyword[i] = tolower(keyword[i]);
      ++i;
    }
    keyword[i]='\0';

    // assign stepcounter
    if (strcmp(keyword,"n") == 0) 
    {
      control=sscanf(line, "%s   %d   %d", keyword, &nx, &ny);
      continue;
    }
    else if (strcmp(keyword, "zero") == 0)
    {
      control=sscanf(line,"%s   %d   %d", keyword, &zero_x, &zero_y);
      continue;
    }
    else if (strcmp(keyword, "dr") == 0)
    {
      control=sscanf(line,"%s   %lf   %lf", keyword, &dr_x, &dr_y);
      continue;
    } 
    else if (strcmp(keyword, "red_mass") == 0)
    {
      control=sscanf(line,"%s   %lf   %lf", keyword, &dr_x, &dr_y);
      continue;
    } 

    // read entries
    count = 0;
    while (sscanf(line,"%lf", &data[count][n_points]) != -1)
    { 
      // eliminate recent number from line
      position = 0;
      while( (line[position] != ' ') && (line[position] != '\t') )
	position ++;

      j = 0;
      for (i = position; i < strlen(line); i++)
	{
	  line[j] = line[i];
	  j++;
	}
      line[j] = '\0';
    
      // remove additional leading whitespace characters
      i=0;
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

      count ++;

    }
    n_points ++;

  } // end second read to store data

  fclose(infile);

  // get minimum energy
  for (j = 0; j < n_points; j++)
  {

    if (e_min > data[n_entries-1][j])
    {
      e_min = data[n_entries-1][j];
      n_min = j;
    }
  }

  // careful: nx requires zero line of y, ny required zero line of x
  nx_min = n_min%ny;
  ny_min = n_min/ny;

  //extract x-potential
  pot_x = (double *) calloc(nx, sizeof(double) );

  for (i = 0; i < nx; i ++)
  {
    idx_x = nx_min + i * ny;

    pot_x[i] = data[n_entries-1][idx_x];

    // printf("\n%d   %12.8lf    %12.8lf    %12.8lf", i, data[0][idx_x], data[1][idx_x], pot_x[i]);
  }

  printf("\n");

  //extract y-potential
  pot_y = (double *) calloc(ny, sizeof(double) );

  for (i = 0; i < ny; i ++)
  {
    idx_y = ny*ny_min  + i;

    pot_y[i] = data[n_entries-1][idx_y];

    //printf("\n%d   %12.8lf    %12.8lf   %12.8lf", i, data[0][idx_y], data[1][idx_y], pot_y[i]);
  }
  printf("\n");

  // compute pair correlation potential
  for (i = 0; i < n_points; i++)
  {

    // careful: order now reversed
      idx_y = i%ny;
      idx_x = i/ny;

      //printf("\n%d  %d  %d", i,idx_x, idx_y);

      data[n_entries][i] = data[n_entries-1][i] - pot_x[idx_x] - pot_y[idx_y];
  
      // printf("\n%12.8lf  %12.8lf   %12.8lf   %12.8lf   %12.8lf", data[0][i], data[1][i], data[n_entries-1][i], data[n_entries][i], pot_x[idx_x] + pot_y[idx_y] );
  }


  // normalise wavefunctions
  for (i = 2; i <= n_entries - 2; i++)
  {

    // generate square of wavefunction
    for (k = 0; k < n_points; k++)
    {
      integrand[k] = data[i][k] * data[i][k];
    }

    integral = integral_2d(integrand, nx, ny, dr_x, dr_y);
    
    // norm is sqrt of integral
    integral = sqrt(integral);
    
    
    // devide by norm
    for (k = 0; k < n_points; k++)
        data[i][k] = data[i][k]/integral;
    
  }


  printf("\n#\n#\n#");

  printf(" --");
  for (i=2; i <= n_entries - 2; i++)
    printf("----------");
  printf("\n#");

  printf("   ");
  for (i=2; i <= n_entries - 2; i++)
    printf("   ");

  printf("--==>> Orthonormality <<==--\n");


  printf("# --");
  for (i=2; i <= n_entries - 2; i++)
    printf("----------");
  printf("\n#");

  printf("   ");

  for (i=2; i <= n_entries - 2; i++)
    printf("%7d   ",i-2);     

  for (i = 2; i <= n_entries - 2; i++)
  {
    printf("\n#%3d",i-2);

      printf("");
      for (j = n_entries - i + 1; j < n_entries-1; j++)
      {
	printf("          ");
      }

      for (j = i; j <= n_entries - 2; j++)
      {
	for (k = 0; k < n_points; k++)
	{
	  // generate integrand
	  integrand[k] = data[i][k] * data[j][k];
	}
	integral = integral_2d(integrand, nx, ny, dr_x, dr_y);

	printf("  %8.5lf", integral);
     }
  }

  printf("\n# --");

for (i=2; i <= n_entries - 2; i++)
  printf("----------");
  printf("\n#");


  printf("\n#\n#");

  printf(" --");

  for (i=2; i <= n_entries - 2; i++)
    printf("----------");
  printf("\n#");


  printf("   ");
  for (i=2; i <= n_entries - 2; i++)
    printf("   ");
  printf("--==>> Potential V <<==--\n");


  printf("# --");
  for (i=2; i <= n_entries - 2; i++)
    printf("----------");
  printf("\n#");

  printf("   ");
  for (i=2; i <= n_entries - 2; i++)
    printf("%7d   ",i-2);     

  for (i = 2; i <= n_entries - 2; i++)
  {
     printf("\n#%3d",i-2);

    printf("");
      for (j = n_entries - i + 1; j < n_entries-1; j++)
      {
	printf("          ");
      }


      for (j = i; j <= n_entries - 2; j++)
      {
	for (k = 0; k < n_points; k++)
	{
	  // generate integrand
	  integrand[k] = data[i][k] * data[j][k] * data[n_entries - 1][k];
	}
	integral = integral_2d(integrand, nx, ny, dr_x, dr_y);

	printf("  %8.5lf", integral);
     }
  }

  printf("\n# --");

for (i=2; i <= n_entries - 2; i++)
  printf("--------------");
  printf("\n#");



for (i=2; i <= n_entries - 2; i++)
  printf("     ");

  printf("--==>> Coupling Potential V <<==--");

  printf("\n# --");

for (i=2; i <= n_entries - 2; i++)
  printf("--------------");
  printf("\n");

  printf("#");
  for (i=2; i <= n_entries - 2; i++)
    printf("%11d   ",i-2);     

  for (i = 2; i <= n_entries - 2; i++)
  {
     printf("\n#%3d",i-2);

      printf("");
      for (j = n_entries - i + 1; j < n_entries-1; j++)
      {
	printf("              ");
      }


      for (j = i; j <= n_entries - 2; j++)
      {
	for (k = 0; k < n_points; k++)
	{
	  // generate integrand
	  integrand[k] = data[i][k] * data[j][k] * data[n_entries][k];
	}
	integral = integral_2d(integrand, nx, ny, dr_x, dr_y);

	printf("  %12.5e", integral);
     }
  }

  printf("\n# --");

for (i=2; i <= n_entries - 2; i++)
  printf("--------------");
  printf("\n");

  printf("#\n#\n#\n");

  for (i=0; i < n_points; i++)
  {
      printf("\n%20.14lf  %20.14lf  %20.14lf  %20.14lf  ", data[0][i], data[1][i], data[n_entries-1][i], data[n_entries][i]);

      idx_y = i%ny;

      if ( (i+1)%ny == 0)
	printf("\n");
  }




  return(1);
}



double integral_2d(double * integrand, int nx, int ny, double dx, double dy)
{
  int i;
  int idx_x, idx_y;

  double weight;

  double integral = 0.0;

  for (i = 0; i < nx *ny; i++)
  {
      idx_x = i/ny;
      idx_y = i%ny;

      if ( (idx_x == 0) || (idx_x == nx-1) ) 
      {
	if ( (idx_y == 0) || (idx_y == ny-1) ) 
	  {
	  weight = 0.25;
	  }
	else
	  weight = 0.5;
      }
      else if ( (idx_y == 0) || (idx_y == ny-1) ) 
	{
	  weight = 0.5;
	}
      else
	weight = 1.0;

      integral = integral + weight * integrand[i];

  }

  return dx*dy*integral;
}
