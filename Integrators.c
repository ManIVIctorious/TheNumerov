
#include <stdio.h>
#include <stdlib.h>

// provided prototypes
double integrate_1d(int n, double dx, double* integrand);
double integrate_2d(int nq1, int nq2, double dx, double* integrand);
double integrate_3d(int nq1, int nq2, int nq3, double dx, double* integrand);
double integrate_4d(int nq1, int nq2, int nq3, int nq4, double dx, double* integrand);

double Integrate(int dimension, int* nq, double dx, double* integrand);


double Integrate(int dimension, int* nq, double dx, double* integrand){

    double integral = 0.0;

    switch(dimension){

        case 1:
            integral = integrate_1d(nq[0], dx, integrand);
            break;

        case 2:
            integral = integrate_2d(nq[0], nq[1], dx, integrand);
            break;

        case 3:
            integral = integrate_3d(nq[0], nq[1], nq[2], dx, integrand);
            break;

        case 4:
            integral = integrate_4d(nq[0], nq[1], nq[2], nq[3], dx, integrand);
            break;

        default:
            fprintf(stderr,
                "\n (-) Error: The requested integration routine is not implemented"
                "\n     at the moment, Aborting...\n\n"
            );
            exit(EXIT_FAILURE);
    }

    return integral;
}


double integrate_1d(int n, double dx, double *integrand){

    int i;
    double integral = 0.0;

// trapezoidal rule for an equispaced grid with grid spacing h:
//  int_a^b f dx = h * ( f(a)/2 + sum_{i = 1}^{n-1} f(a + i*h) + f(b)/2 )

    integral = 0.5 * (integrand[0] + integrand[n-1]);

    for(i = 1; i < (n-1); ++i){
        integral += integrand[i];
    }

    return (integral * dx);
}




double integrate_2d(int nq1, int nq2, double dx, double *integrand){

    int i, j, index;
    double integral = 0.0;

// Analogous to the 1D trapezoidal rule the inner part
//  (without edges and corners) gets summed up completely
    for(i = 1; i < (nq1 - 1); ++i){
        for(j = 1; j < (nq2 - 1); ++j){

            index = i*nq2 + j;
            integral += integrand[index];

        }
    }

// while the edges (without corners) get summed up half as often
// left and right edges
    for(i = 1; i < (nq1 - 1); ++i){

        integral += integrand[i*nq2]/2.0;             // left edge
        integral += integrand[i*nq2 + (nq2-1)]/2.0;   // right edge

    }

// upper and lower edges
    for(j = 1; j < (nq2 - 1); ++j){

        integral += integrand[j]/2.0;                 // upper edge
        integral += integrand[j + (nq1-1)*nq2]/2.0;   // lower edge

    }

// and the corners are even halve of that
    integral += integrand[   0   ]/4.0;     // top, left
    integral += integrand[(nq2-1)]/4.0;     // top, right
    integral += integrand[(nq1-1)*nq2]/4.0; // bottom, left
    integral += integrand[nq1*nq2 - 1]/4.0; // bottom, right

    return (integral * dx * dx);
}





double integrate_3d(int nq1, int nq2, int nq3, double dx, double *integrand){

    int i, j, k, index;
    double integral = 0.0;

//--------------------------------------------------
// inner part without faces, edges and corners
// factor: 1.0
    for(i = 1; i < (nq1 - 1); ++i){
        for(j = 1; j < (nq2 - 1); ++j){
            for(k = 1; k < (nq3 - 1); ++k){

                index = i*nq2*nq3 + j*nq3 + k;
                integral += integrand[index];

            }
        }
    }

//--------------------------------------------------
// six side faces, without edges and corners
//  factor: 1/2
//  front and back face
    for(j = 1; j < (nq2 - 1); ++j){
        for(k = 1; k < (nq3 - 1); ++k){

            index = j*nq3 + k;

            integral += integrand[index]/2.0; // front face

            integral += integrand[index + (nq1-1)*nq2*nq3]/2.0; // back face

        }// endfor k
    }// endfor j

//  top and bottom face
    for(i = 1; i < (nq1 - 1); ++i){
        for(k = 1; k < (nq3 - 1); ++k){

            index = i*nq2*nq3 + k;

            integral += integrand[index]/2.0; // top face

            integral += integrand[index + (nq2-1)*nq3]/2.0; // bottom face

        }// endfor k
    }// endfor i

//  left and right side face
    for(i = 1; i < (nq1 - 1); ++i){
        for(j = 1; j < (nq2 - 1); ++j){

            index = i*nq2*nq3 + j*nq3;

            integral += integrand[index]/2.0;

            integral += integrand[index + (nq3-1)]/2.0;

        }// endfor k
    }// endfor j

//--------------------------------------------------
// twelve edges, without corners
//  factor: 1/4
//  the four left to right edges
    for(k = 1; k < (nq3 - 1); ++k){

        index = k;
        integral += integrand[index]/4.0; // front, top

        index = k + (nq2-1)*nq3;
        integral += integrand[index]/4.0; // front, bottom

        index = k + (nq1-1)*nq2*nq3;
        integral += integrand[index]/4.0; // back, top

        index = k + (nq1-1)*nq2*nq3 + (nq2-1)*nq3;
        integral += integrand[index]/4.0; // back, bottom

    }

//  the four top to bottom edges
    for(j = 1; j < (nq2 - 1); ++j){

        index = j*nq3;
        integral += integrand[index]/4.0; // front, left

        index = j*nq3 + (nq3-1);
        integral += integrand[index]/4.0; // front, right

        index = j*nq3 + (nq1-1)*nq2*nq3;
        integral += integrand[index]/4.0; // back, left

        index = j*nq3 + (nq1-1)*nq2*nq3 + (nq3-1);
        integral += integrand[index]/4.0; // back, right

    }

//  the four front to back facing edges
    for(i = 1; i < (nq1 - 1); ++i){

        index = i*nq2*nq3;
        integral += integrand[index]/4.0; // top, left

        index = i*nq2*nq3 + (nq3-1);
        integral += integrand[index]/4.0; // top, right

        index = i*nq2*nq3 + (nq2-1)*nq3;
        integral += integrand[index]/4.0; // bottom, left

        index = i*nq2*nq3 + (nq2-1)*nq3 + (nq3-1);
        integral += integrand[index]/4.0; // bottom, right

    }

//--------------------------------------------------
// eight corners
//  factor: 1/8

// front
// top, left
    integral += integrand[   0   ]/8.0;
// top, right
    integral += integrand[(nq3-1)]/8.0;
// bottom, left
    integral += integrand[(nq2-1)*nq3]/8.0;
// bottom, right
    integral += integrand[(nq2-1)*nq3 + (nq3-1)]/8.0;

// back
// top, left
    integral += integrand[(nq1-1)*nq2*nq3]/8.0;
// top, right
    integral += integrand[(nq1-1)*nq2*nq3 + (nq3-1)]/8.0;
// bottom, left
    integral += integrand[(nq1-1)*nq2*nq3 + (nq2-1)*nq3]/8.0;
// bottom, right
    integral += integrand[nq1*nq2*nq3 - 1]/8.0;

    return (integral * dx * dx * dx);
}





double integrate_4d(int nq1, int nq2, int nq3, int nq4, double dx, double integrand[])
{
  int i, j, k, l, index;

  double integral = 0.0;
//--------------------------------------------------
//inner part. without edges and boundary areas.
// factor: 1.0
for (i = 1; i < nq1-1; i++)
   {
     for (j = 1; j < nq2-1; j++)
     {
       for ( k=1; k < nq3-1; k++)
       {
        for ( l=1; l < nq4-1; l++)
         {
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index];
         }// endfor l
       }// endfor k
     }// endfor j
   }// endfor i
//--------------------------------------------------
// twelve hypercubes wit factor 8/16
// factor: 4/8
         l=nq4-1;
for (i = 1; i < nq1-1; i++)
   {
     for (j = 1; j < nq2-1; j++)
     {
       for ( k=1; k < nq3-1; k++)
       {

         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*0.5;

         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*0.5;

      } // endfor k
     } // endfor j
    } // endfor i
k=nq3-1;
for (i = 1; i < nq1-1; i++)
   {
     for (j = 1; j < nq2-1; j++)
     {
        for ( l=1; l < nq4-1; l++)
         {
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + l;
         integral = integral + integrand[index]*0.5;

         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*0.5;
       } // endfor l
     } // endfor j
    } // endfor i


j=nq2-1;
for (i = 1; i < nq1-1; i++)
   {
     for ( k=1; k < nq3-1; k++)
      {
        for ( l=1; l < nq4-1; l++)
         {
         index = i*nq2*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*0.5;
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*0.5;
       } // endfor l
      } // endfor k
    } // endfor i

i=nq1-1;
     for (j = 1; j < nq2-1; j++)
     {
       for ( k=1; k < nq3-1; k++)
       {
        for ( l=1; l < nq4-1; l++)
         {
         index = j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*0.5;
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*0.5;
       } // endfor l
      } // endfor k
     } // endfor j

// done
//--------------------------------------------------
// side "walls"
// factor: 4/16
/*
for (i = 1; i < nq1-1; i++)
   {
     for (j = 1; j < nq2-1; j++)
     {
       for ( k=1; k < nq3-1; k++)
       {
        for ( l=1; l < nq4-1; l++)
         {
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*4/16;
         }// endfor l
       }// endfor k
     }// endfor j
   }// endfor i
*/
k=nq3-1;
l=nq4-1;
for (i = 1; i < nq1-1; i++)
   {
     for (j = 1; j < nq2-1; j++)
     {
         //k=l=0
         index = i*nq2*nq3*nq4 + j*nq3*nq4;
         integral = integral + integrand[index]*4/16;
         //k=0,l=nq-1
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + l;
         integral = integral + integrand[index]*4/16;
         //k=nq-1,l=0
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*4/16;
         //k=nq-1,l=nq-1
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*4/16;
     }// endfor j
   }// endfor i

j=nq2-1;
l=nq4-1;
for (i = 1; i < nq1-1; i++)
   {
       for ( k=1; k < nq3-1; k++)
       {
         //j=l=0;
         index = i*nq2*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*4/16;
         //j=0,l=nq-1
         index = i*nq2*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*4/16;
         //j=nq-1,l=0
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*4/16;
         //j&l=nq-1
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*4/16;
       }// endfor k
   }// endfor i

//j=nq2-1;
k=nq3-1;
for (i = 1; i < nq1-1; i++)
   {
        for ( l=1; l < nq4-1; l++)
         {
         //j=k=0;
         index = i*nq2*nq3*nq4 + l;
         integral = integral + integrand[index]*4/16;
         //j=0,k=nq-1
         index = i*nq2*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*4/16;
         //j=nq-1,k=0
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + l;
         integral = integral + integrand[index]*4/16;
         //j&k=nq-1
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*4/16;
         }// endfor l
   }// endfor i

i=nq1-1;
l=nq4-1;
     for (j = 1; j < nq2-1; j++)
     {
       for ( k=1; k < nq3-1; k++)
       {
         //i=l=0;
         index = j*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*4/16;
         //i=0,l=nq-1
         index = j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*4/16;
         //i=nq-1,l=0
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*4/16;
         //i&l=nq-1
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*4/16;
       }// endfor k
     }// endfor j

//i=nq1-1;
k=nq3-1;
     for (j = 1; j < nq2-1; j++)
     {
        for ( l=1; l < nq4-1; l++)
         {
         //i=k=0;
         index = j*nq3*nq4 + l;
         integral = integral + integrand[index]*4/16;
         //i=0,k=nq-1
         index = j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*4/16;
         //i=nq-1,k=0
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + l;
         integral = integral + integrand[index]*4/16;
         //i&k=nq-1
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*4/16;
         }// endfor l
     }// endfor j

//i=nq1-1;
j=nq2-1;
       for ( k=1; k < nq3-1; k++)
       {
        for ( l=1; l < nq4-1; l++)
         {
         //i=j=0;
         index = k*nq4 + l;
         integral = integral + integrand[index]*4/16;
         //i=0,j=nq-1
         index = j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*4/16;
         //i=nq-1,j=0
         index = i*nq2*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*4/16;
         //i&j=nq-1
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*4/16;
         }// endfor l
       }// endfor k
//done

//--------------------------------------------------
// side "lines"
// factor: 2/16
/*
for (i = 1; i < nq1-1; i++)
   {
     for (j = 1; j < nq2-1; j++)
     {
       for ( k=1; k < nq3-1; k++)
       {
        for ( l=1; l < nq4-1; l++)
         {
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*2/16;
         }// endfor l
       }// endfor k
     }// endfor j
   }// endfor i
*/
j=nq2-1;
k=nq3-1;
l=nq4-1;
for (i = 1; i < nq1-1; i++)
   {
         //j,k,l 0
         index = i*nq2*nq3*nq4;
         integral = integral + integrand[index]*2/16;
         //j,k 0 | l nq
         index = i*nq2*nq3*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //j,l 0 | k nq
         index = i*nq2*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*2/16;
         //j 0 | k,l nq
         index = i*nq2*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //l,k 0 | j nq
         index = i*nq2*nq3*nq4 + j*nq3*nq4;
         integral = integral + integrand[index]*2/16;
         //k 0 | j l nq
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //l 0 | j,k nq
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*2/16;
         //0 | j,k,l nq
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*2/16;

   }// endfor i

i=nq1-1;
k=nq3-1;
l=nq4-1;
for (j = 1; j < nq2-1; j++)
   {
         //i,k,l 0
         index = j*nq3*nq4;
         integral = integral + integrand[index]*2/16;
         //i,k 0 | l nq
         index = j*nq3*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //i,l 0 | k nq
         index = j*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*2/16;
         //i 0 | k,l nq
         index = j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //l,k 0 | i nq
         index = i*nq2*nq3*nq4 + j*nq3*nq4;
         integral = integral + integrand[index]*2/16;
         //k 0 | i l nq
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //l 0 | i,k nq
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*2/16;
         //0 | i,k,l nq
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*2/16;

   }// endfor j

i=nq1-1;
j=nq2-1;
l=nq4-1;
for (k = 1; k < nq3-1; k++)
   {
         //i,j,l 0
         index = k*nq4;
         integral = integral + integrand[index]*2/16;
         //i,j 0 | l nq
         index = k*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //i,l 0 | j nq
         index = j*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*2/16;
         //i 0 | j,l nq
         index = j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //l,j 0 | i nq
         index = i*nq2*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*2/16;
         //j 0 | i l nq
         index = i*nq2*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //l 0 | i,j nq
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*2/16;
         //0 | i,j,l nq
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*2/16;

   }// endfor k

i=nq1-1;
j=nq2-1;
k=nq3-1;
for (l = 1; l < nq4-1; l++)
   {
         //i,j,k 0
         index = l;
         integral = integral + integrand[index]*2/16;
         //i,j 0 | k nq
         index = k*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //i,k 0 | j nq
         index = j*nq3*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //i 0 | j,k nq
         index = j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //k,j 0 | i nq
         index = i*nq2*nq3*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //j 0 | i k nq
         index = i*nq2*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //k 0 | i,j nq
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + l;
         integral = integral + integrand[index]*2/16;
         //0 | i,j,k nq
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*2/16;

   }// endfor l
//--------------------------------------------------
// sixteen edges
//factor 1/16
         i=nq1-1;
         j=nq2-1;
         k=nq3-1;
         l=nq4-1;
         integral = integral + integrand[  0  ]*1/16;
         index = l;
         integral = integral + integrand[index]*1/16;
         index = k*nq4;
         integral = integral + integrand[index]*1/16;
         index = k*nq4 + l;
         integral = integral + integrand[index]*1/16;
         index = j*nq3*nq4;
         integral = integral + integrand[index]*1/16;
         index = j*nq3*nq4 + l;
         integral = integral + integrand[index]*1/16;
         index = j*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*1/16;
         index = j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*1/16;
         index = i*nq2*nq3*nq4 ;
         integral = integral + integrand[index]*1/16;
         index = i*nq2*nq3*nq4 + l;
         integral = integral + integrand[index]*1/16;
         index = i*nq2*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*1/16;
         index = i*nq2*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*1/16;
         index = i*nq2*nq3*nq4 + j*nq3*nq4;
         integral = integral + integrand[index]*1/16;
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + l;
         integral = integral + integrand[index]*1/16;
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4;
         integral = integral + integrand[index]*1/16;
         index = i*nq2*nq3*nq4 + j*nq3*nq4 + k*nq4 + l;
         integral = integral + integrand[index]*1/16;
  return integral * dx * dx * dx * dx; //
}
