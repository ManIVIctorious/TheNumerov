
// provided prototypes
double integrate_1d(int n, double dx, double integrand[]);
double integrate_2d(int nq1, int nq2, double dx, double integrand[]);
double integrate_3d(int nq1, int nq2, int nq3, double dx, double integrand[]);
double integrate_4d(int nq1, int nq2, int nq3, int nq4, double dx, double integrand[]);




double integrate_1d(int n, double dx, double integrand[])
{
  int i;

  double integral = 0.0;

  integral = 0.5*(integrand[0] + integrand[n-1]);

  for (i = 1; i < (n-1); i++)
  {

    integral = integral + integrand[i];
  }

  return integral * dx;
}




double integrate_2d(int nq1, int nq2, double dx, double integrand[])
{
  int i, j,index;

  double integral = 0.0;

// inner part, without edges and "RAND"
   for (i = 1; i < nq1-1; i++)
   {
     for (j = 1; j < nq2-1; j++)
     {
       index = i*nq2 + j;
       integral = integral + integrand[index];
     }
   }

//randpunkte links rechts ohne ecken
   for ( j = 1; j < nq1-1; j++)
   {
     index = j * nq2;
     integral = integral + 1/2*(integrand[index]+integrand[index + nq2 -1]);
   }
// randpunkte oben unten
   for ( i = 1; i < nq2-1; i++)
   {
     index = (nq1-1)*nq2+i;
     integral = integral + 1/2*(integrand[i]+integrand[index]);
   }


   integral = integral + 1/4*(integrand[0]+integrand[nq1-1]+integrand[nq1*nq2-1]+ integrand[(nq1-1)*nq2]);

  return integral * dx * dx; //
}





double integrate_3d(int nq1, int nq2, int nq3, double dx, double integrand[])
{
  int i, j, k,index;

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
         index = i*nq2*nq3 + j*nq3 + k;
         integral = integral + integrand[index];
       }
     }
   }
//--------------------------------------------------
// six side faces
// factor: 4/8
for (j = 1; j < nq2-1; j++)
{
     for (k = 1; k < nq3-1; k++)
     {
       index = j*nq3 + k;
       integral = integral + integrand[index]*1/2;
       index = (nq1-1)*nq2*nq3 + j*nq3 + k;
       integral = integral + integrand[index]*1/2;
     }// endfor k
} // endfor j

for (i = 1; i < nq1-1; i++)
{
     for (k = 1; k < nq3-1; k++)
     {
       index = i*nq2*nq3 + k;
       integral = integral + integrand[index]*1/2;
       index = i*nq2*nq3 + (nq2-1)*nq3 + k;
       integral = integral + integrand[index]*1/2;
     }// endfor k
} // endfor i

for (i = 1; i < nq1-1; i++)
{
     for (j = 1; j < nq2-1; j++)
     {
       index = i*nq2*nq3 + j*nq3;
       integral = integral + integrand[index]*1/2;
       index = i*nq2*nq3 + j*nq3 + (nq3-1);
       integral = integral + integrand[index]*1/2;
     }// endfor k
} // endfor j
//--------------------------------------------------
// twelve side edges
// factor: 2/8
for ( k = 1; k< nq3-1; k++)
{
  index = k;
  integral = integral + integrand[index]*2/8;
  index = (nq2-1)*nq3 + k;
  integral = integral + integrand[index]*2/8;
  index = (nq1-1)*nq2*nq3 + k;
  integral = integral + integrand[index]*2/8;
  index = (nq1-1)*nq2*nq3 + ( nq2-1)*nq3 + k;
  integral = integral + integrand[index]*2/8;
}

for ( j = 1; j< nq2-1; j++)
{
  index = j*nq3;
  integral = integral + integrand[index]*2/8;
  index = j*nq3 + (nq3-1);
  integral = integral + integrand[index]*2/8;
  index = (nq1-1)*nq2*nq3+j*nq3;
  integral = integral + integrand[index]*2/8;
  index = (nq1-1)*nq2*nq3 + j * nq3 + (nq3 -1);
  integral = integral + integrand[index]*2/8;
}

for ( i = 1; i< nq1-1; i++)
{
  index = i*nq2*nq3;
  integral = integral + integrand[index]*2/8;
  index = i*nq2*nq3 +(nq3-1);
  integral = integral + integrand[index]*2/8;
  index = i*nq2*nq3 +(nq2-1)*nq3;
  integral = integral + integrand[index]*2/8;
  index = i*nq2*nq3  +(nq2-1)*nq3+(nq3-1);
  integral = integral + integrand[index]*2/8;
}

//--------------------------------------------------
// eigth edges
//factor 1/8
index = 0;
integral = integral + integrand[index]*1/8;
index = nq3-1;
integral = integral + integrand[index]*1/8;
index = (nq2-1)*nq3;
integral = integral + integrand[index]*1/8;
index = (nq2-1)*nq3+nq3-1;
integral = integral + integrand[index]*1/8;
index = (nq1-1) * nq2*nq3;
integral = integral + integrand[index]*1/8;
index = (nq1-1) * nq2*nq3+nq3-1;
integral = integral + integrand[index]*1/8;
index = (nq1-1) * nq2*nq3+(nq2-1)*nq3;
integral = integral + integrand[index]*1/8;
index = (nq1-1) * nq2*nq3+(nq2-1)*nq3+nq3-1;
integral = integral + integrand[index]*1/8;

  return integral * dx * dx * dx; //
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
