
// provided prototypes
double integrate_4d(int nq1, int nq2, int nq3, int nq4, double dx, double integrand[]);

double integrate_4d(int nq1, int nq2, int nq3, int nq4, double dx, double integrand[]){

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
