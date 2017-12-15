
#define ARMA_64BIT_WORD

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include <armadillo>
#include <SymEigsSolver.h>  // Also includes <MatOp/DenseGenMatProd.h>
//#include <GenEigsSolver.h>  // Also includes <MatOp/DenseGenMatProd.h>
#include <MatOp/SparseGenMatProd.h>

using namespace arma; // prefix

extern "C" int EigensolverArmadillo_2D(double *v, int *nq, double ekin_param, double *stencil, int n_stencil, int n_out, double *E, double *X);

extern "C"{
    int EigensolverArmadillo_2D(double *v, int *nq, double ekin_param, double *stencil, int n_stencil, int n_out, double *E, double *X){

        int i, j, k;
        int xsh,ysh;
        int n_points = nq[0] * nq[1];
        int max_entries = 0;
        int n_conv = 0;
    
        double stencil_entry;
        double threshold = 1.0E-7;
    
        mat evec;
        vec eval;
        vec values;
    
    // determine maximum number of entries
        for(i = 0; i < nq[0]; i++){
            for(j=0; j < nq[1]; j++){
    
                for(xsh = -n_stencil/2; xsh < n_stencil/2 + 1; xsh++){
                    if( (i+xsh > -1) && (i+xsh < nq[0]) ){ 
                        for(ysh = -n_stencil/2; ysh < n_stencil/2 + 1; ysh++){
                            if( (j+ysh > -1) && (j+ysh < nq[1]) ){
    
                                stencil_entry = stencil[(xsh+n_stencil/2)*n_stencil +(ysh+n_stencil/2)];
                                if(stencil_entry > threshold || stencil_entry < -threshold){                      
                                    max_entries++;
                                }// end if >0
    
                            } // end ysh + y inbound
                        } // end ysh 
                    } // end xsh+x inbound
                } // end xsh
    
            }// end j
        }// end i 
    
    
        umat locations(2, max_entries);
        values = zeros(max_entries);
        
        unsigned int  indexz, indexs;
        unsigned int  filled = 0;
    
    // fill matrix
        for(i = 0; i < nq[0]; ++i){ 
            for(j = 0; j < nq[1]; ++j){
                indexz = i * nq[1] + j;
                for(xsh = -n_stencil/2; xsh < n_stencil/2 + 1; ++xsh){
    
                    if( (i+xsh > -1) && (i+xsh < nq[0]) ){
                        for(ysh = -n_stencil/2; ysh < n_stencil/2 + 1; ++ysh){
    
                            if( (j+ysh > -1) && (j+ysh < nq[1]) ){ 
    
                                indexs = ( i + xsh ) * nq[1] + ( j + ysh );
                                stencil_entry = stencil[(xsh+n_stencil/2)*n_stencil + (ysh+n_stencil/2)];
    
                                if(stencil_entry > threshold || stencil_entry < -threshold){      
                                //++++++++++++++++++++++++++++++++++++++++++++++++++
                                    locations ( 0, filled) = indexz;
                                    locations ( 1, filled) = indexs; 
                                //++++++++++++++++++++++++++++++++++++++++++++++++++
                                    values(filled) = stencil_entry * ekin_param / 2;
                                //++++++++++++++++++++++++++++++++++++++++++++++++++
                                    filled++;
                                }// end if >0
    
                            } // end if j + ysh
                        }// end for ysh
    
                    }// end if i+ xsh
                } // end for xsh
    
            } // end for j
        }// end for i
    
        sp_mat A(locations, values, n_points, n_points, true, true);
        
        for(i = 0; i < n_points; ++i){
            A(i,i) = A(i,i) + v[i];
        }

    
    // Use Armadillo ARPACK for the eigen decomposition of A
        SparseGenMatProd<double> op(A);
     
    // SYMMETRIC
        SymEigsSolver< double, SMALLEST_MAGN, SparseGenMatProd<double> > eigs(&op, n_out, 5*n_out);
    
    // Initialize and compute
        eigs.init();
    
        n_conv = eigs.compute(10000, 1.0e-14);
    
        if(n_conv > 0){
            eval = eigs.eigenvalues();
            evec = eigs.eigenvectors();
        }


    // fill eigenvalue E and eigenvector arrays X
        for(i = n_conv - 1, k = 0; i >= 0; --i, ++k){
            E[k] = eval[i];

            for(j = 0; j < n_points; ++j){
                X[k*n_points + j] = evec(j,i);
            }
        }

        return n_conv;
    }
}
