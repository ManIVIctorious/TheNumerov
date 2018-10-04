
#include <armadillo>
#include <stdio.h>
#include "typedefinitions.h"


extern "C" int SolverARPACK_Armadillo(settings prefs, int* nq, double* v, double ekin_param, double* stencil, double* E, double* X);
arma::sp_mat FillArmadillo(settings prefs, int* nq, int n_points, double* v, double ekin_param, double* stencil);

extern "C"{

int SolverARPACK_Armadillo(settings prefs, int* nq, double* v, double ekin_param, double* stencil, double* E, double* X){

    int i, j;
    int n_points;

// get number of points
    for(i = 0, n_points = 1; i < prefs.dimension; ++i){
        n_points *= nq[i];
    }

    arma::sp_mat A;

// depending on the dimensionality use the appropriate filling routing
    switch(prefs.dimension){

        case 2:
            A = FillArmadillo(prefs, nq, n_points, v, ekin_param, stencil);
            break;

        default:
            fprintf(stderr, 
                "\n (-) Error: The requested ARPACK Armadillo filling routine"
                "\n     for a dimensionality of \"%d\" is not implemented."
                "\n     Please check your input. Aborting...\n\n"    
                ,prefs.dimension
            );
            exit(EXIT_FAILURE);
    }

// start eigenstate calculation
    arma::mat evec;
    arma::vec eval;

int n_conv = 0;
// In recent versions of armadillo the ARPACK is directly included
//  The eigendecomposition can either be invoked by the armadillo "high-level" eigs routine
//  or by using the underlying ARPACK implementation directly
//  To use the underlying ARPACK implementation some objects/classes have to be put into scope:
    using arma::newarp::SparseGenMatProd;
    using arma::newarp::SymEigsSolver;
    using arma::newarp::EigsSelect;

// Use Armadillo ARPACK for the eigen decomposition of symmetric A
    SparseGenMatProd<double> op(A);
    SymEigsSolver< double, EigsSelect::SMALLEST_MAGN, SparseGenMatProd<double> > eigs(op, prefs.n_out, 5*prefs.n_out);

// Initialize and compute
    eigs.init();
    n_conv = eigs.compute(10000, 1.0e-14);

    if(n_conv > 0){
        eval = eigs.eigenvalues();
        evec = eigs.eigenvectors();
    }else{
        return -1;
    }

// fill eigenvalue E and eigenvector arrays X
    for(i = 0; i < n_conv; ++i){
        E[i] = eval[i];

        for(j = 0; j < n_points; ++j){
            X[i*n_points + j] = evec(j,i);
        }
    }

    return n_conv;
}
}


arma::sp_mat FillArmadillo(settings prefs, int* nq, int n_points, double* v, double ekin_param, double* stencil){

    int i, j;
    int xsh, ysh;
    int max_entries = 0;

    double stencil_entry;
    double threshold = 1.0E-7;

// determine maximum number of entries
    for(i = 0; i < nq[0]; ++i){
        for(j=0; j < nq[1]; ++j){

            for(xsh = -prefs.n_stencil/2; xsh < prefs.n_stencil/2 + 1; ++xsh){
                if( (i+xsh > -1) && (i+xsh < nq[0]) ){

                    for(ysh = -prefs.n_stencil/2; ysh < prefs.n_stencil/2 + 1; ++ysh){
                        if( (j+ysh > -1) && (j+ysh < nq[1]) ){

                            stencil_entry = stencil[(xsh + prefs.n_stencil/2)*prefs.n_stencil + (ysh + prefs.n_stencil/2)];
                            if(stencil_entry > threshold || stencil_entry < -threshold){ max_entries++; }

                        }
                    }

                }
            }

        }
    }


// determine positions and values
    arma::umat locations(2, max_entries);
    arma::vec  values = arma::zeros(max_entries);
    unsigned int filled = 0;

    for(i = 0; i < nq[0]; ++i){
        for(j = 0; j < nq[1]; ++j){

            for(xsh = -prefs.n_stencil/2; xsh < prefs.n_stencil/2 + 1; ++xsh){
                if( (i+xsh > -1) && (i+xsh < nq[0]) ){

                    for(ysh = -prefs.n_stencil/2; ysh < prefs.n_stencil/2 + 1; ++ysh){
                        if( (j+ysh > -1) && (j+ysh < nq[1]) ){
    
                            stencil_entry = stencil[(xsh+ prefs.n_stencil/2)*prefs.n_stencil + (ysh+prefs.n_stencil/2)];
    
                            if(stencil_entry > threshold || stencil_entry < -threshold){
                            //++++++++++++++++++++++++++++++++++++++++++++++++++
                                locations( 0, filled) =     i      *nq[1] +     j;
                                locations( 1, filled) = ( i + xsh )*nq[1] + ( j + ysh );
                            //++++++++++++++++++++++++++++++++++++++++++++++++++
                                values(filled) = stencil_entry * ekin_param / 2;
                            //++++++++++++++++++++++++++++++++++++++++++++++++++
                                filled++;
                            }
    
                        }
                    }
    
                }
            }
    
        }
    }

// actual sparse matrix fill
    arma::sp_mat A(locations, values, n_points, n_points, true, true);

// add potential value
    for(i = 0; i < n_points; ++i){
        A(i,i) = A(i,i) + v[i];
    }

    return A;
}
