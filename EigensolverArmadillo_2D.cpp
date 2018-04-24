
#define ARMA_64BIT_WORD

#include <armadillo>


extern "C" int EigensolverArmadillo_2D(double *v, int *nq, double ekin_param, double *stencil, int n_stencil, int n_out, double *E, double *X);

extern "C"{
    int EigensolverArmadillo_2D(double *v, int *nq, double ekin_param, double *stencil, int n_stencil, int n_out, double *E, double *X){

        int i, j;
        int xsh,ysh;
        int n_points = nq[0] * nq[1];
        int max_entries = 0;
        int n_conv = 0;

        double stencil_entry;
        double threshold = 1.0E-7;

        arma::mat evec;
        arma::vec eval;
        arma::vec values;

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


        arma::umat locations(2, max_entries);
        values = arma::zeros(max_entries);

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

        arma::sp_mat A(locations, values, n_points, n_points, true, true);

        for(i = 0; i < n_points; ++i){
            A(i,i) = A(i,i) + v[i];
        }



    // In recent versions of armadillo the ARPACK is directly included
    //  The eigendecomposition can either be invoked by the armadillo "high-level" eigs routine
    //  or by using the underlying ARPACK implementation directly
    //  To use the underlying ARPACK implementation some objects/classes have to be put into scope:
        using arma::newarp::SparseGenMatProd;
        using arma::newarp::SymEigsSolver;
        using arma::newarp::EigsSelect;

    // Use Armadillo ARPACK for the eigen decomposition of symmetric A
        SparseGenMatProd<double> op(A);
        SymEigsSolver< double, EigsSelect::SMALLEST_MAGN, SparseGenMatProd<double> > eigs(op, n_out, 5*n_out);

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
