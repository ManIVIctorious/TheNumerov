
// provided prototypes
double Integrate4D(int* nq, double dx, double* integrand);

double Integrate4D(int* nq, double dx, double* integrand){

// shift inidces to reach the end of a given direction
    int ishift =  nq[3]  *  nq[2]  *  nq[1]  * (nq[0]-1);
    int jshift =  nq[3]  *  nq[2]  * (nq[1]-1);
    int kshift =  nq[3]  * (nq[2]-1);
    int lshift = (nq[3]-1);

    double integral = 0.0;

//--------------------------------------------------
//inner part without edges and boundary areas
// factor: 1.0
    for(int i = 1; i < (nq[0]-1); ++i){
        for(int j = 1; j < (nq[1]-1); ++j){
            for(int k = 1; k < (nq[2]-1); ++k){
                for(int l = 1; l < (nq[3]-1); ++l){

                    int index = ( (i*nq[1] + j)*nq[2] + k )*nq[3] + l;
                    integral += integrand[index];
                }
            }
        }
    }

//--------------------------------------------------
// twelve hypercubes with factor 8/16
// factor: 1/2
// l = 0; l = lshift
    for(int i = 1; i < (nq[0]-1); ++i){
        for(int j = 1; j < (nq[1]-1); ++j){
            for(int k = 1; k < (nq[2]-1); ++k){

                int index = ((i*nq[1] + j)*nq[2] + k)*nq[3];

                integral += integrand[index         ]*0.5;
                integral += integrand[index + lshift]*0.5;
            }
        }
    }

// k = 0; k = kshift
    for(int i = 1; i < (nq[0]-1); ++i){
        for(int j = 1; j < (nq[1]-1); ++j){
            for(int l = 1; l < (nq[3]-1); ++l){

                int index = (i*nq[1] + j)*nq[2]*nq[3] + l;

                integral += integrand[index         ]*0.5;
                integral += integrand[index + kshift]*0.5;
            }
        }
    }

// j = 0; j = jshift
    for(int i = 1; i < (nq[0]-1); ++i){
        for(int k = 1; k < (nq[2]-1); ++k){
            for(int l = 1; l < (nq[3]-1); ++l){

                int index = (i*nq[1]*nq[2] + k)*nq[3] + l;

                integral += integrand[index         ]*0.5;
                integral += integrand[index + jshift]*0.5;
            }
        }
    }

// i = 0; i = ishift
     for(int j = 1; j < (nq[1]-1); ++j){
        for(int k = 1; k < (nq[2]-1); ++k){
            for(int l = 1; l < (nq[3]-1); ++l){

                int index = (j*nq[2] + k)*nq[3] + l;

                integral += integrand[index         ]*0.5;
                integral += integrand[index + ishift]*0.5;
            }
        }
     }


//--------------------------------------------------
// side "walls" with factor 4/16
// factor: 1/4

// all permutations of
// k = 0; k = kshift;
// l = 0; l = lshift;
    for(int i = 1; i < (nq[0]-1); ++i){
        for(int j = 1; j < (nq[1]-1); ++j){

            int index = (i*nq[1] + j)*nq[2]*nq[3];

            integral += integrand[index                  ]*0.25;
            integral += integrand[index + kshift         ]*0.25;
            integral += integrand[index + lshift         ]*0.25;
            integral += integrand[index + kshift + lshift]*0.25;
        }
    }

// all permutations of
// j = 0; j = jshift;
// l = 0; l = lshift;
    for(int i = 1; i < (nq[0]-1); ++i){
        for(int k = 1; k < (nq[2]-1); ++k){

            int index = (i*nq[1]*nq[2] + k)*nq[3];

            integral += integrand[index                  ]*0.25;
            integral += integrand[index + jshift         ]*0.25;
            integral += integrand[index + lshift         ]*0.25;
            integral += integrand[index + jshift + lshift]*0.25;
        }
    }

// all permutations of
// j = 0; j = jshift;
// k = 0; k = kshift;
    for(int i = 1; i < (nq[0]-1); ++i){
        for(int l = 1; l < (nq[3]-1); ++l){

            int index = i*nq[1]*nq[2]*nq[3] + l;

            integral += integrand[index                  ]*0.25;
            integral += integrand[index + jshift         ]*0.25;
            integral += integrand[index + kshift         ]*0.25;
            integral += integrand[index + jshift + kshift]*0.25;
        }
    }

// all permutations of
// i = 0; i = ishift;
// l = 0; l = lshift;
    for(int j = 1; j < (nq[1]-1); ++j){
        for(int k = 1; k < (nq[2]-1); ++k){

            int index = (j*nq[2] + k)*nq[3];

            integral += integrand[index                  ]*0.25;
            integral += integrand[index + ishift         ]*0.25;
            integral += integrand[index + lshift         ]*0.25;
            integral += integrand[index + ishift + lshift]*0.25;
        }
    }

// all permutations of
// i = 0; i = ishift;
// k = 0; k = kshift;
    for(int j = 1; j < (nq[1]-1); ++j){
        for(int l = 1; l < (nq[3]-1); ++l){

            int index = j*nq[2]*nq[3] + l;

            integral += integrand[index                  ]*0.25;
            integral += integrand[index + ishift         ]*0.25;
            integral += integrand[index + kshift         ]*0.25;
            integral += integrand[index + ishift + kshift]*0.25;
        }
    }

// all permutations of
// i = 0; i = ishift;
// j = 0; j = jshift;
    for(int k = 1; k < (nq[2]-1); ++k){
        for(int l = 1; l < (nq[3]-1); ++l){

            int index = k*nq[3] + l;

            integral += integrand[index                  ]*0.25;
            integral += integrand[index + ishift         ]*0.25;
            integral += integrand[index + jshift         ]*0.25;
            integral += integrand[index + ishift + jshift]*0.25;
        }
    }

//--------------------------------------------------
// side "lines" with a factor of 2/16
// factor: 1/8

// j = 0; j = jshift;
// k = 0; k = kshift;
// l = 0; l = lshift;
    for(int i = 1; i < (nq[0]-1); ++i){

        int index = i*nq[1]*nq[2]*nq[3];

        integral += integrand[index                           ]*0.125;
        integral += integrand[index + jshift                  ]*0.125;
        integral += integrand[index + kshift                  ]*0.125;
        integral += integrand[index + lshift                  ]*0.125;
        integral += integrand[index + jshift + kshift         ]*0.125;
        integral += integrand[index + jshift + lshift         ]*0.125;
        integral += integrand[index + kshift + lshift         ]*0.125;
        integral += integrand[index + jshift + kshift + lshift]*0.125;
   }

// i = 0; i = ishift;
// k = 0; k = kshift;
// l = 0; l = lshift;
    for(int j = 1; j < (nq[1]-1); ++j){

        int index = j*nq[2]*nq[3];

        integral += integrand[index                           ]*0.125;
        integral += integrand[index + ishift                  ]*0.125;
        integral += integrand[index + kshift                  ]*0.125;
        integral += integrand[index + lshift                  ]*0.125;
        integral += integrand[index + ishift + kshift         ]*0.125;
        integral += integrand[index + ishift + lshift         ]*0.125;
        integral += integrand[index + kshift + lshift         ]*0.125;
        integral += integrand[index + ishift + kshift + lshift]*0.125;
    }

// i = 0; i = ishift;
// j = 0; j = jshift;
// l = 0; l = lshift;
    for(int k = 1; k < (nq[2]-1); ++k){

        int index = k*nq[3];

        integral += integrand[index                           ]*0.125;
        integral += integrand[index + ishift                  ]*0.125;
        integral += integrand[index + jshift                  ]*0.125;
        integral += integrand[index + lshift                  ]*0.125;
        integral += integrand[index + ishift + jshift         ]*0.125;
        integral += integrand[index + ishift + lshift         ]*0.125;
        integral += integrand[index + jshift + lshift         ]*0.125;
        integral += integrand[index + ishift + jshift + lshift]*0.125;
    }

// i = 0; i = ishift;
// j = 0; j = jshift;
// k = 0; k = kshift;
    for(int l = 1; l < (nq[3]-1); ++l){

        int index = l;

        integral += integrand[index                           ]*0.125;
        integral += integrand[index + ishift                  ]*0.125;
        integral += integrand[index + jshift                  ]*0.125;
        integral += integrand[index + kshift                  ]*0.125;
        integral += integrand[index + ishift + jshift         ]*0.125;
        integral += integrand[index + ishift + kshift         ]*0.125;
        integral += integrand[index + jshift + kshift         ]*0.125;
        integral += integrand[index + ishift + jshift + kshift]*0.125;
    }


//--------------------------------------------------
// sixteen "corners"
// factor: 1/16
    integral += integrand[   0  ]*0.0625;

    integral += integrand[ishift]*0.0625;
    integral += integrand[jshift]*0.0625;
    integral += integrand[kshift]*0.0625;
    integral += integrand[lshift]*0.0625;

    integral += integrand[ishift + jshift]*0.0625;
    integral += integrand[ishift + kshift]*0.0625;
    integral += integrand[ishift + lshift]*0.0625;
    integral += integrand[jshift + kshift]*0.0625;
    integral += integrand[jshift + lshift]*0.0625;
    integral += integrand[kshift + lshift]*0.0625;

    integral += integrand[ishift + jshift + kshift]*0.0625;
    integral += integrand[ishift + jshift + lshift]*0.0625;
    integral += integrand[ishift + kshift + lshift]*0.0625;
    integral += integrand[jshift + kshift + lshift]*0.0625;

    integral += integrand[ishift + jshift + kshift + lshift]*0.0625;

    return integral * dx * dx * dx * dx;
}
