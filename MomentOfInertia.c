
#ifdef debug_moment_of_inertia
  #include <stdio.h>
#endif


// provided prototypes
void MomentOfInertia(int dim, int n_atoms, double tot_mass, double* mass, double* x, double* y, double* z, double** zeta, double** Q, double* I, int state_index);


void MomentOfInertia(int dim, int n_atoms, double tot_mass, double* mass, double* x, double* y, double* z, double** zeta, double** Q, double* I, int state_index){

    int i;
    double x0, y0, z0;

// make origin match system's center of mass
//--------------------------------------------------
// calculate center of mass
    x0 = y0 = z0 = 0.0;
    for(i = 0; i < n_atoms; ++i){
        x0 += x[i] * mass[i];
        y0 += y[i] * mass[i];
        z0 += z[i] * mass[i];
    }
    x0 /= tot_mass;
    y0 /= tot_mass;
    z0 /= tot_mass;

// translate system
    for(i = 0; i < n_atoms; ++i){
        x[i] -= x0;
        y[i] -= y0;
        z[i] -= z0;
    }

// determine moment of inertia tensor I
//--------------------------------------------------
// initialize tensor to zero
    for(i = 0; i < 9; ++i){ I[i] = 0.0; }

// calculate tensor values
    for(i = 0; i < n_atoms; ++i){
    // main diagonal:
        I[0] += mass[i] * (y[i]*y[i] + z[i]*z[i]);  // 11 xx
        I[4] += mass[i] * (x[i]*x[i] + z[i]*z[i]);  // 22 yy
        I[8] += mass[i] * (x[i]*x[i] + y[i]*y[i]);  // 33 zz

    // upper triangle:
        I[1] -= mass[i] * x[i] * y[i];              // 12 xy
        I[2] -= mass[i] * x[i] * z[i];              // 13 xz
        I[5] -= mass[i] * y[i] * z[i];              // 23 yz
    }
    // lower triangle
    I[3] = I[1];    // 21 yx
    I[6] = I[2];    // 31 zx
    I[7] = I[5];    // 32 zy


    #ifdef debug_moment_of_inertia
    // output center of mass
        fprintf(stderr,
                "\nCenter of mass\n"
                "\t% .12le\t% .12le\t% .12le\n"
                ,x0, y0, z0
            );

    // output moment of inertia tensor
        int j;
        fprintf(stderr, "\nMoment of inertia tensor I\n");
        for(i = 0; i < 3; ++i){
            for(j = 0; j < 3; ++j){
                fprintf(stderr, "\t% .12le", I[i*3 + j]);
            }
            fprintf(stderr, "\n");
        }
    #endif


/* Corrected moment of inertia I' according to Watson:
//{{{

   I'_{a,b} = I_{a,b} - sum_{k,l,m} zeta_{k,m}^a zeta_{l,m}^b Q_k Q_l

   with a,b in {x,y,z} and k,l,m in mode_{1,...,n}

//}}}*/

    int a,b;    // in {x,y,z}
    int k,l,m;  // in mode_{1,...,n}

    for(a = 0; a < 3; ++a){
        for(b = 0; b < 3; ++b){

            for(k = 0; k < dim; ++k){
                for(l = 0; l < dim; ++l){

                    for(m = 0; m < dim; ++m){
                        I[a*3 + b] -= zeta[a][k*dim + m] * zeta[b][l*dim + m] * Q[k][state_index] * Q[l][state_index];
                    }

                }
            }

        }
    }

    #ifdef debug_moment_of_inertia
    // output corrected moment of inertia tensor
        fprintf(stderr, "\nCorrected moment of inertia\n");
        for(a = 0; a < 3; ++a){
            for(b = 0; b < 3; ++b){
                fprintf(stderr, "\t% .12le", I[a*3 + b]);
            }
            fprintf(stderr, "\n");
        }
    #endif
}
