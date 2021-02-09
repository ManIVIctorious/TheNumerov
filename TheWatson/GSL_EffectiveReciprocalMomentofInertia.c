
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include "settings.h"

// dependencies
int  InputComFile(char* inputfile, double* x, double* y, double* z, int max_lines);
void InvertMatrix(gsl_matrix* Matrix, gsl_matrix* InvMatrix, int dimension);

// provided prototypes
void EffectiveReciprocalMomentofInertia(settings prefs, double* q, char* coordsfile);

void EffectiveReciprocalMomentofInertia(settings prefs, double* q, char* coordsfile){

    int i, j;
    int control;

// allocate memory for x, y and z coordinate vectors
    double * x = malloc(prefs.n_atoms * sizeof(double));
    double * y = malloc(prefs.n_atoms * sizeof(double));
    double * z = malloc(prefs.n_atoms * sizeof(double));
    if( x == NULL ){ perror("x"); exit(errno); }
    if( y == NULL ){ perror("y"); exit(errno); }
    if( z == NULL ){ perror("z"); exit(errno); }


// read coordinate file
//--------------------------------------------------------------------------------
    control = InputComFile(coordsfile, x, y, z, prefs.n_atoms);
    if(control != prefs.n_atoms){
        fprintf(stderr,
                "\n (-) Error in reading coordinate file \"%s\""
                "\n     Number of coordinates (%d) does not match number of atoms (%d)"
                "\n     Aborting...\n\n"
                , coordsfile, control, prefs.n_atoms
            );
        exit(EXIT_FAILURE);
    }


#ifdef debug_coords
//{{{
    for(i = 0; i < prefs.n_atoms; ++i){
        fprintf(stderr,
            "%2d\t% .12le\t% .12le\t% .12le\t% .12le\n"
            , i+1, x[i], y[i], z[i], prefs.atom_masses[i]
        );
    }
    fprintf(stderr, "\n");
//}}}
#endif


// calculate moment of inertia tensor and directly apply the Watson corrections
//--------------------------------------------------------------------------------
// make origin match system's center of mass com
//----------------------------------------------------------------------
    double com[3] = {0.0, 0.0, 0.0};

// calculate center of mass
    for(i = 0; i < prefs.n_atoms; ++i){
        com[0] += x[i] * prefs.atom_masses[i];
        com[1] += y[i] * prefs.atom_masses[i];
        com[2] += z[i] * prefs.atom_masses[i];
    }
    com[0] /= prefs.tot_mass;
    com[1] /= prefs.tot_mass;
    com[2] /= prefs.tot_mass;

// translate system to center of mass
    for(i = 0; i < prefs.n_atoms; ++i){
        x[i] -= com[0];
        y[i] -= com[1];
        z[i] -= com[2];
    }


// determine moment of inertia tensor I
//----------------------------------------------------------------------
    double I[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

// calculate tensor values
    for(i = 0; i < prefs.n_atoms; ++i){
    // main diagonal:
        I[0] += prefs.atom_masses[i] * (y[i]*y[i] + z[i]*z[i]);  // 11 xx
        I[4] += prefs.atom_masses[i] * (x[i]*x[i] + z[i]*z[i]);  // 22 yy
        I[8] += prefs.atom_masses[i] * (x[i]*x[i] + y[i]*y[i]);  // 33 zz

    // upper triangle:
        I[1] -= prefs.atom_masses[i] * x[i] * y[i];              // 12 xy
        I[2] -= prefs.atom_masses[i] * x[i] * z[i];              // 13 xz
        I[5] -= prefs.atom_masses[i] * y[i] * z[i];              // 23 yz
    }
    // lower triangle
    I[3] = I[1];    // 21 yx
    I[6] = I[2];    // 31 zx
    I[7] = I[5];    // 32 zy


#ifdef debug_moment_of_inertia
//{{{
    fprintf(stderr, "\n\n# Filename:\t%s\n", coordsfile);
// output center of mass
    fprintf(stderr,
            "\nCenter of mass\n"
            "\t% .12le\t% .12le\t% .12le\n"
            ,com[0], com[1], com[2]
        );

// output moment of inertia tensor
    fprintf(stderr, "\nMoment of inertia tensor I\n");
    for(i = 0; i < 3; ++i){
        for(j = 0; j < 3; ++j){
            fprintf(stderr, "\t% .12le", I[i*3 + j]);
        }
        fprintf(stderr, "\n");
    }
//}}}
#endif


/* Corrected moment of inertia I' according to Watson:
//----------------------------------------------------------------------
//{{{

   I'_{a,b} = I_{a,b} - sum_{k,l,m} zeta_{k,m}^a zeta_{l,m}^b Q_k Q_l

   with a,b in {x,y,z} and k,l,m in mode_{1,...,n}

//}}}*/

    int a, b;       // in {x,y,z}
    int k, l, m;    // in mode_{1,...,n}

    for(a = 0; a < 3; ++a){
      for(b = 0; b < 3; ++b){

        for(k = 0; k < prefs.dimension; ++k){
          for(l = 0; l < prefs.dimension; ++l){
            for(m = 0; m < prefs.dimension; ++m){

                I[a*3 + b] -= prefs.zeta[a][k*prefs.dimension + m] * prefs.zeta[b][l*prefs.dimension + m] * q[k] * q[l];

            }
          }
        }

      }
    }

// x, y and z are not required anymore
    free(x); x = NULL;
    free(y); y = NULL;
    free(z); z = NULL;


#ifdef debug_moment_of_inertia
//{{{
// output corrected moment of inertia tensor
    fprintf(stderr, "\nCorrected moment of inertia\n");
    for(a = 0; a < 3; ++a){
        for(b = 0; b < 3; ++b){
            fprintf(stderr, "\t% .12le", I[a*3 + b]);
        }
        fprintf(stderr, "\n");
    }
//}}}
#endif


// calculate Effective Reciprocal Inertia Tensor mu by inverting I'
//----------------------------------------------------------------------
    gsl_matrix * I_corr = gsl_matrix_calloc(3, 3);
    gsl_matrix * mu     = gsl_matrix_calloc(3, 3);

    for(i = 0; i < 3; ++i){
        for(j = 0; j < 3; ++j){
            gsl_matrix_set(I_corr, i, j, I[i*3 + j]);
        }
    }

    InvertMatrix(I_corr, mu, 3);
    gsl_matrix_free(I_corr); I_corr = NULL;


#ifdef debug_mu
//{{{
    fprintf(stderr, "\nEffective Reciprocal Inertia Tensor (mu)\n");
    for(i = 0; i < 3; ++i){
        for(j = 0; j < 3; ++j){
            fprintf(stderr, "\t% .12le", gsl_matrix_get(mu, i, j));
        }
        fprintf(stderr, "\n");
    }
//}}}
#endif


/* Output unique values of mu (upper triangle with main diagonal)
//----------------------------------------------------------------------
//{{{

    Since I is symmetric (I = Iᵀ) and its correction is symmetric => I' = (I')ᵀ
    its inverse mu = (I')⁻¹ must be symmetric as well:

        A = Aᵀ      A * A⁻¹ = 𝟙     𝟙 = 𝟙ᵀ      (AB)ᵀ = BᵀAᵀ

        A * A⁻¹ = 𝟙 = 𝟙ᵀ = (A * A⁻¹)ᵀ = (A⁻¹)ᵀ * Aᵀ

              A * A⁻¹ = (A⁻¹)ᵀ * Aᵀ       | A = Aᵀ
              A * A⁻¹ = (A⁻¹)ᵀ * A        | A * A⁻¹ = 𝟙 = A⁻¹ * A
              A⁻¹ * A = (A⁻¹)ᵀ * A        | * A⁻¹
        A⁻¹ * A * A⁻¹ = (A⁻¹)ᵀ * A * A⁻¹

                  A⁻¹ = (A⁻¹)ᵀ

    If A symmetric its inverse must be symmetric as well, q.e.d

//}}}*/

// output q
    for(i = 0; i < prefs.dimension; ++i){
        fprintf(prefs.fdout, "\t% .12le", q[i]);
    }
// output mu
    for(i = 0; i < 3; ++i){
        for(j = i; j < 3; ++j){
            fprintf(prefs.fdout, "\t% .12le", gsl_matrix_get(mu, i, j));
        }
    }
    fprintf(prefs.fdout, "\n");

    gsl_matrix_free(mu); mu = NULL;
}
