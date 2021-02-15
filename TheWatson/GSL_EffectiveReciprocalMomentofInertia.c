
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include "settings.h"

// dependencies
int  InputComFile(char* inputfile, double* x, double* y, double* z, int max_lines);
void InvertMatrix(gsl_matrix* Matrix, gsl_matrix* InvMatrix, int dimension);

// provided prototypes
void EffectiveReciprocalMomentofInertia(settings *set, double* q, char* coordsfile);

void EffectiveReciprocalMomentofInertia(settings *set, double* q, char* coordsfile){

// allocate memory for x, y and z coordinate vectors
    double * x = malloc(set->n_atoms * sizeof(double));
    double * y = malloc(set->n_atoms * sizeof(double));
    double * z = malloc(set->n_atoms * sizeof(double));
    if( x == NULL ){ perror("x"); exit(errno); }
    if( y == NULL ){ perror("y"); exit(errno); }
    if( z == NULL ){ perror("z"); exit(errno); }


// read coordinate file
//--------------------------------------------------------------------------------
    int control = InputComFile(coordsfile, x, y, z, set->n_atoms);
    if(control != set->n_atoms){
        fprintf(stderr,
                "\n (-) Error in reading geometry file \"%s\""
                "\n     Number of atoms (%d) different from mode files (%d)"
                "\n     Aborting...\n\n"
                , coordsfile, control, set->n_atoms
            );
        exit(EXIT_FAILURE);
    }
    free(coordsfile); coordsfile = NULL;


#ifdef debug_coords
//{{{
    for(int i = 0; i < set->n_atoms; ++i){
        fprintf(stderr,
            "%2d\t% .12le\t% .12le\t% .12le\t% .12le\n"
            , i+1, x[i], y[i], z[i], set->atom_masses[i]
        );
    }
    fprintf(stderr, "\n");
//}}}
#endif


// calculate moment of inertia tensor and directly apply the Watson corrections
//--------------------------------------------------------------------------------
// make origin match system's centre of mass COM
//----------------------------------------------------------------------
    double COM[3] = {0.0, 0.0, 0.0};

// calculate centre of mass
    for(int i = 0; i < set->n_atoms; ++i){
        COM[0] += x[i] * set->atom_masses[i];
        COM[1] += y[i] * set->atom_masses[i];
        COM[2] += z[i] * set->atom_masses[i];
    }
    COM[0] /= set->tot_mass;
    COM[1] /= set->tot_mass;
    COM[2] /= set->tot_mass;

// translate system to centre of mass
    for(int i = 0; i < set->n_atoms; ++i){
        x[i] -= COM[0];
        y[i] -= COM[1];
        z[i] -= COM[2];
    }


// determine moment of inertia tensor I
//----------------------------------------------------------------------
    double I[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

// calculate tensor values
    for(int i = 0; i < set->n_atoms; ++i){
    // main diagonal:
        I[0] += set->atom_masses[i] * (y[i]*y[i] + z[i]*z[i]);  // 11 xx
        I[4] += set->atom_masses[i] * (x[i]*x[i] + z[i]*z[i]);  // 22 yy
        I[8] += set->atom_masses[i] * (x[i]*x[i] + y[i]*y[i]);  // 33 zz

    // upper triangle:
        I[1] -= set->atom_masses[i] * x[i] * y[i];              // 12 xy
        I[2] -= set->atom_masses[i] * x[i] * z[i];              // 13 xz
        I[5] -= set->atom_masses[i] * y[i] * z[i];              // 23 yz
    }
    // lower triangle
    I[3] = I[1];    // 21 yx
    I[6] = I[2];    // 31 zx
    I[7] = I[5];    // 32 zy


#ifdef debug_moment_of_inertia
//{{{
    fprintf(stderr, "\n\n# Filename:\t%s\n", coordsfile);
// output centre of mass
    fprintf(stderr,
            "\nCentre of mass\n"
            "\t% .12le\t% .12le\t% .12le\n"
            ,COM[0], COM[1], COM[2]
        );

// output moment of inertia tensor
    fprintf(stderr, "\nMoment of inertia tensor I\n");
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
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

// a,b in {x,y,z}
    for(int a = 0; a < 3; ++a){
    for(int b = 0; b < 3; ++b){
    // k,l,m in mode_{1,...,n}
        for(int k = 0; k < set->dimension; ++k){
        for(int l = 0; l < set->dimension; ++l){
        for(int m = 0; m < set->dimension; ++m){

            I[a*3 + b] -= set->zeta[a][k*set->dimension + m]*set->zeta[b][l*set->dimension + m] * q[k]*q[l];

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
    for(int a = 0; a < 3; ++a){
        for(int b = 0; b < 3; ++b){
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

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            gsl_matrix_set(I_corr, i, j, I[i*3 + j]);
        }
    }

    InvertMatrix(I_corr, mu, 3);
    gsl_matrix_free(I_corr); I_corr = NULL;


#ifdef debug_mu
//{{{
    fprintf(stderr, "\nEffective Reciprocal Inertia Tensor (mu)\n");
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
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

// output q and free its memory afterwards
    for(int i = 0; i < set->dimension; ++i){
        fprintf(set->fdout, "\t% .12le", q[i]);
    }
    free(q); q = NULL;

// output mu
    for(int i = 0; i < 3; ++i){
        for(int j = i; j < 3; ++j){
            fprintf(set->fdout, "\t% .12le", gsl_matrix_get(mu, i, j));
        }
    }
    fprintf(set->fdout, "\n");

    gsl_matrix_free(mu); mu = NULL;
}
