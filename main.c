
#include <stdio.h>
#include <stdlib.h>

int InputComFile(char *inputfile, double **x, double **y, double **z);
int InputNormalMode(char *inputfile, double **dx, double **dy, double **dz, double **mass);

int main(int argc, char **argv){

    int i, j, k;
    int control, control2;
    double threshold = 1E-10;

// input files
//    char * comfile   = "2D_scan_resorcinol-A_b3lyp_ccl4_modes_35_36_dr=+0.00000000000000_+0.00000000000000.com";
    char * comfile   = "2D_scan_resorcinol-A_b3lyp_ccl4_modes_35_36_dr=-0.05162332776183_-0.15486998328547.com";
    char * modefile1 = "mode_35";
    char * modefile2 = "mode_36";
// coordinates
    int n_atoms;
    double * x;
    double * y;
    double * z;
// first mode
    double * dx1;
    double * dy1;
    double * dz1;
    double * mass;
// second mode
    double * dx2;
    double * dy2;
    double * dz2;
    double * mass2; // freed

// center of mass coordinates and moment of inertia tensor
    double x0, y0, z0, tot_mass;
    double MomentOfInertia[9];


//------------------------------------------------------------------------------------------------------------------
// Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input Input
//------------------------------------------------------------------------------------------------------------------

// Error code cheat sheet:
//      1: Wrong input (program handling)
//      2: Memory allocation problem
//      3: Problem with input function
//      4: File inconsistencies

// check input argument
    if(comfile == NULL || modefile1 == NULL || modefile2 == NULL){
        fprintf(stderr, "\n (-) Please specify valid coordinates and mode files");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(1);
    }

// input of com file
    x = malloc(sizeof(double));
    y = malloc(sizeof(double));
    z = malloc(sizeof(double));
    if(x == NULL || y == NULL || z == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for x, y or z");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(2);
    }
    n_atoms = InputComFile(comfile, &x, &y, &z);
    if(n_atoms < 0){
        fprintf(stderr, "\n (-) Error in reading input file \"%s\"", comfile);
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(3);
    }

// input of first mode
    dx1   = malloc(sizeof(double));
    dy1   = malloc(sizeof(double));
    dz1   = malloc(sizeof(double));
    mass  = malloc(sizeof(double));
    if(dx1 == NULL || dy1 == NULL || dz1 == NULL || mass == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for dx1, dy1, dz1 or mass");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(2);
    }
    control  = InputNormalMode(modefile1, &dx1, &dy1, &dz1, &mass);

// input of second mode
    dx2   = malloc(sizeof(double));
    dy2   = malloc(sizeof(double));
    dz2   = malloc(sizeof(double));
    mass2 = malloc(sizeof(double));
    if(dx2 == NULL || dy2 == NULL || dz2 == NULL || mass2 == NULL){
        fprintf(stderr, "\n (-) Error in memory allocation for dx2, dy2, dz2 or mass");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(2);
    }
    control2 = InputNormalMode(modefile2, &dx2, &dy2, &dz2, &mass2);

// check if all files contain the same number of atoms
    if(control != n_atoms || control2 != n_atoms){
        fprintf(stderr, "\n (-) Error in reading input files.");
        fprintf(stderr, "\n     Number of atoms doesn't match.");
        fprintf(stderr, "\n     Aborting...\n\n");
        exit(4);
    }

// check for same masses in mode files (atom order/type)
    for(i = 0; i < n_atoms; ++i){
        if( labs(mass[i] - mass2[i]) > threshold ){
            fprintf(stderr, "\n (-) Error in reading input files.");
            fprintf(stderr, "\n     Atom masses differ between modes.\n");
            fprintf(stderr, "\n     Aborting...\n\n");
            exit(4);
        }
    }
    free(mass2); mass2 = NULL;

//// output input for control
//    fprintf(stderr,"\nNumber of Atoms: %d\n", n_atoms);
//    fprintf(stderr,"\nInput coordinates\n", n_atoms);
//    fprintf(stderr,"\t x             ");
//    fprintf(stderr,"\t y             ");
//    fprintf(stderr,"\t z             ");
//    fprintf(stderr,"\t dx1           ");
//    fprintf(stderr,"\t dy1           ");
//    fprintf(stderr,"\t dz1           ");
//    fprintf(stderr,"\t dx2           ");
//    fprintf(stderr,"\t dy2           ");
//    fprintf(stderr,"\t dz2           ");
//    fprintf(stderr,"\t mass          ");
//    fprintf(stderr,"\n");
//    for(i = 0; i < n_atoms; ++i){
//        fprintf(stderr,"\t% .8le\t% .8le\t% .8le", x[i], y[i], z[i]);   
//        fprintf(stderr,"\t% .8le\t% .8le\t% .8le", dx1[i], dy1[i], dz1[i]);   
//        fprintf(stderr,"\t% .8le\t% .8le\t% .8le", dx2[i], dy2[i], dz2[i]);   
//        fprintf(stderr,"\t% .8le", mass[i]);   
//        fprintf(stderr,"\n");
//    }


//------------------------------------------------------------------------------------------------------------------
// Moment of inertia  Moment of inertia  Moment of inertia  Moment of inertia  Moment of inertia  Moment of inertia
//------------------------------------------------------------------------------------------------------------------

// calculate center of mass
    for(i = 0; i < n_atoms; ++i){
        x0 += x[i] * mass[i];
        y0 += y[i] * mass[i];
        z0 += z[i] * mass[i];

        tot_mass += mass[i];
    }
    x0 /= tot_mass;
    y0 /= tot_mass;
    z0 /= tot_mass;

//// output center of mass
//    fprintf(stderr,"\nCenter of mass\n");
//    fprintf(stderr,"\t% 15.8lf", x0);
//    fprintf(stderr,"\t% 15.8lf", y0);
//    fprintf(stderr,"\t% 15.8lf", z0);
//    fprintf(stderr,"\n");

// translate molecule center of mass to origin
    for(i = 0; i < n_atoms; ++i){
        x[i] -= x0;
        y[i] -= y0;
        z[i] -= z0;
    }
    x0 = y0 = z0 = 0.0;

// calculate moment of inertia tensor
    for(i = 0; i < n_atoms; ++i){
    // main diagonal
        MomentOfInertia[0] += mass[i] * (y[i]*y[i] + z[i]*z[i]);  // 11 xx
        MomentOfInertia[4] += mass[i] * (x[i]*x[i] + z[i]*z[i]);  // 22 yy
        MomentOfInertia[8] += mass[i] * (x[i]*x[i] + y[i]*y[i]);  // 33 zz
    // upper triangle
        MomentOfInertia[1] -= mass[i] * x[i] * y[i];              // 12 xy
        MomentOfInertia[2] -= mass[i] * x[i] * z[i];              // 13 xz
        MomentOfInertia[5] -= mass[i] * y[i] * z[i];              // 23 yz
    }
    // lower triangle
    MomentOfInertia[3] = MomentOfInertia[1];    // 21 yx
    MomentOfInertia[6] = MomentOfInertia[2];    // 31 zx
    MomentOfInertia[7] = MomentOfInertia[5];    // 32 zy

// output moment of inertia tensor
    fprintf(stderr,"\nMoment of inertia tensor\n");
    for(i = 0; i < 3; ++i){
        for(j = 0; j < 3; ++j){
            fprintf(stderr,"\t% 15.8lf", MomentOfInertia[i*3 + j]);
        }
        fprintf(stderr,"\n");
    }

//------------------------------------------------------------------------------------------------------------------
// Coriolis coefficients  Coriolis coefficients  Coriolis coefficients  Coriolis coefficients  Coriolis coefficients
//------------------------------------------------------------------------------------------------------------------
    int dimension;

    double zeta[3];

    for(i = 0; i < n_atoms; ++i){
        zeta[0] += dy1[i]*dz2[i] - dz1[i]*dy2[i];
        zeta[1] += dz1[i]*dx2[i] - dx1[i]*dz2[i];
        zeta[2] += dx1[i]*dy2[i] - dy1[i]*dx2[i];
    } 

// output moment of inertia tensor
    fprintf(stderr,"\nCorriolis coefficients\n");
    fprintf(stderr,"\t zeta_x        ");
    fprintf(stderr,"\t zeta_y        ");
    fprintf(stderr,"\t zeta_z        ");
    fprintf(stderr,"\n");
    fprintf(stderr,"\t% .8le", zeta[0]);
    fprintf(stderr,"\t% .8le", zeta[1]);
    fprintf(stderr,"\t% .8le", zeta[2]);
    fprintf(stderr,"\n");

//------------------------------------------------------------------------------------------------------------------
// Coriolis corrected moment of inertia  Coriolis corrected moment of inertia  Coriolis corrected moment of inertia   
//------------------------------------------------------------------------------------------------------------------

    double deviation_a = -0.05162332776183;
    double deviation_b = -0.15486998328547;


    return 0;
}
