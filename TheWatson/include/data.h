#ifndef _DATA_H
#define _DATA_H

typedef struct data {

    int       dimension;    // number of included modes
    int       n_atoms;      // number of atoms
    double    tot_mass;     // total mass of the system
    double *  atom_masses;  // atomic masses
    double ** zeta;         // {x,y,z}[zeta_00, zeta_01, zeta_12, ...]

} data;

#endif
