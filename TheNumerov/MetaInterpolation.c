
// Dependencies
int nx1dInterpolation(double** v, int* nq, double dq, int dimension, int n_spline);

// provided prototypes
int MetaInterpolation(double* *v, int* nq, double dq, int dimension, int n_spline);


int MetaInterpolation(double* *v, int* nq, double dq, int dimension, int n_spline){

    int control;

    switch(dimension){

        default:
            control = nx1dInterpolation(v, nq, dq, dimension, n_spline);
            break;

    }

    return control;
}
