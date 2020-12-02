
// provided prototypes
double Integrate2D(int* nq, double dx, double* integrand);

double Integrate2D(int* nq, double dx, double* integrand){

    double integral = 0.0;

// Analogous to the 1D trapezoidal rule the inner part
//  (without edges and corners) gets summed up completely
    for(int i = 1; i < (nq[0] - 1); ++i){
        for(int j = 1; j < (nq[1] - 1); ++j){

            int index = i*nq[1] + j;
            integral += integrand[index];

        }
    }

// while the edges (without corners) get summed up half as often
    for(int i = 1; i < (nq[0] - 1); ++i){
        integral += integrand[i*nq[1]            ]*0.5; // left  edge
        integral += integrand[i*nq[1] + (nq[1]-1)]*0.5; // right edge
    }
    for(int i = 1; i < (nq[1] - 1); ++i){
        integral += integrand[i                  ]*0.5; // upper edge
        integral += integrand[i + (nq[0]-1)*nq[1]]*0.5; // lower edge
    }

// and the corners halve of that
    integral += integrand[    0           ]*0.25;   // upper, left
    integral += integrand[(nq[1]-1)       ]*0.25;   // upper, right
    integral += integrand[(nq[0]-1)*nq[1] ]*0.25;   // lower, left
    integral += integrand[ nq[0]*nq[1] - 1]*0.25;   // lower, right

    return (integral * dx * dx);
}
