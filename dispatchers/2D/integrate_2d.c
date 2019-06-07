
// provided prototypes
double integrate_2d(int nq1, int nq2, double dx, double *integrand);

double integrate_2d(int nq1, int nq2, double dx, double *integrand){

    int i, j, index;
    double integral = 0.0;

// Analogous to the 1D trapezoidal rule the inner part
//  (without edges and corners) gets summed up completely
    for(i = 1; i < (nq1 - 1); ++i){
        for(j = 1; j < (nq2 - 1); ++j){

            index = i*nq2 + j;
            integral += integrand[index];

        }
    }

// while the edges (without corners) get summed up half as often
// left and right edges
    for(i = 1; i < (nq1 - 1); ++i){

        integral += integrand[i*nq2]/2.0;             // left edge
        integral += integrand[i*nq2 + (nq2-1)]/2.0;   // right edge

    }

// upper and lower edges
    for(j = 1; j < (nq2 - 1); ++j){

        integral += integrand[j]/2.0;                 // upper edge
        integral += integrand[j + (nq1-1)*nq2]/2.0;   // lower edge

    }

// and the corners are even halve of that
    integral += integrand[   0   ]/4.0;     // top, left
    integral += integrand[(nq2-1)]/4.0;     // top, right
    integral += integrand[(nq1-1)*nq2]/4.0; // bottom, left
    integral += integrand[nq1*nq2 - 1]/4.0; // bottom, right

    return (integral * dx * dx);
}
