
// provided prototypes
double Integrate1D(int* nq, double dx, double *integrand);

double Integrate1D(int* nq, double dx, double *integrand){

// trapezoidal rule for an equispaced grid with grid spacing h:
//  int_a^b f dx = h * ( f(a)/2 + sum_{i = 1}^{n-1} f(a + i*h) + f(b)/2 )

    double integral = 0.5 * (integrand[0] + integrand[nq[0]-1]);

    for(int i = 1; i < (nq[0] - 1); ++i){
        integral += integrand[i];
    }

    return (integral * dx);
}
