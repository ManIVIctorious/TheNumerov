 
// provided prototypes
double integrate_1d(int n, double dx, double *integrand);

double integrate_1d(int n, double dx, double *integrand){

    int i;
    double integral = 0.0;

// trapezoidal rule for an equispaced grid with grid spacing h:
//  int_a^b f dx = h * ( f(a)/2 + sum_{i = 1}^{n-1} f(a + i*h) + f(b)/2 )

    integral = 0.5 * (integrand[0] + integrand[n-1]);

    for(i = 1; i < (n-1); ++i){
        integral += integrand[i];
    }

    return (integral * dx);
}
