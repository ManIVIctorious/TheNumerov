
// provided prototypes
double integrate_3d(int nq1, int nq2, int nq3, double dx, double *integrand);

double integrate_3d(int nq1, int nq2, int nq3, double dx, double *integrand){

    int i, j, k, index;
    double integral = 0.0;

//--------------------------------------------------
// inner part without faces, edges and corners
// factor: 1.0
    for(i = 1; i < (nq1 - 1); ++i){
        for(j = 1; j < (nq2 - 1); ++j){
            for(k = 1; k < (nq3 - 1); ++k){

                index = i*nq2*nq3 + j*nq3 + k;
                integral += integrand[index];

            }
        }
    }

//--------------------------------------------------
// six side faces, without edges and corners
//  factor: 1/2
//  front and back face
    for(j = 1; j < (nq2 - 1); ++j){
        for(k = 1; k < (nq3 - 1); ++k){

            index = j*nq3 + k;

            integral += integrand[index]/2.0; // front face

            integral += integrand[index + (nq1-1)*nq2*nq3]/2.0; // back face

        }// endfor k
    }// endfor j

//  top and bottom face
    for(i = 1; i < (nq1 - 1); ++i){
        for(k = 1; k < (nq3 - 1); ++k){

            index = i*nq2*nq3 + k;

            integral += integrand[index]/2.0; // top face

            integral += integrand[index + (nq2-1)*nq3]/2.0; // bottom face

        }// endfor k
    }// endfor i

//  left and right side face
    for(i = 1; i < (nq1 - 1); ++i){
        for(j = 1; j < (nq2 - 1); ++j){

            index = i*nq2*nq3 + j*nq3;

            integral += integrand[index]/2.0;

            integral += integrand[index + (nq3-1)]/2.0;

        }// endfor k
    }// endfor j

//--------------------------------------------------
// twelve edges, without corners
//  factor: 1/4
//  the four left to right edges
    for(k = 1; k < (nq3 - 1); ++k){

        index = k;
        integral += integrand[index]/4.0; // front, top

        index = k + (nq2-1)*nq3;
        integral += integrand[index]/4.0; // front, bottom

        index = k + (nq1-1)*nq2*nq3;
        integral += integrand[index]/4.0; // back, top

        index = k + (nq1-1)*nq2*nq3 + (nq2-1)*nq3;
        integral += integrand[index]/4.0; // back, bottom

    }

//  the four top to bottom edges
    for(j = 1; j < (nq2 - 1); ++j){

        index = j*nq3;
        integral += integrand[index]/4.0; // front, left

        index = j*nq3 + (nq3-1);
        integral += integrand[index]/4.0; // front, right

        index = j*nq3 + (nq1-1)*nq2*nq3;
        integral += integrand[index]/4.0; // back, left

        index = j*nq3 + (nq1-1)*nq2*nq3 + (nq3-1);
        integral += integrand[index]/4.0; // back, right

    }

//  the four front to back facing edges
    for(i = 1; i < (nq1 - 1); ++i){

        index = i*nq2*nq3;
        integral += integrand[index]/4.0; // top, left

        index = i*nq2*nq3 + (nq3-1);
        integral += integrand[index]/4.0; // top, right

        index = i*nq2*nq3 + (nq2-1)*nq3;
        integral += integrand[index]/4.0; // bottom, left

        index = i*nq2*nq3 + (nq2-1)*nq3 + (nq3-1);
        integral += integrand[index]/4.0; // bottom, right

    }

//--------------------------------------------------
// eight corners
//  factor: 1/8

// front
// top, left
    integral += integrand[   0   ]/8.0;
// top, right
    integral += integrand[(nq3-1)]/8.0;
// bottom, left
    integral += integrand[(nq2-1)*nq3]/8.0;
// bottom, right
    integral += integrand[(nq2-1)*nq3 + (nq3-1)]/8.0;

// back
// top, left
    integral += integrand[(nq1-1)*nq2*nq3]/8.0;
// top, right
    integral += integrand[(nq1-1)*nq2*nq3 + (nq3-1)]/8.0;
// bottom, left
    integral += integrand[(nq1-1)*nq2*nq3 + (nq2-1)*nq3]/8.0;
// bottom, right
    integral += integrand[nq1*nq2*nq3 - 1]/8.0;

    return (integral * dx * dx * dx);
}
