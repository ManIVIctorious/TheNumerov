
// provided prototypes
double Integrate3D(int* nq, double dx, double* integrand);

double Integrate3D(int* nq, double dx, double* integrand){

    double integral = 0.0;

// start indices
    int front  = nq[2]*nq[1]*(nq[0]-1); // front top    left
    int right  = (nq[2]-1);             // back  top    right
    int bottom = nq[2]*(nq[1]-1);       // back  bottom left

//--------------------------------------------------
// inner part without faces, edges and corners
// factor: 1.0
    for(int i = 1; i < (nq[0] - 1); ++i){
        for(int j = 1; j < (nq[1] - 1); ++j){
            for(int k = 1; k < (nq[2] - 1); ++k){

                int index = (i*nq[1] + j)*nq[2] + k;
                integral += integrand[index];

            }
        }
    }

//--------------------------------------------------
// six faces, without edges and corners
//  factor: 1/2
//  front (i = 0) and back (i = nq[0]) face
    for(int j = 1; j < (nq[1] - 1); ++j){
        for(int k = 1; k < (nq[2] - 1); ++k){

            int index = j*nq[2] + k;

            integral += integrand[index        ]*0.5; // back  face
            integral += integrand[index + front]*0.5; // front face
        }
    }

//  top (j = 0) and bottom (j = nq[1]) face
    for(int i = 1; i < (nq[0] - 1); ++i){
        for(int k = 1; k < (nq[2] - 1); ++k){

            int index  = i*nq[1]*nq[2] + k;

            integral += integrand[index         ]*0.5; // top    face
            integral += integrand[index + bottom]*0.5; // bottom face
        }
    }

//  left (k = 0) and right (k = nq[2]) side face
    for(int i = 1; i < (nq[0] - 1); ++i){
        for(int j = 1; j < (nq[1] - 1); ++j){

            int index = (i*nq[1] + j)*nq[2];

            integral += integrand[index        ]*0.5; // left  face
            integral += integrand[index + right]*0.5; // right face
        }
    }

//--------------------------------------------------
// twelve edges, without corners
//  factor: 1/4
    for(int k = 1; k < (nq[2] - 1); ++k){

        int index  = k;

        integral += integrand[index                 ]*0.25; // back   top
        integral += integrand[index + bottom        ]*0.25; // back   bottom
        integral += integrand[index + front         ]*0.25; // front  top
        integral += integrand[index + front + bottom]*0.25; // front  bottom
    }

    for(int j = 1; j < (nq[1] - 1); ++j){

        int index = j*nq[2];

        integral += integrand[index                ]*0.25; // back  left
        integral += integrand[index + right        ]*0.25; // back  right
        integral += integrand[index + front        ]*0.25; // front left
        integral += integrand[index + front + right]*0.25; // front right
    }

//  the four front to back facing edges
    for(int i = 1; i < (nq[0] - 1); ++i){

        int index  = i*nq[2]*nq[1];
        int right  = (nq[2]-1);
        int bottom = nq[2]*(nq[1]-1);

        integral += integrand[index                 ]*0.25; // top    left
        integral += integrand[index + right         ]*0.25; // top    right
        integral += integrand[index + bottom        ]*0.25; // bottom left
        integral += integrand[index + bottom + right]*0.25; // bottom right

    }

//--------------------------------------------------
// eight corners
//  factor: 1/8
    integral += integrand[       0              ]*0.125; // back  top    left
    integral += integrand[right                 ]*0.125; // back  top    right
    integral += integrand[bottom                ]*0.125; // back  bottom left
    integral += integrand[bottom + right        ]*0.125; // back  bottom right

    integral += integrand[front                 ]*0.125; // front top    left
    integral += integrand[front + right         ]*0.125; // front top    right
    integral += integrand[front + bottom        ]*0.125; // front bottom left
    integral += integrand[front + bottom + right]*0.125; // front bottom right

    return (integral * dx * dx * dx);
}
