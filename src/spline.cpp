#include "spline.h"
#include "eigen3/Eigen/Dense"

using namespace std;
using namespace Eigen;

vector<vector<double>> spline(vector<double> pos_x, vector<double> pos_y){
    
    int numSegments = pos_x.size() - 1;


    MatrixXd M = MatrixXd::Constant(3*numSegments, 3*numSegments, 0.0);

    // one equation for each segment between the knots of the spline
    // example function f0x(t) = a0x + b0x.t + c0x.(t^2) + d0x.(t^3)

    // o-----------------o-----------------o-----------------o
    
    //(x0,y0)   f0    (x1,y1)    f1     (x2,y2)    f2     (x3,y3)

    for (int i=0; i<numSegments; i++){
        M(i, 3*i) = 1;
        M(i, 3*i + 1) = 1;
        M(i, 3*i + 2) = 1;
    }

    // equations for the slope (first derivative) and 
    // the curvature (second derivative) at every knot between two segments

    for (int i=0; i<(numSegments-1); i++){
        M(i + numSegments, 3*i) = 1;
        M(i + numSegments, 3*i + 1) = 2;
        M(i + numSegments, 3*i + 2) = 3;
        M(i + numSegments, 3*i + 3) = -1;
        M(i + numSegments*2 - 1, 3*i + 1) = 1;
        M(i + numSegments*2 - 1, 3*i + 2) = 3;
        M(i + numSegments*2 - 1, 3*i + 4) = -1;
    }

    // the boundary conditions at extreme knots (second derivative = 0)

    M(numSegments*3 - 2, 1) = 1;
    M(numSegments*3 - 1, 3*numSegments - 2) = 1;
    M(numSegments*3 - 1, 3*numSegments - 1) = 3;

    MatrixXd Minv = M.inverse();

    VectorXd bx(3*numSegments);
    for (int i=0; i<3*numSegments; i++){
        if (i<numSegments)
        bx(i) = pos_x[i+1] - pos_x[i];
        else bx(i) = 0;
    }

    VectorXd by(3*numSegments);
    for (int i=0; i<3*numSegments; i++){
        if (i<numSegments)
        by(i) = pos_y[i+1] - pos_y[i];
        else by(i) = 0;
    }

    VectorXd ax = Minv * bx;
    VectorXd ay = Minv * by;

    // The new vectors x and y created
    float ss = 0.05; // ss is the stepsize according to which the 't' in f0(t) increments
    vector<double> xnew(numSegments/ss);
    vector<double> ynew(numSegments/ss);
    vector<vector<double>> pos_xy(numSegments/ss, vector<double>(2));

    int k = 0;
    int steps = 1/ss;

    for (int i=0; i<numSegments; i++){
        float t = 0.0;
        for (int j=0; j<steps; j++){
            xnew[k] = pos_x[i] + ax[3*i]*t + ax[3*i + 1]*t*t + ax[3*i + 2]*t*t*t;
            ynew[k] = pos_y[i] + ay[3*i]*t + ay[3*i + 1]*t*t + ay[3*i + 2]*t*t*t;
            pos_xy[k][0] = pos_x[i] + ax[3*i]*t + ax[3*i + 1]*t*t + ax[3*i + 2]*t*t*t;
            pos_xy[k][1] = pos_y[i] + ay[3*i]*t + ay[3*i + 1]*t*t + ay[3*i + 2]*t*t*t;
            k++;
            t = t + ss;
        }
    }
    return pos_xy;
}