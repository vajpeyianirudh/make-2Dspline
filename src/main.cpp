#include <iostream>
#include <vector>
#include <fstream>
#include "eigen3/Eigen/Dense"
#include "spline.h"

using namespace std;
using namespace Eigen;

int main(){
    cout << "Hi! Let us make a spline making project." << endl;

    vector<double> scope_posx{11.08, -68.59, -95.8, -47.34, 68.79, 104.94, 24.66, 1.21};
    vector<double> scope_posy{-504.95, -40.55, -279.81, -163.26, -129.54, -235.75, -333.73, 441.4};

    int numSegments = scope_posx.size() - 1;
    cout << "Number of segments = " << numSegments << endl;

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
        bx(i) = scope_posx[i+1] - scope_posx[i];
        else bx(i) = 0;
    }

    VectorXd by(3*numSegments);
    for (int i=0; i<3*numSegments; i++){
        if (i<numSegments)
        by(i) = scope_posy[i+1] - scope_posy[i];
        else by(i) = 0;
    }

    VectorXd ax = Minv * bx;
    VectorXd ay = Minv * by;

    // The new vectors x and y created
    float ss = 0.05; // ss is the stepsize according to which the 't' in f0(t) increments
    vector<double> xnew(numSegments/ss);
    vector<double> ynew(numSegments/ss);

    int k = 0;
    int steps = 1/ss;

    for (int i=0; i<numSegments; i++){
        float t = 0.0;
        for (int j=0; j<steps; j++){
            xnew[k] = scope_posx[i] + ax[3*i]*t + ax[3*i + 1]*t*t + ax[3*i + 2]*t*t*t;
            ynew[k] = scope_posy[i] + ay[3*i]*t + ay[3*i + 1]*t*t + ay[3*i + 2]*t*t*t;
            k++;
            t = t + ss;
        }
    }

    ofstream XYnew;
    XYnew.open("XYnew.csv");
    for (int i=0; i<xnew.size(); i++){
        XYnew << to_string(xnew[i]) << "," << to_string(ynew[i]) << endl;
    }
    XYnew.close();

    // VALUES FROM FUNCTION

    vector<vector<double>> pos_xynew;
    pos_xynew = spline(scope_posx, scope_posy);
    ofstream funcXYnew;
    funcXYnew.open("funcXYnew.csv");
    for (int i=0; i<xnew.size(); i++){
        funcXYnew << to_string(pos_xynew[i][0]) << "," << to_string(pos_xynew[i][1]) << endl;
    }
    funcXYnew.close();

    cout << "spline generation successful" << endl;

    return 0;
}