#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector <VectorXd> &estimations,
                              const vector <VectorXd> &ground_truth) {
    /**
    TODO:
      * Calculate the RMSE here.
      * note:
      *   code taken from answer at 5.23
    */

    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    if (estimations.size() == 0) {
        cout << "Error: estimations size 0\n";
        return rmse;
    }
    //  * the estimation vector size should equal ground truth vector size
    if (estimations.size() != ground_truth.size()) {
        cout << "Error: estimations size not equal to ground truth size\n";
    }

    //accumulate squared residuals
    for (int i = 0; i < estimations.size(); ++i) {
        for (int j = 0; j < 4; ++j) {
            rmse(j) += pow(estimations[i](j) - ground_truth[i](j), 2);
        }
    }

    //calculate the mean
    for (int j = 0; j < 4; ++j) {
        rmse(j) /= estimations.size();
    }

    //calculate the squared root
    for (int j = 0; j < 4; ++j) {
        rmse(j) = sqrt(rmse(j));
    }

    //return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {
    /**
    TODO:
      * Calculate a Jacobian here.
      * note:
      *   code taken from answer at 5.19
    */

    MatrixXd Hj(3, 4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //check division by zero
    if (px == 0 and py == 0) {
        cout << "Error: ZeroDivisionError\n";
        return Hj;
    }

    float pxpxpypy = px * px + py * py;
    float sqrt_pxpxpypy = sqrt(pxpxpypy);

    //compute the Jacobian matrix
    Hj(0, 0) = px / sqrt_pxpxpypy;
    Hj(0, 1) = py / sqrt_pxpxpypy;
    Hj(0, 2) = 0;
    Hj(0, 3) = 0;

    Hj(1, 0) = -py / (pxpxpypy);
    Hj(1, 1) = px / (pxpxpypy);
    Hj(1, 2) = 0;
    Hj(1, 3) = 0;

    Hj(2, 0) = py * (vx * py - vy * px) / pow((pxpxpypy), 1.5);
    Hj(2, 1) = px * (vy * px - vx * py) / pow((pxpxpypy), 1.5);
    Hj(2, 2) = px / sqrt_pxpxpypy;
    Hj(2, 3) = py / sqrt_pxpxpypy;

    return Hj;
}
