#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
    /**
    TODO:
      * predict the state
    */

    std::cout << "predicting...\n";

    // x' = Fx + v
    // from 5.13
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_; // uncertainty increase

}

void KalmanFilter::Update(const VectorXd &z) {
    /**
    TODO:
      * update the state by using Kalman Filter equations
    */

    // this is called for lidar measurements, EKF will be used for radar
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /**
    TODO:
      * update the state by using Extended Kalman Filter equations
      *
      * now has to use:
      *   x' = f(x) +v
      * and
      *   z = h(x') + w
    */
    // get state values
    double px = x_(0);
    double py = x_(1);
    double vx = x_(2);
    double vy = x_(3);

    // build three equations for h(x) function
    VectorXd H_func(3);

    // too close to zero check and rho
    double eq1 = sqrt(px * px + py * py);
    if (eq1 < 0.00001) {
        (px > 0) ? px += 0.001 : px -= 0.001;
        (py > 0) ? py += 0.001 : py -= 0.001;
        eq1 = sqrt(px * px + py * py);
    }

    // phi
    double eq2 = atan2(py, px);

    // rhodot
    double eq3 = (px * vx + py * vy) / sqrt(px * px + py * py);
    H_func << eq1, eq2, eq3;

    VectorXd y = z - H_func;  // y = z - h(x')

    // normalize resulting phi angle to range -pi to +pi
    if (y(1) > M_PI) {
        y(1) -= 2 * M_PI;
    } else if (y(1) < -M_PI) {
        y(1) += 2 * M_PI;
    }

    // used for radar measurements because H cannot be used and h function has to be used instead

    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
