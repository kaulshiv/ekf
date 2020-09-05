#include "kalman_filter.h"
#include <math.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
  cout << "PREDICT" << endl;
  x_ = F_*x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  cout << "UPDATE" << endl;
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S_ = H_ * P_ * Ht + R_;
  MatrixXd Si = S_.inverse();
  MatrixXd K =  P_ * Ht * Si;

  // new state
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(4, 4);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  cout << "UPDATE EKF" << endl;
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  float rho = pow(pow(px, 2) + pow(py, 2),0.5);

  VectorXd h_of_x = VectorXd(3);

  h_of_x << rho, atan2(py,px), (px*vx + py*vy)/rho;

  VectorXd y = z - h_of_x;

  MatrixXd Hjt = H_.transpose();
  MatrixXd S_ = H_ * P_ * Hjt + R_;
  MatrixXd Si = S_.inverse();
  MatrixXd K =  P_ * Hjt * Si;

  // new state
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(4, 4);
  P_ = (I - K * H_) * P_;

  cout << "y = " << y << endl;
  // cout << "Hjt = " << Hjt << endl;
  // cout << "S = " << S_ << endl;



}
