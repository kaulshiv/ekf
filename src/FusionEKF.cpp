#include "FusionEKF.h"
#include <iostream>
#include <math.h>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // set  Hj_ when calculating the Jacobian in the Update Step for Radar
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    
    MatrixXd H_in;
    MatrixXd R_in;
    VectorXd x_in(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float ro = measurement_pack.raw_measurements_(0);
      float theta = measurement_pack.raw_measurements_(1);

      x_in(0) = ro*cos(theta);
      x_in(1) = ro*sin(theta);
      x_in(2) = 0;
      x_in(3) = 0;
      H_in = Hj_;
      R_in = R_radar_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x_in(0) = measurement_pack.raw_measurements_(0);
      x_in(1) = measurement_pack.raw_measurements_(1);
      x_in(2) = 0;
      x_in(3) = 0;
      H_in = H_laser_;
      R_in = R_laser_;
    }


    MatrixXd P_in = MatrixXd(4, 4);
    P_in << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 100, 0,
            0, 0, 0, 100;

    MatrixXd F_in = MatrixXd::Identity(4, 4);
    MatrixXd Q_in = MatrixXd::Identity(4, 4);

    ekf_.Init(x_in, P_in, F_in, H_in, R_in, Q_in);


    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  /**
   * Prediction
   */
  
  float dt =  (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;

  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;
  ekf_.Q_ << pow(dt, 4)/4*9, 0, pow(dt, 3)/2*9, 0,
             0, pow(dt, 4)/4*9, 0, pow(dt, 3)/2*9,
             pow(dt, 3)/2*9, 0, pow(dt, 2)*9, 0,
             0, pow(dt, 3)/2*9, 0, pow(dt, 2)*9; 

  ekf_.Predict();

  /**
   * Update
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    float ro = measurement_pack.raw_measurements_(0);
    float theta = measurement_pack.raw_measurements_(1);
    float ro_dot = measurement_pack.raw_measurements_(2);
    VectorXd z = VectorXd(3);
    z << ro, theta, ro_dot;
    ekf_.UpdateEKF(z);
  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    VectorXd z = VectorXd(2);
    float px = measurement_pack.raw_measurements_(0);
    float py = measurement_pack.raw_measurements_(1);
    z << px, py;
    ekf_.Update(z);
  }

  previous_timestamp_ = measurement_pack.timestamp_;

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
