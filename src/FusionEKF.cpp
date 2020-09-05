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

#define sind(x) (sin(fmod((x),360) * M_PI / 180))
#define cosd(x) (cos(fmod((x),360) * M_PI / 180))

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

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

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
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    
    MatrixXd H_in;
    MatrixXd R_in;
    VectorXd x_in(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
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
      // TODO: Initialize state.

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
            0, 0, 1000, 0,
            0, 0, 0, 1000;

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

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  float dt =  (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  // cout << "dt : " << dt << endl; 
  // cout << "ts : " << measurement_pack.timestamp_ << endl; 


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;
  ekf_.Q_ << pow(dt, 4)/4*9, 0, pow(dt, 3)/2*9, 0,
             0, pow(dt, 4)/4*9, 0, pow(dt, 3)/2*9,
             pow(dt, 3)/2*9, 0, pow(dt, 2)*9, 0,
             0, pow(dt, 3)/2*9, 0, pow(dt, 2)*9;

  cout << "F : " << ekf_.F_<< endl; 
  cout << "Q : " << ekf_.Q_  << endl; 

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    cout << "RADAR" << endl;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    float ro = measurement_pack.raw_measurements_(0);
    float theta = measurement_pack.raw_measurements_(1);
    float ro_dot = measurement_pack.raw_measurements_(2);
    // cout << "ro = " << ro << endl;
    // cout << "theta = " << theta << endl;
    // cout << "ro_dot = " << ro_dot << endl;
    VectorXd z = VectorXd(3);
    z << ro, theta, ro_dot;
    ekf_.UpdateEKF(z);

  } else {
    // TODO: Laser updates
    cout << "LIDAR" << endl;
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
