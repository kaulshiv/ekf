#include "tools.h"
#include <iostream>
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth){ 

   VectorXd rmse(4);
   rmse << 0,0,0,0;

   // check the validity of the following inputs:
   //  * the estimation vector size should not be zero
   //  * the estimation vector size should equal ground truth vector size
   if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
      std::cout << "Invalid estimation or ground_truth data" << std::endl;
      return rmse;
   }

   // accumulate squared residuals
   for (unsigned int i=0; i < estimations.size(); ++i) {

      VectorXd residual = estimations[i] - ground_truth[i];

      // coefficient-wise multiplication
      residual = residual.array()*residual.array();
      rmse += residual;
   }

   // calculate the mean
   rmse = rmse/estimations.size();

   // calculate the squared root
   rmse = rmse.array().sqrt();
   
   std::cout << "rmse" << rmse << std::endl;
   // return the result
   return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
   MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // pre-compute a set of terms to avoid repeated calculation
  float c = px*px+py*py;

  // check division by zero
  if (fabs(c) < 0.0001) {
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return Hj;
  }

  Hj << px/sqrt(c), py/sqrt(c), 0, 0,
        -py/c, px/c, 0, 0,
        py*(vx*py - vy*px)/pow(c, 1.5), px*(vy*px - vx*py)/pow(c, 1.5), px/sqrt(c), py/sqrt(c);

  return Hj;

}
