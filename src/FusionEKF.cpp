#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
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

  //measurement covariance matrix - laser is 2*2 as there are 2 inputs
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar is 3*3 as there are 3 inputs and its covariance matrix should be 3*3
  R_radar_ << 0.09, 0, 0,
        0, 0.0006, 0,
        0, 0, 0.09;
  
  //Initializing H vector for LIDAR, For Radar it will be a 3*4 matrix
  
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;


  //Initializing The Covariance matrix of X			  
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // Initializing the X vector for first measurement
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "EKF : First measurement RADAR" << endl;
      double rho = measurement_pack.raw_measurements_[0]; // range
  	  double phi = measurement_pack.raw_measurements_[1]; // bearing
  	  double rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho
      
	  //Converting polar to cartesian coordinates
  	  double x = rho * cos(phi);
  	  double vx = rho_dot * cos(phi);
  	  double vy = rho_dot * sin(phi);
  	  double y = rho * sin(phi);
	  
      if ( x < 0.0001 ) {
        x = 0.0001;
      }

      if ( y < 0.0001 ) {
        y = 0.0001;
      }
	  
      ekf_.x_ << x, y, 0 , 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      cout << "EKF : First measurement LASER" << endl;
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // Saving first timestamp in seconds
    previous_timestamp_ = measurement_pack.timestamp_ ;
    
    //The first step is done, So no need to call this function again
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
   previous_timestamp_ = measurement_pack.timestamp_;

  ///The F vector with Time stamp
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;

  double noise_ax = 9.0;
  double noise_ay = 9.0;

  double dt_2 = dt * dt; 
  double dt_3 = dt_2 * dt; 
  double dt_4 = dt_3 * dt; 

  //Calculating the value fo Q matrix
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << (dt_4 / 4) * noise_ax, 0, (dt_3 / 2) * noise_ax, 0,
	         0, (dt_4 / 4) * noise_ay, 0, (dt_3 / 2) * noise_ay,
	         (dt_3 / 2) * noise_ax, 0, dt_2 * noise_ax, 0,
 	         0, (dt_3 / 2) * noise_ay, 0, dt_2 * noise_ay;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  
  //Selecting the vector H and R based on radar or lisar  

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
  	ekf_.R_ = R_radar_;
  	ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    ekf_.H_ = H_laser_;
  	ekf_.R_ = R_laser_;
  	ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
