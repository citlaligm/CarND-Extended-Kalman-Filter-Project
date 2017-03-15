#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>       /* cos */
#define PI 3.14159265

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


  /**
    * Finish initializing the FusionEKF.
  */
  //sigma2=0.0225 from lectures
  R_laser_<< 0.0225, 0,
		     0, 0.0225;

  //Measurement Function Matrix
  H_laser_<< 1, 0, 0, 0,
  			 0, 1, 0, 0;

  //Given for the project
  R_radar_<< 0.09, 0, 0,
 		     0, 0.0009, 0,
			 0, 0, 0.09;
  //Measurement Function Matrix
  //TODO: Check if this is a good default.
  Hj_<< 1, 1, 0, 0,
        1, 1, 0, 0,
        1, 1, 1, 1;

  //State Transition Matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
  			 0, 1, 0, 1,
  			 0, 0, 1, 0,
  			 0, 0, 0, 1;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  int noise_ax = 7;
  int noise_ay = 7;
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    ekf_.x_ = VectorXd(4);
    previous_timestamp_ = measurement_pack.timestamp_;

    //Initialize the state ekf_.x_ with the first measurement.
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float theta = measurement_pack.raw_measurements_[1];

      float px = measurement_pack.raw_measurements_[0]*cos(theta);
      float py = measurement_pack.raw_measurements_[0]*sin(theta);
      cout<<"Px: "<<px<<"Py: "<<py<<endl;

      if(px==0 or py ==0){
        cout<<"Error!! "<<py<<endl;
        return;
      }

      ekf_.x_ << px, py, 0, 0;


    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      float px_laser = measurement_pack.raw_measurements_[0];
      float py_laser = measurement_pack.raw_measurements_[0];

      if(px_laser==0 or py_laser ==0) return;
      ekf_.x_ << px_laser, py_laser, 0, 0;

    }

    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
        		 0, 1, 0, 0,
        		 0, 0, 1000, 0,
        		 0, 0, 0, 1000;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;


  if(dt > 0.001){
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;



    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
                 0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
                 dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
                 0, dt_3/2*noise_ay, 0, dt_2*noise_ay;


    ekf_.Predict();




  }


  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	  // Radar updates

    MatrixXd z_predicted (3,1);
	  Hj_ << tools.CalculateJacobian(ekf_.x_);
	  ekf_.H_ = Hj_;
	  ekf_.R_ = R_radar_;

	  //Convert from cartesian to polar notation prior calling the updating function
	  float px = ekf_.x_[0];
    float py = ekf_.x_[1];
    float vx = ekf_.x_[2];
    float vy = ekf_.x_[3];

    //h(x')
	  float range = sqrt(px*px + py*py);
	  float bearing = atan(py/px);
	  float radial_velocity = (px*vx + py*vy)/range;
	  z_predicted << range, bearing, radial_velocity;


	  ekf_.UpdateEKF(measurement_pack.raw_measurements_, z_predicted);

  }

  else {
    // Laser updates
	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;

	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = \n" << ekf_.x_ << endl;
  cout << "P_ = \n" << ekf_.P_ << endl;
}
