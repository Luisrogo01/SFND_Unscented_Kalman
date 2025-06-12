#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "ukf.h"
#define _USE_MATH_DEFINES
#include <cmath>


using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
*/

UKF::UKF() {
  
 
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //Nice work. Your std_a_ should be between 0 and 3 which is the case.
  //std_a_ = 30;
  std_a_ = 0.85; //.85

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30; Nice value for your std_yawdd_. Your stdyawdd ` should be between 0 and 1. which is the case currently.
  // It is low enough not to affect the results negatively.
  std_yawdd_ = 0.3; //.3
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
  */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
  */
  
  /*
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;

  n_x = 5; // State dimension
  n_aug = 7; // Augmented state dimension
  lambda = 3 - n_x; // Design parameter LAMBDA

  // The P matrix is the corresponding covariance matrix. Its diagonal elements tell us how certain we are in the state variables (means).
  // When we initialize P, we would like to describe our certainty on the state variables. We are relatively certain in the x and y values, 
  //because both measurement types give us some relatively accurate value (so we can set 1 for the first two diagonal elements).
  // Because the first measurement can be RADAR measurement and in case of radar measurement we do not really know the exact velocity vector (vx,vy), 
  // we can set a higher value (like 100 or 1000) for the 3rd and 4th diagonal elements.
  // However, it is often used that we set P as the identity matrix, and we hope that the Kalman Filter will converge relatively fast to the correct values.

  x_ = VectorXd::Zero(n_x); // INitial State vector
  //P = MatrixXd::Identity(n_x, n_x); // Initial of State COvariance matrix
  // The initial P value can also affect results. The values you have initialized for your P_ matrix are too high which causes the results to go out of range. 
  // change your values as below
  P = MatrixXd::Zero(n_x, n_x);
  /*P << 1,0,0,0,0,
        0,1,0,0,0,
        0,0,1.35,0,0,
        0,0,0,1.1,0,
        0,0,0,0,1.1; */
  P << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 0.0225, 0,
        0, 0, 0, 0, 0.0225;

        
  Xsig_pred = MatrixXd::Zero(n_x, 2 * n_aug + 1); // Initialization of Predicted sigma points matrix
  weights = VectorXd::Zero(2 * n_aug + 1);// VEctor of Sigma points weights
  previous_timestamp_ = 0;
}

//UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
  */
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_) {
    return;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_) {
    return;
  }

  // STEP 0: SET INITIAL VALUES - State and Covariance
  if (!is_initialized_) {
    /**
    * Initialize the state with the first measurement.
    */
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state
      auto ro = meas_package.raw_measurements_(0);
      auto phi = meas_package.raw_measurements_(1);
      x_ << ro * cos(phi), ro * sin(phi), 0, 0, 0; //  x_ << ro * cos(phi), ro * sin(phi), 1, 1, 1; 

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state
      // Set the state with the initial location and zero velocity
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0; 
      //x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 1, 1, 1;
    }
    /*
    if (measurement_package.sensor_type_ == MeasurementPackage::RADAR) {
            double rho = measurement_package.raw_measurements_(0);
            double phi = measurement_package.raw_measurements_(1);
            double rho_d = measurement_package.raw_measurements_(2);
            double p_x = rho * cos(phi);
            double p_y = rho * sin(phi);
            // To note: this is inaccurate value, aim to initialize velocity which's magnitude/order is close to real value
            double vx = rho_d * cos(phi);
            double vy = rho_d * sin(phi);
            double v = sqrt(vx * vx + vy * vy);
            x_ << p_x, p_y, v,0, 0;

        } else {
            x_ << measurement_package.raw_measurements_(0), measurement_package.raw_measurements_(1), 0., 0, 0;

        } */

    previous_timestamp_ = meas_package.timestamp_;
    // Done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  // TIMESTAMP
  //compute the time elapsed between the current and previous measurements, dt - expressed in seconds
  auto dt = static_cast<double>((meas_package.timestamp_ - previous_timestamp_) * 1e-6);

  previous_timestamp_ = meas_package.timestamp_;

  while (dt > 0.1) {
    constexpr double delta_t = 0.05;
    Prediction(delta_t);
    dt -= delta_t;
  }

  // STEP 1 and 2 : COMPUTE SIGMA POINTS and PREDICT STATE AND ERROR COVARIANCE
  Prediction(dt);

  // STEP 3 - PREDICT MEASUREMENT AND COVARIANCE
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    //set measurement dimension, radar can measure r, phi, and r_dot
    n_z = 3;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
   //set measurement dimension, laser can measure px, py
    n_z = 2;
  }

  //create matrix for sigma points in measurement space
  Zsig = MatrixXd::Zero(n_z, 2 * n_aug + 1);

  //mean predicted measurement
  z_pred = VectorXd::Zero(n_z);

  //measurement covariance matrix S
  S = MatrixXd::Zero(n_z, n_z);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    PredictRadarMeasurement();
    // update NIS
    //NIS_radar_ = (meas_package.raw_measurements_-z_pred_).transpose()*S_.inverse()*(meas_package.raw_measurements_-z_pred_);
    //std::cout << "NIS Radar = " << NIS_radar_ << std::endl;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    PredictLaserMeasurement();
    // update NIS
    //NIS_laser_ = (meas_package.raw_measurements_-z_pred_).transpose()*S_.inverse()*(meas_package.raw_measurements_-z_pred_);
    //std::cout << "NIS LiDAR = " << NIS_laser_ << std::endl;
  }
  
  // STEP 4, 5 and 6: COMPUTE KALMAN GAIN, THE STIMATE, AND ERRROR COVARIANCE 
  Update(meas_package.raw_measurements_);
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
  */

  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug, 2 * n_aug + 1);

  // STEP 1: Compute Sigma Points and Weights - In this case with Agumentation because noise is also not linear

  GenerateAugmentedSigmaPoints(Xsig_aug);
  SigmaPointPrediction(Xsig_aug, delta_t);

  // STEP 2: Predict State and error Covariance
  PredictMeanAndCovariance();
}

void UKF::Update(const VectorXd& z) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar / radar NIS, if desired.
   */

  // STEP 4 : Computing Kalman Gain.
  // create matrix for cross correlation Tc - On Book called as Pxz
  MatrixXd Tc = MatrixXd(n_x, n_z);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; ++i) {  // 2n+1 simga points

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  // STEP 5: Compute the stimate
  x_ = x_ + K * z_diff;
  // STEP 6: Compute the error covariance
  P = P - K*S*K.transpose();

  // print result
  //std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P << std::endl;
}


void UKF::GenerateAugmentedSigmaPoints(MatrixXd& Xsig_aug) {
 // std::cout << "Xsig_aug PRIMERO= " << std::endl << Xsig_aug << std::endl;
  // create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(n_aug);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug, n_aug);

  // create augmented mean state
  x_aug.head(n_x) = x_;

  //std::cout << "x_aug = " << std::endl << x_aug << std::endl;

  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x,n_x) = P;
  P_aug(n_x,n_x) = std_a_*std_a_;
  P_aug(n_x+1,n_x+1) = std_yawdd_*std_yawdd_;

  //std::cout << "P_aug = " << std::endl << P_aug << std::endl;

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();


  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i< n_aug; ++i) {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug) * L.col(i);
    Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * L.col(i);
  }

  // print result
  //std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
}


void UKF::SigmaPointPrediction(MatrixXd& Xsig_aug, double &delta_t) {

  //double delta_t = 0.1; // time diff in sec

  // predict sigma points
  for (int i = 0; i< 2*n_aug+1; ++i) {
    // extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  // print result
  //std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;
}

void UKF::PredictMeanAndCovariance() {

  // set weights
  weights(0) = lambda/(lambda+n_aug);

  for (int i=1; i<2*n_aug+1; ++i) {  // 2n+1 weights
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }

  // predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; ++i) {  // iterate over sigma points
    x_ = x_ + weights(i) * Xsig_pred.col(i);
  }

  // predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }

  // print result
  //std::cout << "Predicted state" << std::endl;
  //std::cout << x_ << std::endl;
  //std::cout << "Predicted covariance matrix" << std::endl;
  //std::cout << P << std::endl;
}

void UKF::PredictRadarMeasurement() {

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug + 1; ++i) {  // 2n+1 simga points
    // extract values for better readability
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       // r
    Zsig(1,i) = atan2(p_y,p_x);                                // phi
    Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   // r_dot
  }

  // mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug+1; ++i) {
    z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //  covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug+ 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;

  S = S + R;

  // print result
  //std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  //std::cout << "S: " << std::endl << S << std::endl;
}

void UKF::PredictLaserMeasurement() {
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 sigma points
        // measurement model
        Zsig(0,i) = Xsig_pred(0,i);           //px
        Zsig(1,i) = Xsig_pred(1,i);           //py
    }
  
    // measurement model
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {
        z_pred = z_pred + weights(i) * Zsig.col(i);
    }

    //  covariance matrix S
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 sigma points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
        while (z_diff(1) <-M_PI) z_diff(1) += 2.*M_PI;
        
        S = S + weights(i) * z_diff * z_diff.transpose();
    }

    // add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R = MatrixXd::Zero(2, 2);
    R<< std_laspx_ * std_laspx_,0,
        0, std_laspy_ * std_laspy_ ;

    S = S + R;
  
  // print result
  //std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  //std::cout << "S: " << std::endl << S << std::endl;
}
