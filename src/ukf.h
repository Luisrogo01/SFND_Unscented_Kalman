#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"
#include <iostream>


class UKF {
 public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF() = default;

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser / radar measurement
   * @param meas_package The measurement at k+1
   */
  void Update(const Eigen::VectorXd& z);


  void GenerateAugmentedSigmaPoints(Eigen::MatrixXd& Xsig_out);

  void SigmaPointPrediction(Eigen::MatrixXd& Xsig_aug, double& delta_t);

  void PredictMeanAndCovariance();

  void PredictRadarMeasurement();

  void PredictLaserMeasurement();

  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Weights of sigma points
  Eigen::VectorXd weights;

  // State dimension
  int n_x;

  // Augmented state dimension
  int n_aug;

  // Sigma point spreading parameter
  double lambda;

  long previous_timestamp_;
  
  // create matrix for sigma points in measurement space
  Eigen::MatrixXd Zsig ;

  // mean predicted measurement
  Eigen::VectorXd z_pred ;
  
  // measurement covariance matrix S
  Eigen::MatrixXd S;

  int n_z; // measurment dimenssion

};

#endif  // UKF_H