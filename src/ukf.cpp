#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

///* row index of CTRV Model State
#define STATE_PX_IDX        0
#define STATE_PY_IDX        1
#define STATE_V_IDX         2
#define STATE_YAW_IDX       3


///* row index of sigma points matrix
#define AUG_SIGMA_V_IDX     2
#define AUG_SIGMA_YAW_IDX   3
#define AUG_SIGMA_YAW_DOT_IDX 4
#define AUG_SIGMA_MU_A_IDX  5
#define AUG_SIGMA_MU_YAW_DD_IDX 6


///* row index of Radar measurement vector
#define RADAR_RHO_IDX       0
#define RADAR_PHI_IDX       1
#define RADAR_RHO_DOT_IDX   2

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_ << 0, 0, 5.19, 0.001, 0.01;

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_.setZero();

  // create sigma point matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.setZero();

  weights_ = VectorXd(2 * n_aug_ + 1);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    weights_(i) = 0.5 / (n_aug_ + lambda_);
  }

  time_us_ = 0;
  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_) {
    time_us_ = meas_package.timestamp_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      x_[0] = rho * cos(phi);
      x_[1] = rho * sin(phi);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_[0] = meas_package.raw_measurements_[0];
      x_[1] = meas_package.raw_measurements_[1];
    }

    is_initialized_ = true;
    return;
  }

  if ((!use_laser_) &&
      (meas_package.sensor_type_ == MeasurementPackage::LASER)) {
    return;
  }
  if ((!use_radar_) &&
      (meas_package.sensor_type_ == MeasurementPackage::RADAR)) {
    return;
  }

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  double delta_t_2 = delta_t * delta_t;
  int aug_cols = 2 * n_aug_ + 1;

  /**
   * Step 1. generateAugmentedSigmaPoints
   */
  MatrixXd Xsig_aug = MatrixXd(n_aug_, aug_cols);

  ///* augmented state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate noise1 noise2]
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.setZero();
  x_aug.head(n_x_) = x_;

  ///* create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.setZero();
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_aug_-1, n_aug_-1) = std_yawdd_ * std_yawdd_;
  P_aug.topLeftCorner(n_x_, n_x_) = P_;

  ///* create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();
  ///* create augmented sigma points
  double factor = sqrt(lambda_ + n_aug_);
  MatrixXd P_aug_sqr = factor * A_aug;

  MatrixXd second = P_aug_sqr.colwise() + x_aug;
  MatrixXd third = (-P_aug_sqr).colwise() + x_aug;

  Xsig_aug << x_aug, second, third;


  /**
   * Step 2. SigmaPointPrediction
   */
  MatrixXd v_k = Xsig_aug.row(AUG_SIGMA_V_IDX);
  MatrixXd yaw_k = Xsig_aug.row(AUG_SIGMA_YAW_IDX);
  MatrixXd yaw_dot_k = Xsig_aug.row(AUG_SIGMA_YAW_DOT_IDX);
  MatrixXd mu_a_k = Xsig_aug.row(AUG_SIGMA_MU_A_IDX);
  MatrixXd mu_yaw_dot_dot_k = Xsig_aug.row(AUG_SIGMA_MU_YAW_DD_IDX);

  MatrixXd noise = MatrixXd(n_x_, aug_cols);
  noise << (mu_a_k.array() * yaw_k.array().cos()) * delta_t_2 / 2,
      (mu_a_k.array() * yaw_k.array().sin()) * delta_t_2 / 2,
      mu_a_k.array() * delta_t,
      mu_yaw_dot_dot_k.array() * delta_t_2 / 2,
      mu_yaw_dot_dot_k.array() * delta_t;

  MatrixXd xk_1 = MatrixXd(n_x_, aug_cols);
  ///* avoid division by zero
  for (int i = 0; i < aug_cols; i++) {
    if (fabs(yaw_dot_k(0, i)) <= 0.001) {
      xk_1.col(i) << v_k(0, i) * cos(yaw_k(0, i)) * delta_t,
          v_k(0, i) * sin(yaw_k(0, i)) * delta_t,
          0,
          yaw_dot_k(0, i) * delta_t,
          0;
    } else {
      xk_1.col(i) << (sin(yaw_k(0, i) + yaw_dot_k(0, i) * delta_t) - sin(yaw_k(0, i))) * v_k(0, i) / yaw_dot_k(0, i),
          (-cos(yaw_k(0, i) + yaw_dot_k(0, i) * delta_t) + cos(yaw_k(0, i))) * v_k(0, i) / yaw_dot_k(0, i),
          0,
          yaw_dot_k(0, i) * delta_t,
          0;
    }
  }

  ///* write predicted sigma points into right column
  Xsig_pred_ = Xsig_aug.block(0, 0, n_x_, aug_cols) + xk_1 + noise;


  /**
   * Step 3. PredictMeanAndCovariance
   */
  x_ = Xsig_pred_ * weights_;

  MatrixXd x_diff = Xsig_pred_.colwise() - x_;

  ///* Normalize yaw angle
  for (int i = 0; i < aug_cols; i++) {
    while (x_diff(STATE_YAW_IDX, i) > M_PI) {
      x_diff(STATE_YAW_IDX, i) -= 2. * M_PI;
    }
    while (x_diff(STATE_YAW_IDX, i) < -M_PI) {
      x_diff(STATE_YAW_IDX, i) += 2. * M_PI;
    }
  }

  MatrixXd x_diff_t = x_diff.transpose();
  x_diff = x_diff.array().rowwise() * weights_.transpose().array();
  ///* Predicted covariance matrix at k+1.
  P_ = x_diff * x_diff_t;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  int n_meas_laser = 2;

  MatrixXd H = MatrixXd(n_meas_laser, n_x_);
  H << 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0;
  MatrixXd Ht = H.transpose();
  MatrixXd PHt = P_ * Ht;

  ///* Predict Laser measurement
  VectorXd z_pred = H * x_;

  MatrixXd R = MatrixXd(n_meas_laser, n_meas_laser);
  R << std_laspx_ * std_laspx_, 0,
      0, std_laspy_ * std_laspy_;
  MatrixXd S = H * PHt + R;
  MatrixXd Si = S.inverse();

  ///* calculate Kalman gain K;
  MatrixXd K = PHt * Si;

  ///* Incomming Laser measurement
  VectorXd z = meas_package.raw_measurements_;

  ///* update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  x_ = x_ + (K * z_diff);
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K * H) * P_;

  double NIS = z_diff.transpose() * Si * z_diff;
  cout << NIS << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  int n_meas_radar = 3;
  int aug_cols = 2 * n_aug_ + 1;
  /**
   * Step 1. PredictRadarMeasurement
   */
  ///* create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_meas_radar, aug_cols);

  ///* mean predicted measurement
  VectorXd z_pred = VectorXd(n_meas_radar);

  ///* transform sigma points into measurement space
  VectorXd px = Xsig_pred_.row(STATE_PX_IDX);
  VectorXd py = Xsig_pred_.row(STATE_PY_IDX);
  VectorXd v = Xsig_pred_.row(STATE_V_IDX);
  VectorXd yaw = Xsig_pred_.row(STATE_YAW_IDX);
  VectorXd c1 = px.array() * px.array() + py.array() * py.array();
  VectorXd rho = c1.array().sqrt();
  VectorXd v1 = px.array() * yaw.array().cos() * v.array();
  VectorXd v2 = py.array() * yaw.array().sin() * v.array();
  VectorXd rho_dot = v1.array() + v2.array();
  rho_dot = rho_dot.array() / rho.array();

  Zsig.row(RADAR_RHO_IDX) = rho.transpose();
  Zsig.row(RADAR_RHO_DOT_IDX) = rho_dot.transpose();

  for (int i = 0; i < aug_cols; i++) {
    Zsig(RADAR_PHI_IDX, i) = atan2(py(i), px(i));
  }

  ///* calculate mean predicted measurement
  z_pred = Zsig * weights_;

  ///* calculate innovation covariance matrix S
  MatrixXd Zsig_diff = Zsig.colwise() - z_pred;
  ///* Normalize phi angle
  for (int i = 0; i < aug_cols; i++) {
    while (Zsig_diff(RADAR_PHI_IDX, i) > M_PI) {
      Zsig_diff(RADAR_PHI_IDX, i) -= 2. * M_PI;
    }
    while (Zsig_diff(RADAR_PHI_IDX, i) < -M_PI) {
      Zsig_diff(RADAR_PHI_IDX, i) += 2. * M_PI;
    }
  }

  MatrixXd Zsig_diff_t = Zsig_diff.transpose();
  Zsig_diff = Zsig_diff.array().rowwise() * weights_.transpose().array();
  ///* measurement covariance matrix S
  MatrixXd S = Zsig_diff * Zsig_diff_t;

  MatrixXd R = MatrixXd(n_meas_radar, n_meas_radar);
  R << std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;
  S += R;
  MatrixXd Si = S.inverse();

  /**
   * Step 2. UpdateState
   */
  ///* Incomming Radar measurement
  VectorXd z = meas_package.raw_measurements_;

  ///* calculate cross correlation matrix
  MatrixXd Tx = Xsig_pred_.colwise() - x_;
  MatrixXd Tz = Zsig.colwise() - z_pred;
  ///* angle normalization
  for (int i = 0; i < aug_cols; i++) {
    while (Tx(STATE_YAW_IDX, i) > M_PI) {
      Tx(STATE_YAW_IDX, i) -= 2. * M_PI;
    }
    while (Tx(STATE_YAW_IDX, i) < -M_PI) {
      Tx(STATE_YAW_IDX, i) += 2. * M_PI;
    }

    while (Tz(RADAR_PHI_IDX, i) > M_PI) {
      Tz(RADAR_PHI_IDX, i) -= 2. * M_PI;
    }
    while (Tz(RADAR_PHI_IDX, i) < -M_PI) {
      Tz(RADAR_PHI_IDX, i) += 2. * M_PI;
    }
  }

  Tx = Tx.array().rowwise() * weights_.transpose().array();
  ///* create matrix for cross correlation
  MatrixXd Tk_1 = Tx * Tz.transpose();

  ///* calculate Kalman gain K;
  MatrixXd K = Tk_1 * Si;

  VectorXd z_diff = z - z_pred;
  while (z_diff[RADAR_PHI_IDX] > M_PI) {
    z_diff[RADAR_PHI_IDX] -= 2. * M_PI;
  }
  while (z_diff[RADAR_PHI_IDX] < -M_PI) {
    z_diff[RADAR_PHI_IDX] += 2. * M_PI;
  }

  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  double NIS = z_diff.transpose() * Si * z_diff;
  cout << NIS << ",";
}
