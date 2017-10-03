#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
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
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;
  Xsig_pred_ = MatrixXd(5, 2*n_aug_+1);

  //P_ = MatrixXd::Identity(5, 5);
   P_ << 0.15,    0, 0, 0, 0,
               0, 0.15, 0, 0, 0,
               0,    0, 1, 0, 0,
               0,    0, 0, 1, 0,
               0,    0, 0, 0, 1;

  //create vector for weights
  weights_ = VectorXd(2*n_aug_+1);

  //set weights
  weights_(0) = lambda_/(lambda_+n_aug_);
  for(int i=1;i<2*n_aug_+1;i++){
      weights_(i) = 0.5/(lambda_ + n_aug_);
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    x_ = VectorXd(5);
    //ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float ro_dot = measurement_pack.raw_measurements_[2];

      x_ << ro*cos(phi), ro*sin(phi), 0, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0, 0;
    }

    if (fabs(x_(0)) < 0.0001 and fabs(x_(1)) < 0.0001){
      x_(0) = 0.0001;
      x_(1) = 0.0001;
    }



    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

 float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  
  Prediction(dt);

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(measurement_pack);
  }else{
    UpdateLidar(measurement_pack);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /*MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();
  MatrixXd first = (sqrt(lambda_+n_x_)*A).colwise() + x_;
  MatrixXd second = (-sqrt(lambda_+n_x_)*A).colwise() + x_;

  Xsig << x_, first, second;

*/

  cout << "Previous x_ ";
  cout << x_ << endl;
  cout << "Previous P_ " << P_ << endl << endl << endl;

  MatrixXd Xsig_aug = GenerateSigmaPoints();
 

  cout << "Xsig_aug " << Xsig_aug << endl; 

  //predict sigma points
  PredictSigmaPoints(Xsig_aug, delta_t);
  

   
  cout << "Xsig_pred_: " << Xsig_pred_ << endl;
  

 
  //std::cout << weights << std::endl;
  
  //predict state mean
  PredictMeanAndCovariance();
  
  

}


MatrixXd UKF::GenerateSigmaPoints(){
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  //x_aug << x, std_a, std_yawdd;
  x_aug << x_, 0, 0;
  //std::cout << x_aug << endl;
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  MatrixXd tmp(2,2);
  tmp << std_a_*std_a_, 0, 
         0, std_yawdd_*std_yawdd_;
  P_aug.bottomRightCorner(2,2) = tmp;
  //cout << P_aug << endl;
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  //create augmented sigma points
 
  // lambda = 3-nx 
  MatrixXd first = (sqrt(3)*A).colwise() + x_aug;
  MatrixXd second = (-sqrt(3)*A).colwise() + x_aug;
  
  Xsig_aug << x_aug, first, second;

  /*Xsig_aug.fill(0.0);
  Xsig_aug.col(0) = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    cout << "In for loop: " << endl;
    Xsig_aug.col(i + 1) = x_ + sqrt(3) * A.col(i);
    Xsig_aug.col(i + 1 + n_x_) = x_ - sqrt(3) * A.col(i);
  }*/


  return Xsig_aug;
}

void UKF::PredictSigmaPoints(MatrixXd& Xsig_aug, float delta_t){
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

void normalize(VectorXd& v){
  for(int i=0;i<v.size();i++){
    if(fabs(v(i)) < 0.0001){
      v(i) = 0.0001;
    }
  }
}

void UKF::PredictMeanAndCovariance(){
  //cout << "Weights: " << weights_ << endl; 
  x_.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
    VectorXd weighted_x = weights_(i)*Xsig_pred_.col(i);
    //cout << "Weighted x: " << weighted_x << endl;
    x_ = x_ + weighted_x;
    //normalize(x_);
    //cout << "Updated x: " << x_ << endl; 
  }


  cout << "After normalizing predicted x_:  " << x_ << endl;

  //predict state covariance matrix
  P_.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
    VectorXd x_diff = (Xsig_pred_.col(i) - x_);
    //normalize(x_diff);
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
    P_ += weights_(i)*x_diff* (x_diff.transpose()); 
  }

}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //extract measurement as VectorXd
  //extract measurement as VectorXd
  int n_z = 3;

  float rho, theta, rhod;
  float x1 = x_(0);
  float y1 = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  VectorXd y;
  VectorXd Hx;
  MatrixXd K, S;

  VectorXd z;
  MatrixXd R = MatrixXd(5, 5);
  // identity 4 matrix
  MatrixXd I = MatrixXd::Identity(5, 5);
  MatrixXd H = MatrixXd(2, 5);
    H << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0;
  Hx = H*x_;
  z = meas_package.raw_measurements_;

    R = MatrixXd(2, 2);
     R << 0.022, 0,
   0, 0.022;
  // updating the state
  y = z - Hx;

  // normalizing the angle
  y(1) = std::atan2( sin( y(1)), cos(y(1)));

  S = H * P_ * H.transpose() + R;
  K = P_ * H.transpose() * S.inverse();

  //new state
  x_ = x_ + (K * y);
  P_ = (I - K * H)*P_;


}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  int n_z = 3;
   //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  

  //transform sigma points into measurement space
  for(int i=0;i<2*n_aug_+1;i++){
    VectorXd curr = Xsig_pred_.col(i);
    float px = curr(0);
    float py = curr(1);
    float v = curr(2);
    float yaw = curr(3);
    float yaw_dot = curr(4);
    
    VectorXd z(3);
    float ro = pow(px*px + py*py, 0.5);
    z << ro, atan2(py, px), (px*cos(yaw)*v + py*sin(yaw)*v)/ro;
      
    Zsig.col(i) = z;
     
    z_pred += weights_(i)*z;
  }

  cout << "Predicting z done" << endl;
  //calculate mean predicted measurement

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  //calculate measurement covariance matrix S
  for(int i=0;i<2*n_aug_+1;i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    S += weights_(i)*z_diff*(z_diff.transpose());
  }

  cout << "Calculated S " << endl;

  MatrixXd R(3,3);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
       
  cout << "Calculated R " << endl;
  S = S + R;

  float ro = measurement_pack.raw_measurements_[0];
  float phi = measurement_pack.raw_measurements_[1];
  float ro_dot = measurement_pack.raw_measurements_[2];

  VectorXd z(3);
  z << ro, phi, ro_dot;

  //create matrix for cross correlation Tc
  
  MatrixXd Tc(n_x_, n_z);

  Tc.fill(0.0);
  //calculate cross correlation matrix
  for(int i=0;i<2 * n_aug_ + 1;i++){
    VectorXd x_diff_mat(5);
    x_diff_mat << (Xsig_pred_.col(i) - x_);
    cout << "Xsig_pred_: " << Xsig_pred_.col(i) <<endl;
    cout <<"x_:  " << x_ << endl;
    //normalize(x_diff_mat);
    //angle normalization
    while (x_diff_mat(3)> M_PI) x_diff_mat(3)-=2.*M_PI;
    while (x_diff_mat(3)<-M_PI) x_diff_mat(3)+=2.*M_PI;
    cout << "x_diff: " << x_diff_mat << endl;

    VectorXd z_diff_mat(3);
    z_diff_mat << (Zsig.col(i)-z_pred);
    cout << "z_diff: " << z_diff_mat << endl;
    //normalize(z_diff_mat);
     //angle normalization
    while (z_diff_mat(1)> M_PI) z_diff_mat(1)-=2.*M_PI;
    while (z_diff_mat(1)<-M_PI) z_diff_mat(1)+=2.*M_PI;

    //Tc += weights_(i)*x_diff_mat*(z_diff_mat.transpose());
    Tc += weights_(i) * x_diff_mat*(z_diff_mat.transpose());
    //cout << "Tc: " << Tc << endl;
  }

  cout << "Calculated Tc " << endl;
  cout << "Tc: " << Tc << endl;
  cout << "S: " << S << endl;
  cout << "S.inverse() " << S.inverse() << endl;
 
  //calculate Kalman gain K;
  MatrixXd K = Tc * (S.inverse());
  //update state mean and covariance matrix

  cout << "Calculated kalman gain " << endl;
  cout << "x_ before: " << x_ << endl;
  cout << "z: " << z << endl;
  cout << "z_pred: " << z_pred << endl;
  cout << "K: " << K << endl;
  x_ = x_ + K*(z - z_pred);

  cout << "x_ after: " << x_  << endl;

  P_ = P_ - K*S*(K.transpose());

  cout << "P_ after: " << P_ << endl;

  cout << "Updated parameters " << endl;

}
