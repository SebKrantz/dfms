#include <RcppArmadillo.h>
#include "helper.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// TODO: Ensure Symmetry???

//' Implementation of a Kalman filter
//' @param X Data matrix (T x n)
//' @param A Transition matrix (rp x rp)
//' @param C Observation matrix (n x rp)
//' @param Q State covariance (rp x rp)
//' @param R Observation covariance (n x n)
//' @param F0 Initial state vector (rp x 1)
//' @param P0 Initial state covariance (rp x rp)
// [[Rcpp::export]]
Rcpp::List KalmanFilter(arma::mat X, arma::mat A, arma::mat C, arma::mat Q,
                        arma::mat R, arma::colvec F0, arma::mat P0) {

  const int T = X.n_rows;
  const int n = X.n_cols;
  const int rp = A.n_rows;

  // In internal code factors are Z (instead of F) and factor covariance V (instead of P),
  // to avoid confusion between the matrices and their predicted (p) and filtered (f) states.
  // Additionally the results matrices for all time periods have a T in the name.

  double loglik = 0, dn = double(n), detS;
  colvec Zp = F0, Zf, et;
  mat K, Vp = P0, Vf, S, VCt;

  // Predicted state mean and covariance
  mat ZTp(T+1, rp, fill::zeros);
  cube VTp(rp, rp, T+1, fill::zeros);

  // Filtered state mean and covariance
  mat ZTf(T, rp, fill::zeros);
  cube VTf(rp, rp, T, fill::zeros);

  // Handling missing values in the filter
  mat Ci, Ri;
  uvec miss, nmiss = find_finite(A.row(0)), a(1);

  for (int i = 0; i < T; ++i) {

    // If missing observations are present at some timepoints, exclude the
    // appropriate matrix slices from the filtering procedure.
    miss = find_finite(X.row(i));
    Ci = C.submat(miss, nmiss);
    Ri = R.submat(miss, miss);
    a[0] = i;

    // Intermediate results
    VCt = Vp * Ci.t();
    S = (Ci * VCt + Ri).i();

    // Prediction error
    et = X.submat(a, miss).t() - Ci * Zp;
    // Kalman gain
    K = VCt * S;
    // Updated state estimate
    Zf = Zp + K * et;
    // Updated state covariance estimate
    Vf = Vp - K * Ci * Vp;

    // Compute likelihood. Skip this part if S is not positive definite.
    detS = det(S);
    if(detS > 0) {
      loglik += -0.5 * (dn * log(2.0 * datum::pi) - log(detS) +
        conv_to<double>::from(et.t() * S * et));
    }

    // Store predicted and filtered data needed for smoothing
    ZTp.row(i) = Zp.t();
    VTp.slice(i) = Vp;
    ZTf.row(i) = Zf.t();
    VTf.slice(i) = Vf;

    // Run a prediction
    Zp = A * ZTf.row(i).t();
    Vp = A * VTf.slice(i) * A.t() + Q;

  }

  return Rcpp::List::create(Rcpp::Named("F") = ZTf,
                            Rcpp::Named("P") = VTf,
                            Rcpp::Named("F_pred") = ZTp,
                            Rcpp::Named("P_pred") = VTp,
                            Rcpp::Named("loglik") = loglik);
}


// TODO: Futher Optimize Smoother

//' Runs a Kalman smoother
//' @param A Transition matrix (rp x rp)
//' @param C Observation matrix (n x rp)
//' @param R Observation covariance (n x n)
//' @param ZTf State estimates
//' @param ZTp State predicted estimates
//' @param VTf_v Variance estimates
//' @param VTp_v Predicted variance estimates
//' @return List of smoothed estimates
// [[Rcpp::export]]
Rcpp::List KalmanSmoother(arma::mat A, arma::mat C, arma::mat R,
                          arma::mat ZTf, arma::mat ZTp,
                          Rcpp::NumericVector VTf_v,
                          Rcpp::NumericVector VTp_v) {

  const int T = ZTf.n_rows;
  const int rp = A.n_rows;
  const int n = C.n_rows;

  cube VTf = array2cube(VTf_v);
  cube VTp = array2cube(VTp_v);

  cube J(rp, rp, T, fill::zeros);
  cube L(n, n, T, fill::zeros);
  cube Ks(rp, n, T, fill::zeros);

  // Smoothed state mean and covariance
  mat ZsT(T, rp, fill::zeros);
  cube VsT(rp, rp, T, fill::zeros);
  cube VVsT(rp, rp, T, fill::zeros);

  // Initialize smoothed data with last observation of filtered data
  ZsT.row(T-1) = ZTf.row(T-1);
  VsT.slice(T-1) = VTf.slice(T-1);

  for (int i = 0; i < T-1; ++i) {
    J.slice(i) = VTf.slice(i) * A.t() * VTp.slice(i+1).i();
  }

  // Smoothed state variable and covariance
  for (int j = 2; j < T+1; ++j) {

    ZsT.row(T-j) = ZTf.row(T-j) +
      (J.slice(T-j) * (ZsT.row(T-j+1) - ZTp.row(T-j+1)).t()).t();

    VsT.slice(T-j) = VTf.slice(T-j) +
      J.slice(T-j) * (VsT.slice(T-j+1) - VTp.slice(T-j+1)) * J.slice(T-j).t();

  }

  // Additional variables used in EM-algorithm
  for (int i = 0; i < T; ++i) {
    L.slice(i) = (C * VTp.slice(i) * C.t() + R).i();
    Ks.slice(i) = VTp.slice(i) * C.t() * L.slice(i);
  }

  VVsT.slice(T-1) = (eye(rp,rp) - Ks.slice(T-1) * C) * A * VTf.slice(T-2);

  for (int j = 2; j < T-1; ++j) {
    VVsT.slice(T-j) = VTf.slice(T-j) * J.slice(T-j-1).t() + J.slice(T-j)
    * (VVsT.slice(T-j+1) - A * VTf.slice(T-j))
    * J.slice(T-j-1).t();
  }

  return Rcpp::List::create(Rcpp::Named("Fs") = ZsT,
                            Rcpp::Named("Ps") = VsT,
                            Rcpp::Named("PPs") = VVsT);
}



//' Kalman Filter and Smoother
//' @param X Data matrix (T x n)
//' @param A Transition matrix (rp x rp)
//' @param C Observation matrix (n x rp)
//' @param Q State covariance (rp x rp)
//' @param R Observation covariance (n x n)
//' @param F0 Initial state vector (rp x 1)
//' @param P0 Initial state covariance (rp x rp)
// [[Rcpp::export]]
Rcpp::List KalmanFilterSmoother(arma::mat X, arma::mat A, arma::mat C, arma::mat Q,
                                arma::mat R, arma::colvec F0, arma::mat P0) {

  const int T = X.n_rows;
  const int n = X.n_cols;
  const int rp = A.n_rows;

  // In internal code factors are Z (instead of F) and factor covariance V (instead of P),
  // to avoid confusion between the matrices and their predicted (p) and filtered (f) states.
  // Additionally the results matrices for all time periods have a T in the name.

  double loglik = 0, dn = double(n), detS;
  colvec Zp = F0, Zf, et;
  mat K, Vp = P0, Vf, S, VCt;

  // Predicted state mean and covariance
  mat ZTp(T+1, rp, fill::zeros);
  cube VTp(rp, rp, T+1, fill::zeros);

  // Filtered state mean and covariance
  mat ZTf(T, rp, fill::zeros);
  cube VTf(rp, rp, T, fill::zeros);

  // Handling missing values in the filter
  mat Ci, Ri;
  uvec miss, nmiss = find_finite(A.row(0)), a(1);

  for (int i = 0; i < T; ++i) {

    // If missing observations are present at some timepoints, exclude the
    // appropriate matrix slices from the filtering procedure.
    miss = find_finite(X.row(i));
    Ci = C.submat(miss, nmiss);
    Ri = R.submat(miss, miss);
    a[0] = i;

    // Intermediate results
    VCt = Vp * Ci.t();
    S = (Ci * VCt + Ri).i();

    // Prediction error
    et = X.submat(a, miss).t() - Ci * Zp;
    // Kalman gain
    K = VCt * S;
    // Updated state estimate
    Zf = Zp + K * et;
    // Updated state covariance estimate
    Vf = Vp - K * Ci * Vp;

    // Compute likelihood. Skip this part if S is not positive definite.
    detS = det(S);
    if(detS > 0) {
      loglik += -0.5 * (dn * log(2.0 * datum::pi) - log(detS) +
        conv_to<double>::from(et.t() * S * et));
    }

    // Store predicted and filtered data needed for smoothing
    ZTp.row(i) = Zp.t();
    VTp.slice(i) = Vp;
    ZTf.row(i) = Zf.t();
    VTf.slice(i) = Vf;

    // Run a prediction
    Zp = A * ZTf.row(i).t();
    Vp = A * VTf.slice(i) * A.t() + Q;

  }

  // Kalman Smoother
  cube J(rp, rp, T, fill::zeros);
  cube L(n, n, T, fill::zeros);
  cube Ks(rp, n, T, fill::zeros);

  // Smoothed state mean and covariance
  mat ZsT(T, rp, fill::zeros);
  cube VsT(rp, rp, T, fill::zeros);
  cube VVsT(rp, rp, T, fill::zeros);

  // Initialize smoothed data with last observation of filtered data
  ZsT.row(T-1) = ZTf.row(T-1);
  VsT.slice(T-1) = VTf.slice(T-1);

  for (int i = 0; i < T-1; ++i) {
    J.slice(i) = VTf.slice(i) * A.t() * VTp.slice(i+1).i();
  }

  // Smoothed state variable and covariance
  for (int j = 2; j < T+1; ++j) {

    ZsT.row(T-j) = ZTf.row(T-j) +
      (J.slice(T-j) * (ZsT.row(T-j+1) - ZTp.row(T-j+1)).t()).t();

    VsT.slice(T-j) = VTf.slice(T-j) +
      J.slice(T-j) * (VsT.slice(T-j+1) - VTp.slice(T-j+1)) * J.slice(T-j).t();

  }

  // Additional variables used in EM-algorithm
  for (int i = 0; i < T; ++i) {
    L.slice(i) = (C * VTp.slice(i) * C.t() + R).i();
    Ks.slice(i) = VTp.slice(i) * C.t() * L.slice(i);
  }

  VVsT.slice(T-1) = (eye(rp,rp) - Ks.slice(T-1) * C) * A * VTf.slice(T-2);

  for (int j = 2; j < T-1; ++j) {
    VVsT.slice(T-j) = VTf.slice(T-j) * J.slice(T-j-1).t() + J.slice(T-j)
    * (VVsT.slice(T-j+1) - A * VTf.slice(T-j))
    * J.slice(T-j-1).t();
  }

  return Rcpp::List::create(Rcpp::Named("Fs") = ZsT,
                            Rcpp::Named("Ps") = VsT,
                            Rcpp::Named("PPs") = VVsT,
                            Rcpp::Named("loglik") = loglik);
}
