#include <RcppArmadillo.h>
#include "helper.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


//' Implementation of a Kalman filter
//' @param y Data matrix (T x n)
//' @param H Observation matrix
//' @param Q State covariance
//' @param R Observation covariance
//' @param F Transition matrix
//' @param x0 Initial state vector
//' @param P0 Initial state covariance
// [[Rcpp::export]]
Rcpp::List KalmanFilter(arma::mat y, arma::mat H, arma::mat Q, arma::mat R,
                        arma::mat F, arma::colvec x0, arma::mat P0) {

  const int T = y.n_rows;
  const int n = y.n_cols;
  const int r = F.n_rows;

  double loglik = 0;
  mat K, Pf, Pp;
  colvec xf, xp, xe;
  // Predicted state mean and covariance
  mat xpT(T+1, r, fill::zeros);
  cube PpT(r, r, T+1, fill::zeros);

  // Filtered state mean and covariance
  mat xfT(T, r, fill::zeros);
  cube PfT(r, r, T, fill::zeros);

  mat tH = H;
  mat tR = R;
  mat S;
  uvec miss;
  uvec nmiss = find_finite(F.row(0));
  uvec a(1);

  xp = x0;
  Pp = P0;

  for (int t=0; t < T; t++) {

    // If missing observations are present at some timepoints, exclude the
    // appropriate matrix slices from the filtering procedure.
    miss = find_finite(y.row(t));
    H = tH.submat(miss, nmiss);
    R = tR.submat(miss, miss);
    a[0] = t;

    S = (H * Pp * H.t() + R).i();

    // Prediction error
    xe = y.submat(a, miss).t() - H * xp;
    // Kalman gain
    K = Pp * H.t() * S;
    // Updated state estimate
    xf = xp + K * xe;
    // Updated state covariance estimate
    Pf = Pp - K * H * Pp;

    // Compute likelihood. Skip this part if S is not positive definite.
    if (det(S) > 0) {
      loglik += -0.5 * (double(n) * log(2.0 * datum::pi) - log(det(S)) +
        conv_to<double>::from(xe.t() * S * xe));
    }

    // Store predicted and filtered data needed for smoothing
    xpT.row(t) = xp.t();
    PpT.slice(t) = Pp;
    xfT.row(t) = xf.t();
    PfT.slice(t) = Pf;

    // Run a prediction
    xp = F * xfT.row(t).t();
    Pp = F * PfT.slice(t) * F.t() + Q;

  }

  return Rcpp::List::create(Rcpp::Named("xF") = xfT,
                            Rcpp::Named("Pf") = PfT,
                            Rcpp::Named("xP") = xpT,
                            Rcpp::Named("Pp") = PpT,
                            Rcpp::Named("loglik") = loglik);
}


//' Runs a Kalman smoother
//' @param F transition matrix
//' @param H observation matrix
//' @param R Observation covariance
//' @param xfT State estimates
//' @param xpTm State predicted estimates
//' @param PfT_v Variance estimates
//' @param PpT_v Predicted variance estimates
//' @return List of smoothed estimates
// [[Rcpp::export]]
Rcpp::List KalmanSmoother(arma::mat F, arma::mat H, arma::mat R,
                          arma::mat xfT, arma::mat xpT,
                          Rcpp::NumericVector PfT_v, Rcpp::NumericVector PpT_v) {

  const int T = xfT.n_rows;
  const int r = F.n_rows;
  const int n = H.n_rows;

  cube PfT = array2cube(PfT_v);
  cube PpT = array2cube(PpT_v);

  cube J(r, r, T, fill::zeros);
  cube L(n, n, T, fill::zeros);
  cube K(r, n, T, fill::zeros);

  cube PsTm(r, r, T, fill::zeros);

  // Smoothed state mean and covariance
  mat xsT(T, r, fill::zeros);
  cube PsT(r,r,T, fill::zeros);
  // Initialize smoothed data with last observation of filtered data
  xsT.row(T-1) = xfT.row(T-1);
  PsT.slice(T-1) = PfT.slice(T-1);

  // cube PsTm(r,r,T, fill::zeros);
  for (int t=0; t < T-1; t++) {
    J.slice(t) = PfT.slice(t) * F.t() * PpT.slice(t+1).i();
  }

  // Smoothed state variable and covariance
  for (int j=2; j < T+1; j++) {

    xsT.row(T-j) = xfT.row(T-j) +
      (J.slice(T-j) * (xsT.row(T-j+1) - xpT.row(T-j+1)).t()).t();

    PsT.slice(T-j) = PfT.slice(T-j) +
      J.slice(T-j) * (PsT.slice(T-j+1) - PpT.slice(T-j+1)) * J.slice(T-j).t();

  }

  // Additional variables used in EM-algorithm
  for (int i=0; i < T; i++) {
    L.slice(i) = (H * PpT.slice(i) * H.t() + R).i();
    K.slice(i) = PpT.slice(i) * H.t() * L.slice(i);
  }

  PsTm.slice(T-1) = (eye(r,r) - K.slice(T-1) * H) * F * PfT.slice(T-2);

  for (int j=2; j < T-1; j++) {
    PsTm.slice(T-j) = PfT.slice(T-j) * J.slice(T-j-1).t() + J.slice(T-j)
    * (PsTm.slice(T-j+1) - F * PfT.slice(T-j))
    * J.slice(T-j-1).t();
  }

  return Rcpp::List::create(Rcpp::Named("xS") = xsT,
                            Rcpp::Named("Ps") = PsT,
                            Rcpp::Named("PsTm") = PsTm);
}
