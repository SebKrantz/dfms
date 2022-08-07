#include <RcppArmadillo.h>
#include "helper.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// TODO: what about likelihood with BM??

// Implementation of a Kalman filter
// X Data matrix (T x n)
// A Transition matrix (rp x rp)
// C Observation matrix (n x rp)
// Q State covariance (rp x rp)
// R Observation covariance (n x n)
// F0 Initial state vector (rp x 1)
// P0 Initial state covariance (rp x rp)
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
  if(nmiss.n_elem == 0) Rcpp::stop("Missing first row of transition matrix");

  for (int i = 0; i < T; ++i) {

    // If missing observations are present at some timepoints, exclude the
    // appropriate matrix slices from the filtering procedure.
    miss = find_finite(X.row(i));
    if(miss.n_elem > 0) {

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

    } else { // If all missing: just prediction.
      Zf = Zp;
      Vf = Vp;
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

  // Store final prediction for time T+1 (prediction for time 1 was initialization values)
  ZTp.row(T) = Zp.t();
  VTp.slice(T) = Vp;

  return Rcpp::List::create(Rcpp::Named("F") = ZTf,
                            Rcpp::Named("P") = VTf,
                            Rcpp::Named("F_pred") = ZTp,
                            Rcpp::Named("P_pred") = VTp,
                            Rcpp::Named("loglik") = loglik);
}

// Runs a Kalman smoother
// A Transition matrix (rp x rp)
// ZTf State estimates
// ZTp State predicted estimates
// VTf_v Variance estimates
// VTp_v Predicted variance estimates
// [[Rcpp::export]]
Rcpp::List KalmanSmoother(arma::mat A,
                          arma::mat ZTf, arma::mat ZTp,
                          Rcpp::NumericVector VTf_v,
                          Rcpp::NumericVector VTp_v) {

  const int T = ZTf.n_rows;
  const int rp = A.n_rows;

  cube VTf = array2cube(VTf_v);
  cube VTp = array2cube(VTp_v);

  // Smoothed state mean and covariance
  mat ZsT(T, rp, fill::zeros);
  cube VsT(rp, rp, T, fill::zeros);

  // Initialize smoothed data with last observation of filtered data
  ZsT.row(T-1) = ZTf.row(T-1);
  VsT.slice(T-1) = VTf.slice(T-1);

  mat At = A.t(), Ji, Vfi;

  // Smoothed state variable and covariance
  // See e.g. Shumway and Stoffer (2002) p297, or astsa::Ksmooth0
  for (int i = T-1; i--; ) {
    Vfi = VTf.slice(i);
    Ji = Vfi * At * VTp.slice(i+1).i();
    ZsT.row(i) = ZTf.row(i) + (Ji * (ZsT.row(i+1) - ZTp.row(i+1)).t()).t();
    VsT.slice(i) = Vfi + Ji * (VsT.slice(i+1) - VTp.slice(i+1)) * Ji.t();
  }

  // Smoothed value for t = 0
  Vfi = VTp.slice(0);
  Ji = Vfi * At * VTp.slice(0).i();

  return Rcpp::List::create(Rcpp::Named("F_smooth") = ZsT,
                            Rcpp::Named("P_smooth") = VsT,
                            Rcpp::Named("F_smooth_0") = ZTp.row(0) + (Ji * (ZsT.row(0) - ZTp.row(0)).t()).t(),
                            Rcpp::Named("P_smooth_0") = Vfi + Ji * (VsT.slice(0) - Vfi) * Ji.t());
}



// Kalman Filter and Smoother
// X Data matrix (T x n)
// A Transition matrix (rp x rp)
// C Observation matrix (n x rp)
// Q State covariance (rp x rp)
// R Observation covariance (n x n)
// F0 Initial state vector (rp x 1)
// P0 Initial state covariance (rp x rp)
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
  if(nmiss.n_elem == 0) Rcpp::stop("Missing first row of transition matrix");

  for (int i = 0; i < T; ++i) {

    // If missing observations are present at some timepoints, exclude the
    // appropriate matrix slices from the filtering procedure.
    miss = find_finite(X.row(i));
    if(miss.n_elem > 0) {

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
      Vf =  0.5 * (Vf + Vf.t()); // Ensure symmetry

      // Compute likelihood. Skip this part if S is not positive definite.
      detS = det(S);
      if(detS > 0) {
        loglik += -0.5 * (dn * log(2.0 * datum::pi) - log(detS) +
          conv_to<double>::from(et.t() * S * et));
      }

    } else { // If all missing: just prediction.
      Zf = Zp;
      Vf = Vp;
    }

    // Store predicted and filtered data needed for smoothing
    ZTp.row(i) = Zp.t();
    VTp.slice(i) = Vp;
    ZTf.row(i) = Zf.t();
    VTf.slice(i) = Vf;

    // Run a prediction
    Zp = A * ZTf.row(i).t();
    Vp = A * VTf.slice(i) * A.t() + Q;
    Vp =  0.5 * (Vp + Vp.t()); // Ensure symmetry

  }

  // Store final prediction for time T+1 (prediction for time 1 was initialization values)
  ZTp.row(T) = Zp.t();
  VTp.slice(T) = Vp;

  // Smoothed state mean and covariance
  mat ZsT(T, rp, fill::zeros);
  cube VsT(rp, rp, T, fill::zeros);
  cube VVsT(rp, rp, T, fill::zeros); // Cov(Z_t, Z_t-1), used in EM

  // Initialize smoothed data with last observation of filtered data
  ZsT.row(T-1) = Zf.t();
  VsT.slice(T-1) = Vf;
  K = (miss.n_elem == 0) ? mat(rp, rp, fill::zeros) : K * Ci;
  VVsT.slice(T-1) = (eye(rp,rp) - K) * A * VTf.slice(T-2);

  mat At = A.t(), Jimt = (VTf.slice(T-2) * At * VTp.slice(T-1).i()).t(), Ji;

  // Smoothed state variable and covariance
  for (int i = T-1; i--; ) {
    Vf = VTf.slice(i);
    Ji = Jimt.t();
    ZsT.row(i) = ZTf.row(i) + (Ji * (ZsT.row(i+1) - ZTp.row(i+1)).t()).t();
    VsT.slice(i) = Vf + Ji * (VsT.slice(i+1) - VTp.slice(i+1)) * Jimt;
    // Cov(Z_t, Z_t-1): Needed for EM
    if(i > 0) Jimt = (VTf.slice(i-1) * At * VTp.slice(i).i()).t();
    VVsT.slice(i) = Vf * Jimt + Ji * (VVsT.slice(i+1) - A * Vf) * Jimt;
  }

  // Smoothed value for t = 0
  Vf = VTp.slice(0);
  Ji = Vf * At * VTp.slice(0).i();

  return Rcpp::List::create(Rcpp::Named("F") = ZTf,
                            Rcpp::Named("P") = VTf,
                            Rcpp::Named("F_pred") = ZTp,
                            Rcpp::Named("P_pred") = VTp,
                            Rcpp::Named("F_smooth") = ZsT,
                            Rcpp::Named("P_smooth") = VsT,
                            Rcpp::Named("PPm_smooth") = VVsT,
                            Rcpp::Named("F_smooth_0") = ZTp.row(0) + (Ji * (ZsT.row(0) - ZTp.row(0)).t()).t(),
                            Rcpp::Named("P_smooth_0") = Vf + Ji * (VsT.slice(0) - Vf) * Ji.t(),
                            Rcpp::Named("loglik") = loglik);
}






