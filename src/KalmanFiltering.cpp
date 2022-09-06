#include <RcppArmadillo.h>
#include "helper.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


// Implementation of a Kalman filter
// X Data matrix (T x n)
// A Transition matrix (rp x rp)
// C Observation matrix (n x rp)
// Q State covariance (rp x rp)
// R Observation covariance (n x n)
// F_0 Initial state vector (rp x 1)
// P_0 Initial state covariance (rp x rp)
// retLL Return log-likelihood.
// [[Rcpp::export]]
Rcpp::List SKF(arma::mat X, arma::mat A, arma::mat C, arma::mat Q,
               arma::mat R, arma::colvec F_0, arma::mat P_0, bool retLL = false) {

  const int T = X.n_rows;
  const int n = X.n_cols;
  const int rp = A.n_rows;
  int n_c;

  // In internal code factors are Z (instead of F) and factor covariance V (instead of P),
  // to avoid confusion between the matrices and their predicted (p) and filtered (f) states.
  // Additionally the results matrices for all time periods have a T in the name.

  double loglik = retLL ? 0 : NA_REAL, dn = 0, detS;
  if(retLL) dn = n * log(2.0 * datum::pi);
  colvec Zp, Zf = F_0, et, xt;
  mat K, Vp, Vf = P_0, S, VCt;

  // Predicted state mean and covariance
  mat ZTp(T, rp, fill::zeros);
  cube VTp(rp, rp, T, fill::zeros);

  // Filtered state mean and covariance
  mat ZTf(T, rp, fill::zeros);
  cube VTf(rp, rp, T, fill::zeros);

  // Handling missing values in the filter
  mat Ci, Ri;
  uvec nmiss, arow = find_finite(A.row(0));
  if(arow.n_elem == 0) Rcpp::stop("Missing first row of transition matrix");

  for (int i = 0; i < T; ++i) {

    // Run a prediction
    Zp = A * Zf;
    Vp = A * Vf * A.t() + Q;
    Vp += Vp.t(); // Ensure symmetry
    Vp *= 0.5;

    // If missing observations are present at some timepoints, exclude the
    // appropriate matrix slices from the filtering procedure.
    xt = X.row(i).t();
    nmiss = find_finite(xt);
    n_c = nmiss.n_elem;
    if(n_c > 0) {
      if(n_c == n) {
        Ci = C;
        Ri = R;
      } else {
        Ci = C.submat(nmiss, arow);
        Ri = R.submat(nmiss, nmiss);
        xt = xt.elem(nmiss);
      }

      // Intermediate results
      VCt = Vp * Ci.t();
      S = inv_sympd(Ci * VCt + Ri); // .i();

      // Prediction error
      et = xt - Ci * Zp;
      // Kalman gain
      K = VCt * S;
      // Updated state estimate
      Zf = Zp + K * et;
      // Updated state covariance estimate
      Vf = Vp - K * Ci * Vp;
      Vf += Vf.t(); // Ensure symmetry
      Vf *= 0.5;

      // Compute likelihood. Skip this part if S is not positive definite.
      if(retLL) {
        detS = det(S);
        if(detS > 0) loglik += log(detS) - conv_to<double>::from(et.t() * S * et) - dn;
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

  }

  if(retLL) loglik *= 0.5;

  return Rcpp::List::create(Rcpp::Named("F") = ZTf,
                            Rcpp::Named("P") = VTf,
                            Rcpp::Named("F_pred") = ZTp,
                            Rcpp::Named("P_pred") = VTp,
                            // Rcpp::Named("F_0") = F_0,
                            // Rcpp::Named("P_0") = P_0,
                            Rcpp::Named("loglik") = loglik);
}

// Runs a Kalman smoother
// A Transition matrix (rp x rp)
// ZTf State estimates
// ZTp State predicted estimates
// VTf_v Variance estimates
// VTp_v Predicted variance estimates
// F_0 Initial state vector (rp x 1)
// P_0 Initial state covariance (rp x rp)
// [[Rcpp::export]]
Rcpp::List FIS(arma::mat A,
               arma::mat ZTf, arma::mat ZTp,
               Rcpp::NumericVector VTf_v,
               Rcpp::NumericVector VTp_v,
               SEXP F_0SEXP, SEXP P_0SEXP) {

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

  mat At = A.t(), Ji, Vfi, Vpi;

  // Smoothed state variable and covariance
  // See e.g. Shumway and Stoffer (2002) p297, or astsa::Ksmooth0
  for (int i = T-1; i--; ) {
    Vfi = VTf.slice(i);
    Vpi = VTp.slice(i+1);
    Ji = Vfi * At * inv_sympd(Vpi); // .i();
    ZsT.row(i) = ZTf.row(i) + (Ji * (ZsT.row(i+1) - ZTp.row(i+1)).t()).t();
    VsT.slice(i) = Vfi + Ji * (VsT.slice(i+1) - Vpi) * Ji.t();
  }

  if(Rf_isNull(F_0SEXP) || Rf_isNull(P_0SEXP))
    return Rcpp::List::create(Rcpp::Named("F_smooth") = ZsT,
                              Rcpp::Named("P_smooth") = VsT);


  // From https://dirk.eddelbuettel.com/papers/rcpp_ku_nov2013-part2.pdf
  Rcpp::NumericVector F_0_v(F_0SEXP);
  Rcpp::NumericMatrix P_0_v(P_0SEXP);
  int n = P_0_v.nrow(), k = P_0_v.ncol();
  arma::mat P_0(P_0_v.begin(), n, k, false);
  arma::rowvec F_0(F_0_v.begin(), n, false);

  // Smoothed value for t = 0
  Vpi = VTp.slice(0);
  Ji = P_0 * At * Vpi.i();

  return Rcpp::List::create(Rcpp::Named("F_smooth") = ZsT,
                            Rcpp::Named("P_smooth") = VsT,
                            Rcpp::Named("F_smooth_0") = F_0 + (Ji * (ZsT.row(0) -  ZTp.row(0)).t()).t(),
                            Rcpp::Named("P_smooth_0") = P_0 + Ji * (VsT.slice(0) - Vpi) * Ji.t());
}



// Kalman Filter and Smoother
// X Data matrix (T x n)
// A Transition matrix (rp x rp)
// C Observation matrix (n x rp)
// Q State covariance (rp x rp)
// R Observation covariance (n x n)
// F_0 Initial state vector (rp x 1)
// P_0 Initial state covariance (rp x rp)
// retLL 0-no likelihood, 1-standard Kalman Filter, 2-BM14
// [[Rcpp::export]]
Rcpp::List SKFS(arma::mat X, arma::mat A, arma::mat C, arma::mat Q,
                arma::mat R, arma::colvec F_0, arma::mat P_0, bool retLL = false) {

  const int T = X.n_rows;
  const int n = X.n_cols;
  const int rp = A.n_rows;
  int n_c;

  // In internal code factors are Z (instead of F) and factor covariance V (instead of P),
  // to avoid confusion between the matrices and their predicted (p) and filtered (f) states.
  // Additionally the results matrices for all time periods have a T in the name.

  double loglik = retLL ? 0 : NA_REAL, dn = 0, detS;
  if(retLL) dn = n * log(2.0 * datum::pi);
  colvec Zp, Zf = F_0, et, xt;
  mat K, Vp, Vf = P_0, S, VCt;

  // Predicted state mean and covariance
  mat ZTp(T, rp, fill::zeros);
  cube VTp(rp, rp, T, fill::zeros);

  // Filtered state mean and covariance
  mat ZTf(T, rp, fill::zeros);
  cube VTf(rp, rp, T, fill::zeros);

  // Handling missing values in the filter
  mat Ci, Ri;
  uvec nmiss, arow = find_finite(A.row(0));
  if(arow.n_elem == 0) Rcpp::stop("Missing first row of transition matrix");

  for (int i = 0; i < T; ++i) {

    // Run a prediction
    Zp = A * Zf;
    Vp = A * Vf * A.t() + Q;
    Vp += Vp.t(); // Ensure symmetry
    Vp *= 0.5;

    // If missing observations are present at some timepoints, exclude the
    // appropriate matrix slices from the filtering procedure.
    xt = X.row(i).t();
    nmiss = find_finite(xt);
    n_c = nmiss.n_elem;
    if(n_c > 0) {
      if(n_c == n) {
        Ci = C;
        Ri = R;
      } else {
        Ci = C.submat(nmiss, arow);
        Ri = R.submat(nmiss, nmiss);
        xt = xt.elem(nmiss);
      }

      // Intermediate results
      VCt = Vp * Ci.t();
      S = inv_sympd(Ci * VCt + Ri); //.i();

      // Prediction error
      et = xt - Ci * Zp;
      // Kalman gain
      K = VCt * S;
      // Updated state estimate
      Zf = Zp + K * et;
      // Updated state covariance estimate
      Vf = Vp - K * Ci * Vp;
      Vf += Vf.t(); // Ensure symmetry
      Vf *= 0.5;

      // Compute likelihood. Skip this part if S is not positive definite.
      if(retLL) {
        detS = det(S);
        if(detS > 0) loglik += log(detS) - conv_to<double>::from(et.t() * S * et) - dn; // Standard Kalman Filter Likelihood
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

  }

  if(retLL) loglik *= 0.5;

  // Smoothed state mean and covariance
  mat ZsT(T, rp, fill::zeros);
  cube VsT(rp, rp, T, fill::zeros);
  cube VVsT(rp, rp, T, fill::zeros); // Cov(Z_t, Z_t-1), used in EM

  // Initialize smoothed data with last observation of filtered data
  ZsT.row(T-1) = Zf.t();
  VsT.slice(T-1) = Vf;
  K = (n_c == 0) ? mat(rp, rp, fill::zeros) : K * Ci;
  VVsT.slice(T-1) = (eye(rp,rp) - K) * A * VTf.slice(T-2);

  mat At = A.t(), Jimt = (VTf.slice(T-2) * At * VTp.slice(T-1).i()).t(), Ji;

  // Smoothed state variable and covariance
  for (int i = T-1; i--; ) {
    Vf = VTf.slice(i);
    Ji = Jimt.t();
    ZsT.row(i) = ZTf.row(i) + (Ji * (ZsT.row(i+1) - ZTp.row(i+1)).t()).t();
    VsT.slice(i) = Vf + Ji * (VsT.slice(i+1) - VTp.slice(i+1)) * Jimt;
    // Cov(Z_t, Z_t-1): Needed for EM
    if(i > 0) {
      Jimt = (VTf.slice(i-1) * At * inv_sympd(VTp.slice(i))).t(); // .i()
      VVsT.slice(i) = Vf * Jimt + Ji * (VVsT.slice(i+1) - A * Vf) * Jimt;
    }
  }

  // Smoothed value for t = 0
  Vp = VTp.slice(0);
  Jimt = (P_0 * At * Vp.i()).t();
  VVsT.slice(0) = Vf * Jimt + Ji * (VVsT.slice(1) - A * Vf) * Jimt;
  Ji = Jimt.t();


  return Rcpp::List::create(Rcpp::Named("F") = ZTf,
                            Rcpp::Named("P") = VTf,
                            Rcpp::Named("F_pred") = ZTp,
                            Rcpp::Named("P_pred") = VTp,
                            Rcpp::Named("F_smooth") = ZsT,
                            Rcpp::Named("P_smooth") = VsT,
                            Rcpp::Named("PPm_smooth") = VVsT,
                            Rcpp::Named("F_smooth_0") = F_0.t() + (Ji * (ZsT.row(0) - ZTp.row(0)).t()).t(),
                            Rcpp::Named("P_smooth_0") = P_0 + Ji * (VsT.slice(0) - Vp) * Jimt,
                            Rcpp::Named("loglik") = loglik);
}






