#include <RcppArmadillo.h>
#include "KalmanFiltering.h"
#include "helper.h"
#include <stdio.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::List Estep(arma::mat y, arma::mat H, arma::mat Q, arma::mat R,
                    arma::mat F, arma::colvec x0, arma::mat P0) {

  const unsigned int T = y.n_rows;
  const unsigned int n = y.n_cols;
  const unsigned int r = F.n_rows;

  // Run Kalman filter and save its output for smoothing
  List kf = KalmanFilter(y, H, Q, R, F, x0, P0);

  mat xF = as<mat>(kf["xF"]);
  mat xP = as<mat>(kf["xP"]);
  NumericVector Pf = as<NumericVector>(kf["Pf"]);
  NumericVector Pp = as<NumericVector>(kf["Pp"]);

  double loglik = as<double>(kf["loglik"]);

  // Run Kalman smoother and save its output for further computations
  List ks = KalmanSmoother(F, H, R, xF, xP, Pf, Pp);

  mat xS = as<mat>(ks["xS"]);
  cube Vsmooth = array2cube(as<NumericVector>(ks["Ps"]));
  cube Wsmooth = array2cube(as<NumericVector>(ks["PsTm"]));

  // Run computations and return all estimates
  mat delta(n, r); delta.zeros();
  mat gamma(r, r); gamma.zeros();
  mat beta(r, r); beta.zeros();

  // For E-step purposes it is sufficient to set missing observations
  // to being 0.
  y(find_nonfinite(y)).zeros();

  for (unsigned int t=0; t<T; t++) {
    delta += y.row(t).t() * xS.row(t);
    gamma += xS.row(t).t() * xS.row(t) + Vsmooth.slice(t);
    if (t > 0) {
      beta += xS.row(t).t() * xS.row(t-1) + Wsmooth.slice(t);
    }
  }

  mat gamma1 = gamma - xS.row(T-1).t() * xS.row(T-1) - Vsmooth.slice(T-1);
  mat gamma2 = gamma - xS.row(0).t() * xS.row(0) - Vsmooth.slice(0);
  colvec x1 = xS.row(0).t();
  mat V1 = Vsmooth.slice(0);

  return Rcpp::List::create(Rcpp::Named("beta_t") = beta,
                            Rcpp::Named("gamma_t") = gamma,
                            Rcpp::Named("delta_t") = delta,
                            Rcpp::Named("gamma1_t") = gamma1,
                            Rcpp::Named("gamma2_t") = gamma2,
                            Rcpp::Named("x1") = x1,
                            Rcpp::Named("V1") = V1,
                            Rcpp::Named("loglik_t") = loglik,
                            Rcpp::Named("xS") = ks["xS"]);
}
