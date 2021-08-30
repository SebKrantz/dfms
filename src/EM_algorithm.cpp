#include <RcppArmadillo.h>
#include "KalmanFiltering.h"
#include "helper.h"
#include <stdio.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::List Estep(arma::mat y, arma::mat C, arma::mat Q, arma::mat R,
                    arma::mat A, arma::colvec x0, arma::mat P0) {

  const unsigned int T = y.n_rows;
  const unsigned int n = y.n_cols;
  const unsigned int rp = A.n_rows;

  // Run Kalman filter and Smoother
  List ks = KalmanFilterSmoother(y, C, Q, R, A, x0, P0);
  double loglik = as<double>(ks["loglik"]);
  mat xS = as<mat>(ks["xS"]);
  cube Vsmooth = array2cube(as<NumericVector>(ks["Ps"]));
  cube Wsmooth = array2cube(as<NumericVector>(ks["PsTm"]));

  // Run computations and return all estimates
  mat delta(n, rp); delta.zeros();
  mat gamma(rp, rp); gamma.zeros();
  mat beta(rp, rp); beta.zeros();

  // For E-step purposes it is sufficient to set missing observations
  // to being 0.
  y(find_nonfinite(y)).zeros();

  for (unsigned int t=0; t<T; ++t) {
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

/*
// [[Rcpp::export]]
Rcpp::List EMstep(arma::mat y, arma::mat C, arma::mat Q, arma::mat R,
                  arma::mat A, arma::colvec x0, arma::mat P0, arma::mat xx, int r) {

  const unsigned int T = y.n_rows;
  const unsigned int n = y.n_cols;
  const unsigned int rp = A.n_rows;

  // Run Kalman filter and Smoother
  List ks = KalmanFilterSmoother(y, C, Q, R, A, x0, P0);
  double loglik = as<double>(ks["loglik"]);
  mat xS = as<mat>(ks["xS"]);
  cube Vsmooth = array2cube(as<NumericVector>(ks["Ps"]));
  cube Wsmooth = array2cube(as<NumericVector>(ks["PsTm"]));

  // Run computations and return all estimates
  mat delta(n, rp); delta.zeros();
  mat gamma(rp, rp); gamma.zeros();
  mat beta(rp, rp); beta.zeros();

  // For E-step purposes it is sufficient to set missing observations
  // to being 0.
  y(find_nonfinite(y)).zeros();

  for (unsigned int t=0; t<T; ++t) {
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

  em_res <- Estep(X, C, Q, R, A, x0, P0)
    beta <- em_res$beta_t
    gamma <- em_res$gamma_t
    delta <- em_res$delta_t
    gamma1 <- em_res$gamma1_t
    gamma2 <- em_res$gamma2_t
    P1sum <- em_res$V1 + tcrossprod(em_res$x1)
    x1sum <- em_res$x1
    loglik <- em_res$loglik_t

    num_iter <- num_iter + 1L

  // M-step computes model parameters as a function of the sufficient
  // statistics that were computed with the E-step. Iterate the procedure
  // until convergence. Due to the model specification, likelihood maximiation
  // in the M-step is just an OLS estimation. In particular, X_t = C*F_t and
  // F_t = A*F_(t-1).

   C.head_cols(r) = delta.head_cols(r) * pinv(gamma.submat(1, r, 1, r));
   mat A_update = beta.head_rows(r) * inv(gamma1);
   A.head_rows(r) = A_update;
   Q.submat(1, r, 1, r) = (gamma2.submat(1, r, 1, r) - A_update * beta.head_rows(r).t()) / (T-1);
   R = (xx.t() * xx - C * delta.t()) / T;
    RR <- diag(R); RR[RR < 1e-7] <- 1e-7; R <- diag(RR)
      R <- diag(diag(R))

## Assign new initial values for next EM-algorithm step
      x0 <- x1sum
        P0 <- P1sum - tcrossprod(x0)

        converged <- em_converged(loglik, previous_loglik, threshold = tol)
        previous_loglik <- loglik

## Iterate at least 25 times
        if(num_iter < 25L) converged <- FALSE
}

 */
