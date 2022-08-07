Rcpp::List KalmanFilter(arma::mat X, arma::mat A, arma::mat C, arma::mat Q,
                        arma::mat R, arma::colvec F0, arma::mat P0, bool retLL = true);

Rcpp::List KalmanSmoother(arma::mat A,
                          arma::mat ZTf, arma::mat ZTp,
                          Rcpp::NumericVector VTf_v,
                          Rcpp::NumericVector VTp_v);

Rcpp::List KalmanFilterSmoother(arma::mat X, arma::mat A, arma::mat C, arma::mat Q,
                                arma::mat R, arma::colvec F0, arma::mat P0, int retLL = 1);
