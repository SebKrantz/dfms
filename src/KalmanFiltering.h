Rcpp::List fKF(arma::mat X, arma::mat A, arma::mat C, arma::mat Q,
                        arma::mat R, arma::colvec F0, arma::mat P0, bool retLL = false);

Rcpp::List fKS(arma::mat A,
                          arma::mat ZTf, arma::mat ZTp,
                          Rcpp::NumericVector VTf_v,
                          Rcpp::NumericVector VTp_v);

Rcpp::List fKFS(arma::mat X, arma::mat A, arma::mat C, arma::mat Q,
                                arma::mat R, arma::colvec F0, arma::mat P0, int retLL = 0);
