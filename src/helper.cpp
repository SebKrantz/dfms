#include <RcppArmadillo.h>
#include "helper.h"
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Since RcppArmadillo can't take cubes as arguments in R code,
// this function takes a 3-dimensional array as a numeric
// vector in R and converts it to an Armadillo cube with appropriate
// dimensions.
arma::cube array2cube( Rcpp::NumericVector myArray ) {

  Rcpp::IntegerVector arrayDims = myArray.attr("dim");

  arma::cube cubedArray(myArray.begin(),
                        arrayDims[0], arrayDims[1], arrayDims[2]);

  return(cubedArray);
}

//// Not currently needed...
// arma::field<arma::mat> array2field2mat( Rcpp::NumericVector myArray) {
//
//   Rcpp::IntegerVector arrayDims = myArray.attr("dim");
//   Rcpp::NumericVector::iterator it = myArray.begin();
//
//   arma::field<arma::mat> fieldedArray(arrayDims[2], arrayDims[3]);
//
//   for (int i=0; i < arrayDims[2]; i++) {
//     for (int j=0; j < arrayDims[3]; j++) {
//       fieldedArray(i, j) = arma::mat(it, arrayDims[0], arrayDims[1]);
//       it += (arrayDims[0] * arrayDims[1]);
//       }
//   }
//
//   return fieldedArray;
// }
//
// arma::field<arma::cube> array2field1cube( Rcpp::NumericVector myArray) {
//
//   Rcpp::IntegerVector arrayDims = myArray.attr("dim");
//   Rcpp::NumericVector::iterator it = myArray.begin();
//
//   arma::field<arma::cube> fieldedArray(arrayDims[3]);
//
//   for (int i=0; i < arrayDims[3]; i++) {
//     fieldedArray(i) = arma::cube(it, arrayDims[0], arrayDims[1], arrayDims[2]);
//     it += (arrayDims[0] * arrayDims[1] * arrayDims[2]);
//   }
//
//   return fieldedArray;
// }
//
// arma::field<arma::cube> array2field2cube( Rcpp::NumericVector myArray) {
//
//   Rcpp::IntegerVector arrayDims = myArray.attr("dim");
//   Rcpp::NumericVector::iterator it = myArray.begin();
//
//   arma::field<arma::cube> fieldedArray(arrayDims[3], arrayDims[4]);
//
//   for (int i=0; i < arrayDims[3]; i++) {
//     for (int j=0; j < arrayDims[4]; j++) {
//       fieldedArray(i, j) = arma::cube(it, arrayDims[0], arrayDims[1], arrayDims[2]);
//       it += (arrayDims[0] * arrayDims[1] * arrayDims[2]);
//     }
//   }
//
//   return fieldedArray;
// }


// [[Rcpp::export]]
SEXP ainv(SEXP x) {
  SEXP res = Rcpp::wrap(arma::inv(Rcpp::as<arma::mat>(x)));
  DUPLICATE_ATTRIB(res, x);
  return res;
}

// [[Rcpp::export]]
SEXP apinv(SEXP x) {
  SEXP res = Rcpp::wrap(arma::pinv(Rcpp::as<arma::mat>(x)));
  DUPLICATE_ATTRIB(res, x);
  return res;
}

