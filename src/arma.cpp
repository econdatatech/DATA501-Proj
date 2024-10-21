// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// Making the matrix inversion function from the Armadillo package available
//' Title Armadillo package matrix inversion function
//' Description Takes a matrix and inverts it.
//' @param x matrix to be inverted
//' @return Inverted matrix
//' @noRd
// [[Rcpp::export]]
arma::mat rcpparma_armaInv(const arma::mat & x) { return arma::inv(x); }
