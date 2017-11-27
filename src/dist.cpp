#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

//' @title Times two
//'
//' @description a description about the function
//'              note this is a dummy function for help
//'
//' @param x vector of numeric values
//'
//' @export
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}
