#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

//' @title Mean euclidean distance for a single pair of observations
//'
//' @description computes the mean euclidian distance
//'              for a single pair of samples all features
//'
//' @param y_i NumericVector p data vector for sample i
//' @param y_j NumericVector p data vector for sample j
//'
//' @return d_ij euclidian distance normalized by number of non-missing
//'         features
//'
//' @export
// [[Rcpp::export]]
double mean_dist_pair(NumericVector y_i, NumericVector y_j) {

    // number of features
    int p = y_i.length();

    // number of non-missing features
    double m = 0.0;

    // distance for pair of observations
    double d_ij = 0.0;

    // loop over features
    for (int k = 0; k < p; k++) {
        // check that both features are non-missing
        if (!R_IsNA(y_i(k)) && !R_IsNA(y_j(k))) {

            // add to euclidian distiance
            d_ij = d_ij + pow((y_i(k) - y_j(k)), 2);

            // count non-missing features
            m = m + 1.0;
        }
    }

    // normalize by number of non-missing features
    d_ij = d_ij / m;

    return d_ij;
}

//' @title Mean euclidean distance
//'
//' @description Computes the mean euclidian distance matrix
//'              for all pairs of samples accross all the features.
//'              This only includes non-missing feature for both pairs.
//'
//'
//' @param y NumericMatrix n x p data matrix
//'
//' @return d n x n euclidian distance matrix
//'
//' @export
// [[Rcpp::export]]
NumericMatrix mean_dist(NumericMatrix y) {

    // number of samples
    int n = y.nrow();

    // distance matrix
    NumericMatrix d(n, n);

    // loop over all pairs of samples
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {

            // compute distance for each pair
            d(i, j) = mean_dist_pair(y(_, i), y(_, j));

        }
    }

    return d;

}
