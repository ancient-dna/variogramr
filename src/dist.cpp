#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

//' @title Mean euclidean distance
//'
//' @description Computes the mean euclidian distance matrix
//'              for all pairs of samples accross all the features.
//'              This only includes non-missing feature for both pairs.
//'
//'
//' @param y NumericMatrix n x p data matrix
//'
//' @return res list storing d the n x n euclidian distance matrix and
//'         m the n x n matrix storing the number of non-missing features
//'
//' @export
// [[Rcpp::export]]
List mean_dist(NumericMatrix y) {

    // number of samples
    int n = y.nrow();
    int p = y.ncol();
    NumericVector e;

    // euclidian distance matrix
    NumericMatrix d(n, n);
    
    // matrix storing number of non-missing features
    NumericMatrix m(n, n);

    // return list
    List res;

    // loop over pairs to need to compute full distance matrix
    for (int i = 0; i < (n - 1); i++) {
        for (int j = (i+1); j < n; j++) {

            // compute error for pair i and j
            e = y(i,_) - y(j,_);

            // compute number of non-missing sites
            m(i, j) = p - sum(is_na(e));
            m(j, i) = m(i, j);

            // compute distance for each pair
            d(i, j) = sum(na_omit(e * e)) / m(i, j);;
            d(j, i) = d(i, j);

        }
    }

    res["d"] = d;
    res["m"] = m;

    return res;

}