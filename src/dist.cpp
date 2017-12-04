#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

struct Distance : public Worker {
    
    // p x n matrix of observed data    
    const RMatrix<double> y;
    
    // n x n unnormalized distance matrix
    RMatrix<double> q;
    
    // n x n matrix storing number of non-missing sites
    RMatrix<double> m;
    
    const int n; // number of samples
    const int p; // number of features
    
    Distance(const NumericMatrix _y, NumericMatrix _q, NumericMatrix _m)
        : y(_y), q(_q), m(_m), n(_y.ncol()), p(_y.nrow()) {}
    
    void operator()(std::size_t begin, std::size_t end) {
        
        // loop over pairs
        for (std::size_t i = begin; i < end; i++) {
            for (std::size_t j = 0; j <= i; j++) {
                
                // extract column i
                RMatrix<double>::Column y_i = y.column(i);
                
                // extract column j
                RMatrix<double>::Column y_j = y.column(j);

                // loop over features
                for(int k=0; k < p; k++){
                    
                    // compute error
                    double e = y_i[k] - y_j[k];
                    if(!std::isnan(e)){
                        
                        // count number of non-missing features
                        m(i,j) ++;
                        m(j,i) = m(i,j);
                        
                        // compute distance
                        q(i,j) += e * e;
                        q(j,i) = q(i,j);
                        
                    }
                }
                
            }
        }
    }
};

//' @title Mean euclidean distance in parallel
//'
//' @description Computes the mean euclidian distance matrix
//'              for all pairs of samples accross all the features.
//'              This only includes non-missing feature for both pairs.
//'
//'
//' @param y NumericMatrix p x n data matrix
//'
//' @return res list storing Q the n x n unormalized euclidian distance matrix 
//'         and M the n x n matrix storing the number of non-missing features
//'
//' @export
// [[Rcpp::export]]
List mean_dist(NumericMatrix y) {
    
    // number of samples
    int n = y.ncol();
    
    // n x n unormalized distance matrix
    NumericMatrix q(n, n);
    
    // n x n matrix of non-missing counts
    NumericMatrix m(n, n);
    
    Distance Distance(y, q, m);
    parallelFor(0, n, Distance);
    
    // return list
    return List::create(Named("Q") = q, 
                        Named("M") = m);
    
}