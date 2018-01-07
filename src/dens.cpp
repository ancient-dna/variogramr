#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

/*
struct Likelihood : public Worker {
    
    // n x p matrix of observed data    
    const RMatrix<double> yt;
    
    // covariance matrix
    const RMatrix<double> s;
    
    // n x n unnormalized distance matrix
    double ll;
    
    const int n; // number of samples
    const int p; // number of features
    
    Likelihood(const NumericMatrix _yt, const NumericMatrix _s)
        : yt(_yt), y(_s), ll(1.0), n(_y.nrow()), p(_y.ncol()) {}
    
    void operator()(std::size_t begin, std::size_t end) {
        
        // snps
        for (std::size_t i = begin; i < end; i++) {
            
            // extract column i
            RMatrix<double>::Column yt_i = yt.column(i);
            
            // loop over samples and count number of non-missing sites
            m = 0;
            for(int j=0; j < n; j++){
                if(!std::isnan(yt_i[j])){
                    m = m + 1;
                }
            }
            
            RVector<double> yt_is(m);
            // to be continued
            
            // loop over samples
            m = 0;
            for(int j=0; j < n; j++){
                if(!std::isnan(yt_i[j])){
                    
                }
            }
                
                
            for(int k=0; k < n; k++){
                    
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
*/