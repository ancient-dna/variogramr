context("Distance")

test_that("Correct euclidian distances are computed with no missing data", {
    # p x n data matrix
    Y <- matrix(c(0, 0, 1,
                  2, 1, 1),
                  byrow=TRUE,
                  ncol=3)
    
    # compute distances from data matrix
    res <- mean_dist(Y)
    
    # unnormalized distaces
    Q <- res$Q
    
    # number of non-missing features
    M <- res$M
    
    # distance matrix
    D <- Q / M
    
    # expected distance matrix
    D_exp <- matrix(c(0, .5, 1,
                      .5, 0, .5,
                      1, .5, 0),
                    byrow=TRUE, 
                    ncol=3)
    
    # expected matrix of non-missing features
    M_exp <- matrix(2, ncol=3, nrow=3)
    
    # test!!!
    expect_equal(M, M_exp)
    expect_equal(D, D_exp)
})

test_that("Correct euclidian distances are computed with missing data", {
    # p x n data matrix
    Y <- matrix(c(0, 0, NA,
                  2, NA, 1),
                byrow=TRUE,
                ncol=3)
    
    # compute distances from data matrix
    res <- mean_dist(Y)
    
    # unnormalized distaces
    Q <- res$Q 
    
    # number of non-missing features
    M <- res$M
    
    # distance matrix
    D <- Q / M
    
    # expected distance matrix
    D_exp <- matrix(c(0, 0, 1,
                      0, 0, NA,
                      1, NA, 0),
                    byrow=TRUE, 
                    ncol=3)
    
    # expected matrix of non-missing features
    M_exp <- matrix(c(2, 1, 1,
                      1, 1, 0,
                      1, 0, 1), 
                    byrow=TRUE,
                    nrow=3)
    
    # test!!!
    expect_equal(M, M_exp)
    expect_equal(D, D_exp)
})