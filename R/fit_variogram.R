#' @import MCMCpack stats4

#' @title Compute exponetial covariance kernal
#' 
#' @param D_x n x n input distance matrix 
#' @param gamma2 double variance of the covariance kernel
#' @param l double lengthscale of the covariance kernel
#' 
#' @return S n x n covariance matrix
#' @export
comp_exp_cov <- function(D_x, gamma2, l){
    S <- gamma2 * exp(D_x / l)
    return(S)
}

#' @title Generalized wishart density
#' 
#' @description Computes the density of a distance matrix
#'              using the generalized Wishart distribution
#'              see McCullagh 2009, Hanks and Hooten 2013, and
#'              Petkova et al. 2016
#'              
#' @param D n x n output distance matrix
#' @param nu degrees of freedom paramter
#' @param S n x n covariance matrix
#' @param log boolean to return log probablity
#' 
#' @return probality or log probablity of D | nu, S
#' @export
dgwish <- function(D, nu, S, log=TRUE){
    n <- nrow(S)
    L <- cbind(diag(n - 1), -1)
    D_gen <- L %*% (-D) %*% t(L)
    S_gen <- L %*% (2 * S) %*% t(L)
    lik <- MCMCpack::dwish(D_gen, nu, S_gen)
    if(log){
        lik <- log(lik)
    }
    return(lik)
}

#' @title Fit variogram
#' 
#' @param D n x n output distance matrix
#' @param D_x n x n input distance matrix (i.e. lat lons or time)
#' @param nu degrees of freedom paramter
#' 
#' @return mle_res list of optimzed covariance kernel params
#' @export
fit_variogram <- function(D, D_x, nu){
    
    # wrapper of genaralized wishart log-likelihod
    neg_log_likelihood <- function(z_gamma2, z_l){
        gamma2 <- exp(z) # transfom back to nonnegative space
        l <- exp(z) # transfom back to nonnegative space
        S <- comp_exp_cov(D_x, gamma2, l)
        nll <- -dgwish(D, nu, S, log=TRUE)
        return(nll)
    }
    
    # optimize parameters in unconstrained space
    fit <- stats4::mle(minuslogl=neg_log_likelihood, start=list(z_gamma2=0.0, z_l=0.0))
    
    # gather results in a list
    mle_res <- list(gamma2=as.numeric(exp(fit@coef[1])), l=as.numeric(exp(fit@coef[2])))
    
    return(mle_res)
}

