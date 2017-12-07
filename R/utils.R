#' @import readr

#' @title Read genotype data
#'
#' @param geno_path string of path to eigenstrat geno file
#' @param n int number of columns or samples in the geno file
#'
#' @return Y matrix of genotypes
#' @export
read_geno <- function(geno_path, n){
    
    Y <- as.matrix(read_fwf(geno_path, progress=TRUE, fwf_widths(rep(1, n))))
    colnames(Y) <- NULL
    Y[Y==9] <- NA
    
    return(Y)
}

#' @title Sample heterozygotes
#'
#' @param Y matrix n x p of genotypes
#'
#' @return Y matrix of het sampled genotypes
#' @export
sample_hets <- function(Y){
    
    p <- sum(Y==1, na.rm=TRUE)
    samps <- as.vector(sample(c(0, 2), size=p, replace=TRUE))
    Y[!is.na(Y) & Y == 1] <- samps

    return(Y)
}
