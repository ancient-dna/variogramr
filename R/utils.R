#' @import readr

#' @title Read genotype data
#'
#' @param geno_path string of path to eigenstrat geno file
#' @param n int number of columns or samples in the geno file
#'
#' @return Y matrix of genotypes
#' @export
read_geno <- function(geno_path, n){
    df <- read_fwf(geno_path, progress=TRUE, na='9', fwf_widths(rep(1, n)))
    YT <- as.matrix(df)
    colnames(YT) <- NULL

    return(t(YT))
}

#' @title Sample heterozygotes
#'
#' @param Y matrix n x p of genotypes
#' @param idx row indicies of people you want to het sample
#'
#' @return Y_samp matrix of het sampled genotypes
#' @export
sample_hets <- function(Y, idx){
    NULL
}
