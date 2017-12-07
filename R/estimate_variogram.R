#' @import binr

#' @title Estimates the variogram using the traditional method
#'
#' @param D_x the geographic or temporal distance matrix
#' @param D the genetic distance matrix
#' @param lag the lag time
#' @param range the range (i.e. the value of the distance where the genetic distances stabilize)
#'
#' @return vg a dataframe
#' @export
gen_variogram_classic <- function(D_x, D, lag=mean(D_x, na.rm = TRUE) / sqrt(nrow(D_x)), 
                                  range=max(D_x, na.rm=TRUE)) {

  # some of this code is borrowed from the Phylin package
  # see: https://cran.r-project.org/web/packages/phylin/index.html
  if (!(class(D_x) == "matrix")) {
    warning("D_x is not a distance matrix. Attempting conversion...")
    D_x <- as.matrix(D_x)
  }
  if (!(class(D) == "matrix")) {
    warning("D is not a distance matrix. Attempting conversion...")
    D <- as.matrix(D)
  }
  if (!identical(dim(X), dim(X))) {
    warning("X and Y are of not the same dimensions")
  }
    
  d <- D[which(D_x < range & upper.tri(D) & !is.na(D_x) & !is.na(D))]
  d_x <- D_x[which(D_x < range & upper.tri(D) & !is.na(D_x) & !is.na(D))]

  bins <- rep(0, length(d_x))
  lagv <- seq(0, range, lag)
  tol <- lag / 2
  for (i in 1:length(lagv)) {
    l <- lagv[i]
    il <- which(d_x > l - tol & d_x <= l + tol, arr.ind=TRUE)
    bins[il] <- i
  }

  bins <- as.factor(bins)
  vg <- data.frame(d_x=d_x, d=d, bins=bins)

  return(vg)

}


#' @title Estimates the variogram based on an automatic binning scheme using the binr package
#'
#' @param D_x the geographic or temporal distance matrix
#' @param D the genetic distance matrix
#' @param target_bins the targeted number of bins
#' @param minpts minimum number of points per bin
#' @param range the range (i.e. the value of the distance where the genetic distances stabilize)
#'
#' @return vg a dataframe
#' @export
gen_variogram <- function(D_x, D, target_bins=10, minpts=30, range=max(D_x, na.rm=TRUE)) {
  
  if (!(class(D_x) == "matrix")) {
    warning("D_x is not a distance matrix. Attempting conversion...")
    D_x <- as.matrix(D_x)
  }
  if (!(class(D) == "matrix")) {
    warning("D is not a distance matrix. Attempting conversion...")
    D <- as.matrix(D)
  }
  if (!identical(dim(D_x), dim(D))) {
    warning("D_x and D are of not the same dimensions")
  }
      
  d <- D[which(D_x < range & upper.tri(D) & !is.na(D_x) & !is.na(D))]
  d_x <- D_x[which(D_x < range & upper.tri(D) & !is.na(D_x) & !is.na(D))]

  b <- bins(d_x, target.bins=target_bins, minpts=minpts)
  my_bins <- cut(d_x, bins.getvals(b), labels=names(b$binct), maxpt=range)
  levels(my_bins) <- 1:length(levels(my_bins))
    
  vg <- data.frame(d_x=d_x, d=d, bins=my_bins)

  return(vg)

}
