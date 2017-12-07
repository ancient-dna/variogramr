#' @import readr
#' @import binr


#' @title estimates the variogram using the traditional method
#'
#' @param x the geographic or temporal distance matrix
#' @param y the genetic distance matrix
#' @param lag the lag time
#' @param range the range (i.e. the value of the distance where the genetic distances stabilize)
#'
#' @return vg a dataframe
#' @export
gen_variogram_classic <- function(x, y, lag = mean(x, na.rm = TRUE)/sqrt(nrow(x)), 
                                  range = max(x, na.rm=TRUE)) {

  ## some of this code is borrowed from the Phylin package
  ## see: https://cran.r-project.org/web/packages/phylin/index.html
  if (!(class(x) == "matrix")) {
    warning("X is not a distance matrix. Attempting conversion...")
    x <- as.matrix(x)
  }
  if (!(class(y) == "matrix")) {
    warning("Y is not a distance matrix. Attempting conversion...")
    y <- as.matrix(y)
  }
  if (!identical(dim(x), dim(y))) {
    warning("X and Y are of not the same dimensions")
  }
    
  yy <- y[which(x < range & upper.tri(y) & !is.na(x) & !is.na(y))]
  xx <- x[which(x < range & upper.tri(y) & !is.na(x) & !is.na(y))]
  bins <- rep(0, length(xx))
  lagv <- seq(0, range, lag)
  tol <- lag/2
  for (i in 1:length(lagv)) {
    l <- lagv[i]
    il <- which(xx > l-tol & xx <= l+tol, arr.ind=TRUE)
    bins[il] <- i
  }
  bins <- as.factor(bins)
  vg <- data.frame(x_dist = xx, y_dist=yy, bins=bins)
  return(vg)
}


#' @title estimates the variogram based on an automatic binning scheme using the binr package
#'
#' @param x the geographic or temporal distance matrix
#' @param y the genetic distance matrix
#' @param target.bins the targeted number of bins
#' @param minpts minimum number of points per bin
#' @param range the range (i.e. the value of the distance where the genetic distances stabilize)
#'
#' @return vg a dataframe
#' @export
gen_variogram <- function(x, y, target.bins = 10, minpts = 30, 
                          range = max(x, na.rm=TRUE)) {
  
  if (!(class(x) == "matrix")) {
    warning("X is not a distance matrix. Attempting conversion...")
    x <- as.matrix(x)
  }
  if (!(class(y) == "matrix")) {
    warning("Y is not a distance matrix. Attempting conversion...")
    y <- as.matrix(y)
  }
  if (!identical(dim(x), dim(y))) {
    warning("X and Y are of not the same dimensions")
  }
      
 yy <- y[which(x < range & upper.tri(y) & !is.na(x) & !is.na(y))]
 xx <- x[which(x < range & upper.tri(y) & !is.na(x) & !is.na(y))]

 b <- bins(xx, target.bins = target.bins, minpts = minpts)
 my_bins <- cut(xx, bins.getvals(b), labels = names(b$binct), maxpt = range)
 levels(my_bins) <- 1:length(levels(my_bins))
    
 vg <- data.frame(x_dist = xx, y_dist=yy, bins=my_bins)
 return(vg)
}
