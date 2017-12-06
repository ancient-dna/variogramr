#' @import readr
#' @import binr


#' @title estimates the variogram based on the phylin package
#'
#' @param x the geographic or temporal distance matrix
#' @param y the genetic distance matrix
#' @param lag the lag time
#' @param range the range (i.e. the value of the distance where the genetic distances stabilize)
#'
#' @return vg a list
#' @export
gen.variogram.phylin <-
  function(x, y, lag = mean(x, na.rm = TRUE)/sqrt(nrow(x)), range = max(x, na.rm=TRUE)) {

    # some variable checking
    if (!(class(x) == "matrix")) {
      #warning("X is not a distance matrix. Attempting conversion...")
      x <- as.matrix(x)
    }
    if (!(class(y) == "matrix")) {
      #warning("Y is not a distance matrix. Attempting conversion...")
      y <- as.matrix(y)
    }
    
    if (!identical(dim(x), dim(y))){
      warnings("X and Y are of not the same dimensions")
    }
    
    yy <- y[which(x < range & upper.tri(y) & !is.na(x) & !is.na(y))]
    xx <- x[which(x < range & upper.tri(y) & !is.na(x) & !is.na(y))]

    lagv <- seq(0, range, lag)
    gamma <- n <- rep(NA, length(lagv))
    tol = lag/2

    for (i in 1:length(lagv) )
    {
      l <- lagv[i]
      il <- which(xx > l-tol & xx <= l+tol, arr.ind=TRUE)
      n[i] <- length(il)
      lagv[i] <- mean(xx[il])
      if (n[i] != 0) {
        gamma[i] <- sum(yy[il])/n[i]
      } else {
        gamma[i] <- lagv[i] <- NA
        
      }
    }

    vg <- list(x=xx, y=yy, lag=lagv, gamma=gamma, n=n)
    vg
  }


#' @title estimates the variogram based on an automatic binning scheme using the binr package
#'
#' @param x the geographic or temporal distance matrix
#' @param y the genetic distance matrix
#' @param target.bins the targeted number of bins
#' @param minpts minimum number of points per bin
#' @param range the range (i.e. the value of the distance where the genetic distances stabilize)
#'
#' @return vg a list
#' @export
gen.variogram <-
    function(x, y, target.bins = 10, minpts = 30, range = max(x, na.rm=TRUE)) {

        # some variable checking
        if (!(class(x) == "matrix")) {
            #warning("X is not a distance matrix. Attempting conversion...")
            x <- as.matrix(x)
        }
        if (!(class(y) == "matrix")) {
            #warning("Y is not a distance matrix. Attempting conversion...")
            y <- as.matrix(y)
        }
      
        if (!identical(dim(x), dim(y))){
          warnings("X and Y are of not the same dimensions")
        }
      
        yy <- y[which(x < range & upper.tri(y) & !is.na(x) & !is.na(y))]
        xx <- x[which(x < range & upper.tri(y) & !is.na(x) & !is.na(y))]

        b=bins(xx, target.bins = target.bins, minpts = minpts)
        nf = cut(xx, bins.getvals(b), labels = names(b$binct), maxpt = range)
        gamma = tapply(yy, nf, mean)
        lag=tapply(xx, nf, mean)

        vg <- list(x=xx, y=yy, lag=lag, gamma=gamma, n=b$binct)
        vg
    }
