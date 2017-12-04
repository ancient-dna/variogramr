#' @import readr
#' @import binr


#' @title estimates the variogram based on the phylin package
#'
#' @param x the geographic or temporal distance matrix
#' @param y the genetic distance matrix
#' @param lag the lag time
#' @param tol the bin-size
#' @param range the range (i.e. the value of the distance where the genetic distances stabilize)
#'
#' @return vg a list
#' @export
gen.variogram.phylin <-
  function(x, y, lag = mean(x)/sqrt(nrow(x)), tol=lag/2, range = max(x)) {

    # some variable checking
    if (!(class(x) == "matrix")) {
      #warning("X is not a distance matrix. Attempting conversion...")
      x <- as.matrix(x)
    }
    if (!(class(y) == "matrix")) {
      #warning("Y is not a distance matrix. Attempting conversion...")
      y <- as.matrix(y)
    }

    lagv <- seq(0, range, lag)
    gamma <- n <- rep(NA, length(lagv))

    for (i in 1:length(lagv) )
    {
      l <- lagv[i]

      #remove duplicates from distance matrix
      il <- which(x > l-tol & x <= l+tol & y != 0, arr.ind=TRUE)
      il <- unique(t(apply(il, 1, sort)))

      n[i] <- nrow(il)
      if (n[i] != 0) {
        gamma[i] <- sum(y[il])/n[i]
        lagv[i] <- mean(x[il])
      } else {
        gamma[i] <- lagv[i] <- NA
      }
    }

    # vg
    vg <- list(x=x, y=y, lag=lagv, gamma=gamma, n=n)
    vg
  }


#' @title estimates the variogram based on binning data based on number per bin
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
    function(x, y, target.bins = 10, minpts = 50, range = max(x)) {

        # some variable checking
        if (!(class(x) == "matrix")) {
            #warning("X is not a distance matrix. Attempting conversion...")
            x <- as.matrix(x)
        }
        if (!(class(y) == "matrix")) {
            #warning("Y is not a distance matrix. Attempting conversion...")
            y <- as.matrix(y)
        }

        y <- y[which(x < range & x > 0)]
        x <- x[which(x < range & x > 0)]

        b=bins(x, target.bins = target.bins, minpts = minpts)
        nf = cut(x, bins.getvals(b), labels = names(b$binct), maxpt = range)
        gamma = tapply(y, nf, mean)
        lag=tapply(x, nf, mean)

        vg <- list(x=x, y=y, lag=lag, gamma=gamma, n=b$binct)
        vg
    }
