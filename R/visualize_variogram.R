#' @import dplyr

#' @title plots the semivariogram
#'
#' @param vg a list that contains the estimated empirical semivariogram
#'
#' @return nothing
#' @export
plot_variogram <- function(vg) {
  lag   <- vg %>% group_by(bins) %>% summarise(avg=mean(x_dist))
  gamma <- vg %>% group_by(bins) %>% summarise(avg=mean(y_dist))
  plot(lag$avg, gamma$avg, xlab = "lag (h)", ylab = "E[D(h)^2]", main = "semivariogram",
       pch =20, type = "l")
  n <- vg %>% group_by(bins) %>% count() 
  text(lag$avg, gamma$avg, labels = n$n)
}


## from georob package
# Cressie Hawkins robust estimator
est_ch <- function( x ) {
  0.5 * mean( sqrt(abs(x)) )^4 / ( 0.457+0.494/length(x) )
}

#' @title uses Cressie-Hawkins (1980) robust estimator
#'
#' @param vg a list that contains the estimated empirical semivariogram
#'
#' @return nothing
#' @export
plot_variogram_ch <- function(vg) {
  lag   <- vg %>% group_by(bins) %>% summarise(avg=mean(x_dist))
  n     <- (vg %>% group_by(bins) %>% count())$n
  vg    <- vg %>% mutate(abs_diff = sqrt(y_dist))
  est <- vg %>% group_by(bins) %>% summarize(ch = est_ch(y_dist))
  plot(lag$avg, est$ch, xlab = "lag(h)", ylab = "E[D(h)^2]", main = "semivariogram",
       pch =20, type = "l")
  text(lag$avg, est$ch, labels = n)
}

#' @title uses the median instead of the mean
#'
#' @param vg a list that contains the estimated empirical semivariogram
#'
#' @return nothing
#' @export
plot_variogram_median <- function(vg) {
  lag   <- vg %>% group_by(bins) %>% summarise(avg=mean(x_dist))
  n     <- (vg %>% group_by(bins) %>% count())$n
  med   <- (vg %>% group_by(bins) %>% summarize(md = median(y_dist)))$md
  plot(lag$avg, med, xlab = "lag(h)", ylab = "E[D(h)^2]", main = "semivariogram",
       pch =20, type = "l")
  text(lag$avg, med, labels = n)
}

