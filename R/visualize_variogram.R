#' @title plots the semivariogram
#'
#' @param vg a list that contains the estimated empirical semivariogram
#'
#' @return nothing
#' @export
plot.variogram <- function(vg){
  plot(vg$lag, vg$gamma, xlab = "lag (h)", ylab = "E[D(h)^2]", main = "semivariogram",
       pch =20, type = "l")
  text(vg$lag, vg$gamma, labels = vg$n)
}
