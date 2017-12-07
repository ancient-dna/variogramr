#' @import dplyr

#' @title plots the semivariogram
#'
#' @param vg a list that contains the estimated empirical semivariogram
#'
#' @return nothing
#' @export
plot.variogram <- function(vg){
  lag   <- vg %>% group_by(bins) %>% summarise(avg=mean(x_dist))
  gamma <- vg %>% group_by(bins) %>% summarise(avg=mean(y_dist))
  plot(lag$avg, gamma$avg, xlab = "lag (h)", ylab = "E[D(h)^2]", main = "semivariogram",
       pch =20, type = "l")
  n = vg %>% group_by(bins) %>% count() 
  text(lag$avg, gamma$avg, labels = n$n)
}



