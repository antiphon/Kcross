#' Guan's KPSS test of stationarity
#'
#' Test of stationarity in a rectangle.
#'
#' @param x ppp
#' @param h variance estimation range parameter
#' @param nx x-step count in Rieman sums
#' @param alpha alpha level for test rejection. Only 0.01, 0.05, 0.1 work.
#'
#'
#' @import spatstat
#' @export

kpss.test <- function(x, h, nx = 100, alpha = 0.05) {

  l <- intensity(x)

  if(missing(h)) h <- 0.15/sqrt(l)

  bbox <- with(x$window, cbind(xrange, yrange))

  nn <- apply(bbox, 2, diff)
  ny <- round( nx * nn[2]/diff(bbox[,1]) )

  x0 <- bbox[1,1]
  y0 <- bbox[1,2]
  dxy <- nn/c(nx, ny)
  stepsx <- seq(bbox[1,1] + dxy[1], bbox[2,1], l = nx)
  stepsy <- seq(bbox[1,2] + dxy[2], bbox[2,2], l = ny)

  dv <- prod( nn / c(nx,ny))

  xy <- cbind(x$x, x$y)
  #
  #
  # The integral
  S <- matrix(0, nx, ny)
  for(i in 1:nx){
    for(j in 1:ny){
        inside <- xy[,1] < stepsx[i] & xy[,2] < stepsy[j]
        S[i,j] <- sum(inside) - (stepsx[i]-x0)*(stepsy[j]-y0)*l
    }
  }
  #
  # The scaling variance:
  #  doublesum <- 2 * c_guan_doublesum_2d_box(xy, bbox, bw = h)
  #  sigma2 <- doublesum - l^2 * pi * h^2 + l
  sigma2 <- l^2 * ( K_box(x, h) - pi * h^2 ) + l
  #
  #Compile
  T1 <- sum(S^2) * dv / prod(nn^2)
  Test <- T1 / sigma2
  #
  # return the image of the S(x,y)
  Sim <- im(t(S), stepsx-dxy[1]/2, stepsy-dxy[2]/2)
  #
  # quantiles from the paper
  quant <- c(0.5013, 0.3244, 0.2542)[alpha==c(0.01, 0.05, 0.1)]
  reject <- quant < Test
  # done
  list(statistic = Test, reject=reject, alpha = alpha,
       quant=quant, S=Sim, sigma2=sigma2, bw = h)
}
