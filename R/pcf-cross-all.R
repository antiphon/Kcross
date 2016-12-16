#' Quick cross 2nd order product density, rectangle translation border correction
#'
#'
#' @export
#' @import spatstat Rcpp
#' @useDynLib Kcross
rho_cross_all_box <- function(x, r, adjust = 1, correction = TRUE) {
  xy <- as.matrix(coords(x))
  types <-  as.integer(marks(x))
  ntypes <- length(unique(types))
  bb <- as.matrix( with(x$window, cbind(xrange, yrange)) )
  intensities <- intensity(x)
  # go
  V <- c_rho_cross_2d_box(xy, bb, types, intensities, r, adjust, correction)
  #
  # compile the pcf's
  G <- array(dim = c(ntypes, ntypes, length(r)))
  for(i in 1:length(r)){
    G[,,i] <- matrix(V[[i]], ncol = ntypes, nrow = ntypes) / 2
  }
  # done
  G
}

#' Quick cross 2nd pcf, rectangle translation border correction
#'
#'
#' @export
#' @import spatstat Rcpp
#' @useDynLib Kcross
pcf_cross_all_box <- function(x, r, ...) {
  G <- rho_cross_all_box(x, r, ...)
  intensities <- intensity(x)
  int2 <- outer(intensities, intensities , "*")
  G <- simplify2array(lapply(1:length(r), function(z) G[,,z]/(int2 * pi)))
  G
}

#' Quick univ 2nd order product density, rectangle translation border correction
#'
#' @param kern 1=epa, 2=box
#'
#' @export
#' @import spatstat Rcpp
#' @useDynLib Kcross
rho_box <- function(x, r, bw, adjust = 1, correction = TRUE, kern = 1){
  xy <- as.matrix(coords(x))
  bb <- as.matrix( with(x$window, cbind(xrange, yrange)) )
  int <- intensity(unmark(x))
  if(missing(bw)) bw <- adjust * 0.15 / sqrt(int)

  # go
  V <- c_rho_2d_box(xy, bb, int, r, bw, correction, kern)
  #
  V
}

#' Quick univ pcf, rectangle translation border correction
#'
#' @param x point pattern ppp
#' @param r range vector
#' @param ... passed on to rho_box
#' @param int_type type of intensity estimator to use.
#'
#' @details See int_type in intensity_adapted
#'
#' @export
#' @import spatstat Rcpp
#' @useDynLib Kcross

pcf_box <- function(x, r, ..., int_type = "n") {
  rho_box(x,r, ...)/(pi*intensity_adapted(unmark(x), r, type=int_type)^2)
}

####################################################
#' Quick mark connect, rectangle translation border correction
#'
#'
#' @export
#' @import spatstat Rcpp
#' @useDynLib Kcross
markconnect_all_box <- function(x, r, ...) {
  G <- rho_cross_all_box(x, r, ...)
  g <- rho_box(x, r, ...)
  ntypes <- length(unique(x$marks))
  for(i in 1:ntypes){
    for(j in i:ntypes){
      G[i,j,] <- G[i,j,]/g
      G[j,i,] <- G[j,i,]/g
    }
  }
  G
}

####################################################
# deltaij, 5.4.14 p 330 ipss
delta_all_box <- function(x, r, ...) {
  D <- G <- pcf_cross_all_box(x, r, ...)
  ints <- intensity(x)
  ntypes <- length(ints)
  for(i in 1:ntypes){
    for(j in i:ntypes){
      D[i,j,] <- ints[j] * G[i,j,] - ints[i]*G[i,i,]
      D[j,i,] <- ints[i] * G[j,i,] - ints[j]*G[j,j,]
    }
  }
  D
}

####################################################

# Check
if(0){
  library(spatstat)
  set.seed(1)
  N <- 1400
  M <- 5
  x <- setmarks(rpoint(N), factor(sample(1:M, N, T)))#
  r <- seq(0, .3, l = 20)
  t0 <- system.time( g <- pcfcross(x, i=1, r=r, j = 2, divisor ="d", correction ="trans") )
  t1 <- system.time( gg <- pcf_cross_all_box(x, r = r) )
  par(mfrow=c(1,1))
  plot(r, g$trans, "l", ylim =c(0,2), lwd=3)
  lines(r, gg[1,2,], col=3, lty=2)
  abline(h = 1, col=2)
  print(range(gg))
  print(rbind(t0,t1))
}

# univariate g
if(0){
  library(spatstat)
  set.seed(12)
  x <- rMatClust(100,.01, 10)
  r <- seq(0, .3, l = 20)
  t0 <- system.time( g <- pcf(x, r=r, divisor ="d", correction ="trans") )
  t1 <- system.time( gg <- pcf_box(x, r = r) )
  par(mfrow=c(1,1))
  plot(r, g$trans, "l", lwd=3)
  lines(r, gg, col=3, lty=2)
  abline(h = 1, col=2)
  print(rbind(t0,t1))
  # close enough
}



# delta. no reference
if(0){
  library(spatstat)
  set.seed(1)
  N <- 1300
  M <- 120
  x <- setmarks(rpoint(N), factor(sample(1:M, N, T)))#
  r <- seq(0, .3, l = 20)
  t1 <- system.time( gg <- delta_all_box(x, r = r) )
  par(mfrow=c(1,1))
  plot(r, gg[1,1,], "l", ylim = c(-2,2)*50)
  for(i in 2:M) lines(r, gg[1,i,])
  for(i in 2:M) lines(r, gg[i,1,], col=2)

}

